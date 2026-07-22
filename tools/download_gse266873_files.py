#!/usr/bin/env python3
"""Download the nine GSE266873 raw 10x matrices directly from GEO.

Each file is written to a temporary .part path and atomically renamed after the
remote Content-Length has been reached. Existing partial files are resumed.
"""
from __future__ import annotations

import concurrent.futures
import re
import time
import urllib.request
from pathlib import Path

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8255nnn/{gsm}/suppl/"
OUT = Path(__file__).resolve().parents[1] / "data" / "raw" / "GSE266873" / "10x"


def listing(gsm: str) -> list[str]:
    url = BASE.format(gsm=gsm)
    html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8")
    return [x for x in re.findall(r'href="([^"]+)"', html)
            if re.search(r"_(barcodes\.tsv|features\.tsv|matrix\.mtx)\.gz$", x)]


def download(item: tuple[str, str]) -> str:
    gsm, name = item
    url = BASE.format(gsm=gsm) + name
    target = OUT / gsm / name
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(target.suffix + ".part")
    for attempt in range(20):
        try:
            start = part.stat().st_size if part.exists() else 0
            request = urllib.request.Request(url, headers={"Range": f"bytes={start}-"})
            with urllib.request.urlopen(request, timeout=180) as response:
                mode = "ab" if start and response.status == 206 else "wb"
                with part.open(mode) as handle:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
            remote = int(urllib.request.urlopen(
                urllib.request.Request(url, method="HEAD"), timeout=60
            ).headers["Content-Length"])
            if part.stat().st_size == remote:
                part.replace(target)
                return f"OK {gsm}/{name} {remote}"
            raise IOError(f"size {part.stat().st_size} != {remote}")
        except Exception as exc:
            if attempt == 19:
                raise
            time.sleep(min(60, 2 ** min(attempt, 5)))
    raise RuntimeError("unreachable")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    samples = [f"GSM{x}" for x in range(8255340, 8255349)]
    jobs = [(gsm, name) for gsm in samples for name in listing(gsm)]
    print(f"Downloading {len(jobs)} files to {OUT}", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
        for message in pool.map(download, jobs):
            print(message, flush=True)


if __name__ == "__main__":
    main()
