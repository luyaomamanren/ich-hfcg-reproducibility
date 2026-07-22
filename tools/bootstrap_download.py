"""Download immutable public GEO inputs when R/libcurl lacks a usable TLS backend.

The analytical workflow remains R-based. This helper only transfers bytes, records
checksums, and never transforms scientific data.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import re
import time
import urllib.request
from pathlib import Path


BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
SERIES = {
    "GSE24265": {
        "prefix": "GSE24nnn",
        "files": ["GSE24265_RAW.tar"],
    },
    "GSE163256": {
        "prefix": "GSE163nnn",
        "files": [
            "GSE163256_gene_lengths.txt.gz",
            "GSE163256_monos_counts.csv.gz",
            "GSE163256_neuts_counts.csv.gz",
        ],
    },
    "GSE266873": {"prefix": "GSE266nnn", "files": ["GSE266873_RAW.tar"]},
    "GSE166638": {"prefix": "GSE166nnn", "files": ["GSE166638_RAW.tar"]},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def list_supplementary(accession: str, prefix: str) -> list[str]:
    url = f"{BASE}/{prefix}/{accession}/suppl/"
    html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", errors="replace")
    return [x for x in re.findall(r'href="([^"]+)"', html) if not x.startswith("?")]


def remote_size(url: str) -> int:
    request = urllib.request.Request(url, headers={"Range": "bytes=0-0", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=180) as response:
        content_range = response.headers.get("Content-Range")
        if content_range:
            return int(content_range.rsplit("/", 1)[1])
        return int(response.headers["Content-Length"])


def download_range(url: str, start: int, end: int, target: Path) -> None:
    if target.exists() and target.stat().st_size == end - start + 1:
        return
    request = urllib.request.Request(
        url,
        headers={"Range": f"bytes={start}-{end}", "User-Agent": "Mozilla/5.0"},
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=240) as response, target.open("wb") as out:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    out.write(block)
            if target.stat().st_size != end - start + 1:
                raise IOError(f"short range {start}-{end}: {target.stat().st_size} bytes")
            return
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)


def download(url: str, destination: Path, workers: int = 24, chunk_mb: int = 4) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return
    size = remote_size(url)
    chunk_size = chunk_mb * 1024 * 1024
    parts_dir = destination.parent / (destination.name + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    ranges = [(start, min(start + chunk_size - 1, size - 1)) for start in range(0, size, chunk_size)]
    jobs = [(start, end, parts_dir / f"{start:012d}-{end:012d}.part") for start, end in ranges]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(download_range, url, start, end, target) for start, end, target in jobs]
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            future.result()
            print(f"  completed chunk {index}/{len(futures)}", flush=True)
    assembled = destination.with_suffix(destination.suffix + ".part")
    with assembled.open("wb") as out:
        for _, _, target in jobs:
            with target.open("rb") as handle:
                while True:
                    block = handle.read(1024 * 1024)
                    if not block:
                        break
                    out.write(block)
    if assembled.stat().st_size != size:
        raise IOError(f"assembled size mismatch for {destination}")
    assembled.replace(destination)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--list-only", action="store_true")
    parser.add_argument("--series", nargs="*", default=["GSE24265", "GSE163256", "GSE266873"])
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--chunk-mb", type=int, default=1)
    args = parser.parse_args()

    selected = {k: SERIES[k] for k in args.series}
    manifest = []
    for accession, spec in selected.items():
        names = spec.get("files") or [
            name for name in list_supplementary(accession, spec["prefix"])
            if re.search(spec["pattern"], name)
        ]
        if args.list_only:
            print(accession, *names, sep="\n  ")
            continue
        for name in names:
            url = f"{BASE}/{spec['prefix']}/{accession}/suppl/{name}"
            target = args.out / accession / name
            print(f"Downloading {url}", flush=True)
            download(url, target, workers=args.workers, chunk_mb=args.chunk_mb)
            manifest.append({
                "accession": accession,
                "file": target.relative_to(args.out).as_posix(),
                "url": url,
                "bytes": target.stat().st_size,
                "sha256": sha256(target),
            })
        soft_name = f"{accession}_family.soft.gz"
        soft_url = f"{BASE}/{spec['prefix']}/{accession}/soft/{soft_name}"
        soft_target = args.out / accession / soft_name
        print(f"Downloading {soft_url}", flush=True)
        download(soft_url, soft_target, workers=args.workers, chunk_mb=args.chunk_mb)
        manifest.append({
            "accession": accession,
            "file": soft_target.relative_to(args.out).as_posix(),
            "url": soft_url,
            "bytes": soft_target.stat().st_size,
            "sha256": sha256(soft_target),
        })
    if not args.list_only:
        (args.out / "download_manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )


if __name__ == "__main__":
    main()
