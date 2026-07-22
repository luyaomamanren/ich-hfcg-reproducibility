"""Fetch public eQTL and regional ICH summary statistics without analysis.

The eQTL Catalogue API and OSF file both support HTTPS, while this Windows R
installation does not currently have a working TLS backend.  This helper is
therefore deliberately limited to byte transfer and TSV/JSON parsing.  All
filtering, harmonisation, MR and colocalisation are performed in R.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


EQTL_API = "https://www.ebi.ac.uk/eqtl/api/v3/associations"
ICH_URL = "https://osf.io/download/axy7g/"


def get(url: str, start: int | None = None, end: int | None = None) -> bytes:
    headers = {"User-Agent": "Mozilla/5.0 (reproducible academic workflow)"}
    if start is not None:
        headers["Range"] = f"bytes={start}-{end}"
    request = urllib.request.Request(url, headers=headers)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=240) as response:
                return response.read()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def remote_size(url: str) -> int:
    request = urllib.request.Request(
        url, headers={"Range": "bytes=0-0", "User-Agent": "Mozilla/5.0"}
    )
    with urllib.request.urlopen(request, timeout=240) as response:
        return int(response.headers["Content-Range"].rsplit("/", 1)[1])


def line_at_or_after(url: str, offset: int, size: int, block: int = 262144):
    start = max(0, min(offset, size - 1))
    raw = get(url, start, min(size - 1, start + block - 1))
    if start == 0:
        line_start = 0
    else:
        first_nl = raw.find(b"\n")
        if first_nl < 0:
            return None
        line_start = start + first_nl + 1
        raw = raw[first_nl + 1 :]
    end = raw.find(b"\n")
    if end < 0:
        return None
    line = raw[:end].decode("utf-8", errors="replace")
    return line_start, line


def lower_bound(url: str, key: str, size: int) -> int:
    low, high = 0, size - 1
    while high - low > 262144:
        mid = (low + high) // 2
        found = line_at_or_after(url, mid, size)
        if found is None:
            high = mid
            continue
        position, line = found
        observed = line.split("\t", 1)[0]
        if observed < key:
            low = position + len(line.encode("utf-8")) + 1
        else:
            high = position
    scan_start = max(0, low - 262144)
    raw = get(url, scan_start, min(size - 1, high + 524288))
    cursor = scan_start
    for part in raw.splitlines(keepends=True):
        line = part.decode("utf-8", errors="replace").rstrip("\r\n")
        if cursor > 0 and cursor == scan_start:
            cursor += len(part)
            continue
        if line.split("\t", 1)[0] >= key:
            return cursor
        cursor += len(part)
    return high


def fetch_eqtl(gene_id: str, target: Path) -> None:
    records: list[dict] = []
    start, page_size = 0, 10000
    while True:
        query = urllib.parse.urlencode(
            {"size": page_size, "start": start, "dataset_id": "QTD000356", "gene_id": gene_id}
        )
        payload = json.loads(get(f"{EQTL_API}?{query}"))
        if isinstance(payload, list):
            page = payload
        else:
            page = payload.get("_embedded", {}).get("associations", [])
        records.extend(page)
        if len(page) < page_size:
            break
        start += page_size
    target.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in records for key in row})
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def fetch_ich_region(chromosome: str, start: int, end: int, target: Path) -> None:
    size = remote_size(ICH_URL)
    lower = f"{chromosome}:{start:09d}"
    upper = f"{chromosome}:{end + 1:09d}"
    byte_start = lower_bound(ICH_URL, lower, size)
    byte_end = lower_bound(ICH_URL, upper, size)
    raw = get(ICH_URL, byte_start, min(size - 1, byte_end + 262144))
    header = get(ICH_URL, 0, 65535).splitlines()[0].decode("utf-8")
    selected = []
    for item in raw.splitlines():
        line = item.decode("utf-8", errors="replace")
        key = line.split("\t", 1)[0]
        if not key.startswith(f"{chromosome}:"):
            continue
        try:
            position = int(key.split(":", 1)[1].split("_", 1)[0])
        except ValueError:
            continue
        if start <= position <= end:
            selected.append(line)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(header + "\n" + "\n".join(selected) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--window", type=int, default=1_000_000)
    parser.add_argument("--genes", nargs="*")
    parser.add_argument("--eqtl-only", action="store_true")
    parser.add_argument("--ich-only", action="store_true")
    args = parser.parse_args()
    with args.config.open(encoding="utf-8") as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))
    if args.genes:
        genes = [row for row in genes if row["gene"] in set(args.genes)]
    manifest = []
    for gene in genes:
        symbol = gene["gene"]
        eqtl_target = args.out / "eqtl" / f"QTD000356_{symbol}.tsv"
        ich_target = args.out / "ich" / f"Meta_ICH_{symbol}_region.tsv"
        if not args.ich_only and not eqtl_target.exists():
            fetch_eqtl(gene["ensembl_gene_id"], eqtl_target)
        start = max(1, int(gene["start_hg19"]) - args.window)
        end = int(gene["end_hg19"]) + args.window
        if not args.eqtl_only and not ich_target.exists():
            fetch_ich_region(gene["chr_hg19"], start, end, ich_target)
        manifest.append({"gene": symbol, "eqtl": str(eqtl_target), "ich": str(ich_target)})
        print(f"{symbol}: eQTL and ICH region downloaded", flush=True)
    (args.out / "genetic_input_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
