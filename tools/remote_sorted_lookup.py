"""Look up exact first-column keys in a large remote, lexicographically sorted TSV.

Uses HTTP byte ranges so a multi-gigabyte public summary-statistics file need not
be downloaded in full. This is data acquisition only; statistical analysis is R.
"""
from __future__ import annotations

import argparse
import urllib.request


def request_range(url: str, start: int, end: int) -> bytes:
    req = urllib.request.Request(
        url,
        headers={"Range": f"bytes={start}-{end}", "User-Agent": "Mozilla/5.0"},
    )
    with urllib.request.urlopen(req, timeout=120) as response:
        return response.read()


def size_of(url: str) -> int:
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as response:
        return int(response.headers["Content-Range"].rsplit("/", 1)[1])


def first_complete_line(data: bytes) -> str:
    parts = data.splitlines()
    if len(parts) < 2:
        return ""
    return parts[1].decode("utf-8", errors="replace")


def lookup(url: str, key: str, block: int = 131072) -> str | None:
    low, high = 0, size_of(url) - 1
    while high - low > block:
        mid = (low + high) // 2
        line = first_complete_line(request_range(url, mid, min(mid + block, high)))
        if not line:
            high = mid
            continue
        observed = line.split("\t", 1)[0]
        if observed < key:
            low = mid
        else:
            high = mid
    start = max(0, low - block)
    data = request_range(url, start, min(size_of(url) - 1, high + 2 * block))
    for raw in data.splitlines():
        line = raw.decode("utf-8", errors="replace")
        if line.split("\t", 1)[0] == key:
            return line
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("keys", nargs="+")
    args = parser.parse_args()
    header = request_range(args.url, 0, 65535).splitlines()[0].decode()
    print(header)
    for key in args.keys:
        line = lookup(args.url, key)
        if line:
            print(line)


if __name__ == "__main__":
    main()
