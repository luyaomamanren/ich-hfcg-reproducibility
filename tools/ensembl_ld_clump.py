"""LD-clump candidate instruments against 1000 Genomes phase 3 Europeans.

Input is a TSV prepared by R with gene, rsid and pvalue columns.  The output
contains every candidate, its retained lead SNP (if excluded), and the observed
pairwise r2. Ensembl's public endpoint permits a maximum 500-kb window; the R
workflow additionally enforces the prespecified cis-region boundary.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


BASE = "https://grch37.rest.ensembl.org"
POPULATION = "1000GENOMES:phase_3:EUR"


def ld_for(rsid: str, threshold: float) -> dict[str, float]:
    population = urllib.parse.quote(POPULATION, safe="")
    url = (
        f"{BASE}/ld/human/{urllib.parse.quote(rsid)}/{population}"
        f"?r2={threshold};window_size=500"
    )
    request = urllib.request.Request(
        url, headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"}
    )
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                rows = json.loads(response.read())
            result = {}
            for row in rows:
                other = row["variation2"] if row["variation1"] == rsid else row["variation1"]
                result[other] = float(row["r2"])
            return result
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)
    return {}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--r2", type=float, default=0.001)
    args = parser.parse_args()
    with args.input.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output = []
    for gene in sorted({row["gene"] for row in rows}):
        candidates = sorted(
            [row for row in rows if row["gene"] == gene], key=lambda x: float(x["pvalue"])
        )
        excluded: dict[str, tuple[str, float]] = {}
        for lead in candidates:
            rsid = lead["rsid"]
            if rsid in excluded:
                continue
            output.append({**lead, "keep": "TRUE", "lead_rsid": rsid, "r2": "1"})
            linked = ld_for(rsid, args.r2)
            for candidate in candidates:
                other = candidate["rsid"]
                if other != rsid and other not in excluded and other in linked:
                    excluded[other] = (rsid, linked[other])
        retained = {row["rsid"] for row in output if row["gene"] == gene and row["keep"] == "TRUE"}
        for row in candidates:
            if row["rsid"] in retained:
                continue
            lead, r2 = excluded.get(row["rsid"], ("not_in_1000G_panel", float("nan")))
            output.append({**row, "keep": "FALSE", "lead_rsid": lead, "r2": r2})
        print(f"{gene}: retained {len(retained)}/{len(candidates)}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) + ["keep", "lead_rsid", "r2"] if rows else []
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(output)


if __name__ == "__main__":
    main()
