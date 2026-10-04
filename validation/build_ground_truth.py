#!/usr/bin/env python3
"""
build_ground_truth.py — KEGG-derived ground truth for the methane reference panel.

For every genome of validation/panel.tsv that KEGG holds (column `kegg`), and every
target of config/targets.yaml that is anchored on KO numbers, the cell is `present`
when KEGG assigns any of the target's KOs to a gene of that genome (REST
`link/ko/<org>`), else `absent`. Nothing here comes from the pipeline's own output.

What this file is and is not:
  * It is KEGG's view. Where a KO does not separate two functions the cell follows
    the KO: an ammonia oxidizer is `present` for pmoA / pmoB / pmoC (K10944 / K10945 /
    K10946 are the shared methane / ammonia monooxygenase KOs). Those cells are
    corrected, with evidence, in curated_function_gt.tsv
    (build_curated_function_gt.py) — the reference the headline metrics use.
  * `mcrA_anme` has no KO of its own (it is K00399 called against ANME-clade seeds),
    so it is absent from this file and fully curated.
  * Genomes that are not in KEGG (11 of 49) have no rows here.

Output: validation/ground_truth.tsv  (genome, target, expected, source)
KEGG answers are cached in validation/.kegg_cache/ (delete to refresh).
"""
from __future__ import annotations

import csv
import sys
import time
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
ROSTER = ROOT / "validation" / "panel.tsv"
TARGETS = ROOT / "config" / "targets.yaml"
OUT = ROOT / "validation" / "ground_truth.tsv"
CACHE = ROOT / "validation" / ".kegg_cache"
SOURCE = "kegg_v1"

# Targets that share their KO with another target and carry no KO-level information
# of their own.
KO_LESS = {"mcrA_anme"}


def read_roster() -> list[dict]:
    with open(ROSTER) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def kegg_ko_set(org: str) -> set[str]:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{org}.ko"
    if not f.exists():
        print(f"  fetching KEGG KO set: {org}", file=sys.stderr)
        with urllib.request.urlopen(f"https://rest.kegg.jp/link/ko/{org}",
                                    timeout=90) as r:
            f.write_bytes(r.read())
        time.sleep(0.4)
    kos = {p[1][3:] for p in (l.split("\t") for l in f.read_text().splitlines())
           if len(p) == 2 and p[1].startswith("ko:")}
    if len(kos) < 500:
        sys.exit(f"error: KEGG returned only {len(kos)} KOs for '{org}' — wrong "
                 f"organism code or a failed download ({f})")
    return kos


def main() -> None:
    targets = yaml.safe_load(open(TARGETS))["targets"]
    genomes = [(r["name"], r["kegg"]) for r in read_roster() if r["kegg"] != "-"]

    rows, n_present = [], 0
    for genome, org in genomes:
        kos = kegg_ko_set(org)
        for t in targets:
            if t["id"] in KO_LESS or not t.get("ko"):
                continue
            present = any(k in kos for k in t["ko"])
            n_present += present
            rows.append((genome, t["id"], "present" if present else "absent", SOURCE))

    with open(OUT, "w") as fh:
        fh.write("genome\ttarget\texpected\tsource\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    print(f"[ground_truth] {len(rows)} cells, {len(genomes)} genomes "
          f"({n_present} present, {len(rows) - n_present} absent) -> {OUT}")


if __name__ == "__main__":
    main()
