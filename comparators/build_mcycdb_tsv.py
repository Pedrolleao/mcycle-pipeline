#!/usr/bin/env python3
"""
build_mcycdb_tsv.py — run MCycDB (Qian et al. 2022) on the panel proteomes and export
the normalized table the benchmark reads (columns: genome, family, count).

MCycDB is the methane-cycle domain database of the NCycDB / SCycDB series
(github.com/qichao1984/MCycDB; 2021 release, repository commit ceba218). Faithful to
its own profiler (MCycDB_FunctionProfiler.PL), protein mode, default settings:
    diamond makedb --in MCycDB_2021.faa --db MCycDB_2021                 (once)
    diamond blastp -k 1 -e 1e-4 -d MCycDB_2021 -q <proteome> -o <hits>
The best-hit subject is mapped to a gene family through id2gene.map; subjects that
are not in the map are the homologue decoys of the database and count for nothing.
`count` = number of query proteins whose best hit is in the family.

The database is NOT redistributed: stage the repository files under comparators/MCyc/
(MCycDB_2021.faa is shipped as a split zip; the stream was inflated and its CRC
checked against the archive's central directory, 753c621e, 380,449,265 bytes,
923,871 sequences). Per-genome hits are cached under comparators/mcycdb_out/.

    python comparators/build_mcycdb_tsv.py [--panel comparators/panel_faa] [--out ...]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MCYC = HERE / "MCyc"
EVALUE = "0.0001"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, default=HERE / "panel_faa")
    ap.add_argument("--hits", type=Path, default=HERE / "mcycdb_out")
    ap.add_argument("--out", type=Path, default=ROOT / "validation" / "benchmark" / "mcycdb.tsv")
    ap.add_argument("--threads", default="8")
    args = ap.parse_args()

    db = MCYC / "MCycDB_2021"
    if not db.with_suffix(".dmnd").exists():
        subprocess.run(["diamond", "makedb", "--in", str(MCYC / "MCycDB_2021.faa"),
                        "--db", str(db), "--quiet"], check=True)
    id2gene = {}
    for line in open(MCYC / "id2gene.map"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            id2gene[p[0]] = p[1]

    args.hits.mkdir(parents=True, exist_ok=True)
    rows = []
    for faa in sorted(args.panel.glob("*.faa")):
        g = faa.stem
        hits = args.hits / f"{g}.diamond"
        if not hits.exists():
            subprocess.run(["diamond", "blastp", "-k", "1", "-e", EVALUE, "-p", args.threads,
                            "-d", str(db), "-q", str(faa), "-o", str(hits), "--quiet"],
                           check=True)
        counts: dict[str, int] = defaultdict(int)
        n_hit = n_family = 0
        for line in open(hits):
            n_hit += 1
            fam = id2gene.get(line.split("\t")[1])
            if fam:
                counts[fam] += 1; n_family += 1
        for fam, n in sorted(counts.items()):
            rows.append({"genome": g, "family": fam, "count": n})
        print(f"[mcycdb] {g}: {n_hit} proteins with a hit, {n_family} in a gene family, "
              f"{len(counts)} families", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "family", "count"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"[mcycdb] {len(rows)} rows, {len({r['genome'] for r in rows})} genomes -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
