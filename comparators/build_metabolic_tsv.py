#!/usr/bin/env python3
"""
build_metabolic_tsv.py — normalize METABOLIC v4.0 output to the table the benchmark
reads (columns: genome, kind, id), as fixed in validation/benchmark/prereg.md.

  kind = ko    a KO of KEGG_identifier_result/<genome>.result.txt with a hit
               (METABOLIC's own KOfam scan) — mapped to targets by KO;
  kind = gene  a row of worksheet 1 (METABOLIC's curated function table) reported
               `Present`, for the gene-named rows that decide a target on their own
               (adapters.METABOLIC_GENE_ROWS: pmoA, pmoB, pmoC, mcrA, ... — METABOLIC
               ships separate pmoA and amoA models under one KO).

    python comparators/build_metabolic_tsv.py \
        [--kegg-dir comparators/metabolic_out/kegg_all] \
        [--worksheets comparators/metabolic_out/worksheet1/*.tsv] [--out ...]
"""
from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "validation" / "benchmark"))
from adapters import METABOLIC_GENE_ROWS  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kegg-dir", type=Path, default=HERE / "metabolic_out" / "kegg_all")
    ap.add_argument("--worksheets", nargs="*",
                    default=sorted(glob.glob(str(HERE / "metabolic_out" / "worksheet1" / "*.tsv"))))
    ap.add_argument("--out", type=Path, default=ROOT / "validation" / "benchmark" / "metabolic.tsv")
    args = ap.parse_args()

    rows, genomes = [], set()
    for fp in sorted(glob.glob(str(args.kegg_dir / "*.result.txt"))):
        g = Path(fp).name[: -len(".result.txt")]
        genomes.add(g)
        for line in open(fp):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].startswith("K") and p[1].strip():
                rows.append((g, "ko", p[0]))
    n_gene = 0
    for ws in args.worksheets:
        with open(ws) as fh:
            rd = csv.reader(fh, delimiter="\t")
            head = next(rd)
            cols = {h[: -len(" Hmm presence")]: i for i, h in enumerate(head)
                    if h.endswith(" Hmm presence")}
            gene_col = head.index("Gene abbreviation")
            for r in rd:
                if r[gene_col].strip() not in METABOLIC_GENE_ROWS:
                    continue
                for g, i in cols.items():
                    if i < len(r) and r[i].strip() == "Present":
                        rows.append((g, "gene", r[gene_col].strip())); n_gene += 1
    rows = sorted(set(rows))
    with open(args.out, "w") as fh:
        fh.write("genome\tkind\tid\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    print(f"[metabolic] {len(genomes)} genomes, {sum(r[1] == 'ko' for r in rows)} KO rows, "
          f"{sum(r[1] == 'gene' for r in rows)} gene rows -> {args.out}")
    return 0 if genomes else 1


if __name__ == "__main__":
    sys.exit(main())
