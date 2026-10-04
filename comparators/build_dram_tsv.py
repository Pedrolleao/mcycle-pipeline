#!/usr/bin/env python3
"""
build_dram_tsv.py — normalize the DRAM v1.4.6 distillate to the table the benchmark
reads (columns: genome, function, present), as fixed in validation/benchmark/prereg.md.

Source: dram_out/b??/distilled/product.tsv (one per batch), the columns of the category
"Methanogenesis and methanotrophy". `function` is the column name without the
category prefix; `present` is DRAM's own True / False.

    python comparators/build_dram_tsv.py [--product ...] [--out ...]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PREFIX = "Methanogenesis and methanotrophy: "


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--product", type=Path, nargs="*",
                    default=sorted((HERE / "dram_out").glob("b??/distilled/product.tsv")))
    ap.add_argument("--out", type=Path, default=ROOT / "validation" / "benchmark" / "dram.tsv")
    args = ap.parse_args()
    rows = []
    for product in args.product:
        with open(product) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                for col, val in r.items():
                    if col.startswith(PREFIX):
                        rows.append((r["genome"], col[len(PREFIX):],
                                     "1" if str(val).strip().lower() == "true" else "0"))
    with open(args.out, "w") as fh:
        fh.write("genome\tfunction\tpresent\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    print(f"[dram] {len({r[0] for r in rows})} genomes, {len({r[1] for r in rows})} functions "
          f"-> {args.out}")
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
