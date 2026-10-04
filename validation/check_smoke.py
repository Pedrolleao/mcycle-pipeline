#!/usr/bin/env python3
"""
check_smoke.py — compare the pipeline's calls on the smoke panel with
validation/smoke_expectations.tsv and exit non-zero on any mismatch.

This is a SMOKE regression (the pipeline still gets the textbook organisms
right), not an accuracy estimate: 18 genomes, expectations for the diagnostic
calls only, thresholds tuned while looking at these same genomes.

Usage:  python validation/check_smoke.py [--results results]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRESENT = {"confirmed", "domain-only"}


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=HERE.parent / "results")
    ap.add_argument("--expect", type=Path, default=HERE / "smoke_expectations.tsv")
    args = ap.parse_args()

    failures, n = [], 0
    cache: dict[str, tuple[dict, dict]] = {}
    for e in read_tsv(args.expect):
        g = e["genome"]
        if g not in cache:
            cdir = args.results / g / "calls"
            if not (cdir / "mcycle_calls.tsv").exists():
                sys.exit(f"FAIL: no calls for {g} under {args.results} — run the "
                         "pipeline on ../test_panel first (make smoke)")
            cache[g] = (
                {r["target_id"]: r for r in read_tsv(cdir / "mcycle_calls.tsv")},
                {r["synergy_id"]: r for r in read_tsv(cdir / "synergy_completeness.tsv")})
        calls, mods = cache[g]
        kind, ident, want = e["kind"], e["id"], e["expected"]
        if kind == "direction":
            m = re.search(r"mcr_(methanogenic|reverse)", calls["mcrA"]["evidence_source"])
            checks = [("mcr direction", m.group(1) if m else "none")]
        elif kind == "module":
            ids = list(mods) if ident == "*" else [ident]
            checks = [(i, "complete" if mods[i]["status"] == "complete"
                       else "not_complete") for i in ids]
        else:
            got = calls[ident]["status"]
            checks = [(ident, "present" if want == "present" and got in PRESENT else got)]
        for label, got in checks:
            n += 1
            if got != want:
                failures.append(f"  {g}: {kind} {label} — expected {want}, got {got}")

    if failures:
        print(f"FAIL: {len(failures)} of {n} smoke expectations not met")
        print("\n".join(failures))
        sys.exit(1)
    print(f"OK: {n} smoke expectations met on {len(cache)} genomes")


if __name__ == "__main__":
    main()
