#!/usr/bin/env python3
"""
trap_independence.py — audit the homology-trap claim for seed / HMM circularity.

The confirmatory endpoint of the benchmark is trap PRECISION. A panel genome whose
genus supplied a BLAST seed or an HMM training sequence for a trap target is not an
independent test of that target: the detector was built to find that protein.
Independence is a property of the (genome, target) cell, not of the genome —
*Methylococcus capsulatus* Bath is a pmoA and mmoX seed strain but an independent
test of xoxF.

Seed-sourced cells come from validation/seed_leakage.tsv (detect_seed_leakage.py:
same genus as a seed of that target, in any detector source). This script lists
the ground-truth-present trap cells with their tag and recomputes the trap metrics
on the independent cells only — the non-circular number, and row 6 of the gate.

    python validation/trap_independence.py        # after detect_seed_leakage.py
"""
from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEAK = HERE / "seed_leakage.tsv"

_cache: set[tuple[str, str]] | None = None


def seed_cells() -> set[tuple[str, str]]:
    global _cache
    if _cache is None:
        if not LEAK.exists():
            raise SystemExit(f"error: {LEAK} missing — run detect_seed_leakage.py")
        with open(LEAK) as fh:
            _cache = {(r["genome"], r["target"])
                      for r in csv.DictReader(fh, delimiter="\t")}
    return _cache


def is_seed(genome: str, target: str) -> bool:
    return (genome, target) in seed_cells()


def main() -> None:
    import score_mcycle as S
    truth = S.load_truth(); calls = S.load_calls({g for g, _ in truth})
    cells = [(g, t) for (g, t) in truth if t in S.TRAP and g in calls]
    present = sorted(c for c in cells if truth[c] == "present")
    print(f"=== {len(present)} PRESENT trap cells (ground truth), tagged ===")
    for g, t in present:
        tag = "SEED-SOURCED" if is_seed(g, t) else "independent"
        pred = "present" if calls[g].get(t) in S.PRESENT else "ABSENT"
        print(f"  {t:<10}{g:<28}pred={pred:<9}[{tag}]")
    print("\n=== trap metrics ===")
    for label, sub in (("ALL trap", cells),
                       ("INDEPENDENT-only trap", [c for c in cells if not is_seed(*c)])):
        tp = fp = fn = 0
        for g, t in sub:
            exp = truth[(g, t)] == "present"; pred = calls[g].get(t) in S.PRESENT
            tp += pred and exp; fp += pred and not exp; fn += (not pred) and exp
        p, r, f = S.prf(tp, fp, fn)
        print(f"  {label:<24}cells={len(sub)} present="
              f"{sum(truth[c] == 'present' for c in sub)}  TP={tp} FP={fp} FN={fn}  "
              f"P={p:.3f} R={r:.3f} F1={f:.3f}")


if __name__ == "__main__":
    main()
