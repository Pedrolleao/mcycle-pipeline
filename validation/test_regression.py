#!/usr/bin/env python3
"""
test_regression.py — accuracy regression gate for mcycle-pipeline.

Fails (exit 1) when the accuracy on the 49-genome reference panel drops below the
floors locked in here. Same eight check rows as the nitrogen and sulfur gates:
scored cells, ALL-F1 CI-lo, ALL-precision CI-lo, ALL-FP cap, trap-precision CI-lo,
trap-independence-precision CI-lo, hold-out-F1 CI-lo, per-pathway-F1 CI-lo (7 pathways).
CI = genome-cluster percentile bootstrap, B = 10,000, seed 1234 (score_mcycle.py).

Run after the pipeline has processed ../ref_panel into results_ref/ (`make
regression` does both; `make regression-score` re-scores existing results).
Default ground truth: curated_function_gt.tsv (GT_FILE=ground_truth.tsv gates on the
KEGG contrast; its trap rows are expected to fail — KEGG assigns the pmo KOs to
ammonia oxidizers).

Floors = CI lower bound observed on 2026-10-04 (first scoring of the hold-out, tool
frozen at commit 9c39500) minus a slack, set once:

  check                         observed               floor   slack
  scored cells                  3173                   2800    guards a partial run
  ALL micro-F1                  0.973 [0.963, 0.981]   0.94    0.02
  ALL precision                 0.976 [0.965, 0.985]   0.94    0.025
  ALL false positives           21                     <= 26   point check
  trap precision                1.000 [1.000, 1.000]   0.95    0.05
  trap-independence precision   1.000 [1.000, 1.000]   0.95    0.05
  hold-out micro-F1             0.961 [0.940, 0.976]   0.91    0.03
  per-pathway F1 (lowest)       0.920 [0.857, 0.972]   0.80    0.06  methylotrophic_methanogenesis

Usage:
  python validation/test_regression.py     # PASS / FAIL table, exit 0 / 1
  pytest validation/test_regression.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_mcycle import compute_metrics  # noqa: E402

MIN_SCORED_CELLS = 2800
MIN_ALL_F1_CI_LO = 0.94
MIN_ALL_PRECISION_CI_LO = 0.94
MAX_ALL_FP = 26
MIN_TRAP_PRECISION_CI_LO = 0.95
MIN_INDEP_TRAP_PRECISION_CI_LO = 0.95
MIN_HOLDOUT_F1_CI_LO = 0.91
MIN_PATHWAY_F1_CI_LO = 0.80


def _ci_lo(metric: dict, key: str) -> float:
    ci = metric.get(f"{key}_ci")
    return ci[0] if ci is not None else float("nan")


def _ci_str(metric: dict, key: str) -> str:
    p, ci = metric[key], metric.get(f"{key}_ci")
    if ci is None or math.isnan(ci[0]) or math.isnan(ci[1]):
        return f"{p:.3f}"
    return f"{p:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]"


def _checks(M: dict) -> list[tuple[str, str, bool]]:
    a, trap = M["aggregates"]["ALL"], M["aggregates"]["trap"]
    itrap, hold = M["indep_trap"], M["split"]["holdout"]
    n = M["scored_cells"]
    checks = [
        (f"scored cells >= {MIN_SCORED_CELLS}", str(n), n >= MIN_SCORED_CELLS),
        (f"ALL micro-F1 CI lo >= {MIN_ALL_F1_CI_LO}",
         _ci_str(a, "f1"), _ci_lo(a, "f1") >= MIN_ALL_F1_CI_LO),
        (f"ALL precision CI lo >= {MIN_ALL_PRECISION_CI_LO}",
         _ci_str(a, "precision"), _ci_lo(a, "precision") >= MIN_ALL_PRECISION_CI_LO),
        (f"ALL false positives <= {MAX_ALL_FP}", str(a["FP"]), a["FP"] <= MAX_ALL_FP),
        (f"homology-trap precision CI lo >= {MIN_TRAP_PRECISION_CI_LO}",
         _ci_str(trap, "precision"), _ci_lo(trap, "precision") >= MIN_TRAP_PRECISION_CI_LO),
        (f"trap-independence precision CI lo >= {MIN_INDEP_TRAP_PRECISION_CI_LO}",
         _ci_str(itrap, "precision"),
         _ci_lo(itrap, "precision") >= MIN_INDEP_TRAP_PRECISION_CI_LO),
        (f"hold-out micro-F1 CI lo >= {MIN_HOLDOUT_F1_CI_LO}",
         _ci_str(hold, "f1"), _ci_lo(hold, "f1") >= MIN_HOLDOUT_F1_CI_LO),
    ]
    for cat, pp in sorted(M["per_pathway"].items()):
        checks.append((f"pathway {cat} F1 CI lo >= {MIN_PATHWAY_F1_CI_LO}",
                       _ci_str(pp, "f1"), _ci_lo(pp, "f1") >= MIN_PATHWAY_F1_CI_LO))
    return checks


def test_regression():
    failures = [name for name, _, ok in _checks(compute_metrics(require_all=True)) if not ok]
    assert not failures, "accuracy regression: " + "; ".join(failures)


def main() -> int:
    results = _checks(compute_metrics(require_all=True))
    width = max(len(name) for name, _, _ in results)
    print("mcycle-pipeline accuracy regression gate (49-genome reference panel)\n")
    n_fail = 0
    for name, value, ok in results:
        n_fail += not ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:<{width}}  observed={value}")
    print()
    if n_fail:
        print(f"REGRESSION: {n_fail} check(s) failed — accuracy dropped below its floor.")
        return 1
    print(f"OK: all {len(results)} checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
