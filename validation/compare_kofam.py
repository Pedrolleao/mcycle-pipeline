#!/usr/bin/env python3
"""
compare_kofam.py — the pipeline against RAW KofamScan-style KO assignment, per target.

"Raw KofamScan": a target is present in a genome when any of its KOs has an hmmscan
hit whose full-sequence score reaches the stock KOfam `ko_list` threshold — no clade
models, no rules, no threshold overrides, and a shared KO counts for every target
that lists it (K00399 calls mcrA AND mcrA_anme; K10944 calls pmoA in an ammonia
oxidizer). That is what mapping KofamScan output onto a pathway does.

Both are scored on the ground-truth cells (GT_FILE, default curated_function_gt.tsv).
The statistics-backed comparison is validation/benchmark/ (M5); this is the quick
per-target table.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_mcycle as S  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KO_LIST = ROOT / "resources" / ".cache" / "ko_list"


def stock_thresholds() -> dict[str, float]:
    out = {}
    for line in open(KO_LIST):
        p = line.rstrip("\n").split("\t")
        try:
            out[p[0]] = float(p[1])
        except (ValueError, IndexError):
            pass
    return out


def raw_kofam(genomes, targets, thr) -> dict[tuple[str, str], bool]:
    calls = {}
    for g in genomes:
        best: dict[str, float] = {}
        f = S.RESULTS / g / "hmm" / f"{g}.hmmscan.tsv"
        if not f.exists():
            continue
        for line in open(f):
            if not line.startswith("#"):
                p = line.split()
                best[p[0]] = max(best.get(p[0], 0.0), float(p[7]))
        for t in targets:
            calls[(g, t["id"])] = any(ko in thr and best.get(ko, 0.0) >= thr[ko]
                                      for ko in t.get("ko") or [])
    return calls


def tally(truth, pred, keep) -> tuple[int, int, int]:
    tp = fp = fn = 0
    for (g, t), exp in truth.items():
        if (g, t) not in pred or not keep(g, t):
            continue
        e, p = exp == "present", pred[(g, t)]
        tp += p and e; fp += p and not e; fn += (not p) and e
    return tp, fp, fn


def main() -> None:
    truth = S.load_truth()
    genomes = sorted({g for g, _ in truth})
    targets = yaml.safe_load(open(ROOT / "config" / "targets.yaml"))["targets"]
    calls = S.load_calls(set(genomes))
    tool = {(g, t["id"]): calls[g].get(t["id"]) in S.PRESENT
            for g in calls for t in targets}
    raw = raw_kofam(genomes, targets, stock_thresholds())
    print(f"ground truth: {S.GT.name}   scope: {S.SCOPE}\n")
    print(f"{'scope':<26}{'mcycle  TP  FP  FN    P     R     F1':<42}raw KOfam  TP  FP  FN    P     R     F1")
    scopes = [("ALL targets", lambda g, t: True),
              ("homology-trap targets", lambda g, t: t in S.TRAP),
              ("non-trap targets", lambda g, t: t not in S.TRAP),
              ("hold-out, ALL", lambda g, t: g in S.HOLDOUT),
              ("hold-out, trap", lambda g, t: g in S.HOLDOUT and t in S.TRAP)]
    for label, keep in scopes:
        a, b = tally(truth, tool, keep), tally(truth, raw, keep)
        pa, pb = S.prf(*a), S.prf(*b)
        print(f"{label:<26}       {a[0]:>4}{a[1]:>4}{a[2]:>4}  {S.fmt(pa[0])} {S.fmt(pa[1])} {S.fmt(pa[2])}"
              f"            {b[0]:>4}{b[1]:>4}{b[2]:>4}  {S.fmt(pb[0])} {S.fmt(pb[1])} {S.fmt(pb[2])}")
    print("\nper target, where the two differ (errors = FP + FN):")
    for t in targets:
        a = tally(truth, tool, lambda g, x, i=t["id"]: x == i)
        b = tally(truth, raw, lambda g, x, i=t["id"]: x == i)
        if a != b:
            print(f"  {t['id']:<11} mcycle FP={a[1]:<3} FN={a[2]:<3}   raw KOfam FP={b[1]:<3} FN={b[2]:<3}"
                  f"{'   *trap' if t['id'] in S.TRAP else ''}")


if __name__ == "__main__":
    main()
