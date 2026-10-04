#!/usr/bin/env python3
"""
explain_cell.py — show the evidence behind one (genome, target) call: the best
hmmscan hits of the genome on each of the target's profiles with their thresholds,
and the call row. Diagnostic only.

    python validation/explain_cell.py <genome> <target> [<target> ...]
    MCYCLE_RESULTS=results_ref (default)
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / os.environ.get("MCYCLE_RESULTS", "results_ref")


def main() -> None:
    genome, wanted = sys.argv[1], sys.argv[2:]
    targets = {t["id"]: t for t in yaml.safe_load(open(ROOT / "config/targets.yaml"))["targets"]}
    tc = {}
    for r in csv.DictReader(open(ROOT / "resources/hmm/tc_cutoffs.tsv"), delimiter="\t"):
        tc[r["profile_id"]] = float(r["tc"]) if r["tc"] else None
    hits: dict[str, list[tuple[float, str, int, int, int]]] = {}
    for line in open(RESULTS / genome / "hmm" / f"{genome}.hmmscan.tsv"):
        if line.startswith("#"):
            continue
        p = line.split()
        hits.setdefault(p[0], []).append((float(p[7]), p[3], int(p[5]), int(p[15]), int(p[16])))
    calls = {r["target_id"]: r for r in csv.DictReader(
        open(RESULTS / genome / "calls" / "mcycle_calls.tsv"), delimiter="\t")}
    for tid in wanted:
        t = targets[tid]; c = calls[tid]
        print(f"{genome}  {tid}: status={c['status']}  evidence={c['evidence_source']}  "
              f"protein={c['protein_id']}  blast={c['blast_acc']} {c['blast_pident']}")
        profiles = list(t.get("ko") or []) + [p for p in hits if p.startswith(tid + "_")]
        for ko in profiles:
            thr = tc.get(ko)
            best = sorted(set(hits.get(ko, [])), reverse=True)[:4]
            shown = ", ".join(f"{s:.0f} ({q}, {ql} aa, hmm {a}-{b})" for s, q, ql, a, b in best)
            print(f"    {ko}  threshold {thr}:  {shown or 'no hit'}")


if __name__ == "__main__":
    main()
