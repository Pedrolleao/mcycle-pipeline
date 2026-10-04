#!/usr/bin/env python3
"""
ko_margins.py — per target, how far apart are true and false hits on its KO
profiles in the TRAINING genomes?

For every KO-anchored target and every training genome with a ground-truth cell,
take the genome's best full-sequence score on the target's profiles as a fraction
of the profile's threshold. Report, per target, the lowest ratio among
ground-truth-present genomes and the highest among ground-truth-absent genomes.
A target is separable when min(present) > max(absent); the threshold sits well
when 1.0 lies between the two.

Calibration aid for M3 — reads hold-out genomes never (panel.tsv split).

    python validation/ko_margins.py            # table of every target
    python validation/ko_margins.py -v mvhD    # per-genome ratios of one target
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / os.environ.get("MCYCLE_RESULTS", "results_ref")


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def main() -> None:
    verbose = sys.argv[2:] if len(sys.argv) > 2 and sys.argv[1] == "-v" else []
    train = [r["name"] for r in read_tsv(ROOT / "validation/panel.tsv") if r["split"] == "train"]
    truth = {(r["genome"], r["target"]): r["expected"]
             for r in read_tsv(ROOT / "validation/curated_function_gt.tsv")}
    targets = yaml.safe_load(open(ROOT / "config/targets.yaml"))["targets"]
    tc = {r["profile_id"]: float(r["tc"]) for r in read_tsv(ROOT / "resources/hmm/tc_cutoffs.tsv") if r["tc"]}
    best: dict[str, dict[str, tuple[float, str, int]]] = {}
    for g in train:
        f = RESULTS / g / "hmm" / f"{g}.hmmscan.tsv"
        if not f.exists():
            continue
        d = best.setdefault(g, {})
        for line in open(f):
            if line.startswith("#"):
                continue
            p = line.split()
            if p[0] not in d or float(p[7]) > d[p[0]][0]:
                d[p[0]] = (float(p[7]), p[3], int(p[5]))
    print(f"{'target':<11}{'KO':<8}{'thr':>7}  {'min present':>12} {'max absent':>11}  verdict")
    for t in targets:
        if t["id"] == "mcrA_anme":
            continue
        rows = []
        for g in best:
            exp = truth.get((g, t["id"]))
            if exp is None:
                continue
            ratio, who = 0.0, ""
            for ko in t.get("ko") or []:
                thr = tc.get(ko)
                if thr and ko in best[g]:
                    s, q, ql = best[g][ko]
                    if s / thr > ratio:
                        ratio, who = s / thr, f"{ko} {s:.0f} {q} ({ql} aa)"
            rows.append((ratio, exp, g, who))
        pres = sorted(r for r in rows if r[1] == "present")
        absn = sorted((r for r in rows if r[1] == "absent"), reverse=True)
        lo = pres[0][0] if pres else float("nan")
        hi = absn[0][0] if absn else float("nan")
        verdict = ("ok" if lo > 1.0 > hi else "ok (no positives)" if not pres and hi < 1.0
                   else "separable, threshold off" if lo > hi else "OVERLAP")
        thr0 = tc.get((t.get("ko") or ["?"])[0], float("nan"))
        print(f"{t['id']:<11}{(t.get('ko') or ['-'])[0]:<8}{thr0:>7.0f}  {lo:>12.2f} {hi:>11.2f}  {verdict}")
        if t["id"] in verbose:
            for ratio, exp, g, who in sorted(rows, reverse=True):
                if ratio >= 0.4:
                    print(f"      {ratio:5.2f}  {exp:<8}{g:<28}{who}")


if __name__ == "__main__":
    main()
