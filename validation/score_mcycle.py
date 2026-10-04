#!/usr/bin/env python3
"""
score_mcycle.py — score the pipeline's calls against the methane ground truth.

Ported from scycle-pipeline/validation/score_scycle.py; the metric battery and the
bootstrap are computed the same way in the three sister tools (genome-cluster
percentile bootstrap, B = 10,000, seed 1234, one RNG threaded through every
interval in the order below).

Only the (genome, target) cells present in the ground truth are scored. A call is
`present` when its status in <results>/<genome>/calls/mcycle_calls.tsv is
`confirmed` or `domain-only`; `absent` and `disqualified` are `absent`.

Frames reported:
  A  whole panel — ALL targets, the homology-trap subset, the rest; per pathway
  B  independent-only trap precision — trap cells minus the seed-sourced ones
     (validation/seed_leakage.tsv; see trap_independence.py)
  C  train / hold-out split — the hold-out genomes of validation/panel.tsv, frozen
     2026-10-04; the hold-out micro-F1 is the generalization headline. Seed-sourced
     hold-out cells are removed ("de-leaked"); by design there are none.

Environment:
  GT_FILE         ground-truth file under validation/ (default
                  curated_function_gt.tsv; ground_truth.tsv = the KEGG contrast)
  MCYCLE_RESULTS  results directory (default results_ref)
  MCYCLE_SCOPE    all | train | holdout — restrict scoring to one side of the
                  split (M3 hardening runs with `train`, so that hold-out results
                  are never looked at while anything is being tuned)
"""
from __future__ import annotations

import csv
import math
import os
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

GT = ROOT / "validation" / os.environ.get("GT_FILE", "curated_function_gt.tsv")
RESULTS = ROOT / os.environ.get("MCYCLE_RESULTS", "results_ref")
SCOPE = os.environ.get("MCYCLE_SCOPE", "all")
TARGETS = ROOT / "config" / "targets.yaml"
ROSTER = ROOT / "validation" / "panel.tsv"
OUT_M = ROOT / "validation" / "mcycle_metrics.tsv"
OUT_C = ROOT / "validation" / "mcycle_confusion.tsv"

PRESENT = {"confirmed", "domain-only"}

# Homology-trap targets: a homologue with ANOTHER function passes, or nearly passes,
# the same KO profile, so the call rests on a clade gate. Fixed 2026-10-04, before
# the pipeline was run on the reference panel:
#   mcrA      vs alkyl-coenzyme M reductases of alkane-oxidizing archaea
#   mcrA_anme vs McrA of methanogens (the direction call)
#   pmoA/B/C  vs ammonia monooxygenase and hydrocarbon Cu-monooxygenases (shared KOs)
#   mmoX      vs butane / propane / alkene / toluene di-iron monooxygenases
#   mxaF/xoxF vs the other PQQ alcohol dehydrogenases (PedE / PedH / ExaA type)
TRAP = {"mcrA", "mcrA_anme", "pmoA", "pmoB", "pmoC", "mmoX", "mxaF", "xoxF"}

BOOTSTRAP_B, BOOTSTRAP_SEED = 10_000, 1234


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def holdout_genomes() -> set[str]:
    return {r["name"] for r in read_tsv(ROSTER) if r["split"] == "holdout"}


HOLDOUT = holdout_genomes()


def prf(tp, fp, fn):
    p = tp / (tp + fp) if (tp + fp) else float("nan")
    r = tp / (tp + fn) if (tp + fn) else float("nan")
    f = (2 * p * r / (p + r)) if (not math.isnan(p) and not math.isnan(r) and (p + r)) \
        else float("nan")
    return p, r, f


def boot_ci(gtallies: dict, name: str, *, B: int = BOOTSTRAP_B,
            rng: random.Random | None = None) -> tuple[float, float]:
    """Genome-cluster 95 % percentile bootstrap CI. gtallies: genome -> (TP,FP,FN,TN).
    Collapses to [point, point] when the scope has no FP and no FN."""
    rng = rng or random.Random(BOOTSTRAP_SEED)
    glist = list(gtallies)
    if not glist:
        return (float("nan"), float("nan"))
    vals = []
    for _ in range(B):
        TP = FP = FN = 0
        for _ in glist:
            t = gtallies[rng.choice(glist)]
            TP += t[0]; FP += t[1]; FN += t[2]
        p, r, f = prf(TP, FP, FN)
        v = {"precision": p, "recall": r, "f1": f}[name]
        if not math.isnan(v):
            vals.append(v)
    if not vals:
        return (float("nan"), float("nan"))
    vals.sort()
    return (vals[int(0.025 * len(vals))], vals[min(len(vals) - 1, int(0.975 * len(vals)))])


def load_truth() -> dict[tuple[str, str], str]:
    truth = {(r["genome"], r["target"]): r["expected"] for r in read_tsv(GT)}
    if SCOPE == "train":
        truth = {k: v for k, v in truth.items() if k[0] not in HOLDOUT}
    elif SCOPE == "holdout":
        truth = {k: v for k, v in truth.items() if k[0] in HOLDOUT}
    elif SCOPE != "all":
        sys.exit(f"error: MCYCLE_SCOPE must be all, train or holdout, not '{SCOPE}'")
    return truth


def load_calls(genomes: set[str]) -> dict[str, dict[str, str]]:
    """genome -> {target: status}; genomes without results are left out."""
    out = {}
    for g in sorted(genomes):
        f = RESULTS / g / "calls" / "mcycle_calls.tsv"
        if f.exists():
            out[g] = {r["target_id"]: r["status"] for r in read_tsv(f)}
    return out


def cat_map() -> dict[str, str]:
    import yaml
    return {t["id"]: t["category"] for t in yaml.safe_load(open(TARGETS))["targets"]}


def compute_metrics(require_all: bool = False) -> dict:
    from trap_independence import is_seed

    truth = load_truth(); cats = cat_map()
    genomes = {g for g, _ in truth}
    calls = load_calls(genomes)
    missing = sorted(genomes - set(calls))
    if missing and require_all:
        sys.exit(f"error: no pipeline results under {RESULTS} for: {', '.join(missing)}")
    rng = random.Random(BOOTSTRAP_SEED)   # one RNG threaded through all intervals

    per, confusion = {}, []
    for (g, t), exp in truth.items():
        if g not in calls:
            continue
        status = calls[g].get(t, "absent")
        pred = status in PRESENT
        exp_pos = exp == "present"
        cell = ("TP" if pred and exp_pos else "FP" if pred else
                "FN" if exp_pos else "TN")
        confusion.append((g, t, exp, status, "present" if pred else "absent", cell))
        per.setdefault(t, {"TP": 0, "FP": 0, "FN": 0, "TN": 0})[cell] += 1

    rows = []
    for t, d in per.items():
        p, r, f = prf(d["TP"], d["FP"], d["FN"])
        rows.append({"target": t, "category": cats.get(t, "?"), "n": sum(d.values()),
                     **d, "precision": p, "recall": r, "f1": f})

    _IDX = {"TP": 0, "FP": 1, "FN": 2, "TN": 3}

    def tallies(keep) -> dict:
        """Per-genome (TP,FP,FN,TN) over the cells `keep(genome, target)` accepts —
        the genome is the bootstrap unit (cells within a genome are correlated)."""
        gt: dict[str, list[int]] = {}
        for (g, t, _e, _s, _p, cell) in confusion:
            if keep(g, t):
                gt.setdefault(g, [0, 0, 0, 0])[_IDX[cell]] += 1
        return {g: tuple(v) for g, v in gt.items()}

    def aggregate(gt: dict, tids: set | None = None) -> dict:
        TP = sum(v[0] for v in gt.values()); FP = sum(v[1] for v in gt.values())
        FN = sum(v[2] for v in gt.values()); TN = sum(v[3] for v in gt.values())
        p, r, f = prf(TP, FP, FN)
        out = {"TP": TP, "FP": FP, "FN": FN, "TN": TN, "precision": p, "recall": r,
               "f1": f, "genomes": sorted(gt)}
        if tids is not None:
            f1s = sorted(m["f1"] for m in rows
                         if m["target"] in tids and not math.isnan(m["f1"]))
            out["median_f1"] = f1s[len(f1s) // 2] if f1s else float("nan")
        out["f1_ci"] = boot_ci(gt, "f1", rng=rng)
        out["precision_ci"] = boot_ci(gt, "precision", rng=rng)
        out["recall_ci"] = boot_ci(gt, "recall", rng=rng)
        return out

    all_t = set(per); trap_t = all_t & TRAP
    aggregates = {
        "ALL": aggregate(tallies(lambda g, t: True), all_t),
        "trap": aggregate(tallies(lambda g, t: t in TRAP), trap_t),
        "nontrap": aggregate(tallies(lambda g, t: t not in TRAP), all_t - trap_t)}

    per_pathway = {}
    for cat in sorted({m["category"] for m in rows}):
        tids = {m["target"] for m in rows if m["category"] == cat}
        per_pathway[cat] = aggregate(tallies(lambda g, t, s=tids: t in s), tids)

    # Frame B — independent-only trap cells
    indep_trap = aggregate(tallies(lambda g, t: t in TRAP and not is_seed(g, t)))

    # Frame C — train / hold-out; the hold-out is de-leaked
    split = {
        "train": aggregate(tallies(lambda g, t: g not in HOLDOUT)),
        "holdout_raw": aggregate(tallies(lambda g, t: g in HOLDOUT)),
        "holdout": aggregate(tallies(lambda g, t: g in HOLDOUT and not is_seed(g, t))),
        "holdout_trap": aggregate(tallies(
            lambda g, t: g in HOLDOUT and t in TRAP and not is_seed(g, t)))}
    n_leaked = sum(1 for (g, t, *_r) in confusion if g in HOLDOUT and is_seed(g, t))

    return {"per_target": rows, "per": per, "aggregates": aggregates,
            "per_pathway": per_pathway, "indep_trap": indep_trap, "split": split,
            "holdout_leaked_cells": n_leaked, "confusion": confusion,
            "missing_genomes": missing, "scored_cells": len(confusion)}


def fmt(v) -> str:
    return "NA" if (isinstance(v, float) and math.isnan(v)) else f"{v:.3f}"


def fmt_ci(ci) -> str:
    if not ci or math.isnan(ci[0]) or math.isnan(ci[1]):
        return ""
    return f" [{ci[0]:.3f}, {ci[1]:.3f}]"


def line(label: str, a: dict) -> str:
    return (f"  {label:<34} cells={a['TP'] + a['FP'] + a['FN'] + a['TN']:<5} "
            f"TP={a['TP']} FP={a['FP']} FN={a['FN']} TN={a['TN']}  "
            f"P={fmt(a['precision'])}{fmt_ci(a.get('precision_ci'))}  "
            f"R={fmt(a['recall'])}{fmt_ci(a.get('recall_ci'))}  "
            f"F1={fmt(a['f1'])}{fmt_ci(a.get('f1_ci'))}")


def main() -> None:
    M = compute_metrics()
    rows, confusion = M["per_target"], M["confusion"]

    if SCOPE == "all" and not M["missing_genomes"]:
        with open(OUT_M, "w") as fh:
            fh.write("target\tcategory\tn\tTP\tFP\tFN\tTN\tprecision\trecall\tf1\n")
            for m in sorted(rows, key=lambda x: (x["category"], x["target"])):
                fh.write(f"{m['target']}\t{m['category']}\t{m['n']}\t{m['TP']}\t{m['FP']}\t"
                         f"{m['FN']}\t{m['TN']}\t{fmt(m['precision'])}\t{fmt(m['recall'])}\t"
                         f"{fmt(m['f1'])}\n")
        with open(OUT_C, "w") as fh:
            fh.write("genome\ttarget\texpected\tstatus\tpredicted\tcell\n")
            for c in confusion:
                fh.write("\t".join(str(x) for x in c) + "\n")

    print(f"ground truth: {GT.name}   results: {RESULTS.name}   scope: {SCOPE}   "
          f"scored cells: {M['scored_cells']}")
    if M["missing_genomes"]:
        print(f"NOT SCORED (no results): {', '.join(M['missing_genomes'])}")
    print(f"\n{'target':<12}{'pathway':<32}{'TP':>4}{'FP':>4}{'FN':>4}{'TN':>4}   P     R     F1")
    for m in sorted(rows, key=lambda x: (x["category"], x["target"])):
        flag = " *" if m["target"] in TRAP else ""
        print(f"  {m['target']:<10}{m['category']:<32}{m['TP']:>4}{m['FP']:>4}{m['FN']:>4}"
              f"{m['TN']:>4}  {fmt(m['precision'])} {fmt(m['recall'])} {fmt(m['f1'])}{flag}")
    print("  (* = homology-trap target)")

    print("\n=== Frame A — whole scope ===")
    print(line("ALL targets", M["aggregates"]["ALL"]))
    print(line("homology-trap targets", M["aggregates"]["trap"]))
    print(line("non-trap targets", M["aggregates"]["nontrap"]))
    print("\nper-pathway:")
    for cat, a in M["per_pathway"].items():
        print(line(cat, a))
    print("\n=== Frame B — independent-only trap cells ===")
    print(line("trap, seed-sourced cells removed", M["indep_trap"]))
    if SCOPE == "all":
        print("\n=== Frame C — train / hold-out ===")
        print(line(f"training ({len(M['split']['train']['genomes'])} genomes)",
                   M["split"]["train"]))
        print(line(f"hold-out ({len(M['split']['holdout']['genomes'])} genomes), de-leaked",
                   M["split"]["holdout"]))
        print(line("hold-out, trap targets", M["split"]["holdout_trap"]))
        print(f"  seed-sourced hold-out cells removed: {M['holdout_leaked_cells']}")

    bad = [c for c in confusion if c[5] in ("FP", "FN")]
    print(f"\n=== {len(bad)} disagreements ===")
    for g, t, exp, status, _pred, cell in sorted(bad, key=lambda c: (c[5], c[1], c[0])):
        print(f"  {cell}  {t:<10}{g:<28}expected={exp:<8}status={status}")


if __name__ == "__main__":
    main()
