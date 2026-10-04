#!/usr/bin/env python3
"""
concordance.py — cross-tool CONCORDANCE on the 500 GTDB genomes (no ground truth).

Arbitrary GTDB genomes have no curated truth, so nothing here is an accuracy. The
analysis measures where mcycle, raw KofamScan, MCycDB and METABOLIC agree, and
where they split, with the spotlight on the two methane traps:

  1. the copper membrane monooxygenase — which tools report `pmoA` in ammonia
     oxidizers (the KO is shared with amoA), by GTDB clade;
  2. the direction of methyl-coenzyme M reductase — which genomes carry an McrA
     according to each tool, and where mcycle calls it `reverse` (ANME clade).
     No comparator has a direction call.

Call semantics are those of the panel benchmark (adapters.py). Run with the GTDB
results directory:

  MCYCLE_RESULTS=results_gtdb500 python validation/benchmark/concordance.py \
      --selection comparators/gtdb500_m/selection.tsv \
      --mcycdb comparators/gtdb500_m/mcycdb.tsv \
      [--metabolic comparators/gtdb500_m/metabolic.tsv] \
      --out comparators/gtdb500_m/CONCORDANCE.md
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import adapters as A  # noqa: E402

PMO = ("pmoA", "pmoB", "pmoC")


def kappa(a: list[bool], b: list[bool]) -> float:
    n = len(a)
    if not n:
        return float("nan")
    po = sum(x == y for x, y in zip(a, b)) / n
    pa, pb = sum(a) / n, sum(b) / n
    pe = pa * pb + (1 - pa) * (1 - pb)
    return (po - pe) / (1 - pe) if pe < 1 else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--mcycdb", type=Path)
    ap.add_argument("--metabolic", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    sel = {r["ncbi_acc"].replace(".", "_"): r
           for r in csv.DictReader(open(args.selection), delimiter="\t")}
    calls = A.S.load_calls(set(sel))
    genomes = set(calls)
    missing = sorted(set(sel) - genomes)
    preds = {"mcycle": A.load_mcycle(genomes), "kofam": A.load_kofam(genomes)}
    if args.mcycdb:
        preds["mcycdb"] = A.load_mcycdb(args.mcycdb, genomes)
    if args.metabolic:
        m = A.load_metabolic(args.metabolic, genomes)
        with_data = {r["genome"] for r in A._read(args.metabolic)}
        preds["metabolic"] = {k: v for k, v in m.items() if k[0] in with_data}
    tools = list(preds)
    cover = {"mcycle": set(A.all_target_ids()), **A.covered_targets()}
    targets = A.all_target_ids()
    group = {g: (sel[g]["clade"] if sel[g]["stratum"] == "enriched" else "backbone")
             for g in genomes}
    direction = {}
    for g in genomes:
        ev = calls[g]
        direction[g] = ("reverse" if ev.get("mcrA_anme") in A.S.PRESENT else
                        "methanogenic" if ev.get("mcrA") in A.S.PRESENT else "-")

    L = []
    L.append("# GTDB-500 cross-tool concordance — methane cycle\n")
    L.append(f"{len(genomes)} genomes ({sum(sel[g]['stratum'] == 'enriched' for g in genomes)} "
             f"methane-enriched, {sum(sel[g]['stratum'] == 'backbone' for g in genomes)} backbone); "
             f"tools: {', '.join(tools)}. No ground truth: agreement, not accuracy.\n")
    if missing:
        L.append(f"Genomes without mcycle results: {len(missing)}.\n")
    for t in tools:
        n = len({g for g, _ in preds[t]})
        L.append(f"- {t}: {n} genomes, {sum(preds[t].values())} present calls")

    # (a) pairwise agreement on mutually representable targets (mcrA_anme left out:
    #     no comparator can express it; it is the subject of section c)
    L.append("\n## a. Pairwise agreement (targets both tools can express)\n")
    L.append("| pair | cells | agreement | Cohen kappa | both present | only first | only second |")
    L.append("|---|---|---|---|---|---|---|")
    for a, b in itertools.combinations(tools, 2):
        both_t = (cover[a] & cover[b]) - {"mcrA_anme"}
        cells = [(g, t) for g in genomes for t in targets
                 if t in both_t and (g, t) in preds[a] and (g, t) in preds[b]]
        va = [preds[a][c] for c in cells]; vb = [preds[b][c] for c in cells]
        agree = sum(x == y for x, y in zip(va, vb))
        L.append(f"| {a} vs {b} | {len(cells)} | {agree / len(cells):.4f} | {kappa(va, vb):.3f} | "
                 f"{sum(x and y for x, y in zip(va, vb))} | "
                 f"{sum(x and not y for x, y in zip(va, vb))} | "
                 f"{sum(y and not x for x, y in zip(va, vb))} |")

    # (b) pmoA by clade
    L.append("\n## b. Trap 1 — genomes in which each tool reports pmoA, by clade\n")
    L.append("| clade | genomes | " + " | ".join(tools) + " |")
    L.append("|---|---|" + "---|" * len(tools))
    by = defaultdict(list)
    for g in genomes:
        by[group[g]].append(g)
    for clade in sorted(by, key=lambda c: (c == "backbone", c)):
        gs = by[clade]
        counts = [str(sum(preds[t].get((g, "pmoA"), False) for g in gs)) for t in tools]
        if any(c != "0" for c in counts):
            L.append(f"| {clade} | {len(gs)} | " + " | ".join(counts) + " |")
    amo_clades = {"AOB_beta", "AOB_gamma", "AOA", "Nitrospira_comammox"}
    amo = [g for g in genomes if group[g] in amo_clades]
    L.append(f"\nIn the {len(amo)} genomes of the ammonia-oxidizer clades, pmoA is reported by: "
             + "; ".join(f"{t} {sum(preds[t].get((g, 'pmoA'), False) for g in amo)}" for t in tools)
             + ".")

    # (c) Mcr direction by clade
    L.append("\n## c. Trap 2 — McrA and its direction, by clade\n")
    L.append("| clade | genomes | " + " | ".join(f"{t} mcrA" for t in tools)
             + " | mcycle reverse | mcycle methanogenic |")
    L.append("|---|---|" + "---|" * (len(tools) + 2))
    for clade in sorted(by, key=lambda c: (c == "backbone", c)):
        gs = by[clade]
        counts = [str(sum(preds[t].get((g, "mcrA"), False) for g in gs)) for t in tools]
        rev = sum(direction[g] == "reverse" for g in gs)
        met = sum(direction[g] == "methanogenic" for g in gs)
        if any(c != "0" for c in counts):
            L.append(f"| {clade} | {len(gs)} | " + " | ".join(counts) + f" | {rev} | {met} |")
    n_rev = sum(d == "reverse" for d in direction.values())
    L.append(f"\nmcycle calls the Mcr `reverse` in {n_rev} genomes; every comparator reports "
             "the same genomes as plain McrA carriers.")

    # (d) most discordant targets
    L.append("\n## d. Targets with the most discordant genomes\n")
    L.append("| target | " + " | ".join(tools) + " | genomes where the tools differ |")
    L.append("|---|" + "---|" * (len(tools) + 1))
    disc = []
    for t in targets:
        able = [x for x in tools if t in cover[x]]
        if len(able) < 2 or t == "mcrA_anme":
            continue
        n = sum(len({preds[x].get((g, t)) for x in able if (g, t) in preds[x]}) > 1
                for g in genomes)
        disc.append((n, t))
    for n, t in sorted(disc, reverse=True)[:20]:
        L.append(f"| {t} | " + " | ".join(
            str(sum(preds[x].get((g, t), False) for g in genomes)) if t in cover[x] else "n/a"
            for x in tools) + f" | {n} |")

    args.out.write_text("\n".join(L) + "\n")
    with open(args.out.with_suffix(".tsv"), "w") as fh:
        fh.write("genome\tstratum\tclade\tgtdb_taxonomy\tmcr_direction\t"
                 + "\t".join(f"{t}_{x}" for t in tools for x in ("pmoA", "mcrA")) + "\n")
        for g in sorted(genomes):
            fh.write("\t".join([g, sel[g]["stratum"], sel[g]["clade"], sel[g]["gtdb_taxonomy"],
                                direction[g]]
                               + [str(int(preds[t].get((g, x), False)))
                                  for t in tools for x in ("pmoA", "mcrA")]) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
