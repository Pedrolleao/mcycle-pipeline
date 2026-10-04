#!/usr/bin/env python3
"""
build_curated_function_gt.py — curated-function ground truth for the methane panel.

The KEGG ground truth (ground_truth.tsv) follows KO presence. The headline question
of the validation — does the tool separate functions that share a KO — needs a
reference in which those cells carry the FUNCTION. This file is that reference:

  1. every KEGG cell, verbatim (source `kegg_v1`);
  2. minus the cells listed as `correction` in curated_cells.tsv, which are replaced
     (source `curated_function`): pmoA / pmoB / pmoC in ammonia oxidizers and in a
     hydrocarbon-monooxygenase carrier, where K10944 / K10945 / K10946 are shared;
  3. plus `literature` cells for the genomes KEGG does not hold (source
     `curated_literature`) — only cells a cited paper states for that genome;
  4. plus `mcrA_anme`, which has no KO: `present` where phenotype_gt.tsv gives the
     Mcr direction `reverse`, `absent` everywhere else (source `curated_phenotype`).

No cell comes from the pipeline's own output. curated_cells.tsv and
phenotype_gt.tsv were written before the pipeline was run on the panel; later
changes are logged, with their evidence, in CHANGELOG.md.

Output: validation/curated_function_gt.tsv
        (genome, target, expected, source, rationale)
Run after build_ground_truth.py.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KEGG_GT = HERE / "ground_truth.tsv"
CELLS = HERE / "curated_cells.tsv"
PHENO = HERE / "phenotype_gt.tsv"
ROSTER = HERE / "panel.tsv"
OUT = HERE / "curated_function_gt.tsv"


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def main() -> None:
    roster = [r["name"] for r in read_tsv(ROSTER)]
    gt: dict[tuple[str, str], tuple[str, str, str]] = {
        (r["genome"], r["target"]): (r["expected"], r["source"], "")
        for r in read_tsv(KEGG_GT)}
    in_kegg = {g for g, _ in gt}

    n_corr = n_lit = 0
    for c in read_tsv(CELLS):
        key = (c["genome"], c["target"])
        if c["genome"] not in roster:
            sys.exit(f"error: curated cell for unknown genome {c['genome']}")
        rationale = f"{c['evidence']} [doi:{c['doi']}]"
        if c["kind"] == "correction":
            if key not in gt:
                sys.exit(f"error: correction {key} has no KEGG cell to correct")
            if gt[key][0] == c["expected"]:
                sys.exit(f"error: correction {key} equals the KEGG cell — remove it")
            gt[key] = (c["expected"], "curated_function", rationale); n_corr += 1
        elif c["kind"] == "literature":
            if c["genome"] in in_kegg:
                sys.exit(f"error: literature cell {key} for a genome KEGG holds — "
                         "use kind=correction and say why KEGG is wrong")
            gt[key] = (c["expected"], "curated_literature", rationale); n_lit += 1
        else:
            sys.exit(f"error: unknown kind '{c['kind']}' for {key}")

    for p in read_tsv(PHENO):
        reverse = p["mcr_direction"] == "reverse"
        why = {"reverse": "anaerobic methanotroph (ANME): Mcr runs in reverse",
               "methanogenic": "methanogen: McrA present but not of an ANME clade",
               "none": "no methyl-coenzyme M reductase in this organism"}[p["mcr_direction"]]
        gt[(p["genome"], "mcrA_anme")] = (
            "present" if reverse else "absent", "curated_phenotype",
            f"{why} ({p['citation']}) [doi:{p['doi']}]")

    order = {g: i for i, g in enumerate(roster)}
    with open(OUT, "w") as fh:
        fh.write("genome\ttarget\texpected\tsource\trationale\n")
        for (g, t), (exp, src, why) in sorted(gt.items(),
                                              key=lambda kv: (order[kv[0][0]], kv[0][1])):
            fh.write(f"{g}\t{t}\t{exp}\t{src}\t{why}\n")
    n_p = sum(v[0] == "present" for v in gt.values())
    print(f"[curated_gt] {len(gt)} cells, {len({g for g, _ in gt})} genomes "
          f"({n_p} present, {len(gt) - n_p} absent); {n_corr} corrections, "
          f"{n_lit} literature cells, {len(roster)} mcrA_anme cells -> {OUT}")


if __name__ == "__main__":
    main()
