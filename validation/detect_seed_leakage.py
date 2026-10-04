#!/usr/bin/env python3
"""
detect_seed_leakage.py — flag (genome, target) cells whose detector was built from
the same organism, or a congener, as the genome under test.

Ported from ncycle-pipeline. The hold-out is meant to measure generalization to
genera absent from every seed list; a positive cell whose BLAST-gate seeds or
custom-HMM training set contain the same genus measures memorization instead. The
panel was designed so that no hold-out genus is a seed genus (PANEL_PLAN.md) — this
script is the mechanical check of that claim, re-run whenever seeds or HMM training
sets change, and it supplies the per-cell independence used by
trap_independence.py and score_mcycle.py.

Method:
  1. Seed organisms per target are parsed from the `OS=` field of every active
     detector source: resources/blast_db/{blast_gated,unstable}_refs.fasta and the
     custom-HMM training sets targets/<id>/expanded.fasta (and refs.fasta).
  2. A panel genome matches a seed when the seed's genus is one of the genome's
     `taxa` tokens in panel.tsv (severity `genus`), and also its species
     (severity `species`).
  3. pmoB and pmoC have no seeds of their own but follow the pmoA call
     (gate_pmo_subunits), so pmoA seeds count for them.
  4. Only ground-truth `present` cells are flagged: a leaked seed can manufacture a
     true positive, never a false positive.

Output: validation/seed_leakage.tsv
        (genome, split, target, severity, seed_acc, seed_organism)
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEED_SOURCES = [
    ROOT / "resources" / "blast_db" / "blast_gated_refs.fasta",
    ROOT / "resources" / "blast_db" / "unstable_refs.fasta",
    *sorted((ROOT / "targets").glob("*/expanded.fasta")),
    *sorted((ROOT / "targets").glob("*/refs.fasta")),
]
ROSTER = ROOT / "validation" / "panel.tsv"
GT = ROOT / "validation" / "curated_function_gt.tsv"
OUT = ROOT / "validation" / "seed_leakage.tsv"

# target -> target whose seeds decide its call
FOLLOWS = {"pmoB": "pmoA", "pmoC": "pmoA"}


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def parse_seed_organisms() -> dict[str, list[tuple[str, str, str]]]:
    """target -> [(genus, species, accession)] over every seed source. Headers
    follow `>target||acc ... OS=Genus species ... OX=`; `Candidatus` is dropped."""
    out: dict[str, set[tuple[str, str, str]]] = {}
    os_re = re.compile(r"OS=(.+?)\s+OX=")
    for src in SEED_SOURCES:
        if not src.exists():
            continue
        target_dir = src.parent.name if src.parent.parent.name == "targets" else None
        for line in src.read_text().splitlines():
            if not line.startswith(">"):
                continue
            head = line[1:]
            if "||" in head:
                target, _, rest = head.partition("||")
                acc = rest.split()[0]
            elif target_dir:
                target, acc = target_dir, head.split()[0]
            else:
                continue
            m = os_re.search(head)
            if not m:
                continue
            toks = [t for t in m.group(1).split() if t.lower() != "candidatus"]
            if not toks:
                continue
            genus = toks[0].lower()
            species = toks[1].lower().rstrip(".,") if len(toks) > 1 else ""
            if species in ("sp", "sp.", "bacterium", "archaeon"):
                species = ""
            out.setdefault(target, set()).add((genus, species, acc))
    return {t: sorted(v) for t, v in out.items()}


def main() -> int:
    seeds = parse_seed_organisms()
    roster = {r["name"]: r for r in read_tsv(ROSTER)}
    rows = []
    for cell in read_tsv(GT):
        if cell["expected"] != "present":
            continue
        genome, target = cell["genome"], cell["target"]
        tokens = set(roster[genome]["taxa"].split(";"))
        best = None
        for genus, species, acc in seeds.get(FOLLOWS.get(target, target), []):
            if genus not in tokens:
                continue
            cand = ((2, "species") if species and species in tokens else (1, "genus"),
                    acc, f"{genus} {species or 'sp.'}")
            if best is None or cand[0][0] > best[0][0]:
                best = cand
        if best:
            rows.append({"genome": genome, "split": roster[genome]["split"],
                         "target": target, "severity": best[0][1],
                         "seed_acc": best[1], "seed_organism": best[2]})

    with open(OUT, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["genome", "split", "target", "severity",
                                           "seed_acc", "seed_organism"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    n_hold = [r for r in rows if r["split"] == "holdout"]
    print(f"[leakage] {len(rows)} seed-sourced positive cells "
          f"({sum(r['severity'] == 'species' for r in rows)} species-level, "
          f"{sum(r['severity'] == 'genus' for r in rows)} genus-level); "
          f"{len(n_hold)} in the hold-out -> {OUT}")
    for r in rows:
        print(f"   {r['split']:<8}{r['genome']:<28}{r['target']:<11}{r['severity']:<9}"
              f"{r['seed_acc']:<12}{r['seed_organism']}")
    if n_hold and "--allow-holdout" not in sys.argv:
        print("error: a hold-out genus supplies a seed — the split is broken",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
