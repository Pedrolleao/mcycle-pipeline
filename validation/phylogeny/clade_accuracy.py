#!/usr/bin/env python3
"""
clade_accuracy.py — orthogonal check of the two clade calls on the GTDB-500 genomes.

mcycle decides the Mcr direction and pmoA-vs-amoA from the SEQUENCE of one protein
(clade HMMs). The reference used here is independent of that sequence: the GTDB
taxonomy of the GENOME, which rests on 53 / 120 marker genes. A genome of an ANME
family is expected to carry a `reverse` McrA, a genome of a methanogen lineage a
`methanogenic` one, a genome of an alkane-oxidizer lineage no McrA call; a
methanotroph lineage is expected to carry pmoA, an ammonia-oxidizer lineage is not.
Lineages with no settled physiology are listed and not scored.

Strata, because the clade HMMs were trained on GTDB-typed genomes:
  training genome   the genome supplied a training sequence          (in-sample)
  training genus    another genome of its genus did
  new genus         no genome of its genus did                       (the real test)
and, as in the sister tools, characterized (named genus) vs candidate (placeholder
genus) lineages.

  MCYCLE_RESULTS=results_gtdb500 python validation/phylogeny/clade_accuracy.py \
      --selection comparators/gtdb500_m/selection.tsv --out comparators/gtdb500_m/DIR_ACCURACY.md
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation"))
import score_mcycle as S  # noqa: E402

ANME = ["f__Methanospirareceae", "f__Methanocomedenaceae", "f__Methanogasteraceae",
        "f__Methanoperedenaceae", "g__Methanovorans"]
METHANOGEN = ["c__Methanobacteria", "c__Methanococci", "c__Methanomicrobia", "c__Methanocellia",
              "c__Methanopyri", "f__Methanosarcinaceae", "f__Methanotrichaceae",
              "f__Methermicoccaceae", "o__Methanomassiliicoccales", "o__Methanofastidiosales",
              "f__Methanomethylicaceae", "c__Methanonatronarchaeia"]
ALKANE = ["f__EX4572-44", "f__Syntropharchaeaceae", "g__Alkanophaga"]
MOB = ["f__Methylomonadaceae", "f__Methylococcaceae", "f__Methylothermaceae", "g__Methylocystis",
       "g__Methylosinus", "g__Methylocapsa", "f__Methylacidiphilaceae", "f__Methylomirabilaceae"]
AMO = ["f__Nitrosomonadaceae", "f__Nitrosococcaceae", "f__Nitrosopumilaceae",
       "f__Nitrososphaeraceae", "g__Nitrospira"]


def in_lineage(lin: str, pats: list[str]) -> bool:
    return any(re.search(re.escape(p) + r"(_[A-Z]+)?(;|$)", lin) for p in pats)


def wilson(k: int, n: int) -> str:
    if not n:
        return "n/a"
    z = 1.96; p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return f"{p:.3f} [{max(0, c - h):.3f}, {min(1, c + h):.3f}]"


def training_sets(target: str) -> tuple[set[str], set[str]]:
    """(assemblies, genera) behind the clade models of a target."""
    asm, genera = set(), set()
    f = ROOT / "targets" / target / "refs.fasta"
    for line in open(f):
        if line.startswith(">"):
            asm.add(line.split("||")[1].split("|")[0].replace(".", "_"))
            m = re.search(r"OS=(\S+)", line)
            if m:
                genera.add(m.group(1))
    return asm, genera


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    sel = {r["ncbi_acc"].replace(".", "_"): r
           for r in csv.DictReader(open(args.selection), delimiter="\t")}
    rows = []
    for trap, target in (("mcr", "mcrA_anme"), ("pmo", "pmoA")):
        t_asm, t_gen = training_sets(target)
        for g, r in sel.items():
            f = S.RESULTS / g / "calls" / "mcycle_calls.tsv"
            if not f.exists():
                continue
            calls = {c["target_id"]: c for c in S.read_tsv(f)}
            lin = r["gtdb_taxonomy"]; genus = lin.split(";")[5][3:]
            if trap == "mcr":
                a, m = calls["mcrA_anme"]["status"], calls["mcrA"]["status"]
                call = ("reverse" if a == "confirmed" else
                        "methanogenic" if m in S.PRESENT else "no_call")
                family_hit = a != "absent" or m != "absent" or calls["mcrB"]["status"] != "absent"
                expect = ("reverse" if in_lineage(lin, ANME) else
                          "no_call" if in_lineage(lin, ALKANE) else
                          "methanogenic" if in_lineage(lin, METHANOGEN) else "unknown")
                clade = calls["mcrA_anme"]["pfam_hits"]
            else:
                p = calls["pmoA"]["status"]
                call = "pmoA" if p == "confirmed" else "not_pmoA"
                family_hit = p != "absent"
                expect = ("pmoA" if in_lineage(lin, MOB) else
                          "not_pmoA" if in_lineage(lin, AMO) else "unknown")
                clade = calls["pmoA"]["pfam_hits"]
            if not family_hit:
                continue
            rel = ("training genome" if g in t_asm else
                   "training genus" if genus.split("_")[0] in t_gen else "new genus")
            named = bool(re.fullmatch(r"[A-Z][a-z]+(_[A-Z]+)?", genus))
            rows.append({"trap": trap, "genome": g, "stratum": r["stratum"], "clade": r["clade"],
                         "genus": genus, "named_genus": "yes" if named else "no",
                         "training_relation": rel, "expected": expect, "mcycle": call,
                         "model": clade,
                         "agree": "" if expect == "unknown" else str(int(expect == call)),
                         "gtdb_taxonomy": lin})

    with open(args.out.with_suffix(".tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    L = ["# Clade calls against the genome taxonomy (GTDB-500)\n",
         "Reference = GTDB r232 lineage of the genome (marker-gene phylogeny), independent of "
         "the McrA / PmoA sequence the tool reads. Agreement with 95 % Wilson intervals. "
         "Only genomes carrying a member of the protein family are listed.\n"]
    for trap, title in (("mcr", "Mcr direction (reverse / methanogenic / no call)"),
                        ("pmo", "pmoA vs other copper monooxygenases")):
        sub = [r for r in rows if r["trap"] == trap]
        scored = [r for r in sub if r["agree"] != ""]
        L.append(f"\n## {title}\n")
        L.append(f"Genomes with a family member: {len(sub)}; with a lineage expectation: "
                 f"{len(scored)}; lineage of unsettled physiology: {len(sub) - len(scored)}.\n")
        L.append("| stratum | agree / n | agreement [95 % CI] |")
        L.append("|---|---|---|")

        def line(label, keep):
            s = [r for r in scored if keep(r)]
            k = sum(r["agree"] == "1" for r in s)
            L.append(f"| {label} | {k} / {len(s)} | {wilson(k, len(s))} |")
        line("**all**", lambda r: True)
        for rel in ("new genus", "training genus", "training genome"):
            line(rel, lambda r, x=rel: r["training_relation"] == x)
        line("named genus", lambda r: r["named_genus"] == "yes")
        line("placeholder genus", lambda r: r["named_genus"] == "no")
        line("new genus, named", lambda r: r["training_relation"] == "new genus" and r["named_genus"] == "yes")
        line("new genus, placeholder", lambda r: r["training_relation"] == "new genus" and r["named_genus"] == "no")
        L.append("\nExpected (rows) x mcycle (columns):\n")
        cats = sorted({r["mcycle"] for r in sub})
        L.append("| expected | " + " | ".join(cats) + " |")
        L.append("|---|" + "---|" * len(cats))
        for e in sorted({r["expected"] for r in sub}):
            L.append(f"| {e} | " + " | ".join(
                str(sum(r["expected"] == e and r["mcycle"] == c for r in sub)) for c in cats) + " |")
        bad = [r for r in scored if r["agree"] == "0"]
        if bad:
            L.append(f"\nDisagreements ({len(bad)}):\n")
            for r in bad:
                L.append(f"- {r['genome']} ({r['clade']}; {r['genus']}; {r['training_relation']}): "
                         f"expected {r['expected']}, mcycle {r['mcycle']}")
        unk = [r for r in sub if r["agree"] == ""]
        if unk:
            L.append(f"\nNot scored — lineage of unsettled physiology ({len(unk)}):\n")
            for r in unk:
                L.append(f"- {r['genome']}: {';'.join(r['gtdb_taxonomy'].split(';')[1:6])} — "
                         f"mcycle {r['mcycle']}")
    args.out.write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
