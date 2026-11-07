#!/usr/bin/env python3
"""
make_gap_analysis.py — plain-text per-sample report.

Sections:
  PRESENT          — the targets found (confirmed / domain-only), with their evidence
  COMPLEXES        — partial complexes with the subunits missing; absent and complete
                     complexes as one-line lists
  PROCESS MODULES  — complete; partial with the slots missing (n/N counts the
                     required slots, as in the figures); ruled out with the reason
                     (Mcr direction or a missing key gene); absent as a list

The first two lines (tool, "SAMPLE: <name>   input: …") are the stable header.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PRESENT = {"confirmed", "domain-only"}


def load_calls(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_complexes(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_present(row: dict) -> str:
    bits = [row["target_id"].ljust(14)]
    if row["pfam_hits"]:
        bits.append(row["pfam_hits"].ljust(20))
    bits.append(f"e={row['best_pfam_evalue']}" if row["best_pfam_evalue"]
                else "")
    if row["blast_acc"]:
        bits.append(f"blast {row['blast_acc']} {row['blast_pident']}%")
    bits.append(f"[{row['status']}]")
    return "  " + " ".join(b for b in bits if b)


def _fraction(row: dict) -> str:
    """n/N of the REQUIRED slots of a module (what the figures show); the direction
    and key-gene rules decide the status and are not slots."""
    n = len([x for x in row["requires_present"].split(",") if x])
    n_tot = n + len([x for x in row["requires_missing"].split(",") if x])
    return f"{n}/{n_tot}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--input-kind", required=True, choices=["protein", "nucleotide"],
                    help="proteome (.faa) or nucleotide assembly (.fna, Prodigal)")
    ap.add_argument("--prodigal-mode", default="", help="single|meta for nucleotide input")
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--complexes", required=True, type=Path)
    ap.add_argument("--synergies", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path,
                    help="targets.yaml (the rule depends on it; not read here)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    calls = load_calls(args.calls)
    complexes = load_complexes(args.complexes)
    synergies = load_complexes(args.synergies)

    if args.input_kind == "nucleotide":
        kind = "nucleotide (Prodigal" + (f", {args.prodigal_mode}" if args.prodigal_mode else "") + ")"
    else:
        kind = "proteome"
    lines: list[str] = []
    lines.append("mcycle-pipeline — methane-cycle gap analysis")
    lines.append(f"SAMPLE: {args.sample}   input: {kind}")
    lines.append("=" * 72)
    lines.append("")

    # --- PRESENT ----------------------------------------------------------
    present = [r for r in calls if r["status"] in PRESENT]
    lines.append(f"PRESENT ({len(present)} target(s)):")
    if not present:
        lines.append("  (nothing detected)")
    else:
        by_cat: dict[str, list[dict]] = {}
        for r in present:
            by_cat.setdefault(r["category"], []).append(r)
        for cat in by_cat:
            lines.append(f"  [{cat}]")
            for r in sorted(by_cat[cat], key=lambda x: x["target_id"]):
                lines.append("  " + fmt_present(r))
    lines.append("")

    # --- COMPLEXES --------------------------------------------------------
    complete = [c["complex_id"] for c in complexes if c["status"] == "complete"]
    partial = [c for c in complexes if c["status"] != "complete" and int(c["n_present"]) > 0]
    absent = [c["complex_id"] for c in complexes if int(c["n_present"]) == 0]
    lines.append(f"COMPLEXES ({len(complete)} complete, {len(partial)} partial, "
                 f"{len(absent)} absent of {len(complexes)}):")
    lines.append("  complete: " + (", ".join(complete) or "(none)"))
    if partial:
        lines.append("  partial (subunits missing):")
        for c in partial:
            lines.append(f"    {c['complex_id']}: {c['n_present']}/{c['n_total']} "
                         f"— missing {c['members_missing']} (present: {c['members_present']})")
            lines.append(f"      {c['application']}")
    else:
        lines.append("  partial: (none)")
    lines.append("  absent: " + (", ".join(absent) or "(none)"))
    lines.append("")

    # --- PROCESS MODULES --------------------------------------------------
    # Same reading as the figures (_completeness.load_synergies): a violated rule
    # makes a module "ruled out" only when some of its required genes are there;
    # with none of them it is simply absent.
    def n_required_present(s: dict) -> int:
        return len([x for x in s["requires_present"].split(",") if x])
    m_complete = [s for s in synergies if s["status"] == "complete"]
    m_ruled = [s for s in synergies if s["status"] != "complete" and s.get("forbids_violated")
               and n_required_present(s) > 0]
    m_absent = [s for s in synergies if s["status"] != "complete" and s not in m_ruled
                and n_required_present(s) == 0]
    m_partial = [s for s in synergies if s["status"] != "complete"
                 and s not in m_ruled and s not in m_absent]
    lines.append(f"PROCESS MODULES ({len(m_complete)} complete, {len(m_partial)} partial, "
                 f"{len(m_ruled)} ruled out, {len(m_absent)} absent of {len(synergies)}):")
    lines.append("  complete: " + (", ".join(s["synergy_id"] for s in m_complete) or "(none)"))
    if m_partial:
        lines.append("  partial (required slots missing; a+b = both needed, a|b = either):")
        for s in m_partial:
            lines.append(f"    {s['synergy_id']}: {_fraction(s)} — missing {s['requires_missing']}")
            lines.append(f"      {s['benefit']}")
    else:
        lines.append("  partial: (none)")
    if m_ruled:
        lines.append("  ruled out (a rule of the module is violated; its shared genes do not count):")
        for s in m_ruled:
            lines.append(f"    {s['synergy_id']}: {s['forbids_violated']} "
                         f"(shared genes present: {s['requires_present'] or 'none'})")
    else:
        lines.append("  ruled out: (none)")
    lines.append("  absent (no gene found): "
                 + (", ".join(s["synergy_id"] for s in m_absent) or "(none)"))
    lines.append("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n")
    print(f"[gap_analysis] wrote {args.out}")


if __name__ == "__main__":
    main()
