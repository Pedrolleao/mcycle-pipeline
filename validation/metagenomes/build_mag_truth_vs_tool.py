#!/usr/bin/env python3
"""
build_mag_truth_vs_tool.py — the MAG realism table (M8): published phenotype of each
MAG against what mcycle, METABOLIC and DRAM report, from nucleotide input.

Inputs
  validation/metagenomes/mag_panel.tsv   roster + expectations (fixed before the runs)
  results_mags/<mag>/                    mcycle, run with --prodigal-mode meta
  validation/metagenomes/checkm.tsv      CheckM lineage_wf --reduced_tree
  comparators/mags_metabolic_out/        METABOLIC v4.0 (worksheet 1 + KO lists)
  comparators/mags_dram_out/d?/distilled/product.tsv   DRAM v1.4.6
Comparator columns are left empty when their output is not there yet.

A MAG whose published phenotype needs an Mcr but whose assembly holds no homologue
of any Mcr subunit is marked `not in assembly`: no tool can call a gene that was
not binned.

Output: validation/metagenomes/mag_truth_vs_tool.tsv (+ a Markdown copy on stdout)
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RES = ROOT / "results_mags"
MCR_KOS = ("K00399", "K00401", "K00402")


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def family_scores(mag: str) -> dict[str, float]:
    best: dict[str, float] = {}
    for line in open(RES / mag / "hmm" / f"{mag}.hmmscan.tsv"):
        if not line.startswith("#"):
            p = line.split()
            best[p[0]] = max(best.get(p[0], 0.0), float(p[7]))
    return best


def training_genera() -> dict[str, set[str]]:
    out = {}
    for t in ("mcrA_anme", "pmoA"):
        out[t] = {m.group(1) for line in open(ROOT / "targets" / t / "refs.fasta")
                  if line.startswith(">") and (m := re.search(r"OS=(\S+)", line))}
    return out


def metabolic() -> dict[str, dict[str, bool]]:
    ws = ROOT / "comparators" / "mags_metabolic_out" / "worksheet1" / "mags.tsv"
    out: dict[str, dict[str, bool]] = {}
    if not ws.exists():
        return out
    with open(ws) as fh:
        rd = csv.reader(fh, delimiter="\t"); head = next(rd)
        cols = {h[: -len(" Hmm presence")]: i for i, h in enumerate(head) if h.endswith(" Hmm presence")}
        gi = head.index("Gene abbreviation")
        for r in rd:
            if r[gi] in ("mcrA", "pmoA", "amoA"):
                for g, i in cols.items():
                    out.setdefault(g, {})[r[gi]] = out.get(g, {}).get(r[gi], False) or r[i] == "Present"
    return out


def dram() -> dict[str, dict[str, bool]]:
    out: dict[str, dict[str, bool]] = {}
    for p in sorted((ROOT / "comparators" / "mags_dram_out").glob("d?/distilled/product.tsv")):
        for r in read_tsv(p):
            out[r["genome"]] = {
                "mcr": r.get("Methanogenesis and methanotrophy: Key functional gene", "") == "True",
                "pmo": r.get("Methanogenesis and methanotrophy: methane => methanol, with oxygen (pmo)", "") == "True"}
    return out


def main() -> int:
    roster = read_tsv(HERE / "mag_panel.tsv")
    checkm = {r["Bin Id"]: r for r in read_tsv(HERE / "checkm.tsv")} if (HERE / "checkm.tsv").exists() else {}
    tg = training_genera(); met = metabolic(); drm = dram()
    rows = []
    for r in roster:
        mag = r["name"]
        calls = {c["target_id"]: c for c in read_tsv(RES / mag / "calls" / "mcycle_calls.tsv")}
        fam = family_scores(mag)
        m = re.search(r"mcr_(methanogenic|reverse)", calls["mcrA"]["evidence_source"])
        direction = m.group(1) if m else "none"
        mmo = "pmmo" if calls["pmoA"]["status"] == "confirmed" else "none"
        mcr_in_assembly = any(fam.get(k, 0) > 0 for k in MCR_KOS)
        if direction == r["expected_direction"] and mmo == r["expected_mmo"]:
            verdict = "as published"
        elif r["expected_direction"] != "none" and not mcr_in_assembly:
            verdict = "Mcr not in assembly"
        else:
            verdict = "DIFFERS"
        model = calls["mcrA_anme"]["pfam_hits"] if calls["mcrA_anme"]["status"] == "confirmed" else (
            calls["pmoA"]["pfam_hits"] if mmo == "pmmo" else "")
        genus_tokens = set(re.findall(r"[A-Z][a-z]+", r["lineage"]))
        in_training = sorted(g for t in tg.values() for g in t if g in genus_tokens)
        c = checkm.get(mag, {})
        rows.append({
            "mag": mag, "accession": r["source"], "contigs": r["contigs"],
            "checkm_completeness": c.get("Completeness", ""), "checkm_contamination": c.get("Contamination", ""),
            "published_phenotype": r["phenotype"],
            "expected_direction": r["expected_direction"], "mcycle_direction": direction,
            "expected_mmo": r["expected_mmo"], "mcycle_mmo": mmo, "mcycle_clade_model": model,
            "mcycle_mcrB": calls["mcrB"]["status"] + ("" if "acr_like" not in calls["mcrB"]["evidence_source"] else " (alkyl-CoM reductase)"),
            "mcycle_verdict": verdict,
            "genus_in_hmm_training": ",".join(in_training) or "no",
            "metabolic_mcrA": {True: "Present", False: "Absent"}.get(met.get(mag, {}).get("mcrA"), ""),
            "metabolic_pmoA": {True: "Present", False: "Absent"}.get(met.get(mag, {}).get("pmoA"), ""),
            "metabolic_amoA": {True: "Present", False: "Absent"}.get(met.get(mag, {}).get("amoA"), ""),
            "dram_mcr": {True: "True", False: "False"}.get(drm.get(mag, {}).get("mcr"), ""),
            "dram_pmo": {True: "True", False: "False"}.get(drm.get(mag, {}).get("pmo"), ""),
            "citation": r["citation"], "doi": r["doi"]})
    out = HERE / "mag_truth_vs_tool.tsv"
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    show = ["mag", "checkm_completeness", "checkm_contamination", "expected_direction", "mcycle_direction",
            "expected_mmo", "mcycle_mmo", "mcycle_clade_model", "mcycle_verdict", "genus_in_hmm_training",
            "metabolic_mcrA", "metabolic_pmoA", "dram_mcr", "dram_pmo"]
    print("| " + " | ".join(show) + " |"); print("|" + "---|" * len(show))
    for r in rows:
        print("| " + " | ".join(str(r[k]) for k in show) + " |")
    n = sum(r["mcycle_verdict"] == "as published" for r in rows)
    print(f"\nmcycle: {n} / {len(rows)} as published; "
          f"{sum(r['mcycle_verdict'] == 'Mcr not in assembly' for r in rows)} with the Mcr operon "
          f"missing from the assembly; {sum(r['mcycle_verdict'] == 'DIFFERS' for r in rows)} differ.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
