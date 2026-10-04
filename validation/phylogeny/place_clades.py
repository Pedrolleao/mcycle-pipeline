#!/usr/bin/env python3
"""
place_clades.py — phylogenetic placement of the GTDB-500 McrA and PmoA / AmoA
sequences among typed references, as a check on the clade-HMM calls.

For one trap (`--trap mcr` or `--trap pmo`):
  references  the typed sequences of targets/<id>/tc_calibration.tsv (clade known
              from the GTDB lineage of the genome they come from), de-replicated
              and capped per group;
  queries     the protein behind mcycle's call (or disqualified family hit) in
              every GTDB-500 genome that did not supply a reference;
  tree        MAFFT --auto, ClipKIT (kpic-smart-gap), IQ-TREE 2 (ModelFinder,
              1,000 ultrafast bootstraps), midpoint-rooted;
  placement   each query takes the type of the references in the smallest clade
              that contains it and at least one reference: a single type -> that
              type, several -> `mixed`. Support = UFBoot of that clade.

The gene tree is not independent of the HMM (both read the same protein), so this
does not measure accuracy. It answers two narrower questions: does the HMM call
agree with where the protein sits in its family tree, and which calls does NO
sequence method resolve (query in a mixed or weakly supported clade)?

  MCYCLE_RESULTS=results_gtdb500 python validation/phylogeny/place_clades.py --trap mcr \
      --selection comparators/gtdb500_m/selection.tsv --outdir comparators/gtdb500_m
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation"))
import score_mcycle as S  # noqa: E402

ENVS = Path.home() / "miniconda3" / "envs"
CLIPKIT = ENVS / "Clipkit_env" / "bin" / "clipkit"
IQTREE = ENVS / "Iqtree_env" / "bin" / "iqtree2"
TRAPS = {"mcr": ("mcrA_anme", ["mcrA_anme", "mcrA"], 400),
         "pmo": ("pmoA", ["pmoA"], 180)}
MAX_REFS_PER_GROUP = 14
SUPPORT = 95


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name:
            seqs[name] += line.strip().rstrip("*")
    return seqs


def ref_type(trap: str, group: str, use: str) -> str:
    if trap == "mcr":
        return "ANME" if use == "positive" else ("alkyl-CoM" if group == "acr" else "methanogen")
    if use == "unresolved":
        return "paralogue"
    return "pmoA" if use == "positive" else ("hydrocarbon" if group == "cumo_hc" else "amoA")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--trap", choices=list(TRAPS), required=True)
    ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--threads", default="16")
    args = ap.parse_args()
    target, call_ids, min_len = TRAPS[args.trap]
    tdir = ROOT / "targets" / target
    work = args.outdir / f"placement_{args.trap}"
    work.mkdir(parents=True, exist_ok=True)

    # ---- references: typed, capped per group, longest first
    harvest = read_fasta(tdir / "harvest.faa")
    cal = list(csv.DictReader(open(tdir / "tc_calibration.tsv"), delimiter="\t"))
    by_group = defaultdict(list)
    for r in cal:
        by_group[r["group"] + ("" if r["use"] != "unresolved" else "_paralogue")].append(r)
    refs, types, ref_asm = {}, {}, set()
    for group, rows in sorted(by_group.items()):
        seen_genus: Counter = Counter()
        for r in sorted(rows, key=lambda r: (-int(r["length"]), r["seq_id"])):
            if sum(1 for k in refs if types[k][1] == group) >= MAX_REFS_PER_GROUP:
                break
            if seen_genus[r["genus"]] >= 2:
                continue
            seen_genus[r["genus"]] += 1
            name = f"REF_{len(refs):03d}"
            refs[name] = harvest[r["seq_id"]]
            types[name] = (ref_type(args.trap, r["group"], r["use"]), group, r["genus"])
    for r in cal:
        ref_asm.add(r["seq_id"].split("|")[0].replace(".", "_"))

    # ---- queries: the protein behind the call in genomes that gave no reference
    sel = {r["ncbi_acc"].replace(".", "_"): r
           for r in csv.DictReader(open(args.selection), delimiter="\t")}
    prot_dir = args.selection.parent / "proteomes"
    queries, qmeta = {}, {}
    for g, r in sorted(sel.items()):
        f = S.RESULTS / g / "calls" / "mcycle_calls.tsv"
        if g in ref_asm or not f.exists():
            continue
        calls = {c["target_id"]: c for c in S.read_tsv(f)}
        row = next((calls[t] for t in call_ids
                    if calls[t]["status"] in ("confirmed", "disqualified", "domain-only")
                    and calls[t]["protein_id"]), None)
        if row is None:
            continue
        seq = read_fasta(prot_dir / f"{g}.faa").get(row["protein_id"], "")
        if len(seq) < min_len:
            continue
        if args.trap == "mcr":
            call = ("ANME" if calls["mcrA_anme"]["status"] == "confirmed" else
                    "methanogen" if calls["mcrA"]["status"] in S.PRESENT else "not McrA")
        else:
            call = "pmoA" if calls["pmoA"]["status"] == "confirmed" else "not pmoA"
        name = f"Q_{len(queries):03d}"
        queries[name] = seq
        qmeta[name] = {"genome": g, "clade": r["clade"], "stratum": r["stratum"],
                       "genus": r["gtdb_taxonomy"].split(";")[5][3:], "mcycle": call,
                       "model": row["pfam_hits"]}
    print(f"[place] {args.trap}: {len(refs)} references "
          f"({dict(Counter(t[0] for t in types.values()))}), {len(queries)} queries")

    # ---- alignment, trimming, tree
    faa, aln, trim = work / "all.faa", work / "all.aln", work / "all.trim.aln"
    with open(faa, "w") as fh:
        for k, v in {**refs, **queries}.items():
            fh.write(f">{k}\n{v}\n")
    tree_file = work / "all.treefile"
    if not tree_file.exists():
        with open(aln, "w") as fh:
            subprocess.run(["mafft", "--auto", "--quiet", "--thread", args.threads, str(faa)],
                           check=True, stdout=fh)
        subprocess.run([str(CLIPKIT), str(aln), "-m", "kpic-smart-gap", "-o", str(trim)],
                       check=True, stdout=subprocess.DEVNULL)
        subprocess.run([str(IQTREE), "-s", str(trim), "--prefix", str(work / "all"), "-m", "MFP",
                        "-B", "1000", "-T", args.threads, "--seed", "1234", "--quiet",
                        "-redo"], check=True)

    tree = Phylo.read(str(tree_file), "newick")
    tree.root_at_midpoint()
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parents[child] = clade
    out_rows = []
    for leaf in tree.get_terminals():
        if leaf.name not in queries:
            continue
        node, support = leaf, None
        while node in parents:
            node = parents[node]
            found = {types[t.name][0] for t in node.get_terminals() if t.name in refs}
            if found:
                break
        found = {types[t.name][0] for t in node.get_terminals() if t.name in refs}
        support = node.confidence if node.confidence is not None else (
            float(node.name) if node.name and node.name.replace(".", "").isdigit() else None)
        placed = next(iter(found)) if len(found) == 1 else "mixed"
        m = qmeta[leaf.name]
        confident = support is not None and support >= SUPPORT and placed != "mixed"
        expected_call = {"ANME": "ANME", "methanogen": "methanogen", "alkyl-CoM": "not McrA",
                         "pmoA": "pmoA", "amoA": "not pmoA", "hydrocarbon": "not pmoA",
                         "paralogue": "not pmoA"}.get(placed, "")
        out_rows.append({**m, "placement": placed,
                         "clade_types": ",".join(sorted(found)),
                         "ufboot": "" if support is None else f"{support:.0f}",
                         "confident": "yes" if confident else "no",
                         "agree": "" if placed == "mixed" else str(int(expected_call == m["mcycle"]))})

    tsv = args.outdir / f"DIR_PLACEMENT_{args.trap}.tsv"
    with open(tsv, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(sorted(out_rows, key=lambda r: (r["clade"], r["genome"])))

    conf = [r for r in out_rows if r["confident"] == "yes"]
    weak = [r for r in out_rows if r["confident"] == "no" and r["placement"] != "mixed"]
    mixed = [r for r in out_rows if r["placement"] == "mixed"]
    L = [f"## Placement — {'McrA' if args.trap == 'mcr' else 'PmoA / AmoA'}\n",
         f"{len(refs)} typed references, {len(out_rows)} query sequences from GTDB-500 genomes "
         f"that supplied no reference. IQ-TREE 2 (ModelFinder, 1,000 UFBoot), midpoint root.\n",
         "| placement | queries | mcycle agrees |", "|---|---|---|",
         f"| single-type clade, UFBoot >= {SUPPORT} | {len(conf)} | "
         f"{sum(r['agree'] == '1' for r in conf)} |",
         f"| single-type clade, UFBoot < {SUPPORT} | {len(weak)} | "
         f"{sum(r['agree'] == '1' for r in weak)} |",
         f"| mixed clade (no sequence method resolves it) | {len(mixed)} | n/a |", ""]
    cats = sorted({r["mcycle"] for r in out_rows})
    L += ["Placement (rows) x mcycle call (columns):\n",
          "| placement | " + " | ".join(cats) + " |", "|---|" + "---|" * len(cats)]
    for p in sorted({r["placement"] for r in out_rows}):
        L.append(f"| {p} | " + " | ".join(
            str(sum(r["placement"] == p and r["mcycle"] == c for r in out_rows)) for c in cats) + " |")
    bad = [r for r in out_rows if r["agree"] == "0"]
    if bad:
        L.append(f"\nDisagreements ({len(bad)}):\n")
        for r in bad:
            L.append(f"- {r['genome']} ({r['clade']}; {r['genus']}): placed {r['placement']} "
                     f"(UFBoot {r['ufboot'] or 'n/a'}), mcycle {r['mcycle']}")
    if mixed:
        L.append(f"\nIn a mixed clade ({len(mixed)}):\n")
        for r in mixed:
            L.append(f"- {r['genome']} ({r['clade']}; {r['genus']}): clade holds "
                     f"{r['clade_types']}; mcycle {r['mcycle']}")
    (args.outdir / f"DIR_PLACEMENT_{args.trap}.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
