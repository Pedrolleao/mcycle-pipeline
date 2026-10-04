#!/usr/bin/env python3
"""
build_clade_hmms.py — build, cross-validate and calibrate the clade HMMs of a
clade-gated target (pmoA, mcrA_anme) from the sequences harvest_clade_refs.py
collected.

A single profile over every methanotroph PmoA is the family profile again (it is
what K10944 is) and cannot exclude the AmoA of a gammaproteobacterial ammonia
oxidizer, which sits inside the family. One profile PER CLADE can: each is narrow,
and each gets its own threshold. A protein is called when it passes any of them.

For every positive group of targets/<id>/clades.yaml:
  1. training set — the harvested sequences of that group, minus the ones the
     `positive_rule` of clades.yaml rejects (pmoA: paralogues under 60 % identity to
     every curated PmoA seed — pxmA, pmoA3 — whose function is unresolved; they are
     kept out of positives AND negatives), de-replicated with CD-HIT at 97 %;
  2. model — MAFFT alignment, hmmbuild (-n `<id>__<clade>`);
  3. leave-one-genus-out — for every genus of the clade (every species when the
     clade is one genus) the model is rebuilt without it and the left-out
     sequences are scored: the score a NEW genus of the clade can be expected to
     reach;
  4. threshold — halfway between the best-scoring negative and the worst
     left-out positive. If they overlap the threshold goes just above the best
     negative (precision first) and the manifest says which left-out genera fall
     under it.

Writes, under targets/<id>/:
  <id>__<clade>.hmm        the models (spliced into the scan database by build_hmm_db.py)
  refs.fasta               every training sequence, `>id||seq clade=… OS=Genus species OX=gtdb`
                           (read by validation/detect_seed_leakage.py)
  tc_calibration.tsv       every harvested sequence under every model, with the
                           leave-one-out score where it has one
  manifest.yaml            provenance, per-clade numbers, thresholds and rationale

    python workflow/scripts/build_clade_hmms.py --target pmoA
"""
from __future__ import annotations

import argparse
import csv
import datetime
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
SEEDS = ROOT / "resources" / "blast_db" / "blast_gated_refs.fasta"


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, check=True, **kw)


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name:
            seqs[name] += line.strip()
    return seqs


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")


def seed_identity(target: str, seqs: dict[str, str], work: Path) -> dict[str, float]:
    """Best identity of each sequence to a curated BLAST seed of the target."""
    seeds = {k: v for k, v in read_fasta(SEEDS).items() if k.startswith(f"{target}||")}
    if not seeds:
        return {}
    write_fasta(work / "seeds.faa", seeds); write_fasta(work / "q.faa", seqs)
    run(["diamond", "makedb", "--in", str(work / "seeds.faa"), "-d", str(work / "seeds"),
         "--quiet"])
    run(["diamond", "blastp", "-q", str(work / "q.faa"), "-d", str(work / "seeds"),
         "-o", str(work / "q.tsv"), "--outfmt", "6", "qseqid", "pident", "--quiet",
         "--query-cover", "50", "--max-target-seqs", "25", "--very-sensitive"])
    best: dict[str, float] = {}
    for line in open(work / "q.tsv"):
        q, pid = line.split("\t")
        best[q] = max(best.get(q, 0.0), float(pid))
    return best


def build_hmm(name: str, seqs: dict[str, str], work: Path, out: Path) -> None:
    write_fasta(work / "in.faa", seqs)
    with open(work / "in.aln", "w") as fh:
        run(["mafft", "--auto", "--quiet", str(work / "in.faa")], stdout=fh)
    run(["hmmbuild", "--amino", "--informat", "afa", "-n", name, str(out), str(work / "in.aln")],
        stdout=subprocess.DEVNULL)


def scores(hmm: Path, faa: Path, work: Path) -> dict[str, float]:
    tbl = work / "s.tbl"
    run(["hmmsearch", "--tblout", str(tbl), "-E", "10", str(hmm), str(faa)],
        stdout=subprocess.DEVNULL)
    out = {}
    for line in open(tbl):
        if not line.startswith("#"):
            f = line.split()
            out[f[0]] = max(out.get(f[0], 0.0), float(f[5]))
    return out


def dereplicate(seqs: dict[str, str], work: Path, ident: float) -> dict[str, str]:
    write_fasta(work / "d.faa", seqs)
    run(["cd-hit", "-i", str(work / "d.faa"), "-o", str(work / "d.out"), "-c", str(ident),
         "-d", "0"], stdout=subprocess.DEVNULL)
    return read_fasta(work / "d.out")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    args = ap.parse_args()
    tid = args.target
    tdir = ROOT / "targets" / tid
    cfg = yaml.safe_load(open(tdir / "clades.yaml"))
    rows = list(csv.DictReader(open(tdir / "harvest.tsv"), delimiter="\t"))
    seqs = read_fasta(tdir / "harvest.faa")
    meta = {r["seq_id"]: r for r in rows}
    for r in rows:
        lin = r["gtdb_lineage"].split(";")
        r["genus"] = lin[5][3:] if len(lin) > 5 else "?"
        r["species"] = lin[6][3:] if len(lin) > 6 else "?"

    work = Path(tempfile.mkdtemp(prefix=f"clade_{tid}_"))
    rule = cfg.get("positive_rule") or {}
    ident = seed_identity(tid, seqs, work) if rule.get("seed_identity_min") else {}
    for r in rows:
        r["seed_identity"] = ident.get(r["seq_id"], 0.0)
        r["use"] = r["role"]
        if r["role"] == "positive" and rule.get("seed_identity_min") and \
                r["seed_identity"] < rule["seed_identity_min"]:
            r["use"] = "unresolved"      # divergent paralogue: neither positive nor negative

    write_fasta(work / "all.faa", seqs)
    for old in tdir.glob(f"{tid}__*.hmm"):
        old.unlink()
    clades = [g["name"] for g in cfg["groups"] if g["role"] == "positive"]
    table: dict[str, dict[str, float]] = {}
    logo: dict[str, dict[str, float]] = {c: {} for c in clades}
    manifest_clades, refs = {}, {}
    for clade in clades:
        pos = {r["seq_id"]: seqs[r["seq_id"]] for r in rows
               if r["group"] == clade and r["use"] == "positive"}
        train = dereplicate(pos, work, cfg.get("cdhit_identity", 0.97))
        if len(train) < 3:
            print(f"  ! {clade}: only {len(train)} training sequences — no model built",
                  file=sys.stderr)
            continue
        hmm = tdir / f"{tid}__{clade}.hmm"
        build_hmm(f"{tid}__{clade}", train, work, hmm)
        table[clade] = scores(hmm, work / "all.faa", work)
        for sid in train:
            m = meta[sid]
            refs[f"{tid}||{sid} clade={clade} OS={m['genus'].split('_')[0]} "
                 f"{m['species'].split(' ')[-1]} OX=gtdb"] = seqs[sid]

        # leave-one-genus-out (leave-one-species-out when the clade is one genus)
        genera = sorted({meta[s]["genus"] for s in train})
        unit = "genus" if len(genera) > 1 else "species"
        for left in sorted({meta[s][unit] for s in train}):
            rest = {s: v for s, v in train.items() if meta[s][unit] != left}
            out = {s: v for s, v in pos.items() if meta[s][unit] == left}
            if len(rest) < 3 or not out:
                continue
            build_hmm("logo", rest, work, work / "logo.hmm")
            write_fasta(work / "out.faa", out)
            sc = scores(work / "logo.hmm", work / "out.faa", work)
            for s in out:
                logo[clade][s] = sc.get(s, 0.0)

        neg = {r["seq_id"]: table[clade].get(r["seq_id"], 0.0) for r in rows
               if r["use"] == "negative"}
        max_neg_id = max(neg, key=neg.get); max_neg = neg[max_neg_id]
        resub = [table[clade].get(s, 0.0) for s in pos]
        lo_pos = min(logo[clade].values()) if logo[clade] else min(resub)
        if lo_pos > max_neg:
            tc = (lo_pos + max_neg) / 2
            why = (f"midpoint of the best negative ({max_neg:.1f}, {meta[max_neg_id]['group']} "
                   f"{meta[max_neg_id]['genus']}) and the worst left-out positive ({lo_pos:.1f})")
            under = []
        else:
            tc = max_neg * 1.02
            under = sorted({meta[s][unit] for s, v in logo[clade].items() if v < tc})
            why = (f"OVERLAP: best negative {max_neg:.1f} ({meta[max_neg_id]['group']} "
                   f"{meta[max_neg_id]['genus']}) >= worst left-out positive {lo_pos:.1f}; "
                   f"threshold set 2 % above the best negative; left-out {unit} "
                   f"under it: {', '.join(under) or 'none'}")
        manifest_clades[clade] = {
            "n_harvested": len(pos), "n_training": len(train),
            "n_genera": len(genera), "leave_out_unit": unit,
            "tc_bitscore": round(tc, 1), "tc_rationale": why,
            "max_negative": round(max_neg, 1),
            "min_training_score": round(min(resub), 1),
            "min_left_out_score": round(lo_pos, 1),
            "left_out_passing": f"{sum(v >= tc for v in logo[clade].values())}/{len(logo[clade])}",
        }
        print(f"  {tid}__{clade:<8} train={len(train):<3} genera={len(genera):<3} "
              f"max_neg={max_neg:7.1f}  min_left_out={lo_pos:7.1f}  TC={tc:7.1f}  "
              f"left-out passing {manifest_clades[clade]['left_out_passing']}")

    write_fasta(tdir / "refs.fasta", refs)
    with open(tdir / "tc_calibration.tsv", "w") as fh:
        built = [c for c in clades if c in table]
        fh.write("seq_id\tgroup\tuse\tgenus\tspecies\tlength\tko\tko_score\tseed_identity\t"
                 + "\t".join(f"{c}\t{c}_left_out" for c in built) + "\n")
        for r in rows:
            fh.write("\t".join([r["seq_id"], r["group"], r["use"], r["genus"], r["species"],
                                r["length"], r["ko"], r["score"], f"{r['seed_identity']:.1f}"]
                               + [x for c in built for x in (
                                   f"{table[c].get(r['seq_id'], 0.0):.1f}",
                                   f"{logo[c][r['seq_id']]:.1f}" if r["seq_id"] in logo[c] else "")])
                     + "\n")
    n_unres = sum(r["use"] == "unresolved" for r in rows)
    manifest = {
        "target_id": tid, "tier": 2,
        "method": "clade HMMs from GTDB-typed genomes (harvest_clade_refs.py, "
                  "build_clade_hmms.py)",
        "source": "GTDB r232 taxonomy; genomes from NCBI; genes called with Prodigal (meta)",
        "holdout_exclusion": "no sequence from a hold-out genus of validation/panel.tsv "
                             "(NCBI or GTDB naming) — enforced by harvest_clade_refs.py",
        "positive_rule": rule or "every harvested homologue of a positive group",
        "n_sequences_harvested": len(rows),
        "n_negatives": sum(r["use"] == "negative" for r in rows),
        "n_unresolved_paralogues": n_unres,
        "clade_models": manifest_clades,
        "curator": "pedro.leao.ru@gmail.com",
        "curation_date": datetime.date.today().isoformat(),
    }
    with open(tdir / "manifest.yaml", "w") as fh:
        yaml.safe_dump(manifest, fh, sort_keys=False, width=100, allow_unicode=True)
    shutil.rmtree(work, ignore_errors=True)
    print(f"[clade_hmms] {tid}: {len(manifest_clades)} models, {len(refs)} training "
          f"sequences, {n_unres} unresolved paralogues -> {tdir}")


if __name__ == "__main__":
    main()
