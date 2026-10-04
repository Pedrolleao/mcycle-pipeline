#!/usr/bin/env python3
"""
harvest_clade_refs.py — collect typed reference sequences for a clade-gated target
from GTDB-classified genomes.

The clade HMMs of pmoA (methanotroph clades vs ammonia / hydrocarbon monooxygenases)
and mcrA_anme (ANME clades vs methanogen McrA) need positives and negatives whose
clade is known independently of any sequence annotation. UniProt names of these
proteins are unreliable (PmoA and AmoA share one name; ANME proteins are filed
under "uncultured archaeon"). The clade of a GENOME is not: GTDB places it from
120 / 53 marker genes. So the sequences are taken from genomes:

  1. select  one genome per GTDB species from the lineages listed in
             targets/<id>/clades.yaml, at most `per_genus` per genus and
             `max_genomes` per group — never from a hold-out genus of
             validation/panel.tsv (NCBI or GTDB naming);
  2. fetch   the genome sequences from NCBI (Datasets API);
  3. extract call genes (Prodigal, meta mode) and keep every protein that hits one
             of the target's KO profiles at >= `min_fraction` of its threshold and
             is >= `min_len_fraction` of the profile length.

Output, under targets/<id>/:
  harvest.tsv   one row per sequence: id, assembly, GTDB lineage, group, role
                (positive | negative), length, KO, score
  harvest.faa   the sequences; header `>id group=<group> role=<role> asm=<accession>`
Genomes and gene calls are cached in ../clade_refs/ (outside the repository).

    python workflow/scripts/harvest_clade_refs.py --target pmoA
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
CACHE = (ROOT / ".." / "clade_refs").resolve()
GTDB = [Path("/home/dmin/Grants/Nitrogen_Cycle/ncycle-pipeline/comparators/gtdb_pilot")
        / f for f in ("ar53_taxonomy.tsv", "bac120_taxonomy.tsv")]
KO_HMMS = ROOT / "resources" / ".cache" / "ko_hmms"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def banned_genera() -> set[str]:
    """Base names (lower case) of every hold-out genus, NCBI and GTDB naming."""
    banned = set()
    hold = {r["source"].split("_")[1].split(".")[0]: r
            for r in read_tsv(ROOT / "validation" / "panel.tsv") if r["split"] == "holdout"}
    for r in hold.values():
        banned.add(r["taxa"].split(";")[0])          # first token = NCBI genus
    for f in GTDB:
        for line in open(f):
            acc, tax = line.rstrip("\n").split("\t")
            if acc.split("_")[2].split(".")[0] in hold:
                banned.add(tax.split(";")[5][3:].split("_")[0].lower())
    banned.discard("")
    # tokens that are not genera (isolate labels of uncultured hold-out genomes)
    return {b for b in banned if b not in {"g37anme1"}} | {"mycolicibacterium"}


def select(groups: list[dict], banned: set[str]) -> list[dict]:
    lineages: dict[str, tuple[str, str]] = {}
    for f in GTDB:
        for line in open(f):
            acc, tax = line.rstrip("\n").split("\t")
            lineages[acc] = tax
    picked = []
    for grp in groups:
        per_species: dict[str, str] = {}
        for acc, tax in lineages.items():
            if not any(f"{lin};" in tax + ";" for lin in grp["lineages"]):
                continue
            if any(f"{lin};" in tax + ";" for lin in grp.get("exclude", [])):
                continue
            genus = tax.split(";")[5][3:]
            if genus.split("_")[0].lower() in banned:
                continue
            sp = tax.split(";")[6]
            # prefer RefSeq assemblies, then the lower accession (stable choice)
            rank = lambda a: (not a.startswith("RS_"), a)
            if sp not in per_species or rank(acc) < rank(per_species[sp]):
                per_species[sp] = acc
        by_genus: dict[str, list[str]] = {}
        for sp, acc in sorted(per_species.items(),
                              key=lambda kv: (not kv[1].startswith("RS_"), kv[0])):
            by_genus.setdefault(lineages[acc].split(";")[5], []).append(acc)
        chosen: list[str] = []
        # round-robin over genera so that the cap does not favour big genera
        for i in range(grp.get("per_genus", 2)):
            for genus in sorted(by_genus):
                if i < len(by_genus[genus]) and len(chosen) < grp.get("max_genomes", 20):
                    chosen.append(by_genus[genus][i])
        for acc in chosen:
            picked.append({"assembly": acc[3:], "group": grp["name"], "role": grp["role"],
                           "lineage": lineages[acc]})
    return picked


def fetch(accessions: list[str]) -> None:
    gdir = CACHE / "genomes"
    gdir.mkdir(parents=True, exist_ok=True)
    todo = [a for a in accessions if not (gdir / f"{a}.fna").exists()]
    for i in range(0, len(todo), 15):
        batch = todo[i:i + 15]
        url = API + ",".join(batch) + "/download?include_annotation_type=GENOME_FASTA"
        for attempt in range(3):
            try:
                data = urllib.request.urlopen(url, timeout=600).read()
                z = zipfile.ZipFile(io.BytesIO(data))
                break
            except Exception as e:                      # noqa: BLE001
                print(f"  retry {attempt + 1} ({e})", file=sys.stderr)
                time.sleep(5)
        else:
            print(f"  ! batch failed: {batch}", file=sys.stderr)
            continue
        for name in z.namelist():
            if name.endswith("_genomic.fna"):
                acc = name.split("/")[-2]
                (gdir / f"{acc}.fna").write_bytes(z.read(name))
        print(f"  fetched {min(i + 15, len(todo))}/{len(todo)}", flush=True)
        time.sleep(1)


def call_genes(acc: str, threads: int) -> Path | None:
    fna = CACHE / "genomes" / f"{acc}.fna"
    faa = CACHE / "proteins" / f"{acc}.faa"
    if not fna.exists():
        return None
    if not faa.exists():
        faa.parent.mkdir(parents=True, exist_ok=True)
        tmp = faa.with_suffix(".faa.tmp")       # never leave a half-written proteome
        subprocess.run(["prodigal", "-p", "meta", "-q", "-i", str(fna), "-a", str(tmp),
                        "-o", "/dev/null"], check=True)
        tmp.rename(faa)
    return faa


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        elif name:
            seqs[name] += line.strip().rstrip("*")
    return seqs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--select-only", action="store_true")
    args = ap.parse_args()

    tdir = ROOT / "targets" / args.target
    cfg = yaml.safe_load(open(tdir / "clades.yaml"))
    banned = banned_genera()
    print(f"[harvest] hold-out genera excluded: {', '.join(sorted(banned))}")
    picked = select(cfg["groups"], banned)
    for grp in cfg["groups"]:
        n = [p for p in picked if p["group"] == grp["name"]]
        print(f"  {grp['name']:<14}{grp['role']:<9}{len(n):>4} genomes, "
              f"{len({p['lineage'].split(';')[5] for p in n})} genera")
    for extra in cfg.get("extra_genomes", []):
        picked.append({"assembly": extra["assembly"], "group": extra["group"],
                       "role": extra["role"], "lineage": extra.get("lineage", "")})
    if args.select_only:
        return
    fetch([p["assembly"] for p in picked])

    tc = {r["profile_id"]: float(r["tc"]) for r in
          read_tsv(ROOT / "resources" / "hmm" / "tc_cutoffs.tsv") if r["tc"]}
    kofam_tc = {}
    for line in open(ROOT / "resources" / ".cache" / "ko_list"):
        p = line.split("\t")
        if p[0] in cfg["ko"] and p[1] not in ("-", "threshold"):
            kofam_tc[p[0]] = float(p[1])
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=args.threads) as pool:   # gene calling, in parallel
        list(pool.map(lambda p: call_genes(p["assembly"], 1), picked))
    rows, seqs_out = [], {}
    for n, p in enumerate(picked, 1):
        faa = call_genes(p["assembly"], args.threads)
        if faa is None:
            print(f"  ! no genome for {p['assembly']}", file=sys.stderr)
            continue
        seqs = read_fasta(faa)
        best: dict[str, tuple[float, str]] = {}
        for ko in cfg["ko"]:
            tbl = CACHE / "hits" / f"{p['assembly']}.{ko}.tbl"
            if not tbl.exists():
                tbl.parent.mkdir(parents=True, exist_ok=True)
                subprocess.run(["hmmsearch", "--cpu", "4", "-E", "1e-5",
                                "--tblout", str(tbl), str(KO_HMMS / f"{ko}.hmm"), str(faa)],
                               check=True, stdout=subprocess.DEVNULL)
            for line in open(tbl):
                if line.startswith("#"):
                    continue
                f = line.split()
                score = float(f[5])
                if score >= cfg["min_fraction"] * kofam_tc.get(ko, tc.get(ko, 1e9)):
                    if f[0] not in best or score / kofam_tc[ko] > best[f[0]][0]:
                        best[f[0]] = (score / kofam_tc[ko], f"{ko}\t{score:.1f}")
        for prot, (_ratio, ko_score) in sorted(best.items()):
            if len(seqs[prot]) < cfg["min_len"]:
                continue
            sid = f"{p['assembly']}|{prot}"
            rows.append([sid, p["assembly"], p["lineage"], p["group"], p["role"],
                         str(len(seqs[prot])), ko_score])
            seqs_out[sid] = (p, seqs[prot])
        if n % 25 == 0:
            print(f"  scanned {n}/{len(picked)}", flush=True)

    with open(tdir / "harvest.tsv", "w") as fh:
        fh.write("seq_id\tassembly\tgtdb_lineage\tgroup\trole\tlength\tko\tscore\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    with open(tdir / "harvest.faa", "w") as fh:
        for sid, (p, seq) in seqs_out.items():
            fh.write(f">{sid} group={p['group']} role={p['role']} asm={p['assembly']}\n{seq}\n")
    print(f"[harvest] {len(rows)} sequences from {len(picked)} genomes -> {tdir}/harvest.*")


if __name__ == "__main__":
    main()
