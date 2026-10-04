#!/usr/bin/env python3
"""
prepare_proteomes.py — proteomes of the 500 GTDB genomes of selection.tsv.

Backbone (380): the Prodigal proteomes ncycle's gtdb500 study already holds are
linked, so that mcycle, MCycDB and the METABOLIC output reused from that study all
refer to the same proteins. Enriched (120): genome sequence from NCBI (Datasets API;
genomes already cached by harvest_clade_refs.py are reused), genes called with
Prodigal (-p meta). File names follow the sister studies: GCA_000403075_1.faa.

Scratch goes to comparators/gtdb500_m/{gtdb_dl,proteomes}/ (not tracked).
"""
from __future__ import annotations

import csv
import io
import subprocess
import sys
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
N_PROT = Path("/home/dmin/Grants/Nitrogen_Cycle/ncycle-pipeline/comparators/gtdb500/proteomes")
CLADE_CACHE = (HERE / ".." / ".." / ".." / "clade_refs" / "genomes").resolve()
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"


def main() -> int:
    rows = list(csv.DictReader(open(HERE / "selection.tsv"), delimiter="\t"))
    prot = HERE / "proteomes"; dl = HERE / "gtdb_dl"
    prot.mkdir(exist_ok=True); dl.mkdir(exist_ok=True)
    missing_backbone, todo = [], []
    for r in rows:
        name = r["ncbi_acc"].replace(".", "_")
        dst = prot / f"{name}.faa"
        if dst.exists():
            continue
        if r["stratum"] == "backbone":
            src = N_PROT / f"{name}.faa"
            if src.exists():
                dst.symlink_to(src)
            else:
                missing_backbone.append(r["ncbi_acc"])
        else:
            todo.append(r["ncbi_acc"])
    todo += missing_backbone
    need = []
    for acc in todo:
        fna = dl / f"{acc}.fna"
        if fna.exists():
            continue
        if (CLADE_CACHE / f"{acc}.fna").exists():
            fna.symlink_to(CLADE_CACHE / f"{acc}.fna")
        else:
            need.append(acc)
    print(f"[prepare] {len(todo)} proteomes to build; {len(need)} genomes to download")
    for i in range(0, len(need), 15):
        batch = need[i:i + 15]
        url = API + ",".join(batch) + "/download?include_annotation_type=GENOME_FASTA"
        for attempt in range(3):
            try:
                z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(url, timeout=600).read()))
                break
            except Exception as e:                       # noqa: BLE001
                print(f"  retry {attempt + 1}: {e}", file=sys.stderr); time.sleep(5)
        else:
            print(f"  ! batch failed: {batch}", file=sys.stderr); continue
        for n in z.namelist():
            if n.endswith("_genomic.fna"):
                (dl / f"{n.split('/')[-2]}.fna").write_bytes(z.read(n))
        print(f"  fetched {min(i + 15, len(need))}/{len(need)}", flush=True); time.sleep(1)

    def call(acc: str) -> None:
        fna = dl / f"{acc}.fna"; out = prot / f"{acc.replace('.', '_')}.faa"
        if out.exists() or not fna.exists():
            return
        tmp = out.with_suffix(".faa.tmp")
        subprocess.run(["prodigal", "-p", "meta", "-q", "-i", str(fna), "-a", str(tmp),
                        "-o", "/dev/null"], check=True)
        tmp.rename(out)
    with ThreadPoolExecutor(max_workers=24) as pool:
        list(pool.map(call, todo))
    have = sum((prot / f"{r['ncbi_acc'].replace('.', '_')}.faa").exists() for r in rows)
    print(f"[prepare] {have} / {len(rows)} proteomes ready in {prot}")
    return 0 if have == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
