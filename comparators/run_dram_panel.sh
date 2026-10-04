#!/usr/bin/env bash
# DRAM v1.4.6 on the 49 panel proteomes: annotate_genes + distill, default settings
# (KOfam + distillation forms; see the comparator install notes of the sulfur campaign).
source "$(conda info --base)/etc/profile.d/conda.sh"; set +u
conda activate DRAM14
HERE="$(cd "$(dirname "$0")" && pwd)"
rm -rf "$HERE/dram_out"
echo "[$(date '+%F %T')] DRAM annotate_genes"
DRAM.py annotate_genes -i "$HERE/dram_in/*.faa" -o "$HERE/dram_out" --threads 10 && echo "[OK] annotate" || echo "[FAIL] annotate"
DRAM.py distill -i "$HERE/dram_out/annotations.tsv" -o "$HERE/dram_out/distilled" && echo "[OK] distill" || echo "[FAIL] distill"
echo "[$(date '+%F %T')] DONE"
