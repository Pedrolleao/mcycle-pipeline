#!/usr/bin/env bash
# METABOLIC v4.0 on the 49 panel proteomes, in sequential batches of 10 (one pristine
# input directory per batch; METABOLIC writes total.faa into its input and must never
# run twice at once). Collects the per-genome KO lists and worksheet 1 of every batch.
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; set +u
conda activate METABOLIC_v4.0
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO=/home/dmin/Grants/Sulfur_Cycle/comparators/METABOLIC
LOG="$HERE/metabolic_panel.log"; : > "$LOG"
mkdir -p "$HERE/metabolic_out/kegg_all" "$HERE/metabolic_out/worksheet1"
cd "$REPO" || exit 2
for bd in "$HERE"/metabolic_in/b??; do
  b=$(basename "$bd"); out="$HERE/metabolic_out/$b"
  rm -rf "$out"; rm -f "$bd"/total.faa
  echo "[$(date '+%F %T')] $b START ($(ls "$bd"/*.faa | wc -l) proteomes)" | tee -a "$LOG"
  perl METABOLIC-G.pl -in "$bd" -o "$out" -t 12 -kofam-db full > "$HERE/metabolic_out/$b.log" 2>&1
  rc=$?
  cp "$out"/KEGG_identifier_result/*.result.txt "$HERE/metabolic_out/kegg_all/" 2>/dev/null
  cp "$out"/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv \
     "$HERE/metabolic_out/worksheet1/$b.tsv" 2>/dev/null
  echo "[$(date '+%F %T')] $b DONE rc=$rc results=$(ls "$out"/KEGG_identifier_result/*.result.txt 2>/dev/null | wc -l)" | tee -a "$LOG"
done
echo "[$(date '+%F %T')] ALLDONE total=$(ls "$HERE"/metabolic_out/kegg_all/*.result.txt 2>/dev/null | wc -l) / 49" | tee -a "$LOG"
