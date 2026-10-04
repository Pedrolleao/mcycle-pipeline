#!/usr/bin/env bash
# Re-run of GTDB batch b01 after the queue: its first attempt never started (the log
# redirect pointed into a directory that did not exist yet).
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; set +u
conda activate METABOLIC_v4.0
HERE="$(cd "$(dirname "$0")" && pwd)"
G="$HERE/gtdb500_m"
until grep -q ALLDONE "$HERE/metabolic_queue.log" 2>/dev/null; do sleep 30; done
cd /home/dmin/Grants/Sulfur_Cycle/comparators/METABOLIC || exit 2
rm -rf "$G/metabolic_out/b01"; rm -f "$G/metabolic_in/b01/total.faa"
echo "[$(date '+%F %T')] START b01" >> "$HERE/metabolic_queue.log"
perl METABOLIC-G.pl -in "$G/metabolic_in/b01" -o "$G/metabolic_out/b01" -t 14 -kofam-db full > "$G/metabolic_out/b01.log" 2>&1
rc=$?
cp "$G"/metabolic_out/b01/KEGG_identifier_result/*.result.txt "$G/metabolic_out/kegg_all/" 2>/dev/null
cp "$G"/metabolic_out/b01/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv "$G/metabolic_out/worksheet1/b01.tsv" 2>/dev/null
echo "[$(date '+%F %T')] B01_DONE rc=$rc total=$(ls "$G"/metabolic_out/kegg_all/*.result.txt 2>/dev/null | wc -l) / 120" >> "$HERE/metabolic_queue.log"
