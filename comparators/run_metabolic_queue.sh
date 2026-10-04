#!/usr/bin/env bash
# Queue for METABOLIC v4.0 after the panel run: (1) the 12 MAGs of the realism study,
# (2) the 120 methane-enriched GTDB genomes in batches of 10. Strictly sequential —
# METABOLIC must never run twice at once. Waits for the panel run to end first.
set -u
source "$(conda info --base)/etc/profile.d/conda.sh"; set +u
conda activate METABOLIC_v4.0
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO=/home/dmin/Grants/Sulfur_Cycle/comparators/METABOLIC
LOG="$HERE/metabolic_queue.log"; : > "$LOG"
until grep -q ALLDONE "$HERE/metabolic_panel.log" 2>/dev/null; do sleep 30; done
cd "$REPO" || exit 2
run_one() {   # input dir, output dir, kegg collection dir, worksheet copy
  rm -rf "$2"; rm -f "$1"/total.faa
  echo "[$(date '+%F %T')] START $1 ($(ls "$1"/*.faa | wc -l) proteomes)" | tee -a "$LOG"
  perl METABOLIC-G.pl -in "$1" -o "$2" -t 14 -kofam-db full > "$2.log" 2>&1
  rc=$?
  mkdir -p "$3" "$(dirname "$4")"
  cp "$2"/KEGG_identifier_result/*.result.txt "$3/" 2>/dev/null
  cp "$2"/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv "$4" 2>/dev/null
  echo "[$(date '+%F %T')] DONE $1 rc=$rc results=$(ls "$2"/KEGG_identifier_result/*.result.txt 2>/dev/null | wc -l)" | tee -a "$LOG"
}
mkdir -p "$HERE/mags_metabolic_out"
run_one "$HERE/mags_in" "$HERE/mags_metabolic_out/run" "$HERE/mags_metabolic_out/kegg_all" \
        "$HERE/mags_metabolic_out/worksheet1/mags.tsv"
echo "[$(date '+%F %T')] MAGS_DONE" | tee -a "$LOG"
G="$HERE/gtdb500_m"
for bd in "$G"/metabolic_in/b??; do
  b=$(basename "$bd")
  run_one "$bd" "$G/metabolic_out/$b" "$G/metabolic_out/kegg_all" "$G/metabolic_out/worksheet1/$b.tsv"
done
echo "[$(date '+%F %T')] ALLDONE gtdb=$(ls "$G"/metabolic_out/kegg_all/*.result.txt 2>/dev/null | wc -l) / 120" | tee -a "$LOG"
