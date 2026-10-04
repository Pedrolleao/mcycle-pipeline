#!/usr/bin/env bash
# DRAM v1.4.6 on the 12 MAGs of the realism study (two parallel batches), after the
# panel DRAM run has ended.
source "$(conda info --base)/etc/profile.d/conda.sh"; set +u
conda activate DRAM14
HERE="$(cd "$(dirname "$0")" && pwd)"
until grep -q "DONE [0-9] / 7 batches" "$HERE/dram_run.log" 2>/dev/null; do sleep 30; done
rm -rf "$HERE/mags_dram_out" "$HERE"/mags_in/d?; mkdir -p "$HERE/mags_dram_out"
i=0
for f in "$HERE"/mags_in/*.faa; do
  b="d$((i % 3 + 1))"; mkdir -p "$HERE/mags_in/$b"; ln -sf "$f" "$HERE/mags_in/$b/"; i=$((i+1))
done
for bd in "$HERE"/mags_in/d?; do
  b=$(basename "$bd")
  ( DRAM.py annotate_genes -i "$bd/*.faa" -o "$HERE/mags_dram_out/$b" --threads 8 > "$HERE/mags_dram_out/$b.log" 2>&1 \
    && DRAM.py distill -i "$HERE/mags_dram_out/$b/annotations.tsv" -o "$HERE/mags_dram_out/$b/distilled" >> "$HERE/mags_dram_out/$b.log" 2>&1 \
    && echo "[$(date '+%F %T')] $b OK" || echo "[$(date '+%F %T')] $b FAIL" ) &
done
wait
echo "[$(date '+%F %T')] MAGS_DRAM_DONE"
