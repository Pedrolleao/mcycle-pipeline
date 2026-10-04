#!/usr/bin/env bash
# DRAM v1.4.6 on the 49 panel proteomes: annotate_genes + distill, default settings
# (KOfam + distillation forms; see the comparator install notes of the sulfur campaign).
# DRAM annotates one genome at a time (~15 min each with the full KOfam set), so the
# panel is split into 7 batches that run side by side; each batch is annotated and
# distilled on its own and the product tables are concatenated by build_dram_tsv.py.
source "$(conda info --base)/etc/profile.d/conda.sh"; set +u
conda activate DRAM14
HERE="$(cd "$(dirname "$0")" && pwd)"
rm -rf "$HERE/dram_out" "$HERE"/dram_in/b??; mkdir -p "$HERE/dram_out"
i=0
for f in "$HERE"/dram_in/*.faa; do
  b=$(printf "b%02d" $((i % 7 + 1))); mkdir -p "$HERE/dram_in/$b"; ln -sf "$f" "$HERE/dram_in/$b/"; i=$((i+1))
done
echo "[$(date '+%F %T')] DRAM annotate_genes, 7 batches"
for bd in "$HERE"/dram_in/b??; do
  b=$(basename "$bd")
  ( DRAM.py annotate_genes -i "$bd/*.faa" -o "$HERE/dram_out/$b" --threads 4 > "$HERE/dram_out/$b.log" 2>&1 \
    && DRAM.py distill -i "$HERE/dram_out/$b/annotations.tsv" -o "$HERE/dram_out/$b/distilled" >> "$HERE/dram_out/$b.log" 2>&1 \
    && echo "[$(date '+%F %T')] $b OK" || echo "[$(date '+%F %T')] $b FAIL" ) &
done
wait
echo "[$(date '+%F %T')] DONE $(ls "$HERE"/dram_out/b??/distilled/product.tsv 2>/dev/null | wc -l) / 7 batches"
