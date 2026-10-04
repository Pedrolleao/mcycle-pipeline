#!/usr/bin/env bash
# Last step of the GTDB-500 study, once METABOLIC has finished the 120 enriched genomes:
# merge its output with the backbone output reused from ncycle's study, then compute the
# four-tool concordance. Run from mcycle-pipeline/ inside the cycle-pipeline env.
set -euo pipefail
G=comparators/gtdb500_m
N=/home/dmin/Grants/Nitrogen_Cycle/ncycle-pipeline/comparators/gtdb500/metabolic_out
mkdir -p "$G/metabolic_all"
ln -sf "$PWD/$G"/metabolic_backbone/kegg_all/*.result.txt "$G/metabolic_all/"
ln -sf "$PWD/$G"/metabolic_out/kegg_all/*.result.txt "$G/metabolic_all/"
echo "KO lists: $(ls "$G"/metabolic_all/*.result.txt | wc -l) / 500"
python comparators/build_metabolic_tsv.py --kegg-dir "$G/metabolic_all" \
    --worksheets "$N/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv" \
                 "$G"/metabolic_out/worksheet1/*.tsv \
    --out "$G/metabolic.tsv"
MCYCLE_RESULTS=results_gtdb500 python validation/benchmark/concordance.py \
    --selection "$G/selection.tsv" --mcycdb "$G/mcycdb.tsv" --metabolic "$G/metabolic.tsv" \
    --out "$G/CONCORDANCE.md"
