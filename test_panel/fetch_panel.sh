#!/usr/bin/env bash
# Populate test_panel/ from panel.tsv: download each genome from NCBI (datasets API;
# genome FASTA or protein FASTA according to the `input` column), or copy a local file.
set -euo pipefail
cd "$(dirname "$0")"
API=https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession
grep -v '^#' panel.tsv | tail -n +2 | while IFS=$'\t' read -r name source input _; do
  out="$name.$input"
  [ -s "$out" ] && { echo "have  $out"; continue; }
  if [[ "$source" == GC[AF]_* ]]; then
    kind=$([ "$input" = faa ] && echo PROT_FASTA || echo GENOME_FASTA)
    tmp=$(mktemp -d)
    curl -fsSL --retry 3 -o "$tmp/d.zip" "$API/$source/download?include_annotation_type=$kind"
    unzip -q -o "$tmp/d.zip" -d "$tmp"
    if [ "$input" = faa ]; then f=$(find "$tmp" -name protein.faa | head -1)
    else f=$(find "$tmp" -name '*_genomic.fna' | head -1); fi
    [ -n "$f" ] && [ -s "$f" ] || { echo "FAIL  $name ($source): no $input in the NCBI package" >&2; rm -rf "$tmp"; continue; }
    cp "$f" "$out"; rm -rf "$tmp"; echo "fetch $out  ($source)"
  else
    [ -s "$source" ] || { echo "FAIL  $name: $source not found" >&2; continue; }
    cp "$source" "$out"; echo "copy  $out"
  fi
done
