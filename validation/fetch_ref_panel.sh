#!/usr/bin/env bash
# Populate ../ref_panel/ from validation/panel.tsv: one genome FASTA per roster row,
# downloaded from NCBI by assembly accession (datasets API). A genome the smoke panel
# already holds under the same name is copied from ../test_panel/ instead (same
# accession in both rosters). Network access is needed for the downloads.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here/../../ref_panel}"
smoke="$here/../../test_panel"
mkdir -p "$out"
API=https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession
fail=0
while IFS=$'\t' read -r name source input _; do
  f="$out/$name.$input"
  [ -s "$f" ] && { echo "have  $name.$input"; continue; }
  if [ -s "$smoke/$name.$input" ] && grep -q "^$name	$source	" "$smoke/panel.tsv"; then
    cp "$smoke/$name.$input" "$f"; echo "copy  $name.$input  (smoke panel)"; continue
  fi
  kind=$([ "$input" = faa ] && echo PROT_FASTA || echo GENOME_FASTA)
  tmp=$(mktemp -d)
  if curl -fsSL --retry 3 -o "$tmp/d.zip" "$API/$source/download?include_annotation_type=$kind" \
     && unzip -q -o "$tmp/d.zip" -d "$tmp"; then
    if [ "$input" = faa ]; then g=$(find "$tmp" -name protein.faa | head -1)
    else g=$(find "$tmp" -name '*_genomic.fna' | head -1); fi
    if [ -n "$g" ] && [ -s "$g" ]; then cp "$g" "$f"; echo "fetch $name.$input  ($source)"
    else echo "FAIL  $name ($source): no $input in the NCBI package" >&2; fail=1; fi
  else echo "FAIL  $name ($source): download failed" >&2; fail=1; fi
  rm -rf "$tmp"; sleep 0.5
done < <(grep -v '^#' "$here/panel.tsv" | tail -n +2)
exit $fail
