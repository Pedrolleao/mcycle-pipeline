# Validation changelog

Dated record of every change to the reference panel, the ground truth and the
pre-registration after they were frozen. Newest last.

## 2026-10-04 — panel frozen (M1)
- `panel.tsv`: 49 genomes, 27 training / 22 hold-out. Rationale in `PANEL_PLAN.md`.
- Identity QC 49 / 49 pass (`panel_qc.tsv`).
- No ground truth existed and the pipeline had not been run on any added genome when
  the split was fixed.

## 2026-10-04 — ground truth built (M2)
- `panel.tsv`, column `kegg` only: `Mfumariolicum_SolV` tbd -> `mfh` (KEGG holds assembly
  GCA_949774925.1, the same sequence as the panel's GCA_000953475.1 after rotation and
  reverse complement — checked base by base); `Mluminyensis_B10` tbd -> `-`.
- `ground_truth.tsv` (KEGG, 37 genomes x 82 KO-anchored targets = 3,034 cells).
- `curated_cells.tsv`: 21 corrections (pmoA / pmoB / pmoC in six ammonia oxidizers and
  the hydrocarbon-monooxygenase carrier *M. chubuense* NBB4) and 90 literature cells for
  11 genomes KEGG does not hold. `phenotype_gt.tsv`: phenotype of all 49 genomes.
- `curated_function_gt.tsv`: 3,173 cells over 49 genomes.
- All of it was written before the pipeline was run on any genome added in M1.
- Looked at and deliberately NOT corrected, for lack of verified evidence:
  *M. burtonii* `mtmB` / `mtbB` (KEGG and UniProt list only `mttB`, `mtbC`, `mtbA`; the
  organism grows on methylamines, so the genes may be unannotated pyrrolysine genes);
  *M. album* BG8 `mmoX` (no sMMO is the textbook statement, but the genome paper does
  not say it — the cell is unscored); *N. maritimus* `mtrA` (KEGG assigns K00577).
