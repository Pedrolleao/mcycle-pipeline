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

## 2026-10-05 — tool amendment after the hold-out was scored: truncated McrA (M10)
Decided by the user on 2026-10-04 after the GTDB-500 check (REPORT section 4) found
the defect. Applied only after M6 was closed and committed with the frozen tool.
- **Tool states.** Frozen: commit `9c39500` — every number of REPORT sections 1-5.
  Amended: commit `c7bf634`.
- **Change.** `gate_mcr_subunits` (`workflow/scripts/apply_rules.py`): a sub-threshold
  McrA homologue is evidence of an alkyl-CoM reductase — and mcrB / mcrG are then
  `disqualified` — only if its score per aligned profile position is below threshold /
  profile length (775.5 / 556 = 1.39 bits). A gene fragment of a canonical McrA scores
  above that (1.62-1.86 measured) and no longer triggers the rule. No fitted cut-off.
  The four genomes that exposed the defect are GTDB-500 genomes; no panel genome was
  used. Patch: `amendments/2026-10-04_truncated_mcra.patch`.
- **Effect** (`amendments/2026-10-05_truncated_mcra_changes.tsv`: all 48,057 calls of the
  smoke panel, the reference panel, the MAG study and the GTDB-500 set, before vs after):
  8 calls change — mcrB and mcrG `disqualified -> confirmed` in GCA_009780795,
  GCA_012799835, GCA_015662225, GCA_035391185.
  - Reference panel: 0 calls change. Hold-out micro-F1, full-panel micro-F1, trap
    precision, the gate (14 / 14) and the comparator benchmark are identical before and
    after (`mcycle_confusion.tsv` unchanged).
  - MAG study: 0 calls change. Smoke: 168 / 168.
  - GTDB-500: genomes outside alkane-oxidizer lineages with mcrB / mcrG labelled as
    alkyl-CoM reductase subunits: 5 -> 1. Alkane-oxidizer genomes so labelled: 3 -> 3.
    Agreement with raw KofamScan 0.9926 -> 0.9928 (`CONCORDANCE_amended.md`).
    `DIR_ACCURACY_amended.md`: unchanged, because the McrA itself is still not called.
- **Not changed, by decision — a separate cycle with its own validation** (`ROADMAP.md`):
  calling the truncated McrA itself, and non-euryarchaeal McrA that fall under the KO
  threshold at full length (the one genome still mislabelled, *Ca.* Methanomethylicus
  GCA_024464205).
- The post-amendment GTDB counts are in-sample for this rule: the genomes that showed
  the defect are the ones on which the repair is counted.
