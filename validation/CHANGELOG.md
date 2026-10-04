# Validation changelog

Dated record of every change to the reference panel, the ground truth and the
pre-registration after they were frozen. Newest last.

## 2026-10-04 — panel frozen (M1)
- `panel.tsv`: 49 genomes, 27 training / 22 hold-out. Rationale in `PANEL_PLAN.md`.
- Identity QC 49 / 49 pass (`panel_qc.tsv`).
- No ground truth existed and the pipeline had not been run on any added genome when
  the split was fixed.
