# Validation report — mcycle-pipeline

Living report; sections are added as the campaign phases close (see `../WORKPLAN.md`).
All intervals are genome-cluster percentile bootstrap 95 % CIs (B = 10,000, seed 1234).

## 1. Reference-panel accuracy (M4, 2026-10-04)

Panel: 49 genomes, 27 training / 22 hold-out, frozen before any ground truth
(`PANEL_PLAN.md`). Ground truth: `curated_function_gt.tsv` — 3,173 cells = KEGG +
21 corrections + 90 literature cells + 49 `mcrA_anme` cells, written before the
pipeline was run on the panel. Tool state: commit `9c39500` (hardening closed on the
training genomes, `HARDENING.md`); the hold-out genomes were run once, after that
commit, and nothing was changed afterwards.

| metric | value [95 % CI] | cells | TP / FP / FN |
|---|---|---|---|
| **hold-out micro-F1** (22 genomes, de-leaked) | **0.961 [0.940, 0.976]** | 1,470 | 370 / 13 / 17 |
| full-panel micro-F1 | 0.973 [0.963, 0.981] | 3,173 | 851 / 21 / 26 |
| full-panel precision | 0.976 [0.965, 0.985] | | |
| training micro-F1 (in-sample) | 0.983 [0.974, 0.990] | 1,703 | 481 / 8 / 9 |
| trap precision, full panel | 1.000 [1.000, 1.000] | 330 | 71 / 0 / 0 |
| trap precision, independent cells only | 1.000 [1.000, 1.000] | 305 | 46 / 0 / 0 |
| trap precision, hold-out | 1.000 [1.000, 1.000] | 155 | 29 / 0 / 0 |

Seed-sourced hold-out cells removed by de-leaking: 0 (`seed_leakage.tsv`: 25
seed-sourced positive cells, all in training genomes).

Per pathway (full panel):

| pathway | F1 [95 % CI] | TP / FP / FN |
|---|---|---|
| mcr_core | 0.994 [0.983, 1.000] | 166 / 0 / 2 |
| co2_reduction | 0.984 [0.969, 0.995] | 180 / 2 / 4 |
| acetoclastic | 0.974 [0.949, 0.992] | 111 / 1 / 5 |
| methylotrophic_methanogenesis | 0.920 [0.857, 0.972] | 40 / 2 / 5 |
| energy_conservation | 0.952 [0.916, 0.977] | 118 / 9 / 3 |
| aerobic_methane_oxidation | 1.000 [1.000, 1.000] | 87 / 0 / 0 |
| formaldehyde_c1 | 0.955 [0.925, 0.977] | 149 / 7 / 7 |

Regression gate (`make regression-score`): 14 / 14 checks pass; floors and the
observed values they were set from are in the header of `test_regression.py`.

### What the trap result does and does not show
- The hold-out decoys were not called: AmoA of *Nitrosospira*, of the
  gammaproteobacterial *Ca.* Nitrosoglobus and of comammox *Nitrospira*, the
  hydrocarbon monooxygenase of *M. chubuense* NBB4, the alkyl-CoM reductases of
  *Ca.* Syntrophoarchaeum, the PQQ dehydrogenases of *Gluconobacter*. The hold-out
  positives were: PmoA of *Methylomicrobium*, *Methylocaldum*, *Methylacidimicrobium*;
  sMMO of *Methyloferula* and *Methylocaldum*; ANME-1 McrA of a second GTDB genus.
  McrA of *Methanococcoides* — the closest relative of ANME-3 — was not called ANME.
- The intervals are [1.000, 1.000] because there is no error to resample, not because
  the estimate is certain: with 29 hold-out trap positives and 126 negatives, the
  one-sided 95 % bounds are about 0.90 on recall and 0.98 on specificity (rule of three).
- Not tested by the hold-out (no genus available outside the seeds): ANME-2d, ANME-2a/2c/3,
  NC10, alpha-proteobacterial pMMO. They are tested in-sample, by leave-one-genus-out
  (`HARDENING.md`), and at scale in M6-M8.

### Against the KEGG ground truth (no corrections)
`GT_FILE=ground_truth.tsv`, 37 genomes, 3,034 cells: full-panel micro-F1 0.960
[0.943, 0.973], precision 0.977 [0.967, 0.986]; hold-out micro-F1 0.944 [0.905, 0.967].
Trap recall drops to 0.724 by construction — KEGG assigns the pmo KOs to ammonia
oxidizers, the tool does not call them. All 21 corrections are `present -> absent`
and favour a tool that separates pmo from amo, so both numbers are reported.

### Against raw KofamScan thresholds (`compare_kofam.py`)
Same hmmscan tables, stock KOfam thresholds, no rules, no clade models:

| scope | mcycle P / R / F1 | raw KOfam P / R / F1 |
|---|---|---|
| ALL targets | 0.976 / 0.970 / 0.973 | 0.922 / 0.952 / 0.937 |
| trap targets | 1.000 / 1.000 / 1.000 | 0.664 / 1.000 / 0.798 |
| hold-out, ALL | 0.966 / 0.956 / 0.961 | 0.908 / 0.946 / 0.927 |
| hold-out, trap | 1.000 / 1.000 / 1.000 | 0.617 / 1.000 / 0.763 |

The trap false positives of raw KOfam are pmoA / pmoB / pmoC in ammonia and
hydrocarbon oxidizers (22) and `mcrA_anme` in methanogens (14; a KO cannot call the
direction). This is a descriptive table; the tested comparison is M5.

### Disagreements left (47 cells; full list: `mcycle_confusion.tsv`)
Recurring classes: `acs` (5 FN, 1 FP — K01895 against its paralogues), `fdh` / `fdhA`
(selenocysteine splits and K22516 / K00123 adjudication, 10 cells), `frhB` (5 FP —
F420-binding subunits of other complexes), `hdrD` (4 FP), `sgaA` / `hprA` (broad
KOs, 6 cells), corrinoid proteins `mttC` / `mtbC` (5 FN), `mvhA` / `mvhD` (3 FN).

Hold-out cells where the ground truth, not the tool, may be wrong — NOT corrected,
listed for a later evidence-based review (`CHANGELOG.md`):
- `mtmB`, `mtbB` in *Methanococcoides burtonii* (called by the tool; KEGG lists only
  `mttB`; the organism grows on mono- and dimethylamine; flagged in the changelog
  before the hold-out was run);
- `ftr` in *Methanonatronarchaeum thermophilum* (called; Sorokin 2017 reports Ftr lost).

## 2. Comparator benchmark (M5, 2026-10-04)

Pre-registered (`benchmark/prereg.md`, commit `3cc200c`) and run as registered; full
write-up in `benchmark/COMPARISON_REPORT.md`, figure `benchmark/figures/curated_panel_accuracy.*`.

| comparator | trap precision | difference [95 % CI] | ALL micro-F1 | difference [95 % CI] |
|---|---|---|---|---|
| mcycle | 1.000 | — | 0.973 | — |
| raw KofamScan | 0.664 | +0.336 [0.210, 0.488] | 0.937 | +0.037 [0.028, 0.047] |
| METABOLIC v4.0 | 0.689 | +0.311 [0.202, 0.446] | 0.936 | +0.037 [0.028, 0.049] |
| DRAM v1.4.6 | 0.527 | +0.473 [0.333, 0.616] | 0.350 | +0.623 [0.561, 0.675] |
| MCycDB 2021 | 0.626 | +0.374 [0.250, 0.526] | 0.795 | +0.178 [0.143, 0.219] |

All eight contrasts: bootstrap p < 0.0001, BH q < 0.0001. The trap differences hold
with `mcrA_anme` removed from the trap set (+0.190 to +0.378, every interval excludes
0) and on the 22 hold-out genomes alone (+0.356 to +0.524). METABOLIC and MCycDB,
which model pmoA separately from amoA, still report pmoA in the two gammaproteobacterial
ammonia oxidizers and in the hydrocarbon-monooxygenase carriers; METABOLIC, DRAM and
MCycDB report the alkyl-CoM reductases of the alkane oxidizers as McrA.
