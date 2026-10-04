# Comparator benchmark — mcycle-pipeline vs raw KofamScan, METABOLIC, DRAM, MCycDB

Run 2026-10-04 exactly as pre-registered (`prereg.md`, commit `3cc200c`, written
before METABOLIC, DRAM or MCycDB touched the panel). 49 genomes, the same Prodigal
proteins for every tool, ground truth `curated_function_gt.tsv` (3,173 cells),
genome-cluster bootstrap B = 10,000, seed 1234. mcycle is frozen at commit `9c39500`.
Full output: `benchmark_final.txt`; every number: `benchmark_results.tsv`.

## Confirmatory family (8 tests, Benjamini-Hochberg)

| comparator | endpoint | mcycle | comparator | difference [95 % CI] | bootstrap p | BH q | verdict |
|---|---|---|---|---|---|---|---|
| raw KofamScan | trap precision | 1.000 | 0.664 | +0.336 [0.210, 0.488] | < 0.0001 | < 0.0001 | mcycle better |
| METABOLIC v4.0 | trap precision | 1.000 | 0.689 | +0.311 [0.202, 0.446] | < 0.0001 | < 0.0001 | mcycle better |
| DRAM v1.4.6 | trap precision | 1.000 | 0.527 | +0.473 [0.333, 0.616] | < 0.0001 | < 0.0001 | mcycle better |
| MCycDB 2021 | trap precision | 1.000 | 0.626 | +0.374 [0.250, 0.526] | < 0.0001 | < 0.0001 | mcycle better |
| raw KofamScan | ALL micro-F1 | 0.973 | 0.937 | +0.037 [0.028, 0.047] | < 0.0001 | < 0.0001 | mcycle better |
| METABOLIC v4.0 | ALL micro-F1 | 0.973 | 0.936 | +0.037 [0.028, 0.049] | < 0.0001 | < 0.0001 | mcycle better |
| DRAM v1.4.6 | ALL micro-F1 | 0.973 | 0.350 | +0.623 [0.561, 0.675] | < 0.0001 | < 0.0001 | mcycle better |
| MCycDB 2021 | ALL micro-F1 | 0.973 | 0.795 | +0.178 [0.143, 0.219] | < 0.0001 | < 0.0001 | mcycle better |

"< 0.0001" = none of the 10,000 bootstrap differences was <= 0. All eight contrasts
meet the decision rule, including its second half:

**Trap precision without `mcrA_anme`** (the target no comparator can express;
pre-specified sensitivity analysis, 281 cells): KofamScan 0.756, difference +0.244
[0.100, 0.415]; METABOLIC 0.810, +0.190 [0.081, 0.333]; DRAM 0.622, +0.378 [0.200,
0.561]; MCycDB 0.719, +0.281 [0.157, 0.439]. Every interval excludes 0, so the trap
result does not rest on the direction call.

## Per tool (subunit resolution)

| tool | ALL precision | ALL recall | ALL F1 | trap precision | trap recall |
|---|---|---|---|---|---|
| mcycle | 0.976 [0.965, 0.985] | 0.970 [0.956, 0.982] | 0.973 [0.963, 0.981] | 1.000 | 1.000 |
| raw KofamScan | 0.922 [0.897, 0.941] | 0.952 [0.934, 0.968] | 0.937 [0.920, 0.950] | 0.664 [0.513, 0.792] | 1.000 |
| METABOLIC | 0.921 [0.893, 0.941] | 0.952 [0.934, 0.967] | 0.936 [0.918, 0.950] | 0.689 [0.557, 0.797] | 1.000 |
| DRAM | 0.765 [0.673, 0.849] | 0.227 [0.186, 0.279] | 0.350 [0.298, 0.411] | 0.527 [0.380, 0.667] | 0.690 [0.603, 0.782] |
| MCycDB | 0.710 [0.644, 0.763] | 0.904 [0.885, 0.921] | 0.795 [0.749, 0.831] | 0.626 [0.478, 0.750] | 0.944 [0.877, 0.988] |

DRAM's recall is low by construction: its distillate has 11 methane functions, and
the 63 targets without one are scored absent (pre-registered, as in the sister
benchmarks). On the 20 targets DRAM does cover, the F1 difference is +0.112
[0.057, 0.179]. On covered targets only, the other differences are: KofamScan +0.029
[0.021, 0.041], METABOLIC +0.029 [0.020, 0.040], MCycDB +0.151 [0.116, 0.197].

## Where the trap errors are (false positives / false negatives per target)

| tool | mcrA | mcrA_anme | pmoA | pmoB | pmoC | mmoX | mxaF | xoxF |
|---|---|---|---|---|---|---|---|---|
| mcycle | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| raw KofamScan | 0 / 0 | 14 / 0 | 6 / 0 | 8 / 0 | 8 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| METABOLIC | 2 / 0 | 16 / 0 | 4 / 0 | 2 / 0 | 4 / 0 | 0 / 0 | 4 / 0 | 0 / 0 |
| DRAM | 2 / 0 | 16 / 0 | 8 / 0 | 8 / 0 | 8 / 0 | 2 / 0 | 0 / 10 | 0 / 12 |
| MCycDB | 2 / 1 | 15 / 0 | 3 / 0 | 4 / 0 | 3 / 0 | 1 / 0 | 8 / 0 | 4 / 3 |

- **pmoA in non-methanotrophs.** Raw KOfam calls it in all six ammonia oxidizers.
  METABOLIC, which ships a dedicated pmoA model, and MCycDB, which has separate pmoA
  and amoA families, both clear the beta-AOB, the AOA and comammox — and both still
  call pmoA in the two gammaproteobacterial ammonia oxidizers (*Nitrosococcus oceani*,
  *Ca.* Nitrosoglobus terrae) and in the hydrocarbon-monooxygenase carriers
  (*M. chubuense* NBB4; METABOLIC also *Nocardioides* sp. CF8). Those are the cases
  the clade models were built for.
- **mcrA in alkane oxidizers.** METABOLIC, DRAM and MCycDB report the alkyl-CoM
  reductases of *Ca.* Ethanoperedens and *Ca.* Syntrophoarchaeum as McrA.
- **Direction.** Scored as pre-registered (a tool without a direction call gets its
  mcrA call), every comparator is wrong on the 14-16 methanogens. This is a statement
  about what the tools report, not an error of theirs on a question they answer.
- mcycle's own errors are all outside the trap set (21 FP, 26 FN; REPORT section 1).

## Hold-out only (22 genomes; exploratory, not in the family)

| comparator | ALL F1 difference | trap precision difference |
|---|---|---|
| raw KofamScan | +0.034 [0.023, 0.053] | +0.383 [0.193, 0.619] |
| METABOLIC | +0.033 [0.021, 0.054] | +0.356 [0.191, 0.571] |
| DRAM | +0.620 [0.518, 0.693] | +0.524 [0.300, 0.732] |
| MCycDB | +0.168 [0.119, 0.237] | +0.413 [0.221, 0.659] |

The differences on genomes mcycle was never tuned on are the same size as on the
whole panel.

## Step resolution (14 steps, 554 cells; exploratory)

Micro-F1: mcycle 0.994 [0.983, 1.000]; KofamScan 0.933; METABOLIC 0.922; MCycDB 0.716;
DRAM 0.532.

## Exploratory: TN-inclusive metrics (own BH family)

MCC differences, ALL / trap: KofamScan +0.051 / +0.244; METABOLIC +0.052 / +0.223;
DRAM +0.636 / +0.525; MCycDB +0.248 / +0.307; all intervals exclude 0.

## How to read this

- The overall-F1 gap to KofamScan and METABOLIC is small (0.037): mcycle is built on
  the same KOfam profiles, and most of the cycle is not a trap. The gain is where it
  was designed to be — the copper monooxygenases, the direction of Mcr, the subunits
  of look-alike enzymes.
- mcycle's trap precision has no error on this panel, so its interval is degenerate;
  the differences are carried by the comparators' errors (REPORT section 1 gives the
  rule-of-three bounds).
- The raw-KofamScan contrast was not blind (it had been tabulated before the
  pre-registration); the other three were.
- 27 of the 49 genomes informed mcycle's thresholds and models; the hold-out table
  above is the comparison without that advantage.
- Against the KEGG ground truth the pmo corrections disappear and every KO-based tool
  is, by construction, right on those cells; that contrast is in REPORT section 1.

## Reproduce
```bash
bash comparators/run_metabolic_panel.sh && python comparators/build_metabolic_tsv.py
bash comparators/run_dram_panel.sh      && python comparators/build_dram_tsv.py
python comparators/build_mcycdb_tsv.py
cd validation/benchmark && python benchmark_stats.py --metabolic metabolic.tsv \
    --dram dram.tsv --mcycdb mcycdb.tsv > benchmark_final.txt
```
