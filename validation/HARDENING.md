# Hardening log (M3) — training genomes only

Everything below was decided on the 27 TRAINING genomes of `validation/panel.tsv`
(`MCYCLE_SCOPE=train`). The 22 hold-out genomes were not run through the pipeline
until this log was closed (2026-10-04). Tools: `validation/score_mcycle.py`,
`validation/ko_margins.py` (score margins per target), `validation/explain_cell.py`.

## Training-panel score, curated-function ground truth

| stage | cells | TP | FP | FN | micro-F1 [95 % CI] | trap P | trap R |
|---|---|---|---|---|---|---|---|
| baseline (pipeline as built) | 1,703 | 472 | 17 | 18 | 0.964 [0.951, 0.977] | 1.000 | 0.976 |
| + rules and thresholds | 1,703 | 478 | 8 | 12 | 0.980 | 1.000 | 0.976 |
| + clade HMMs (final) | 1,703 | 481 | 8 | 9 | 0.983 [0.974, 0.991] | 1.000 | 1.000 |

These are in-sample numbers: they say the changes do what they were meant to do, not
how well the tool generalizes. That is the hold-out's job (M4).

## Changes

### 1. Subunits follow the subunit that defines the enzyme (`apply_rules.py`)
- **mcrB / mcrG vs alkyl-CoM reductase** (`gate_mcr_subunits`). *Ca.* Ethanoperedens:
  EcrA 642 on K00399 (threshold 775.5, not called), EcrB 586 on K00401 (threshold
  504.3, called), EcrG 288 on K00402 (360.9). mcrB / mcrG are disqualified when the
  genome has no McrA call and does carry an McrA homologue at >= 0.5 of the threshold.
  A genome with no McrA homologue at all keeps mcrB / mcrG (fragmented MAGs).
- **mmoY / Z / B / C / D vs other di-iron monooxygenases** (`gate_mmo_subunits`). A
  lone reductase of *Methylocystis* sp. SC2 scores 410 on K16161 (threshold 407.7;
  true MmoC 514-543). Without mmoX the other sMMO components are disqualified.
- **fdhA vs fdh** (`resolve_fdh`). Archaeal F420-dependent FdhA score 0.72-0.95 of
  the K22516 threshold and pass the generic K00123; bacterial FdhA score 0.76-0.86 on
  K22516. With fdhB (K00125) in the genome, an FdhA-family protein at >= 0.70 of the
  K22516 threshold is fdhA, and is not also `fdh`. Fixed fdhA in *M. kandleri* and
  *M. thermautotrophicus* and `fdh` in *M. maripaludis*, *M. thermautotrophicus*.
- Tried and REVERTED: frhB following frhA / frhG — removed two false positives
  (*A. fulgidus*, *Methylocystis*) and created two false negatives (*Methanothrix*
  x 2, where KEGG assigns K00441 without frhA). No net gain, so not kept.

### 2. Query-coverage floor
- **mvhD** `hmm_min_qcov: 0.5`. HdrA-MvhD fusion proteins of 495-812 residues pass
  K14127 (204-225, threshold 180.1) in *Methanosarcina*, *Methanothrix* x 2 and
  *D. vulgaris*; true MvhD are 135-141 residues. Four false positives removed.

### 3. KO thresholds (`ko_tc`, each with `tc_rationale` in `targets.yaml`)
Only where the gap between the lowest true and the highest false hit is wide:

| target | KO | KOfam | new | lowest true hit | highest false hit |
|---|---|---|---|---|---|
| hmd | K13942 | 575.0 | 540 | 566 (*M. kandleri*) | 72 |
| fwdD | K00203 | 122.9 | 110 | 119 (*M. oxyfera*) | 39 |
| hdrB | K03389 | 333.6 | 260 | 276 (*A. fulgidus*) | 90 |
| hdrC | K03390 | 125.3 | 100 | 108 (*A. fulgidus*) | 64 |
| fdh | K00122 | 238.4 | 400 | 640 (*M. silvestris*) | 240 (*M. kandleri*) |

`ko_tc` now also takes a mapping `{KO: threshold}` (`build_hmm_db.py`).
Left alone because the gap is narrow or the labels overlap: hdrD (false 381 vs true
433), frhB (297 vs 332), hprA (336 vs 433), sgaA, mvhA, acs.

### 4. Clade HMMs replace the two identity gates
New: `workflow/scripts/harvest_clade_refs.py` (typed sequences from GTDB-classified
genomes; hold-out genera excluded by name, NCBI and GTDB) and
`workflow/scripts/build_clade_hmms.py` (one model per clade, leave-one-genus-out,
threshold halfway between the best negative and the worst left-out positive).
Models, training sets, calibration tables and manifests are under `targets/pmoA/`
and `targets/mcrA_anme/`. `clade_hmm_required: true` in `targets.yaml`: a protein
passing a clade model is `confirmed`; a family (KO) hit without one is `disqualified`.

| model | training seqs | genera | best negative | worst left-out positive | threshold | left-out passing |
|---|---|---|---|---|---|---|
| pmoA__gamma | 46 | 43 | 402.5 | 444.5 | 423.5 | 54 / 54 |
| pmoA__alpha | 18 | 3 | 370.2 | 396.0 | 383.1 | 41 / 41 |
| pmoA__verruco | 8 | 2 | 312.5 | 347.5 | 330.0 | 1 / 1 |
| pmoA__nc10 | 4 | 1 | 282.3 | 549.5 | 415.9 | 4 / 4 (species) |
| mcrA_anme__anme1 | 13 | 11 | 794.8 | 938.9 | 866.8 | 14 / 14 |
| mcrA_anme__anme2ab | 6 | 4 | 982.2 | 1235.1 | 1108.7 | 9 / 9 |
| mcrA_anme__anme2c | 4 | 3 | 983.0 | 1046.0 | 1014.5 | 2 / 2 |
| mcrA_anme__anme2d | 9 | 3 | 978.4 | 967.6 | 998.0 | 10 / 11 |
| mcrA_anme__anme3 | 5 | 1 | 1227.2 | 1201.8 | 1251.7 | 3 / 5 (species) |

Negatives: pmoA — AmoA of beta-AOB (13), gamma-AOB (6) and AOA (20), and the butane
monooxygenase of *Nocardioides* sp. CF8. mcrA_anme — McrA of methanogens from 77
genera and three alkyl-CoM reductases.

Limits found here, to be carried into the report:
- **ANME-3 is not cleanly separable** from methylotrophic Methanosarcinaceae
  (*Methanomethylovorans*, *Methanolobus* score 1,202-1,227 on the ANME-3 model). The
  threshold is precision-first; two of five left-out ANME-3 species fall under it.
- **ANME-2d**: one lineage (`g__CAZCRM01`, a 229-scaffold MAG whose McrA also fails
  K00399) falls under the threshold.
- 24 divergent PmoA paralogues (pxmA of gammaproteobacteria and *Methylocystis*,
  pmoA3 of Verrucomicrobia; < 60 % identical to every curated seed) are in neither
  the positives nor the negatives — their function is unresolved and they are not
  called.
- No comammox and no mycobacterial sequence is in the negatives: both are hold-out
  genera. The hold-out is their first test.
- pmoA__verruco and pmoA__nc10 rest on two genera and one genus.

### 5. Decoy margins measured (no change needed)
| target | decoy (training) | decoy score | threshold | lowest true hit |
|---|---|---|---|---|
| mmoX K16157 | sBMO BmoX, *T. butanivorans* | 884 | 1168.7 | 1178 |
| xoxF K23995 | PedH, *P. putida* | 548 | 859.1 | — |
| mxaF K14028 | PedE / PedH, *P. putida* | 397 | 1059.1 | — |
| mcrA K00399 | EcrA, *Ca.* Ethanoperedens | 642 | 775.5 | 834 (ANME-1) |
| nod K27148 | qNor-type, *M. oxyfera* | 609 | 1155.0 | 1351 |
| pmoA K10944 | pBmoA, *Nocardioides* CF8 | 338 | 353.0 | (now decided by the clade HMMs) |

The mmoX margin is thin on the true side (1178-1201 against 1168.7).

### 6. Suspected KO-tier noise, examined
`mtrA` in *N. maritimus*: KEGG assigns K00577 to the same genome — an agreement, not
noise. `mtbC` in *M. maripaludis*: likewise in agreement with KEGG.

## What is left wrong on the training genomes (17 cells)
- fdhA / fdh in *Methanothrix* x 2 (KEGG assigns K22516 without an FdhB; the protein
  scores higher on K00123) — 4 cells.
- Selenocysteine formate dehydrogenases of *E. coli* and *D. vulgaris*: the halves do
  not add up to the K00123 threshold — 2 cells.
- Near-threshold or overlapping KOs: sgaA x 3, hprA, hdrD, frhB x 2, mvhA, mvhD
  (*M. kandleri*), acs x 2.
