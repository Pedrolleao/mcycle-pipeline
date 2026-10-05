# GTDB-500 cross-tool study — methane cycle

**Status: complete (2026-10-05).** 500 GTDB r232 species representatives, four tools:
mcycle (frozen at commit `9c39500`), raw KofamScan, METABOLIC v4.0, MCycDB 2021. DRAM
is left out at this scale, as in the nitrogen and sulfur studies. Arbitrary GTDB
genomes have no curated truth, so this is a **concordance** study — where the tools
agree and where they split — not an accuracy estimate; accuracy is the 49-genome
panel (`../../validation/REPORT.md`, sections 1-2).

## 1. Genomes
- **380 backbone** genomes: the cross-phylum set of the nitrogen study, reused
  verbatim, so that the three sister studies share it.
- **120 methane-enriched** genomes: 5 random species representatives (seed 1234) from
  each of 24 clades — methanogens of nine lineages, ANME-1 / -2ab / -2c / -2d / -3,
  alkane oxidizers with alkyl-CoM reductase, five aerobic methanotroph lineages, and
  four ammonia-oxidizer lineages (`select_gtdb_mcyc.py`, `selection.tsv`).
- Proteins: Prodigal calls, the same for every tool. Backbone proteomes and their
  METABOLIC output are those of the nitrogen study; the enriched genomes were fetched
  from NCBI, called with Prodigal (`-p meta`) and run through METABOLIC here, in 12
  batches of 10 (`prepare_proteomes.py`, `../run_metabolic_queue.sh`).
- 32 of the 120 enriched genomes are also among the 320 genomes harvested for the clade
  HMMs (25 of them supplied a training or calibration sequence); `DIR_ACCURACY.tsv` says,
  per genome, whether it was a training genome, of a training genus, or of a new genus.

## 2. Wall-clock (32-core machine)
| stage | 500 genomes |
|---|---|
| mcycle (hmmscan + rules + reports), 8 cores | ~38 min |
| MCycDB (DIAMOND), 6 threads | ~50 min |
| METABOLIC, 120 enriched genomes, 14 threads | 12 batches x ~17 min = ~3.4 h |

## 3. Agreement (`CONCORDANCE.md`, section a)
| pair | cells | agreement | Cohen kappa |
|---|---|---|---|
| mcycle vs raw KofamScan | 41,000 | 0.993 | 0.964 |
| mcycle vs METABOLIC | 41,000 | 0.992 | 0.962 |
| raw KofamScan vs METABOLIC | 41,000 | 0.994 | 0.969 |
| mcycle vs MCycDB | 39,500 | 0.875 | 0.568 |
| raw KofamScan vs MCycDB | 39,500 | 0.871 | 0.553 |
| METABOLIC vs MCycDB | 39,500 | 0.874 | 0.567 |

The three KOfam-based tools agree on more than 99 % of cells — they share profiles —
and differ where mcycle's rules and clade models act. MCycDB reports almost twice as
many present calls (9,028 against 4,791-4,961): its best-hit search assigns broad
families (frmA in 476 genomes, pta in 441, mtkB in 341, fdhA in 325) that the
threshold-based tools report in a few dozen.

## 4. Trap 1 — pmoA in ammonia oxidizers
Genomes in which each tool reports pmoA:

| clade (5 genomes each) | mcycle | raw KofamScan | METABOLIC | MCycDB |
|---|---|---|---|---|
| beta-AOB | 0 | 4 | 0 | 0 |
| gamma-AOB | 0 | 5 | 5 | 5 |
| AOA | 0 | 1 | 0 | 0 |
| comammox *Nitrospira* | 0 | 3 | 0 | 0 |
| **20 ammonia-oxidizer genomes** | **0** | **13** | **5** | **5** |
| Methylomonadaceae | 5 | 5 | 5 | 5 |
| Methylococcaceae | 2 | 3 | 3 | 3 |
| alpha methanotrophs | 4 | 5 | 5 | 5 |
| Verrucomicrobia | 5 | 5 | 5 | 5 |
| NC10 | 4 | 4 | 4 | 4 |

The panel result repeats at scale. A KO alone calls most ammonia oxidizers
methanotrophs; METABOLIC and MCycDB, which model pmoA apart from amoA, clear the
beta-AOB, the AOA and comammox but report pmoA in **every** gammaproteobacterial
ammonia oxidizer; mcycle in none. In the methanotroph clades mcycle is the more
conservative tool by 2 genomes (one *Methylumidiphilus*, one *Methylocapsa* MAG): their
only copper-monooxygenase subunit A is a pxmA-type paralogue, which mcycle recognizes
as outside the pmoA clades and does not call (`DIR_ACCURACY.md`).

## 5. Trap 2 — McrA and its direction
- mcycle calls the Mcr `reverse` in 14 genomes — all in ANME families — and
  `methanogenic` in 42 (34 of the enriched clades, 8 backbone); no genome of
  an ANME family is called methanogenic and no methanogen reverse. No comparator has a
  direction call: the same 14 genomes are plain McrA carriers to them.
- In the alkane-oxidizer clade METABOLIC and MCycDB report mcrA in 5 of 5 genomes,
  mcycle and raw KofamScan in 2 — the two *Ca.* Methanoliparum, which carry a canonical
  Mcr beside their alkyl-CoM reductase.
- McrA that is not called by mcycle although the lineage expects one: 11 genomes
  (6 without the gene in the assembly, 4 with a truncated gene, 1 full-length
  non-euryarchaeal McrA under the threshold); `DIR_ACCURACY.md`. MCycDB and METABOLIC
  report mcrA in a few more of these fragmented genomes (e.g. Methanofastidiosales
  4 and 3 of 5, against 2).

## 6. Orthogonal checks of the clade calls
`DIR_ACCURACY.{md,tsv}` (against the genome's GTDB lineage) and `DIR_PLACEMENT.md`
(gene trees): direction 53 / 53 where an McrA is called; pmoA 34 / 36, 16 / 16 in
genera absent from the training sets. Details in `../../validation/REPORT.md`, section 4.

## 7. Reproduce
```bash
python comparators/gtdb500_m/select_gtdb_mcyc.py
python comparators/gtdb500_m/prepare_proteomes.py
MCYCLE_CONFIG=config/config_gtdb500.yaml python mcycle.py --input comparators/gtdb500_m/proteomes --cores 8 --skip-db-setup
python comparators/build_mcycdb_tsv.py --panel comparators/gtdb500_m/proteomes \
    --hits comparators/gtdb500_m/mcycdb_out --out comparators/gtdb500_m/mcycdb.tsv
bash comparators/run_metabolic_queue.sh        # then run_metabolic_gtdb_b01.sh
bash comparators/gtdb500_m/finalize.sh         # metabolic.tsv + CONCORDANCE.{md,tsv}
MCYCLE_RESULTS=results_gtdb500 python validation/phylogeny/clade_accuracy.py \
    --selection comparators/gtdb500_m/selection.tsv --out comparators/gtdb500_m/DIR_ACCURACY.md
MCYCLE_RESULTS=results_gtdb500 python validation/phylogeny/place_clades.py --trap mcr \
    --selection comparators/gtdb500_m/selection.tsv --outdir comparators/gtdb500_m
```
