# Clade calls against the genome taxonomy (GTDB-500)

Reference = GTDB r232 lineage of the genome (marker-gene phylogeny), independent of the McrA / PmoA sequence the tool reads. Agreement with 95 % Wilson intervals. Only genomes carrying a member of the protein family are listed.


## Mcr direction (reverse / methanogenic / no call)

Genomes with a family member: 70; with a lineage expectation: 67; lineage of unsettled physiology: 3.

| stratum | agree / n | agreement [95 % CI] |
|---|---|---|
| **all** | 56 / 67 | 0.836 [0.729, 0.906] |
| new genus | 42 / 49 | 0.857 [0.733, 0.929] |
| training genus | 6 / 10 | 0.600 [0.313, 0.832] |
| training genome | 8 / 8 | 1.000 [0.676, 1.000] |
| named genus | 47 / 57 | 0.825 [0.706, 0.902] |
| placeholder genus | 9 / 10 | 0.900 [0.596, 0.982] |
| new genus, named | 35 / 42 | 0.833 [0.694, 0.917] |
| new genus, placeholder | 7 / 7 | 1.000 [0.646, 1.000] |

Expected (rows) x mcycle (columns):

| expected | methanogenic | no_call | reverse |
|---|---|---|---|
| methanogenic | 39 | 7 | 0 |
| no_call | 0 | 3 | 0 |
| reverse | 0 | 4 | 14 |
| unknown | 3 | 0 | 0 |

**Direction, where an McrA was called and the lineage has one: 53 / 53 (1.000 [0.932, 1.000]); wrong-direction calls: 0.** The other disagreements are McrA that was not called:

Disagreements (11), by cause:

- **family subunit A absent from the assembly: 6**
  - GCA_964521515_1 (Methanofastidiosales; Methanofastidiosum; new genus): expected methanogenic, mcycle no_call; family score 0, 0 aa
  - GCA_035530415_1 (ANME-2ab; Kmv04; training genus): expected reverse, mcycle no_call; family score 0, 0 aa
  - GCA_034658995_1 (ANME-2c; Methanogaster_D; training genus): expected reverse, mcycle no_call; family score 0, 0 aa
  - GCA_030640725_1 (ANME-2c; Methanogaster_D; training genus): expected reverse, mcycle no_call; family score 0, 0 aa
  - GCA_902386255_1 (p__Halobacteriota; Methanoperedens; training genus): expected reverse, mcycle no_call; family score 0, 0 aa
  - GCA_002067565_1 (p__Methanobacteriota; Methanobacterium_A; new genus): expected methanogenic, mcycle no_call; family score 0, 0 aa
- **subunit A truncated (gene fragment): 4**
  - GCA_015662225_1 (Methanococci; Methanofervidicoccus; new genus): expected methanogenic, mcycle no_call; family score 764, 411 aa
  - GCA_035391185_1 (Methanotrichaceae; Methanocrinis; new genus): expected methanogenic, mcycle no_call; family score 590, 336 aa
  - GCA_009780795_1 (Methanomassiliicoccales; Methanomicula; new genus): expected methanogenic, mcycle no_call; family score 539, 335 aa
  - GCA_012799835_1 (Methanofastidiosales; Methanofastidiosum; new genus): expected methanogenic, mcycle no_call; family score 430, 243 aa
- **full-length McrA homologue under the KO threshold: 1**
  - GCA_024464205_1 (Methanomethylicaceae; Methanomethylicus; new genus): expected methanogenic, mcycle no_call; family score 578, 564 aa

Not scored — lineage of unsettled physiology (3):

- GCF_020885915_1: p__Halobacteriota;c__Methanoliparia;o__Methanoliparales;f__Methanoliparaceae;g__Methanoliparum — mcycle methanogenic
- GCA_013329575_1: p__Halobacteriota;c__Methanoliparia;o__Methanoliparales;f__Methanoliparaceae;g__Methanoliparum — mcycle methanogenic
- GCA_964398245_1: p__UBP6;c__UBA1177;o__UBA1177;f__JAKLFG01;g__CALITV01 — mcycle methanogenic

## pmoA vs other copper monooxygenases

Genomes with a family member: 36; with a lineage expectation: 36; lineage of unsettled physiology: 0.

| stratum | agree / n | agreement [95 % CI] |
|---|---|---|
| **all** | 34 / 36 | 0.944 [0.819, 0.985] |
| new genus | 16 / 16 | 1.000 [0.806, 1.000] |
| training genus | 11 / 13 | 0.846 [0.578, 0.957] |
| training genome | 7 / 7 | 1.000 [0.646, 1.000] |
| named genus | 31 / 33 | 0.939 [0.804, 0.983] |
| placeholder genus | 3 / 3 | 1.000 [0.438, 1.000] |
| new genus, named | 15 / 15 | 1.000 [0.796, 1.000] |
| new genus, placeholder | 1 / 1 | 1.000 [0.207, 1.000] |

Expected (rows) x mcycle (columns):

| expected | not_pmoA | pmoA |
|---|---|---|
| not_pmoA | 14 | 0 |
| pmoA | 2 | 20 |

Disagreements (2), by cause:

- **full-length Cu-MMO subunit A outside the pmoA clade models: 2**
  - GCA_903819965_1 (MOB_Methylococcaceae; Methylumidiphilus; training genus): expected pmoA, mcycle not_pmoA; family score 366, 252 aa
  - GCA_030948025_1 (MOB_alpha; Methylocapsa; training genus): expected pmoA, mcycle not_pmoA; family score 374, 251 aa
