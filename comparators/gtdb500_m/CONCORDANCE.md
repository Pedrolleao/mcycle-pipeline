# GTDB-500 cross-tool concordance — methane cycle

500 genomes (120 methane-enriched, 380 backbone); tools: mcycle, kofam, mcycdb, metabolic. No ground truth: agreement, not accuracy.

- mcycle: 500 genomes, 4851 present calls
- kofam: 500 genomes, 4791 present calls
- mcycdb: 500 genomes, 9028 present calls
- metabolic: 500 genomes, 4961 present calls

## a. Pairwise agreement (targets both tools can express)

| pair | cells | agreement | Cohen kappa | both present | only first | only second |
|---|---|---|---|---|---|---|
| mcycle vs kofam | 41000 | 0.9926 | 0.964 | 4635 | 202 | 100 |
| mcycle vs mcycdb | 39500 | 0.8746 | 0.568 | 4296 | 291 | 4664 |
| mcycle vs metabolic | 41000 | 0.9920 | 0.962 | 4702 | 135 | 194 |
| kofam vs mcycdb | 39500 | 0.8710 | 0.553 | 4176 | 313 | 4784 |
| kofam vs metabolic | 41000 | 0.9936 | 0.969 | 4685 | 50 | 211 |
| mcycdb vs metabolic | 39500 | 0.8740 | 0.567 | 4317 | 4643 | 333 |

## b. Trap 1 — genomes in which each tool reports pmoA, by clade

| clade | genomes | mcycle | kofam | mcycdb | metabolic |
|---|---|---|---|---|---|
| AOA | 5 | 0 | 1 | 0 | 0 |
| AOB_beta | 5 | 0 | 4 | 0 | 0 |
| AOB_gamma | 5 | 0 | 5 | 5 | 5 |
| MOB_Methylococcaceae | 5 | 2 | 3 | 3 | 3 |
| MOB_Methylomonadaceae | 5 | 5 | 5 | 5 | 5 |
| MOB_NC10 | 5 | 4 | 4 | 4 | 4 |
| MOB_Verrucomicrobia | 5 | 5 | 5 | 5 | 5 |
| MOB_alpha | 5 | 4 | 5 | 5 | 5 |
| Nitrospira_comammox | 5 | 0 | 3 | 0 | 0 |
| backbone | 380 | 0 | 1 | 0 | 0 |

In the 20 genomes of the ammonia-oxidizer clades, pmoA is reported by: mcycle 0; kofam 13; mcycdb 5; metabolic 5.

## c. Trap 2 — McrA and its direction, by clade

| clade | genomes | mcycle mcrA | kofam mcrA | mcycdb mcrA | metabolic mcrA | mcycle reverse | mcycle methanogenic |
|---|---|---|---|---|---|---|---|
| ANME-1 | 5 | 2 | 2 | 2 | 2 | 2 | 0 |
| ANME-2ab | 5 | 3 | 3 | 3 | 3 | 3 | 0 |
| ANME-2c | 5 | 1 | 1 | 1 | 1 | 1 | 0 |
| ANME-2d | 5 | 4 | 4 | 5 | 5 | 4 | 0 |
| ANME-3 | 5 | 4 | 4 | 4 | 4 | 4 | 0 |
| Methanobacteria | 5 | 3 | 3 | 4 | 3 | 0 | 3 |
| Methanocellia | 5 | 5 | 5 | 5 | 5 | 0 | 5 |
| Methanococci | 5 | 4 | 4 | 5 | 5 | 0 | 4 |
| Methanofastidiosales | 5 | 2 | 2 | 4 | 3 | 0 | 2 |
| Methanomassiliicoccales | 5 | 4 | 4 | 5 | 5 | 0 | 4 |
| Methanomethylicaceae | 5 | 3 | 3 | 4 | 4 | 0 | 3 |
| Methanomicrobia | 5 | 3 | 3 | 3 | 3 | 0 | 3 |
| Methanosarcinaceae | 5 | 5 | 5 | 5 | 5 | 0 | 5 |
| Methanotrichaceae | 5 | 3 | 3 | 4 | 4 | 0 | 3 |
| alkane_AcrA | 5 | 2 | 2 | 5 | 5 | 0 | 2 |
| backbone | 380 | 8 | 8 | 9 | 8 | 0 | 8 |

mcycle calls the Mcr `reverse` in 14 genomes; every comparator reports the same genomes as plain McrA carriers.

## d. Targets with the most discordant genomes

| target | mcycle | kofam | mcycdb | metabolic | genomes where the tools differ |
|---|---|---|---|---|---|
| frmA | 26 | 26 | 476 | 26 | 450 |
| pta | 130 | 48 | 441 | 167 | 398 |
| mtkB | 17 | 17 | 341 | 17 | 324 |
| fdhA | 34 | 25 | 325 | 25 | 300 |
| mtkA | 16 | 16 | 295 | 16 | 279 |
| mcl | 21 | 21 | 240 | 21 | 221 |
| acs | 220 | 220 | 435 | 225 | 216 |
| hxlB | 129 | 129 | 320 | 129 | 213 |
| fdhB | 38 | 38 | 241 | 38 | 203 |
| fdh | 102 | 91 | 178 | 91 | 185 |
| mtbC | 15 | 15 | 185 | 15 | 172 |
| mmoC | 3 | 5 | 168 | 5 | 165 |
| hprA | 39 | 39 | 177 | 39 | 152 |
| hdrD | 179 | 179 | 310 | 179 | 133 |
| hxlA | 130 | 130 | 219 | 130 | 119 |
| ackA | 207 | 207 | 325 | 206 | 119 |
| mtaA | 23 | 23 | 126 | 23 | 103 |
| mxaF | 5 | 5 | 103 | 27 | 102 |
| hdrA | 185 | 185 | 278 | 185 | 99 |
| mvhA | 83 | 82 | 137 | 82 | 75 |
