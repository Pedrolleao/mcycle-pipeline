# GTDB-500 cross-tool concordance — methane cycle

500 genomes (120 methane-enriched, 380 backbone); tools: mcycle, kofam, mcycdb. No ground truth: agreement, not accuracy.

- mcycle: 500 genomes, 4851 present calls
- kofam: 500 genomes, 4791 present calls
- mcycdb: 500 genomes, 9028 present calls

## a. Pairwise agreement (targets both tools can express)

| pair | cells | agreement | Cohen kappa | both present | only first | only second |
|---|---|---|---|---|---|---|
| mcycle vs kofam | 41000 | 0.9926 | 0.964 | 4635 | 202 | 100 |
| mcycle vs mcycdb | 39500 | 0.8746 | 0.568 | 4296 | 291 | 4664 |
| kofam vs mcycdb | 39500 | 0.8710 | 0.553 | 4176 | 313 | 4784 |

## b. Trap 1 — genomes in which each tool reports pmoA, by clade

| clade | genomes | mcycle | kofam | mcycdb |
|---|---|---|---|---|
| AOA | 5 | 0 | 1 | 0 |
| AOB_beta | 5 | 0 | 4 | 0 |
| AOB_gamma | 5 | 0 | 5 | 5 |
| MOB_Methylococcaceae | 5 | 2 | 3 | 3 |
| MOB_Methylomonadaceae | 5 | 5 | 5 | 5 |
| MOB_NC10 | 5 | 4 | 4 | 4 |
| MOB_Verrucomicrobia | 5 | 5 | 5 | 5 |
| MOB_alpha | 5 | 4 | 5 | 5 |
| Nitrospira_comammox | 5 | 0 | 3 | 0 |
| backbone | 380 | 0 | 1 | 0 |

In the 20 genomes of the ammonia-oxidizer clades, pmoA is reported by: mcycle 0; kofam 13; mcycdb 5.

## c. Trap 2 — McrA and its direction, by clade

| clade | genomes | mcycle mcrA | kofam mcrA | mcycdb mcrA | mcycle reverse | mcycle methanogenic |
|---|---|---|---|---|---|---|
| ANME-1 | 5 | 2 | 2 | 2 | 2 | 0 |
| ANME-2ab | 5 | 3 | 3 | 3 | 3 | 0 |
| ANME-2c | 5 | 1 | 1 | 1 | 1 | 0 |
| ANME-2d | 5 | 4 | 4 | 5 | 4 | 0 |
| ANME-3 | 5 | 4 | 4 | 4 | 4 | 0 |
| Methanobacteria | 5 | 3 | 3 | 4 | 0 | 3 |
| Methanocellia | 5 | 5 | 5 | 5 | 0 | 5 |
| Methanococci | 5 | 4 | 4 | 5 | 0 | 4 |
| Methanofastidiosales | 5 | 2 | 2 | 4 | 0 | 2 |
| Methanomassiliicoccales | 5 | 4 | 4 | 5 | 0 | 4 |
| Methanomethylicaceae | 5 | 3 | 3 | 4 | 0 | 3 |
| Methanomicrobia | 5 | 3 | 3 | 3 | 0 | 3 |
| Methanosarcinaceae | 5 | 5 | 5 | 5 | 0 | 5 |
| Methanotrichaceae | 5 | 3 | 3 | 4 | 0 | 3 |
| alkane_AcrA | 5 | 2 | 2 | 5 | 0 | 2 |
| backbone | 380 | 8 | 8 | 9 | 0 | 8 |

mcycle calls the Mcr `reverse` in 14 genomes; every comparator reports the same genomes as plain McrA carriers.

## d. Targets with the most discordant genomes

| target | mcycle | kofam | mcycdb | genomes where the tools differ |
|---|---|---|---|---|
| frmA | 26 | 26 | 476 | 450 |
| pta | 130 | 48 | 441 | 394 |
| mtkB | 17 | 17 | 341 | 324 |
| fdhA | 34 | 25 | 325 | 300 |
| mtkA | 16 | 16 | 295 | 279 |
| mcl | 21 | 21 | 240 | 221 |
| acs | 220 | 220 | 435 | 215 |
| hxlB | 129 | 129 | 320 | 213 |
| fdhB | 38 | 38 | 241 | 203 |
| fdh | 102 | 91 | 178 | 185 |
| mtbC | 15 | 15 | 185 | 172 |
| mmoC | 3 | 5 | 168 | 165 |
| hprA | 39 | 39 | 177 | 152 |
| hdrD | 179 | 179 | 310 | 133 |
| hxlA | 130 | 130 | 219 | 119 |
| ackA | 207 | 207 | 325 | 118 |
| mtaA | 23 | 23 | 126 | 103 |
| hdrA | 185 | 185 | 278 | 99 |
| mxaF | 5 | 5 | 103 | 98 |
| mvhA | 83 | 82 | 137 | 75 |
