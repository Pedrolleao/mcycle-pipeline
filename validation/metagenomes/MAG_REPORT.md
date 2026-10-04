# MAG realism study (M8)

Twelve inputs — ten published metagenome-assembled genomes with a stated phenotype and two
negative genomes — run from **nucleotide sequence** with `--prodigal-mode meta`, i.e. through
the code path a user with MAGs takes. Roster and expectations were fixed and committed
(`mag_panel.tsv`, commit `ad5a9ef`) before any of them was run. mcycle is the frozen tool of
commit `9c39500`. METABOLIC v4.0 and DRAM v1.4.6 ran on the proteins mcycle called. Quality:
CheckM `lineage_wf --reduced_tree`. Table: `mag_truth_vs_tool.tsv`
(`build_mag_truth_vs_tool.py`).

| MAG | CheckM compl. / contam. | published | mcycle direction | mcycle MMO | clade model | verdict | METABOLIC mcrA / pmoA | DRAM mcr / pmo |
|---|---|---|---|---|---|---|---|---|
| ANME1_Methanospirare_FWG175 | 64.07 / 37.25 | reverse | reverse | none | anme1 | as published | Present / Absent | True / False |
| ANME2a_Methanocomedens_CONS7142H05b1 | 65.03 / 0.98 | reverse | none | none | — | Mcr not in assembly | Absent / Absent | False / False |
| ANME2c_Methanogaster_AMVER4-21 | 57.03 / 4.25 | reverse | none | none | — | Mcr not in assembly | Absent / Absent | False / False |
| ANME3_Methanovorans_HMMV | 98.04 / 3.27 | reverse | reverse | none | anme3 | as published | Present / Absent | True / False |
| ANME2d_HGW_Methanoperedenaceae_1 | 96.41 / 1.96 | reverse | reverse | none | anme2d | as published | Present / Absent | True / False |
| Mmethylicus_V2 | 99.07 / 0.93 | methanogenic | methanogenic | none | — | as published | Present / Absent | True / False |
| Mflorens_bog38 | 95.92 / 0.65 | methanogenic | methanogenic | none | — | as published | Present / Absent | True / False |
| Mliparum_NM1a | 92.46 / 1.31 | methanogenic | methanogenic | none | — | as published | Present / Absent | True / False |
| Argoarchaeum_EthArch1 | 92.27 / 1.66 | no Mcr, no MMO | none | none | — | as published | Present / Absent | True / False |
| USCalpha_Mlahnbergensis | 86.85 / 2.70 | pMMO | none | pmmo | alpha | as published | Absent / Present | False / True |
| Bfragilis_NCTC9343 | 99.26 / 0.00 | no Mcr, no MMO | none | none | — | as published | Absent / Absent | False / False |
| Pmarinus_MED4 | 99.46 / 0.27 | no Mcr, no MMO | none | none | — | as published | Absent / Absent | False / False |

**mcycle: 10 of 12 as published; 2 where the assembly holds no Mcr gene at all; 0 wrong calls.**

## What the study shows

- **Direction from a MAG.** ANME-1, ANME-3 and the second ANME-2d genus are called `reverse`;
  the three methanogen MAGs — a Verstraetearchaeota Mcr, *Ca.* Methanoflorens, and
  *Ca.* Methanoliparum, which carries an alkyl-CoM reductase beside its Mcr — are called
  `methanogenic`. METABOLIC and DRAM report McrA / "key functional gene" present in all six
  and have no way to tell the two groups apart.
- **The ethane oxidizer.** *Ca.* Argoarchaeum ethanivorans has no McrA call in mcycle and its
  `mcrB` is disqualified as an alkyl-CoM reductase subunit. METABOLIC reports mcrA Present and
  DRAM reports the methanogenesis key gene: both read an ethane oxidizer as a methanogen.
- **Atmospheric methane oxidizer.** The USC-alpha MAG is called pMMO by the alpha clade model
  (and by both comparators).
- **Two misses are missing genes, not missed genes.** The ANME-2a and ANME-2c MAGs (CheckM
  65 % and 57 % complete) have no hit to any Mcr subunit profile (K00399, K00401, K00402):
  the operon was not binned. No tool reports it. A lineage expectation is not evidence
  that a gene is in a given assembly.
- **Negatives** are clean in all three tools.

## Caveats

- **In-sample MAGs.** Two MAGs themselves supplied a training sequence to the model that
  calls them — the ANME-3 MAG (HMMV) and the ANME-2d MAG (HGW-Methanoperedenaceae-1) — and
  the *Ca.* Argoarchaeum MAG was a calibration negative. Four more belong to genera
  represented in the training sets by other genomes (*Methanospirare*, *Methanocomedens*,
  *Methanogaster*, *Methylocapsa*). Columns `genome_in_hmm_sets` and `genus_in_hmm_training`
  of the table. For these, the study tests the nucleotide / MAG code path, not
  generalization; generalization is the hold-out (REPORT section 1) and the new-genus
  stratum of the GTDB check (section 4). The three methanogen MAGs and the negatives are
  independent of every model.
- CheckM's generic euryarchaeal marker set reports 64 % completeness and 37 % contamination
  for the 11-contig ANME-1 assembly; it is given as measured.
- The ANME-2d MAG has no activity data; its expectation is by lineage (stated in the roster).
- This study did not contain a MAG with a truncated McrA; the GTDB-500 set did (4 genomes),
  and there the tool fails — see REPORT section 4 and `ROADMAP.md`.
