# Benchmark pre-registration — mcycle-pipeline vs external tools

Written and committed on **2026-10-04, before METABOLIC, DRAM or MCycDB were run on
any panel genome.** Design copied from the nitrogen and sulfur sister benchmarks so
that the three report identical endpoints. Later changes go under "Amendments" at the
bottom, dated; the confirmatory family is never changed after the fact.

## What was already known when this was written (disclosure)
- mcycle's own results on the full panel, hold-out included (`../REPORT.md` section 1).
  The tool is frozen at commit `9c39500`; nothing in it changes for this benchmark.
- The raw-KofamScan contrast was already seen as a descriptive table
  (`../compare_kofam.py`, REPORT section 1). It is therefore **not blind**; it stays
  in the family for parity with the sister benchmarks and is flagged as such.
- METABOLIC, DRAM and MCycDB results on the panel: **not known**. Only their static
  definition files were read, to write the vocabulary maps below.

## Hypotheses
- **Primary (confirmatory):** on the homology-trap targets, mcycle has higher
  **precision** than each comparator.
- **Secondary (confirmatory family, same correction):** mcycle has higher **overall
  micro-F1** than each comparator.
- **Exploratory (not corrected with the family):** per-pathway F1; step-resolution
  F1; hold-out-only contrasts; phenotype accuracy (`../phenotype_gt.tsv`);
  TN-inclusive specificity and MCC, with their own BH correction, as in the sisters.

## Panel, inputs, ground truth
- The 49 genomes of `../panel.tsv`, whole panel (frozen 2026-10-04).
- Every tool receives the same proteins: the Prodigal gene calls the pipeline made
  from the genome sequence (`results_ref/<genome>/prodigal/<genome>.faa`). mcycle's
  joining of genes split at an in-frame stop is part of mcycle, not of the input.
- Ground truth: `../curated_function_gt.tsv` (3,173 cells). The KEGG ground truth
  (`../ground_truth.tsv`) is reported as a contrast. Neither was built from any tool's
  output, and neither is edited for this benchmark.

## Comparators (default settings, no asymmetric tuning)
| tool | version | what counts as a call |
|---|---|---|
| mcycle-pipeline | commit `9c39500` | status `confirmed` or `domain-only` |
| raw KofamScan | KOfam release of the cache (2026-05-24), stock `ko_list` thresholds | any KO of the target at or above its threshold |
| METABOLIC | v4.0 (`METABOLIC-G.pl`, env `METABOLIC_v4.0`) | see vocabulary |
| DRAM | v1.4.6 (`DRAM.py annotate_genes` + `distill`, env `DRAM14`) | see vocabulary |
| MCycDB | 2021 release, repository commit `ceba218` (Qian et al. 2022) | best hit, profiler defaults `diamond blastp -k 1 -e 1e-4` |

## Trap-target set (fixed in `../score_mcycle.py` before the panel was scored)
`mcrA`, `mcrA_anme`, `pmoA`, `pmoB`, `pmoC`, `mmoX`, `mxaF`, `xoxF`.

## Vocabulary — fixed before the runs, applied identically to every genome
**Subunit resolution: the 83 mcycle targets.**
- *raw KofamScan*: KO -> every target that lists it.
- *METABOLIC*: the per-genome KO list `KEGG_identifier_result/<genome>.result.txt`
  (METABOLIC's KOfam scan) -> targets by KO, as in the nitrogen benchmark — EXCEPT the
  targets for which METABOLIC has a dedicated, gene-named row in its function table
  (`hmm_table_template.txt`, worksheet 1 of `METABOLIC_result.xlsx`): there the
  worksheet row decides, because that is METABOLIC's own call and it can separate
  what the KO cannot (METABOLIC ships distinct `pmoA.hmm` and `amoA.hmm`, both filed
  under K10944). Gene-named rows used: `pmoA`, `pmoB`, `pmoC`, `mcrA`, `mcrB`, `mcrC`,
  `mmoB`, `mmoD`, `mxaF`, `fae`, `frmA`, `cdhD`, `cdhE`. An `amoA` / `amoB` / `amoC`
  row never counts as pmo.
- *DRAM*: the distillate (`product.tsv`), category "Methanogenesis and methanotrophy".
  A function reported present is expanded to the targets of its KOs (credit and blame
  symmetric, as in the sisters):
  `Key functional gene` -> mcrA, mcrB, mcrG; `acetate => methane, pt 1 / 2 / 3` -> acs /
  ackA / pta; `methanol => methane` -> mtaB; `trimethylamine => dimethylamine` -> mttB;
  `dimethylamine => monomethylamine` -> mtbB; `monomethylamine => ammonia` -> mtmB;
  `putative but not defining CO2 => methane` -> fwdA; `methane => methanol, with oxygen
  (pmo)` -> pmoA, pmoB, pmoC; `methane => methanol, with oxygen (mmo)` -> mmoX, mmoY,
  mmoZ, mmoB, mmoC, mmoD. Targets with no distillate function are predicted absent
  (a coverage limit of the tool, as in the sisters).
- *MCycDB*: gene family -> target of the same name, plus: `hdrA1`, `hdrA2` -> hdrA;
  `hdrB1`, `hdrB2` -> hdrB; `hdrC1`, `hdrC2` -> hdrC; `xoxF1`, `xoxF2`, `xoxF4`,
  `xoxF5` -> xoxF; `fmdB` -> fwdB; `fmdC` -> fwdC; `fae-hps` -> fae and hxlA; `fghA` ->
  frmB; `gck` -> gckA; `fdoG`, `fdhF` -> fdh. An `amoA` / `amoB` / `amoC` hit never
  counts as pmo. No family: fwdA, fwdD, sgaA (predicted absent).
- **`mcrA_anme` for tools that have no direction call** (all four comparators): the
  tool's call is its `mcrA` call — the same "a shared marker counts for every target
  that uses it" rule applied to KOs. Because this target can only hurt a tool that
  cannot express it, the primary endpoint is ALSO reported with `mcrA_anme` removed
  from the trap set (pre-specified sensitivity analysis); a win is claimed only if it
  holds in both.
- **Coverage sensitivity** (pre-specified, descriptive): each comparator is also
  scored on the targets it covers only.

**Step resolution** (mcycle mapped up: a step is present when its diagnostic marker
is; ground truth collapsed by the same rule): `mcr` (mcrA), `mtr` (mtrA),
`co2_to_methyl` (fwdA), `acetoclastic` (cdhC), `methylotrophic` (any of mtaB, mtmB,
mtbB, mttB), `heterodisulfide` (any of hdrA, hdrD), `reverse_methanogenesis`
(mcrA_anme), `pmmo` (pmoA), `smmo` (mmoX), `methanol_oxidation` (any of mxaF, xoxF),
`nod` (nod), `formaldehyde_h4mpt` (mtdB), `rump` (hxlA), `serine_cycle` (mtkA).

## Scoring and statistics (as in `benchmark_stats.py` of the sisters)
- A cell counts only if it is in the ground truth. Positive class = present.
- Unit of independence = genome. Genome-cluster bootstrap, B = 10,000, seed 1234;
  percentile 95 % CIs on precision, recall, F1 per tool and subset.
- Primary significance: paired bootstrap of the difference (mcycle minus comparator,
  same resampled genomes); two-sided p from the bootstrap distribution.
- Secondary: exact McNemar on discordant cells (ignores clustering; reported).
- **Multiple testing: Benjamini-Hochberg over 4 comparators x {trap precision, ALL
  micro-F1} = 8 tests.** Specificity / MCC contrasts: separate BH family.

## Decision rules
- A comparator is "beaten on the trap" iff the 95 % CI of the trap-precision
  difference excludes 0, BH q < 0.05, the point estimate favours mcycle, AND the same
  holds with `mcrA_anme` removed from the trap set.
- An overall-F1 difference whose CI includes 0 is reported as "not resolved by this
  panel", not as parity and not as a win.
- If a comparator beats mcycle on an endpoint, that is reported in the same table.

## Known limitations (written before the results)
- The trap set contains one target no comparator can express (`mcrA_anme`) and three
  (`pmoA/B/C`) whose ground truth was corrected in the direction that favours a tool
  separating pmo from amo. Hence the sensitivity analysis above and the KEGG contrast.
- mcycle's trap precision on this panel has no false positive, so its bootstrap CI
  is degenerate; the differences are driven by the comparators' errors.
- 27 of the 49 genomes informed mcycle's thresholds and models (training). The
  hold-out-only contrast (22 genomes) is the one free of that advantage; it is
  exploratory because it halves the power.
- Comparators may have been built with some panel genomes in their references.
- The vocabulary maps are a bias source; they are fixed here.

## Amendments
(none)
