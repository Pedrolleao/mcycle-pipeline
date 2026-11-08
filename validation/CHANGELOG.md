# Validation changelog

Dated record of every change to the reference panel, the ground truth and the
pre-registration after they were frozen. Newest last.

## 2026-10-04 — panel frozen (M1)
- `panel.tsv`: 49 genomes, 27 training / 22 hold-out. Rationale in `PANEL_PLAN.md`.
- Identity QC 49 / 49 pass (`panel_qc.tsv`).
- No ground truth existed and the pipeline had not been run on any added genome when
  the split was fixed.

## 2026-10-04 — ground truth built (M2)
- `panel.tsv`, column `kegg` only: `Mfumariolicum_SolV` tbd -> `mfh` (KEGG holds assembly
  GCA_949774925.1, the same sequence as the panel's GCA_000953475.1 after rotation and
  reverse complement — checked base by base); `Mluminyensis_B10` tbd -> `-`.
- `ground_truth.tsv` (KEGG, 37 genomes x 82 KO-anchored targets = 3,034 cells).
- `curated_cells.tsv`: 21 corrections (pmoA / pmoB / pmoC in six ammonia oxidizers and
  the hydrocarbon-monooxygenase carrier *M. chubuense* NBB4) and 90 literature cells for
  11 genomes KEGG does not hold. `phenotype_gt.tsv`: phenotype of all 49 genomes.
- `curated_function_gt.tsv`: 3,173 cells over 49 genomes.
- All of it was written before the pipeline was run on any genome added in M1.
- Looked at and deliberately NOT corrected, for lack of verified evidence:
  *M. burtonii* `mtmB` / `mtbB` (KEGG and UniProt list only `mttB`, `mtbC`, `mtbA`; the
  organism grows on methylamines, so the genes may be unannotated pyrrolysine genes);
  *M. album* BG8 `mmoX` (no sMMO is the textbook statement, but the genome paper does
  not say it — the cell is unscored); *N. maritimus* `mtrA` (KEGG assigns K00577).

## 2026-10-05 — tool amendment after the hold-out was scored: truncated McrA (M10)
Decided by the user on 2026-10-04 after the GTDB-500 check (REPORT section 4) found
the defect. Applied only after M6 was closed and committed with the frozen tool.
- **Tool states.** Frozen: commit `9c39500` — every number of REPORT sections 1-5.
  Amended: commit `c7bf634`.
- **Change.** `gate_mcr_subunits` (`workflow/scripts/apply_rules.py`): a sub-threshold
  McrA homologue is evidence of an alkyl-CoM reductase — and mcrB / mcrG are then
  `disqualified` — only if its score per aligned profile position is below threshold /
  profile length (775.5 / 556 = 1.39 bits). A gene fragment of a canonical McrA scores
  above that (1.62-1.86 measured) and no longer triggers the rule. No fitted cut-off.
  The four genomes that exposed the defect are GTDB-500 genomes; no panel genome was
  used. Patch: `amendments/2026-10-04_truncated_mcra.patch`.
- **Effect** (`amendments/2026-10-05_truncated_mcra_changes.tsv`: all 48,057 calls of the
  smoke panel, the reference panel, the MAG study and the GTDB-500 set, before vs after):
  8 calls change — mcrB and mcrG `disqualified -> confirmed` in GCA_009780795,
  GCA_012799835, GCA_015662225, GCA_035391185.
  - Reference panel: 0 calls change. Hold-out micro-F1, full-panel micro-F1, trap
    precision, the gate (14 / 14) and the comparator benchmark are identical before and
    after (`mcycle_confusion.tsv` unchanged).
  - MAG study: 0 calls change. Smoke: 168 / 168.
  - GTDB-500: genomes outside alkane-oxidizer lineages with mcrB / mcrG labelled as
    alkyl-CoM reductase subunits: 5 -> 1. Alkane-oxidizer genomes so labelled: 3 -> 3.
    Agreement with raw KofamScan 0.9926 -> 0.9928 (`CONCORDANCE_amended.md`).
    `DIR_ACCURACY_amended.md`: unchanged, because the McrA itself is still not called.
- **Not changed, by decision — a separate cycle with its own validation** (`ROADMAP.md`):
  calling the truncated McrA itself, and non-euryarchaeal McrA that fall under the KO
  threshold at full length (the one genome still mislabelled, *Ca.* Methanomethylicus
  GCA_024464205).
- The post-amendment GTDB counts are in-sample for this rule: the genomes that showed
  the defect are the ones on which the repair is counted.

## 2026-10-05 — KOfam built from a snapshot in the repository (no change to any call)
Made when the repository was prepared for cloning (commit `b7ea805`).
- **Why.** genome.jp serves KOfam from one rolling URL. The release now there
  (`profiles.tar.gz` of 2026-09-29) has a different threshold for 85 of the tool's 90
  KOs (K00399: 775.53 -> 741.47; K00192: 612.00 -> 423.30), so a database built from a
  fresh download is not the validated one.
- **Change.** `resources/kofam_pinned/` holds the 90 profiles and their `ko_list` rows
  from the validated release (2026-05-24, SHA-256 `b03d20b9…` of the full tarball);
  `build_hmm_db.py` builds from it and keeps the download as `--upstream` / fallback.
  `apply_rules.py` and `config/targets.yaml` are untouched.
- **Check, from a fresh clone.** `mcycle_targets.hmm`, `tc_cutoffs.tsv` and both seed
  FASTAs (re-fetched from UniProt) are byte-identical to those of the validation; the
  smoke panel meets 168 / 168 expectations and its 18 call tables are byte-identical
  to the previous run; `make regression` passes its 14 checks with the numbers of the
  report (hold-out micro-F1 0.961 [0.940, 0.976], full panel 0.973 [0.963, 0.981],
  3,173 cells) and rewrites the committed score tables without a difference.

## 2026-11-07 — code-review changes, no call changes (R1–R8 of `WORKPLAN_2026-11_review.md`)
Eight small, separately committed changes after the Fable review of commit `97fd402`;
every one re-ran `make smoke` (168 / 168) and compared the 18 call tables byte for byte
with the snapshot of that commit (`validation/baseline_2026-11.sha256`, 54 tables:
`mcycle_calls.tsv`, `complex_completeness.tsv`, `synergy_completeness.tsv`). Nothing in
the tool's calls, module statuses or figures changed; the only cell that differs is the
`application` text of the `mcr` row of `complex_completeness.tsv` (R6, a developer note
removed).
- **R1** the launcher refuses sample names outside the Snakemake wildcard pattern,
  duplicate stems and malformed FASTA before anything runs (used to die at DAG build, run
  only the last duplicate, or hang in `diamond blastp`).
- **R2** the error panel shows the failed job's log tail and the right Snakemake log path.
- **R3** `hmmscan` on one thread per job: clean smoke run 67 s → 30 s on 8 cores, all 18
  `hmmscan` tables identical.
- **R4** every figure rendered once and saved as SVG + PNG (22 / 22 figures pixel- and
  text-identical, `validation/compare_figures.py`).
- **R5** per-run config `<results_dir>/mcycle_run_config.yaml` (tool version, date, command
  line, samples) instead of rewriting `config/config.yaml`; the checkout stays clean.
- **R6** `gap_analysis.txt` rewritten (present / complexes / process modules, with the
  input kind in the header); the note "⚠️ direction set by resolve_mcr_direction" became a
  YAML comment. The `targets.yaml` edit rebuilt both databases: `mcycle_targets.hmm` and
  `tc_cutoffs.tsv` byte-identical; the two seed FASTAs differ in UniProt's header text
  only — the 36 diamond tables of the smoke panel made with the validated database are
  byte-identical to those made with the rebuilt one.
- **R7** (one module-fraction convention) measured and **not applied**: it would flip
  `acetoclastic_methanogenesis` from partial to absent in 7 genomes of the reference,
  GTDB-500 and MAG sets (details in the work plan) — a change for its own cycle.
- **R8** static grids above 120 genomes become a pointer figure (`figures.max_rows`); the
  HTML report's column labels stay in view while scrolling.
- The shared files were ported to `ncycle-pipeline` and `scycle-pipeline` (R1, R2, R4, R5,
  R8; R3 to scycle, which also runs `hmmscan`); on their 5-genome test data every call
  table and figure is identical before and after the port.

## 2026-11-07 — amendment: module completeness counts the required slots only (R7)
- **What.** `compute_complex_completeness.py` no longer counts a matching Mcr direction as
  an extra slot of the three methanogenesis modules. `completeness`, `n_present` and
  `n_total` in `synergy_completeness.tsv` (and the `synergy__*` columns of
  `mcycle_matrix.tsv`) now count the required slots, which is what the figures, the HTML
  report and `gap_analysis.txt` already showed. The direction and key-gene rules still
  decide the status and are listed under `forbids_satisfied` / `forbids_violated`.
- **Why.** The TSV and the figures disagreed (6 / 13 against 6 / 12 for the same module),
  and the extra slot lifted `acetoclastic_methanogenesis` to exactly 5 / 10 = 0.500
  `partial` in methyl-reducing methanogens that carry Mcr and HdrD but none of the
  acetoclastic genes.
- **Effect** (old against new script on the unchanged calls of every result set;
  `validation/amendments/2026-11-07_module_fraction_changes.tsv` lists each row): 5,549
  module rows change `n_present` / `n_total`, 1,474 of them the `completeness` value,
  7 the status — `acetoclastic_methanogenesis` partial → absent in `Malvi_Mx1201`
  (reference panel), `GCA_001717005.1`, `GCA_023402275.1`, `GCA_029887795.1`,
  `GCA_035404685.1`, `GCA_963586845.1` (GTDB-500) and `Mmethylicus_V2` (MAG study), all
  Methanomassiliicoccales / Methanomethylicia / Methanofastidiosales. No gene call
  changes; no status change in the smoke panel (168 / 168); the reference-panel scores
  (gene calls) are unaffected. The smoke panel's 22 figures are pixel-identical before
  and after; 67 cells of `mcycle_matrix.tsv` (all `synergy__*`) change value.
- **Shared code.** `_completeness.load_synergies` reads the per-tool convention from
  `_domain.MODULE_RULES_COUNT_AS_SLOTS` (False here, True in ncycle and scycle, whose
  TSVs and figures are unchanged).

## 2026-11-08 — BLAST seeds pinned in the repository (no call changes)
- `resources/seeds_pinned/unstable_refs.fasta` (19 records) and `blast_gated_refs.fasta`
  (9 records): the UniProt seed FASTAs behind the DIAMOND databases, copied from the
  rebuild of 2026-11-07 (R6). The files of the 2026-10-04 fetch, used for the validation
  of 2026-10-05, differ from these in UniProt header text only (36 diamond tables and
  54 call tables of the smoke panel byte-identical, reference panel 14/14 with unchanged
  score tables; `WORKPLAN_2026-11_review.md`, R6 notes).
- `workflow/scripts/build_blast_db.py` copies and indexes the pinned files when they
  cover every curated accession of `config/targets.yaml`, and fetches from UniProt only
  with `--upstream` or after a seed is added — the same logic as ncycle / scycle. A
  `targets.yaml` edit therefore no longer needs the network.
- Acceptance: rebuild from the snapshot, `make smoke` 168/168, all 54 call tables
  byte-identical to `validation/baseline_2026-11.sha256`.
