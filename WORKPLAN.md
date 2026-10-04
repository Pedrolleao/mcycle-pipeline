# WORKPLAN — bring mcycle-pipeline to the validation level of ncycle / scycle

**Living document.** A new session reads this first, picks the first phase that is not
DONE, and updates the status table and the phase notes before it ends. Written
2026-10-04; execution started the same day.

## Goal and definition of done

`mcycle-pipeline` is built and smoke-tested (18 genomes, 168 expectations, thresholds
tuned on those same genomes). `ncycle-pipeline` and `scycle-pipeline` went through a
full validation campaign. Done = methane satisfies every row of the sister-tool
contract (`../../Unify_Tools/shared/sister-tool-contract.md`) and reports the same
metric battery, computed the same way (genome-cluster bootstrap, B = 10,000, seed 1234,
95 % CI):

1. de-leaked hold-out micro-F1 — the generalization headline
2. full-panel micro-F1
3. independent-only trap precision
4. per-pathway F1 (7 pathways)
5. GTDB-500 cross-tool concordance, with the trap divergence at scale
6. MAG realism — published phenotype vs tool calls, from nucleotide input
7. comparator benchmark vs KofamScan / METABOLIC / DRAM / a methane domain database,
   BH-FDR, pre-registered endpoints

and ships the identical 8-row regression gate (`make regression`): scored cells ·
ALL-F1 CI-lo · ALL-precision CI-lo · ALL-FP cap · trap-precision CI-lo ·
trap-independence-precision CI-lo · hold-out-F1 CI-lo · per-pathway-F1 CI-lo.

## Status

| Phase | What | State |
|---|---|---|
| M0 | Housekeeping and decisions | DONE 2026-10-04 |
| M1 | Reference panel: roster, identity QC, train / hold-out split (frozen) | DONE 2026-10-04 |
| M2 | Dual ground truth | DONE 2026-10-04 |
| M3 | Hardening on the TRAINING genomes only (seeds, clade HMMs, thresholds) | DONE 2026-10-04 |
| M4 | Scoring, leak detection, trap independence, regression gate | DONE 2026-10-04 |
| M5 | Pre-registration, then comparator benchmark | DONE 2026-10-04 |
| M6 | GTDB-500 concordance | IN PROGRESS (mcycle + MCycDB running; METABOLIC on the 120 enriched pending) |
| M7 | Orthogonal (phylogeny-anchored) validation of the two clade calls | DONE 2026-10-04 |
| M8 | MAG realism study | DONE 2026-10-04 |
| M9 | Report, docs, contract, reproducibility pins | TODO |

Order: M0 → M1 → M2 → M3 → M4 → M5; then M6, M7 (needs M6), M8 in any order; M9 last.
M1 before M3 is not negotiable: the hold-out must be frozen before anything is tuned.

## Where things are

- Tool: `Methane_Cycle/mcycle-pipeline/` (local git repo, branch `main`, no remote). Atlas with every measured
  number so far: `../Info-methane.md`. Config: `config/targets.yaml` (83 targets,
  17 complexes, 11 modules). Smoke panel: `../test_panel/` (`panel.tsv`, `fetch_panel.sh`).
- Run: `python run.py --input <dir> --prodigal-mode single --cores 8`; `make smoke`.
  Conda env `cycle-pipeline` (shared by the three tools). `python` is not on PATH
  outside the env — activate it or call `~/miniconda3/envs/cycle-pipeline/bin/python`.
- Templates to port from — the sulfur repo is the most recent and already harmonized:
  `../../Sulfur_Cycle/scycle-pipeline/validation/` and `comparators/`; nitrogen-only
  pieces in `../../Nitrogen_Cycle/ncycle-pipeline/validation/`
  (`detect_seed_leakage.py`, `validate_panel.py`, `verify_seeds.py`).
  History and lessons of that campaign: `../../Unify_Tools/WORKPLAN.md`.
- Comparator installs already on this machine: conda envs `METABOLIC_v4.0` and `DRAM14`;
  METABOLIC repo and DRAM data under `../../Sulfur_Cycle/comparators/`; GTDB r232
  taxonomy tables and the 380-genome backbone selection under
  `../../Nitrogen_Cycle/ncycle-pipeline/comparators/{gtdb_pilot,gtdb500}/`.
  Also present: `CheckM_env`, `Iqtree_env`, `Clipkit_env`.
- Methane engine differs from the sisters in four places the ported scripts must
  respect: modules use `requires_any` and `key` (see `compute_complex_completeness.py`);
  the direction call is `mcr_methanogenic` / `mcr_reverse` on the mcrA/B/G
  `evidence_source`; pmoB / pmoC follow pmoA (`gate_pmo_subunits`); genes split at an
  in-frame stop are joined (`join_inframe_stops`, tag `joined:`).
- The report template, `make_html_report.py` and `fonts/` are byte-identical in the
  three repos; a change to one must be copied to all three. Browser test for the report:
  `node ~/tools/shotter/report-export-test.mjs <report.html> <outdir>`.

## Ground rules (each one cost the sister campaigns real time)

- **Identifiers are verified, never recalled**: KO against the cached `ko_list`, UniProt
  and NCBI accessions against their APIs, at the moment they are written down.
- **Ground truth comes from the literature and KEGG, never from this tool's output.**
  Every curated correction records its evidence.
- **Freeze before you look**: the hold-out roster (M1) and the pre-registration (M5) are
  written and dated before the corresponding results exist. Later panel changes are
  allowed only as disclosed, dated additions in a CHANGELOG.
- **No seed, HMM training sequence or threshold calibration may come from a hold-out
  genus.** Two smoke-panel genomes are already seed strains (*M. capsulatus* Bath,
  *M. trichosporium*) — they are training genomes by definition.
- **Check that each panel file is the organism it claims** (`validate_panel.py`); a
  mislabelled proteome invalidated three nitrogen ground-truth corrections.
- **Genome sequence is the default input**, so the MAG code paths (Prodigal, in-frame
  stop joining, locus maps) are what gets validated.
- **Comparators**: one pristine input directory per tool (METABOLIC writes `total.faa`
  into its input); never two METABOLIC runs at once; batches of ≤ 12 genomes at GTDB
  scale. After any adapter fix, re-run every consumer (panel benchmark, GTDB
  concordance, figures).
- **Operational**: `run.py` overwrites `results/` and rewrites `config/config.yaml` —
  use a separate config and results dir per study
  (`MCYCLE_CONFIG=config/config_<study>.yaml`, as `config_p3_mags.yaml` does in sulfur).
  `run.py` runs Snakemake without `--keep-going`: one empty proteome aborts the DAG.
  Network calls need the sandbox disabled.

## Phases

### M0 — Housekeeping and decisions
- Ask the user the open decisions at the bottom of this file.
- `git init` + first commit if agreed, so every later change is diffable.
- Give the tool its own KOfam cache copy or keep the symlink (contract: self-contained,
  SHA-pinned; the build downloads its own copy if the links are removed).
- Done when: decisions recorded here; repo state agreed.
- **DONE 2026-10-04.** `git init` (branch `main`, author identity as in the sister
  repos, no remote — the user wants the repo kept private: never push or create a
  remote without asking, and any future remote must be private). KOfam cache stays a
  symlink to the nitrogen cache (gitignored; the build fetches its own SHA-pinned copy
  if the links are removed).

### M1 — Reference panel (target ≈ 45 genomes, ≈ 20 of them hold-out)
- Write `validation/PANEL_PLAN.md`: roster with NCBI accession, guild, published
  physiology (one citation each), train / hold-out, and why it is there. Guilds to cover
  — candidates to verify, not a final list:
  - methanogens of every order: Methanobacteriales, Methanococcales, Methanopyrales,
    Methanomicrobiales, Methanocellales, Methanosarcinales (incl. an obligate
    methylotroph and a second *Methanothrix*), Methanomassiliicoccales,
    Methanonatronarchaeia; one H₂ + methanol specialist (*Methanosphaera*)
  - ANME-2d from a second genus-level lineage; marine ANME only if a closed or
    high-quality genome exists (otherwise they go to M8)
  - alkane oxidizers with alkyl-CoM reductase (decoys for Mcr)
  - aerobic methanotrophs: γ type I and type X, α (pMMO-only, sMMO-only, both),
    Verrucomicrobia, NC10
  - methylotrophs without a methane monooxygenase
  - decoys: β- and γ-AOB, AOA, comammox; a hydrocarbon-oxidizing Cu-monooxygenase
    carrier; butane / propane / toluene di-iron monooxygenase carriers; a lanthanide
    ethanol-dehydrogenase carrier; an acetogen and sulfate reducers (Wood–Ljungdahl,
    Hdr-like); general negatives
- Hold-out = genera absent from every seed list and from every tuning step so far;
  the 18 smoke genomes are training (thresholds were tuned on them).
- Extend `../test_panel/fetch_panel.sh` usage to a `ref_panel/` directory; port
  `validate_panel.py` and run it.
- Done when: roster frozen and dated; all files fetched and identity-checked.
- **DONE 2026-10-04.** 49 genomes, 27 training / 22 hold-out, frozen in
  `validation/panel.tsv` (+ `PANEL_PLAN.md`, `CHANGELOG.md`); FASTA in `../ref_panel/`
  (`make ref-panel`); identity QC 49 / 49 (`make validate-panel`, `panel_qc.tsv`).
  Things later phases must know:
  - "genus" = NCBI genus and GTDB r232 genus; the hold-out genera are listed in
    PANEL_PLAN.md and are banned as seed / HMM-training / calibration sources.
  - No hold-out is possible for ANME-2d, NC10, alpha pMMO and acetoclastic
    methanogenesis (every genus is a seed or smoke genus) — say so in the report.
  - ANME-2d second lineage (`g__Methanoperedens_A`) and marine ANME-2a / -2c / -3 have
    only fragmented MAGs: they are M8 material. ANME-1 is in the panel (G60 training
    and intended McrA seed source, G37 hold-out).
  - 11 genomes are not in KEGG (5 hold-out, 6 training); four smoke genomes have
    `kegg = tbd`. `rest.kegg.jp/list/organism` now answers 400 — use `list/genome`.
  - *Thauera butanivorans* is the type-strain assembly (156 contigs), not the newer
    complete genome of an uncharacterized strain.

### M2 — Dual ground truth
- Port `build_ground_truth.py` → KEGG-KO ground truth (`ground_truth.tsv`: genome,
  target, expected, source) through KEGG REST `link/ko/<org>`. Expect partial coverage:
  uncultured and candidate genomes are not in KEGG.
- Port `build_curated_function_gt.py` → `curated_function_gt.tsv` = KEGG ground truth +
  function-level corrections on the cells where KO presence ≠ function, each with
  evidence: pmoA/B/C in ammonia oxidizers (KEGG assigns the shared KOs), `mcrA_anme`
  (no KO — fully manual), `mmoX` vs other di-iron monooxygenases, `xoxF` vs ethanol
  dehydrogenases, recoded genes missing from KEGG gene lists.
- Methane-specific addition: `phenotype_gt.tsv` — per genome, the published
  methanogenesis type(s), methane-oxidation type and Mcr direction, for scoring the
  module and direction calls (secondary endpoint).
- Done when: both files build reproducibly; every correction has a rationale row.
- **DONE 2026-10-04** (`make ground-truth`). `ground_truth.tsv` = KEGG, 37 genomes,
  3,034 cells. `curated_function_gt.tsv` = 3,173 cells over 49 genomes: KEGG + 21
  corrections + 90 literature cells (`curated_cells.tsv`, evidence and DOI per row) +
  49 `mcrA_anme` cells derived from `phenotype_gt.tsv`. Notes for later phases:
  - The corrections are all `present -> absent` (pmoA/B/C in ammonia oxidizers and
    *M. chubuense*): they favour a tool that separates pmo from amo. Report metrics on
    BOTH ground truths, as sulfur does.
  - KEGG needed no correction for mmoX, mxaF or xoxF in the decoys (it does not assign
    K16157 / K14028 / K23995 to *Thauera*-type, PedH-type or *Gluconobacter* enzymes that
    it holds) — those trap cells are plain KEGG `absent`.
  - KEGG-less genomes carry 4-22 scored cells each; *S. pneumoniae* carries one.
  - Disagreements found once the tool is run may only become corrections with
    independent evidence, logged in `CHANGELOG.md` (three candidates are listed there).

### M3 — Hardening, on training genomes only
From the backlog in `../Info-methane.md`:
- Seeds for marine ANME McrA (ANME-1, -2a, -2c, -3). UniProt had no verified
  full-length sequence; allow `build_blast_db.py` to take a curated local FASTA with a
  provenance file, or use NCBI protein accessions.
- Train clade HMMs where an identity gate is brittle — pmoA vs amoA (gate 70, decoy at
  63.5) and ANME McrA (gate 80, decoy at 70): `build_custom_hmms.py`,
  `calibrate_tc.py`, manifests under `targets/<id>/`, leave-one-clade-out with a ported
  `logo_cv.py`.
- Calibrate `mmoX`, `xoxF`, `nod`, NC10 pmoA and the alkyl-CoM reductase boundary
  against the new training decoys.
- Examine the probable KO-tier noise (`mtrA` in AOA, `mtbC`, `mtsA`).
- Done when: every change is justified by training genomes only and logged; `make
  smoke` still passes.
- **DONE 2026-10-04.** Log with every number: `validation/HARDENING.md`. Training
  micro-F1 0.964 -> 0.983, trap precision and recall 1.000 (in-sample); smoke 168 / 168.
  - Engine: `gate_mcr_subunits`, `gate_mmo_subunits`, `resolve_fdh`, per-KO `ko_tc`,
    query-coverage floor for KO profiles, clade models `<target>__<clade>` with
    `clade_hmm_required`.
  - pmoA and mcrA_anme are now decided by clade HMMs (4 + 5 models), built with
    `harvest_clade_refs.py` + `build_clade_hmms.py` from GTDB-typed genomes cached in
    `../clade_refs/` (317 genomes). These two scripts replace the UniRef90-based
    `build_custom_hmms.py expand` / `calibrate_tc.py` / `logo_cv.py` route of the
    sisters: leave-one-genus-out is built into `build_clade_hmms.py`.
  - Runs: `MCYCLE_CONFIG=config/config_ref.yaml python run.py --input ../ref_panel_train
    ...` writes `results_ref/`; scoring with `MCYCLE_SCOPE=train`. `../ref_panel_holdout/`
    holds the 22 hold-out genomes and had NOT been run when M3 closed.
  - Known limits to report: ANME-3 not cleanly separable; no comammox / mycobacterial
    negatives (hold-out genera); verruco and NC10 pmoA models rest on 2 and 1 genera.
  - `hmmbuild` takes `-n`, not `--name` (the older `build_custom_hmms.py` uses the
    wrong flag and would fail if it were ever run).

### M4 — Scoring and gate
- Port `score_scycle.py` → `score_mcycle.py` (define `TRAP` and `HOLDOUT`), N's
  `detect_seed_leakage.py`, `trap_independence.py`, `compare_kofam.py`,
  `test_regression.py`; add `make regression` and `make regression-score`.
- Floors = observed CI lower bounds minus a small slack, set once and written down with
  the observed values (as in the sulfur gate).
- Done when: the 8-row gate passes and the numbers are in `validation/REPORT.md`.
- **DONE 2026-10-04.** Hold-out micro-F1 0.961 [0.940, 0.976]; full panel 0.973
  [0.963, 0.981]; trap precision 1.000 (0 FP in 259 trap negatives). Gate 14 / 14
  (`make regression-score`). Report: `validation/REPORT.md` section 1. The tool was at
  commit `9c39500` when the hold-out was first run — any later change to rules, seeds,
  thresholds or models must be disclosed as post-hold-out, and the pre-change numbers
  kept. `compare_kofam.py` is the descriptive raw-KOfam table; `mcrA_anme` under a raw
  KO is "McrA present" — M5 must fix, in the pre-registration, how tools that cannot
  call direction are scored on that target.

### M5 — Pre-registration, then benchmark
- Write `validation/benchmark/prereg.md` BEFORE running any comparator: primary
  endpoint trap precision vs each comparator; secondary ALL micro-F1; exploratory
  per-pathway F1 and phenotype accuracy; trap-target set; BH-FDR family.
- Port `adapters.py` (methane vocabulary and step collapse), `benchmark_stats.py`,
  `plot_benchmark.py`. Raw-KofamScan baseline comes free from the pipeline's own
  `*.hmmscan.tsv` — check the column layout the loader assumes.
- Run METABOLIC and DRAM on the Prodigal proteomes; write `build_metabolic_tsv.py` /
  `build_dram_tsv.py` mappings for methane. Domain database: MCycDB (the methane
  counterpart of NCycDB / SCycDB) — confirm it exists and is obtainable, stage it
  untracked, write `build_mcycdb_tsv.py`.
- Done when: `benchmark_results.tsv` and `COMPARISON_REPORT.md` exist for 4 comparators.
- **DONE 2026-10-04.** Pre-registration committed at `3cc200c` before any comparator
  ran; benchmark run as registered (B = 10,000). mcycle better on all 8 contrasts of
  the family (BH q < 0.0001), also without `mcrA_anme` and on the hold-out alone.
  `validation/benchmark/{prereg.md, adapters.py, benchmark_stats.py, benchmark_final.txt,
  benchmark_results.tsv, COMPARISON_REPORT.md, figures/}`; normalized comparator tables
  `metabolic.tsv`, `dram.tsv`, `mcycdb.tsv`; runners and builders under `comparators/`.
  - METABOLIC: 5 batches of 10, ~22 min each. DRAM: 7 parallel batches of 7 with 4
    threads finished in ~70 min (one sequential run would have taken ~14 h).
  - MCycDB lives under `comparators/MCyc/` (untracked); the split zip has to be
    inflated by hand (see `comparators/build_mcycdb_tsv.py`).
  - `plot_benchmark.py` is the sister file, unchanged; its figures 2 and 3 expect the
    sister schema of CONCORDANCE / DIR_ACCURACY and are skipped for methane.
  - Never overwrite a shell script that a running bash is executing (bash reads the
    file as it goes): the DRAM runner was replaced mid-run and the old process carried
    on into the new text. The outcome was the intended batched run, but by accident.

### M6 — GTDB-500 concordance
- `select_gtdb_mcyc.py`: reuse the 380-genome backbone verbatim, add ≈ 120
  methane-enriched species representatives (methanogen orders, ANME, alkane-oxidizing
  archaea, the methanotroph families, ammonia oxidizers as decoys).
- Run mcycle on all; reuse the backbone's existing METABOLIC output; run the enriched
  set in small batches. Port `concordance.py` with the spotlight on the two methane
  traps: amoA called as pmoA, and Mcr direction.
- Done when: `CONCORDANCE.{md,tsv}` and `GTDB500_REPORT.md` exist.
- **IN PROGRESS.** `comparators/gtdb500_m/`: `select_gtdb_mcyc.py` -> `selection.tsv`
  (380 backbone + 120 enriched, 24 clades x 5); `prepare_proteomes.py` (backbone
  proteomes linked from ncycle's gtdb500, enriched called with Prodigal -p meta) — 500 /
  500 ready. Running: mcycle (`MCYCLE_CONFIG=config/config_gtdb500.yaml`, protein input,
  `results_gtdb500/`) and MCycDB (`gtdb500_m/mcycdb.tsv`). To do: METABOLIC on the 120
  enriched (batches of <= 12, only after the panel METABOLIC run has ended), merge with
  the backbone KO lists and worksheet 1 of
  `Nitrogen_Cycle/ncycle-pipeline/comparators/gtdb500/metabolic_out/`, then
  `MCYCLE_RESULTS=results_gtdb500 python validation/benchmark/concordance.py ...`.

### M7 — Orthogonal validation of the clade calls
- Curate typed reference sets with provenance: McrA (methanogen / ANME clades /
  alkyl-CoM reductase) and Cu-monooxygenase subunit A (pmoA clades / amoA clades /
  hydrocarbon monooxygenases).
- Best-hit scorer, then tree placement (MAFFT, ClipKIT, IQ-TREE) for the GTDB calls,
  split by characterized genera vs candidate lineages; report which calls no sequence
  method resolves.
- Done when: `DIR_ACCURACY` and `DIR_PLACEMENT` tables exist for both traps.
- **DONE 2026-10-04.** `comparators/gtdb500_m/DIR_ACCURACY.{md,tsv}` (clade call vs
  the GTDB lineage of the genome — the independent reference, since mcycle's call is
  itself sequence-based) and `DIR_PLACEMENT{,_mcr,_pmo}.{md,tsv}` (gene trees;
  `validation/phylogeny/place_clades.py`, ClipKIT from `Clipkit_env`, IQ-TREE 2 from
  `Iqtree_env`). Direction 53 / 53 where an McrA is called; pmoA 34 / 36, 16 / 16 in
  new genera. **Defect found:** an McrA truncated at a contig end (4 GTDB genomes) or
  a full-length non-euryarchaeal McrA under the threshold (1) is not called, and
  `gate_mcr_subunits` then disqualifies mcrB / mcrG as alkyl-CoM reductase subunits.
  Not fixed — the tool stays at `9c39500` for the whole campaign. Fix sketch in
  `ROADMAP.md` (bits per aligned position separate a fragment, ~1.8, from an alkyl-CoM
  reductase, ~1.1); it must be followed by `make regression` and a changelog entry.

### M8 — MAG realism
- Dossier of ≈ 8–10 published MAGs with a stated phenotype and CheckM quality: marine
  ANME-1 / -2a / -2c, a second ANME-2d, a non-euryarchaeal Mcr lineage, an alkane
  oxidizer, uncultured methanotrophs (e.g. atmospheric-methane clades, NC10), a
  methanogen MAG, two negatives.
- Run from nucleotide with `--prodigal-mode meta`; build the truth-vs-tool table with
  METABOLIC and DRAM alongside. Expect this phase to surface real bugs — in nitrogen it
  found one the isolate panel could not.
- Done when: `mag_truth_vs_tool.tsv` and a short write-up exist.
- **DONE 2026-10-04.** Roster and expectations committed (`ad5a9ef`) before the runs.
  `validation/metagenomes/{mag_panel.tsv, mag_truth_vs_tool.tsv, MAG_REPORT.md,
  checkm.tsv, build_mag_truth_vs_tool.py}`; genomes in `../mag_panel/`; mcycle results
  in `results_mags/` (`config/config_mags.yaml`); METABOLIC in
  `comparators/mags_metabolic_out/`, DRAM in `comparators/mags_dram_out/`. mcycle 10 / 12
  as published, 2 with no Mcr gene in the assembly, 0 wrong. Two ANME MAGs were
  training genomes of the model that calls them (column `genome_in_hmm_sets`).

### M9 — Close out
- `validation/REPORT.md` with the full battery; update `README.md`, `ROADMAP.md`,
  `../Info-methane.md`, `../README.md`.
- Reproducibility: lock file for the env, archive of `resources/blast_db/`, CI workflow
  adapted from the sister template.
- Add a methane column to the sister-tool contract if the user wants it there.

## Decisions (answered by the user 2026-10-04)

1. Git: yes, commit at the end of each phase; repo stays private (local only).
2. Panel: ≈ 45 genomes, ≈ 20 hold-out.
3. Clade HMMs: in scope for M3 (pmoA vs amoA, ANME McrA).
4. MCycDB is the domain-database comparator, staged locally and untracked.
5. GTDB-500: full 4-tool run, METABOLIC included (background batches of ≤ 12).
6. Methane gets a column in the `Unify_Tools` contract at M9; the manuscript is not
   touched without asking again.
7. Engine merge: not asked, out of scope.
8. KOfam cache: keep the symlink.

The questions as they were put:

1. Put `mcycle-pipeline` under git now, and commit at each phase?
2. Panel size and split: ≈ 45 genomes with ≈ 20 hold-out, or larger?
3. Is training clade HMMs in scope (M3), or should the BLAST gates stay as they are and
   be reported as such?
4. MCycDB as the domain-database comparator — acceptable, and can it be staged locally?
5. GTDB-500 with METABOLIC is many hours of compute; run it in full, or the tool +
   KofamScan + domain database first?
6. Should methane join the `Unify_Tools` contract and manuscript, or stay separate?
7. Optional and separate from this plan: merging the three engines into one.
