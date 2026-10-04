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
| M1 | Reference panel: roster, identity QC, train / hold-out split (frozen) | TODO |
| M2 | Dual ground truth | TODO |
| M3 | Hardening on the TRAINING genomes only (seeds, clade HMMs, thresholds) | TODO |
| M4 | Scoring, leak detection, trap independence, regression gate | TODO |
| M5 | Pre-registration, then comparator benchmark | TODO |
| M6 | GTDB-500 concordance | TODO |
| M7 | Orthogonal (phylogeny-anchored) validation of the two clade calls | TODO |
| M8 | MAG realism study | TODO |
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

### M4 — Scoring and gate
- Port `score_scycle.py` → `score_mcycle.py` (define `TRAP` and `HOLDOUT`), N's
  `detect_seed_leakage.py`, `trap_independence.py`, `compare_kofam.py`,
  `test_regression.py`; add `make regression` and `make regression-score`.
- Floors = observed CI lower bounds minus a small slack, set once and written down with
  the observed values (as in the sulfur gate).
- Done when: the 8-row gate passes and the numbers are in `validation/REPORT.md`.

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

### M6 — GTDB-500 concordance
- `select_gtdb_mcyc.py`: reuse the 380-genome backbone verbatim, add ≈ 120
  methane-enriched species representatives (methanogen orders, ANME, alkane-oxidizing
  archaea, the methanotroph families, ammonia oxidizers as decoys).
- Run mcycle on all; reuse the backbone's existing METABOLIC output; run the enriched
  set in small batches. Port `concordance.py` with the spotlight on the two methane
  traps: amoA called as pmoA, and Mcr direction.
- Done when: `CONCORDANCE.{md,tsv}` and `GTDB500_REPORT.md` exist.

### M7 — Orthogonal validation of the clade calls
- Curate typed reference sets with provenance: McrA (methanogen / ANME clades /
  alkyl-CoM reductase) and Cu-monooxygenase subunit A (pmoA clades / amoA clades /
  hydrocarbon monooxygenases).
- Best-hit scorer, then tree placement (MAFFT, ClipKIT, IQ-TREE) for the GTDB calls,
  split by characterized genera vs candidate lineages; report which calls no sequence
  method resolves.
- Done when: `DIR_ACCURACY` and `DIR_PLACEMENT` tables exist for both traps.

### M8 — MAG realism
- Dossier of ≈ 8–10 published MAGs with a stated phenotype and CheckM quality: marine
  ANME-1 / -2a / -2c, a second ANME-2d, a non-euryarchaeal Mcr lineage, an alkane
  oxidizer, uncultured methanotrophs (e.g. atmospheric-methane clades, NC10), a
  methanogen MAG, two negatives.
- Run from nucleotide with `--prodigal-mode meta`; build the truth-vs-tool table with
  METABOLIC and DRAM alongside. Expect this phase to surface real bugs — in nitrogen it
  found one the isolate panel could not.
- Done when: `mag_truth_vs_tool.tsv` and a short write-up exist.

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
