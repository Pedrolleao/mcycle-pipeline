# WORKPLAN — fixes and improvements from the code review of 2026-11-07

**Living document.** Written 2026-11-07 from the review of commit `97fd402`. Nothing here
has been started. The validation campaign (`WORKPLAN.md`) is closed; this plan is a
series of small, separately committed changes that leave every call untouched, plus two
items at the end that do change what the tool reports and therefore wait for their own
validated cycle.

## Principle

The tool is frozen and validated (`validation/REPORT.md`, amendment of 2026-10-05). So:

1. **Output-identical first.** Phases R1–R6 may not change one cell of any
   `mcycle_calls.tsv`. Each phase ends with `make smoke` (168 / 168) and a byte-for-byte
   diff of the 18 call tables against the baseline taken in R0.
2. **One phase, one commit**, with the diff result in the commit message. A phase that
   fails its acceptance test is reverted, not patched forward.
3. **Shared files are ported.** `mcycle.py`, `cycle_ui.py`, the figure and report scripts
   and `report_template.html` are byte-identical in `ncycle-pipeline` and
   `scycle-pipeline`. A change to one of them is finished only when the same patch is
   applied there (the sister tools' launchers read `_domain.py`, so most patches apply
   verbatim). Push only after asking.
4. **Anything that changes a reported value** (R7 fraction convention, R9, R10) is
   written up before it is run, with before / after tables, and gets a dated entry in
   `validation/CHANGELOG.md`.

## Status

| Phase | What | Changes calls? | State |
|---|---|---|---|
| R0 | Baseline: snapshot of the current outputs and timings | no | DONE 2026-11-07 |
| R1 | Launcher refuses what the workflow cannot run (names, duplicates, malformed FASTA) | no | DONE 2026-11-07 |
| R2 | Error panel shows the cause and the right log path | no | DONE 2026-11-07 |
| R3 | hmmscan on one thread per job | no | DONE 2026-11-07 |
| R4 | Every figure rendered once, saved twice | no | DONE 2026-11-07 |
| R5 | Per-run config instead of rewriting `config/config.yaml` | no | DONE 2026-11-07 |
| R6 | Text outputs: gap_analysis wording and relevance, header, dev note | no | DONE 2026-11-07 |
| R7 | Module fraction: one convention in TSV, figures and report | **also 7 module statuses** (measured) | MEASURED 2026-11-07, NOT APPLIED — user's decision |
| R8 | Large runs: cap the static grids, sticky header in the HTML report, column-name note | no | DONE 2026-11-07 |
| R9 | hmmsearch instead of hmmscan | needs the gate | DESIGN ONLY |
| R10 | Port to ncycle / scycle, CHANGELOG, README, push | no | TODO |

Order: R0 → R1 → R2 → R3 → R4 are independent of each other after R0 and can be done in
any order; R5 touches the same launcher code as R1 / R2 and goes after them; R6 and R7
after R5 (they rebuild outputs); R8 last of the output-identical set; R10 closes.

## Acceptance test used by every phase

```bash
conda activate cycle-pipeline
make smoke                                   # ends in: OK: 168 smoke expectations met on 18 genomes
for s in mcycle_results/*/calls; do
  d=$(basename $(dirname $s))
  cmp -s $s/mcycle_calls.tsv validation/baseline_2026-11/$d/mcycle_calls.tsv || echo "DIFF $d"
done                                         # prints nothing
python validation/check_report_consistency.py mcycle_results   # exit 0
```

The baseline directory is created in R0 and is **not committed** (add it to `.gitignore`);
its checksums are.

## Phases

### R0 — Baseline (30 min)

- Copy the 18 current `mcycle_results/<s>/calls/` directories to
  `validation/baseline_2026-11/<s>/` and write `validation/baseline_2026-11.sha256`
  (committed) with the checksums of every `mcycle_calls.tsv`,
  `complex_completeness.tsv`, `synergy_completeness.tsv`.
- Also keep the 18 `hmm/<s>.hmmscan.tsv` files there: R3 is checked at that level.
- Record the wall time of a clean `make smoke` (`rm -rf mcycle_results/*` first) with
  `--cores 8`, and the per-rule profile of `.snakemake/log/` (the parser used in the
  review: timestamps of `jobid:` vs `Finished jobid:`). These are the before numbers for
  R3 and R4.
- Done when: the snapshot exists, `check_smoke.py` passes on it, the timings are in this
  file under "Measurements".

### R1 — Launcher refuses what the workflow cannot run (2 h)

Found: a file name with a space is accepted and the run dies at DAG build
(`common.smk` allows `[A-Za-z0-9_.-]+` only); `a.faa` + `a.fna` write two `a:` keys and
only the last runs; a FASTA with illegal characters makes `diamond blastp` hang forever;
`bootstrap_env_and_reexec` accepts `mamba` but exec's the literal `conda`.

Changes, all in `mcycle.py`, stage "Input":
- `discover_samples`: validate every stem against the same regex as `common.smk`
  (define it once in the launcher; a comment in `common.smk` points to it). On a
  violation exit with the file name and the rule ("letters, digits, `_ . -`"). Do not
  rename silently: the sample name is the row label in every output.
- Duplicate stems: exit with both paths. (An option to suffix them by extension can come
  later; refusing is the safe default on a cluster.)
- Sequence sanity: while `fasta_is_nucleotide` already reads the first 3,000 residues,
  reject a file whose first record has no `>` header or whose sampled residues contain
  characters outside `A–Z a–z * - .` — that is what hung diamond. Message names the
  file and the first bad character.
- `bootstrap_env_and_reexec`: exec the binary that was found.
- `--mode` is dead (only `protein`): remove it or mark it hidden; `modes_present` in
  `build_databases` is unused, drop it.
- Tests (keep as `validation/launcher_cases/` with a 10-line script, no pipeline run
  needed, `--dry-run --plain`): space in name → exit 1 with the name; `a.faa` + `a.fna`
  → exit 1 listing both; `not a fasta` → exit 1; a `?` in a sequence → exit 1; the smoke
  panel → dry run plans 18 samples as before.
- README "Run": one sentence on the naming rule.

### R2 — Error panel shows the cause (1–2 h)

Found: `_RE_LOG` expects `Complete log:` but Snakemake 9 prints `Complete log(s):`, so
the hint is always `.snakemake/log/`; the panel body is the last 2,500 characters of the
error block, i.e. the shell command and `message: None`, while the cause is in the job's
log file.

Changes in `mcycle.py`, `run_snakemake_live`:
- `_RE_LOG = r"Complete log(?:\(s\)| file)?:\s*(\S+)"`.
- From the error block, take the `log:` line (strip the "(check log file(s)…)" suffix)
  and, when the file exists and is non-empty, show its last 15 lines in the panel under
  a "job log" heading, above the Snakemake block. For rules without `log:` (the report
  rules) the traceback arrives on the pipe before "Error in rule"; keep the last 25 pipe
  lines before the block and show them when no job log was found.
- Head line: prefer `Error in rule <name>` over `RuleException:` when both are present.
- Plain mode is unaffected (Snakemake's own output passes through).
- Test: the `not a fasta` case of R1 run with `--skip-validation` (add this hidden flag
  for exactly this purpose, or run the rule by hand) under
  `COLUMNS=110 script -qec "…" /dev/null`; the panel must contain
  `Error detecting input file format` and the real log path.

### R3 — hmmscan on one thread per job (30 min + a smoke run)

Measured on `Msoehngenii_GP6` (3,049 proteins, 99 profiles): `--cpu 4` 5.1 s,
`--cpu 1` 4.0 s, `domtblout` byte-identical; four `--cpu 1` jobs at once 3.3 s against
20.4 s for four `--cpu 4` jobs in a row. On the GTDB-500 run hmmscan was 3,492 of the
summed CPU-seconds and the step that set the wall time.

- `workflow/rules/protein_mode.smk`, rule `hmmscan`: `threads: 1`. Leave diamond at 4
  for now (0.1 s per job, it does not matter) or set it to 2 for symmetry; note the
  choice in the commit.
- Acceptance: all 18 `hmm/<s>.hmmscan.tsv` byte-identical to the R0 snapshot after
  stripping the `# Program / # Date / # Current dir` comment lines (`grep -v '^#'`),
  calls identical, `make smoke` wall time recorded under "Measurements".
- Port to ncycle / scycle (same rule file).

### R4 — Render once, save twice (2 h)

Found: `cross_sample_report`, `pathway_heatmap`, `complex_heatmap_figure`,
`synergy_heatmap_figure` call their script twice, once per extension; on the 500-genome
run that is ~160 CPU-seconds and the single-threaded tail of the DAG.
`make_cycle_map.py` and `make_locus_maps.py` already take `--out A.svg A.png`.

- Give `make_cross_sample_report.py` (`--heatmap-out`), `make_pathway_heatmap.py`,
  `make_complex_heatmap.py`, `make_synergy_heatmap.py` (`--out`) `nargs="+"` and a shared
  `_viz.save_all(fig, outs)` that saves each path by extension and closes the figure
  once. Update the four rules to one invocation each.
- Matplotlib writes a creation date into SVG and PNG metadata; pass
  `metadata={"Date": None}` (SVG) / `pil_kwargs` is not needed — compare pixels, not
  bytes: a 10-line script opening before / after PNGs with Pillow and asserting equal
  arrays, and SVGs equal after deleting the `<dc:date>` line.
- Acceptance: pixel-identical PNGs and SVGs for the smoke run, figure-stage time before
  / after under "Measurements".

### R5 — Per-run config (2–3 h)

Found: the launcher rewrites the tracked `config/config.yaml`, so `git describe --dirty`
reports `-dirty` after the first run, concurrent runs clobber each other, and `git status`
is permanently noisy (the README has to explain it).

- `mcycle.py`: read the base config (`MCYCLE_CONFIG` or `config/config.yaml`), add the
  `samples:` block and `results_dir`, write it to `<results_dir>/mcycle_run_config.yaml`
  (created before Snakemake starts), and pass that with `--configfile`. The base config
  is never written. The run config is also provenance: add `tool_version` and the
  command line as comments at its top.
- `workflow/Snakefile`: the `configfile:` directive stays as a default for people who run
  Snakemake by hand; document that `samples:` must then be filled in manually.
- Check `grep -rn "config_ref\|config_gtdb500\|config_mags\|samples:" validation
  comparators Makefile` for anything that reads the samples block back from the study
  configs; the committed study configs keep their blocks (they are the record of the
  runs) but nothing may depend on the launcher refreshing them.
- `tool_version()`: after this change `-dirty` is meaningful again; keep it.
- README: delete the two sentences about `config.yaml` showing as modified; add one on
  `mcycle_run_config.yaml`.
- Acceptance: `git status` clean after `make smoke`; two launcher runs with different
  inputs and outputs started at the same time both finish; calls identical.

### R6 — Text outputs (2 h)

Found in `gap_analysis.txt`: every absent complex is "MISSING — CRITICAL … blocks
application" (13 lines for *A. fulgidus*, including pMMO and the methanol
methyltransferase); the synergy header promises "application still possible but
degraded" above rows that say "ruled out"; `mode: protein` is printed for nucleotide
input; the `application` text of the `mcr` complex carries the developer note
"⚠️ direction set by resolve_mcr_direction" into `gap_analysis.txt` and
`complex_completeness.tsv`.

- `make_gap_analysis.py`:
  - rename the sections: "PRESENT", "COMPLEXES — partial" (missing subunits, as now),
    "COMPLEXES — absent" as a one-line list of ids, not one block each; "PROCESS
    MODULES": complete / partial (missing slots) / ruled out (reason). Keep the
    file name and the first two header lines so nothing that greps it breaks
    (`grep -rn gap_analysis validation comparators` first).
  - input kind: pass `--input-kind nucleotide|protein` from the rule (the Snakefile knows
    `SAMPLES[s]["kind"]`), print "input: nucleotide (Prodigal, single)" or "input:
    proteome".
- The developer note: strip it from `config/targets.yaml` (`complexes.mcr.application`,
  keep the information as a YAML comment). **Editing `targets.yaml` makes `_stale()`
  rebuild both databases on the next launcher run** (needs network for the UniProt
  seeds); the CHANGELOG of 2026-10-05 showed the rebuild is byte-identical, verify it
  again (`sha256sum resources/hmm/mcycle_targets.hmm resources/hmm/tc_cutoffs.tsv
  resources/blast_db/*.fasta` before and after). If the network is not available, strip
  the note at output time in `compute_complex_completeness.py` instead and leave the
  YAML for later.
- Acceptance: calls identical; `check_smoke.py` passes; the three example files
  (*A. fulgidus*, *M. nitroreducens*, *E. coli*) read in full in the commit message.

### R7 — One module-fraction convention (1 h, changes matrix values)

Found: `compute_complex_completeness.py` counts a violated key-gene rule as a slot, so
hydrogenotrophic methanogenesis in *A. fulgidus* is 6 / 13 = 0.462 in
`synergy_completeness.tsv` and `mcycle_matrix.tsv`, while `_completeness.load_synergies`
(figures, HTML) shows 6 / 12.

- Decide: the TSV follows the figures (n_total = required slots only; the direction and
  key-gene rules decide `status` and are listed under `forbids_*`, not counted). This
  only moves the `completeness` value of modules that are already `absent`; `status`
  does not change, so `check_smoke.py` (module status only) is unaffected.
- Check first that nothing in `validation/` scores the completeness number:
  `grep -rn "completeness\|synergy__" validation/*.py comparators/*.py`.
- Record the before / after values for the smoke panel in `validation/CHANGELOG.md`
  under a dated entry "no call changes; module completeness of ruled-out modules
  re-based".

**Measured 2026-11-07 (dry comparison, old vs new `compute_complex_completeness.py`
over the existing calls of all result sets: smoke 18, reference 49, GTDB-500 500,
MAGs 12).** The premise "status does not change" is wrong: with the direction slot
no longer counted, a module with an odd number of required slots can cross the 0.5
threshold. `acetoclastic_methanogenesis` (9 required slots + the Mcr direction) goes
from 5 / 10 = 0.500 `partial` to 4 / 9 = 0.444 `absent` in 7 genomes — reference
panel `Malvi_Mx1201`; GTDB-500 `GCA_001717005_1`, `GCA_023402275_1`,
`GCA_029887795_1`, `GCA_035404685_1`, `GCA_963586845_1`; MAG `Mmethylicus_V2`.
No status changes in the smoke panel; 1,474 completeness values change without a
status change (e.g. *A. fulgidus* hydrogenotrophic 6 / 13 = 0.462 → 6 / 12 = 0.500).
Because 7 reported module statuses would move, this is a change for its own validated
cycle (amendment with before / after tables in `validation/CHANGELOG.md`), not a
re-basing of values — **not applied**. Alternatives to decide between: (a) apply it
as a disclosed amendment; (b) keep the TSV as is and make the figures show the TSV's
n / N (then the figures would read 5 / 10 for a 9-slot module); (c) leave both as they
are and document the two conventions in the README. The figures and the HTML report
keep reading the required slots from the `requires_present` / `requires_missing`
columns either way.

### R8 — Large runs (2–3 h)

Found: at 500 genomes `mcycle_heatmap.svg` is 33 MB, its PNG 53 MB, `complexes.png`
21 MB; the HTML report handles 500 rows in 80 ms, but its gene grid is 11,000 px tall
and the column labels scroll away; `pfam_hits` / `best_pfam_evalue` hold KO numbers and
clade-model names.

- Static grids: above `--max-rows` genomes (default 120) `make_cross_sample_report.py`,
  `make_pathway_heatmap.py`, `make_complex_heatmap.py`, `make_synergy_heatmap.py` write
  the pointer figure that `make_cycle_map.py` already writes above 48 panels ("N genomes
  — open mcycle_report.html"). The matrix TSV is always written.
- `report_template.html`: make the gene-name header of the gene grid and of the complex
  grid sticky while scrolling the grid (a second small SVG with the labels, `position:
  sticky; top: 0`, same column x positions). Check with
  `node ~/tools/shotter/shot.mjs results_gtdb500/mcycle_report.html --theme both --full`
  and a `click-shot.mjs --eval` that scrolls 3,000 px and screenshots.
- README "Outputs": one sentence that `pfam_hits` names the signature that fired (KO,
  clade model or Pfam) and `best_pfam_evalue` its E-value; renaming the columns is a
  sister-tool-wide change and is not done here.
- Acceptance: calls identical; the smoke figures unchanged (18 < 120); the GTDB-500
  figures regenerated by hand with `--max-rows 120` to see the pointer.

### R9 — hmmsearch instead of hmmscan (design only; own cycle)

Measured: `hmmsearch --cpu 4` 1.07 s against 5.07 s for hmmscan on the same proteome;
with `-Z 99 --domZ 99` the full-sequence scores and E-values of every shared row are
identical, but hmmscan reports 17 more weak secondary domain rows (182 shared of 199)
because it sets domZ per query. `apply_rules.py` reads domain ranges for
`join_inframe_stops` (profile ranges of the two halves) and `gate_mcr_subunits`
(covered profile length), and the per-domain `hmm_min_qcov` filter, so the calls can
differ.

If taken up: write the design (which Z / domZ, whether to drop the per-domain
dependence by taking env coverage over the union of domains), run the dry-run
comparison over every existing result set as in M10 (`validation/amendments/`), then
`make regression` and the GTDB-500 clade checks; CHANGELOG entry with the cell-level
diff. Expected gain after R3 is smaller than it looks (R3 already removes most of the
wall time), so decide after the R3 measurement.

### R10 — Close out (1–2 h)

- Port the shared-file patches (R1, R2, R3, R4, R5, R8) to `ncycle-pipeline` and
  `scycle-pipeline`; run their smoke targets.
- `validation/CHANGELOG.md`: one dated entry listing R1–R8 with "no call changes" and
  the R7 re-basing.
- README: the sentences named above; `ROADMAP.md`: link this plan under "Done after the
  campaign".
- Ask before pushing the three repositories.

## Measurements (filled in as phases close)

| What | Before (R0) | After | Phase |
|---|---|---|---|
| `make smoke` wall time, 8 cores (clean, incl. check_smoke) | 67 s | 30 s | R3 |
| hmmscan summed job seconds, smoke panel (18 jobs, `--cpu 4`) | 101 s (5.6 s / job) | 93 s with `--cpu 1`, 8 jobs at once (5.2 s / job) | R3 |
| figure stage (cross_sample_report + 7 pathway + complexes + synergies), job seconds | 12 s (2 + 7 + 2 + 1); the 10 scripts timed by hand: 14.2 s as two calls each | 8.3 s as one call with two outputs (mean of 3) | R4 |
| summed job seconds, all 175 jobs | 179 s | | R3 / R4 |
| `git status` after a run | modified: config/config.yaml | clean | R5 |

R8 notes (2026-11-07): GTDB-500 figures regenerated by hand with the default cap: the
overview, complex and module grids and the pathway grids with more than 120 active rows
(acetoclastic: 149) become the pointer figure (28 KB SVG instead of 33 MB); `mcr_core`
(99 active rows) is still drawn. The 500-genome HTML report rebuilt with the new template:
scrolled 3,000 px, the header group sits at the viewport top (`translate(0 2548)`), the
export clone carries no transform, no console errors; the smoke report shot in both themes
is clean.

R6 notes (2026-11-07): editing `targets.yaml` rebuilt both databases (network was up).
`mcycle_targets.hmm` and `tc_cutoffs.tsv` came back byte-identical; the two seed FASTAs
did not (`unstable_refs.fasta` f518807e… → 29f733f9…, `blast_gated_refs.fasta`
c0198350… → b46ba66b…), but the difference is UniProt's header text only: the 36 diamond
tables of the smoke panel made with the validated database (results of 2026-10-05) are
byte-identical to those made with the rebuilt one, and so are all 18 call tables. The
build script keeps the full UniProt description line after `>target||acc`, so the file
checksum tracks UniProt's annotation releases, not the sequences; the sister tools pin
their seeds under `resources/seeds_pinned/` with a per-sequence manifest — worth copying
here in a later phase (see "Possible follow-ups").

R0 notes (2026-11-07): the snapshot in `validation/baseline_2026-11/` (18 × calls/*.tsv,
hmm/<s>.hmmscan.tsv, plus figures/, mcycle_heatmap.*, mcycle_matrix.tsv) comes from a
clean `make smoke` of commit `97fd402`; its 54 call tables are byte-identical to the
`mcycle_results/` of 2026-10-05 (the validation state). Per-rule job seconds come from
the `jobid:` / `Finished jobid:` timestamps of `.snakemake/log/2026-11-07T171946…log`
(1 s resolution, so sub-second rules read 0).

## Possible follow-ups (not in this plan)

- Pin the UniProt seed FASTAs as the sister tools do (`resources/seeds_pinned/` plus a
  per-sequence sha256 manifest), so a `targets.yaml` edit no longer needs the network and
  the seeds are the validated ones by construction.

## Decisions for the user

1. R1: refuse duplicate stems, or suffix them (`a_faa`, `a_fna`)?
2. R6: edit `targets.yaml` (triggers a database rebuild, network needed) or strip the
   note at output time?
3. R7: go ahead with re-basing the completeness value of ruled-out modules?
4. R8: 120 as the cap for the static grids?
5. R9: wanted at all after the R3 numbers are in?
