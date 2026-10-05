# mcycle-pipeline

Maps MAGs / isolate proteomes to their participation in the **methane cycle**:
methanogenesis (CO₂-reducing, acetoclastic, methylotrophic), anaerobic methane
oxidation by reverse methanogenesis, aerobic methane oxidation, and the fate of the
C1 units that follow. Third sister of the nitrogen (`ncycle-pipeline`) and sulfur
(`scycle-pipeline`) tools, on the same engine: detection is **KO-primary**
(KOfam HMMs + adaptive per-KO thresholds), with **clade HMMs** for the two homology
traps (pmoA vs amoA; ANME vs methanogen McrA) and curated seeds as corroboration.
Covers **83 targets** (across 7 process modules), **17 obligatory complexes**, and
**11 process-completeness synergies**.

**Status: validated on a 49-genome reference panel (2026-10-04).** The validation
campaign of the sister tools has been run for methane — frozen train / hold-out split,
dual ground truth, regression gate, pre-registered comparator benchmark, GTDB-500
concordance, MAG study. Headline numbers (95 % genome-cluster bootstrap CIs; full
report in [`validation/REPORT.md`](validation/REPORT.md)):

| | |
|---|---|
| hold-out micro-F1 (22 genomes of genera never used for seeds, models or thresholds) | **0.961 [0.940, 0.976]** |
| full-panel micro-F1 (49 genomes, 3,173 cells) | 0.973 [0.963, 0.981] |
| homology-trap precision | 1.000 — 71 true calls, 0 false in 259 trap negatives |
| vs raw KofamScan / METABOLIC / DRAM / MCycDB, trap precision | +0.34 / +0.31 / +0.47 / +0.37, all BH q < 0.0001 |
| Mcr direction on 500 GTDB genomes, where McrA is called | 53 / 53 agree with the genome's lineage |

These numbers are those of the tool frozen for the validation (commit `9c39500`); one
minimal fix was applied afterwards and changes none of them (`validation/REPORT.md`,
section 6). What they do not cover is in *Known limits* below.

## Install

Developed and validated on Linux (x86-64). Needs `git`, conda or mamba, `curl` and
`unzip`; about 3 GB of disk for the environment. Network access is needed for the install,
for the first database build (28 seed sequences from UniProt) and for fetching the
test genomes from NCBI; the pipeline itself runs offline.

```bash
git clone https://github.com/Pedrolleao/mcycle-pipeline.git
cd mcycle-pipeline
conda env create -f envs/mcycle.yaml      # env `cycle-pipeline`; or: mamba env create …
conda activate cycle-pipeline
make smoke          # fetch 18 genomes (59 MB), build the databases, run, check 168 expected calls
```

`make smoke` ending in `OK: 168 smoke expectations met on 18 genomes` means the
install reproduces the reference calls.
`envs/mcycle.lock.yaml` is the exact environment of the validation (linux-64), for
when the open version ranges of `envs/mcycle.yaml` resolve to something that behaves
differently.

## Run

```bash
python mcycle.py --input <dir-of-.faa-or-.fna> --cores 8                 # results in mcycle_results/
python mcycle.py --input <dir> --output <results-dir> --cores 8          # results where you want them
# Outside the conda env, mcycle.py re-runs itself inside `cycle-pipeline` (and creates it
# from envs/mcycle.yaml if it does not exist). To use another env with the same
# dependencies:
#   MCYCLE_ENV=<env-name> python mcycle.py --input <dir>
```

`mcycle.py` auto-detects protein (`.faa`) vs nucleotide (`.fna`, → Prodigal) input,
builds the databases on first run (and again whenever `config/targets.yaml` is newer
than them), then dispatches Snakemake. It asks whether nucleotide input is isolate
genomes or metagenome assemblies unless `--prodigal-mode single|meta` is given, and
writes the samples it found into the `samples:` block of `config/config.yaml` — so
`git status` shows that file as modified after a run.

Results go to `mcycle_results/` inside this directory unless `--output DIR` is given (a
path relative to where you run the command). When the run finishes the launcher prints
the results directory and the path of `mcycle_report.html`, the page to open first.

## How it works

1. **Gene calls** — proteomes used directly; nucleotide assemblies → Prodigal.
2. **HMM scan** — `hmmscan` against `resources/hmm/mcycle_targets.hmm`: the **KOfam
   profile HMM for every KO in `config/targets.yaml`** (90 profiles) plus **nine clade
   HMMs** (`targets/pmoA/`, `targets/mcrA_anme/`). Thresholds (KOfam `ko_list`,
   per-target `ko_tc` overrides with their `tc_rationale`, clade-model thresholds from
   the manifests) live in `resources/hmm/tc_cutoffs.tsv`.
3. **BLAST** — `diamond blastp` against **curated UniProt seeds**
   (`resources/blast_db/`), tagged `>{target_id}||{acc}`; reported as the nearest
   reference, and still the corroboration for `mcrA` and `mmoX`.
4. **Calls** — `apply_rules.py` integrates evidence per target (signature precedence
   clade HMM > KO > Pfam). What is methane-specific:
   - **`pmoA`** is decided by four clade HMMs (gamma, alpha, Verrucomicrobia, NC10),
     each with a threshold set between the best-scoring ammonia / hydrocarbon
     monooxygenase and the worst left-out methanotroph genus. A copper monooxygenase
     subunit A that passes none is `disqualified` — KOfam alone calls the amoA of
     *Nitrosococcus* a methane monooxygenase. **`pmoB` / `pmoC`** share their KOs
     with amoB / amoC and follow the pmoA call.
   - **`mcrA_anme`** is McrA evaluated against five ANME clade HMMs (ANME-1, -2a/b,
     -2c, -2d, -3). If one passes, the genome's Mcr direction is `reverse`, otherwise
     `methanogenic`; the call is tagged on the mcrA / mcrB / mcrG `evidence_source`
     (e.g. `ko|mcr_reverse`) and names the clade model.
   - **Subunits follow the subunit that defines the enzyme.** `mcrB` / `mcrG` are
     `disqualified` when the genome has no McrA call and carries a full McrA homologue
     under the threshold (alkyl-coenzyme M reductases of alkane oxidizers; a gene
     fragment of a canonical McrA does not count); `mmoY/Z/B/C/D` need `mmoX`.
   - **`fdhA` vs `fdh`.** With the F420-binding `fdhB` in the genome, an FdhA-family
     protein is the F420-dependent `fdhA` and not the NAD-linked `fdh`.
   - **`mcrA`** and **`mmoX`** are corroborated by seeds without being gated: a KO
     hit with no seed above the identity floor is reported as `domain-only`.
   - **Genes split at an in-frame stop** *(nucleotide input)*. Methanogens read UGA
     as selenocysteine and UAG as pyrrolysine; Prodigal stops there, so such a gene
     arrives as two neighbouring ORFs that each miss the KOfam threshold. Two ORFs
     are scored as one when they are neighbours on one strand, the second continues
     in the reading frame of the first across at most 150 nt, they hit consecutive,
     non-overlapping parts of the same profile, and their scores together reach the
     threshold. The call is tagged `joined:<other ORF>` in `evidence_source` —
     sequence alone cannot tell recoding from a nonsense mutation.
5. **Reports** — per-sample `calls/mcycle_calls.tsv`, `complex_completeness.tsv`,
   `synergy_completeness.tsv`, `report/gap_analysis.txt`, `report/mcycle_map.*`;
   cross-sample `mcycle_matrix.tsv`, figures, and an interactive `mcycle_report.html`.

## Outputs

All paths are under the results directory (`mcycle_results/` by default, or `--output DIR`). Figures are written
as SVG (vector) and PNG (300 DPI).

**Per sample — `<sample>/`**

| file | what it is |
|---|---|
| `calls/mcycle_calls.tsv` | one row per target: status, evidence source, protein, HMM hit + E-value, BLAST reference + identity. For nucleotide input the columns `contig`, `start`, `end`, `strand` locate the called gene (empty for a pre-called proteome). One protein is reported per target; `other_copies` lists further proteins that reach the same status |
| `calls/complex_completeness.tsv`, `calls/synergy_completeness.tsv` | completeness of the 17 complexes and 11 process modules |
| `calls/mcycle_loci.tsv` | *(nucleotide input only)* called genes grouped into loci: genes on one contig with at most 5 other genes between them |
| `report/gap_analysis.txt` | plain-text summary |
| `report/mcycle_map.svg/.png` | the genome's calls drawn on the methane cycle (see below) |
| `report/loci.svg/.png` | *(nucleotide input only)* gene-arrow maps of every locus, grouped by pathway, to a common bp scale |

**Across samples**

| file | what it is |
|---|---|
| `mcycle_matrix.tsv` | genomes × (targets, complexes, modules) |
| `mcycle_heatmap.svg/.png` | overview dot grid; genomes ordered by gene-content similarity |
| `figures/pathway_<pathway>.svg/.png` | one dot grid per pathway |
| `figures/complexes.svg/.png`, `figures/synergies.svg/.png` | complex / process-module completeness |
| `figures/mcycle_maps.svg/.png` | every genome's methane-cycle map side by side (up to 48 genomes) |
| `mcycle_report.html` | self-contained interactive report (no network needed): the gene grid and the complex / module grid with hover evidence, row search / ordering, and a per-genome panel with the cycle map, locus maps and the full calls table. Light and dark themes. Every figure in it (gene grid, complex / module grid, cycle map, each locus map) has a **Save PNG (300 dpi)** button: it downloads that figure as currently shown — row filter and order, hidden pathways, selected genome, light or dark theme — with its title and legend, rendered at 300 dpi (a grid too large for a browser canvas is saved at the highest resolution that fits, and says so). The page follows the group's *Simple Terminal* design system (`design/Simple`): JetBrains Mono, hairline `[ bracketed ]` frames, its dark palette or its Light variant according to the system theme, with a LIGHT / DARK selector in the top-right corner to pin either. The font is inlined from `workflow/scripts/fonts/` (SIL OFL 1.1, licence alongside), so the report looks the same offline and the PNG export uses it too; pathway colours stay the validated palette of the static figures. |

**Reading the glyphs** (same in every figure and in `mcycle_report.html`): solid disc =
confirmed; half-filled = domain-only (HMM signature, no BLAST support); ring with a
cross = disqualified (failed the homology-trap gate); faint ring = absent. In the
complex / module grids: solid = complete, ring with `n/N` = partial, faint ring with
a cross = ruled out. Colour always means pathway; the palette lives in
`workflow/scripts/_domain.py`.

**The cycle map.** A loop around a ladder of H₄MPT-bound C1 intermediates.
Methanogenesis runs down the right-hand side (CO₂ → CH₄), with acetate and methylated
compounds entering from the right; aerobic methane oxidation climbs the left
(CH₄ → CH₃OH → HCHO → formate → CO₂, or → biomass). Two things depend on the genome.
In a `reverse` (ANME) genome the methanogenesis chain is drawn from CH₄ to CO₂. And
the shared `fwdABC` / `ftr` / `mch` genes are drawn where they act: on the
CO₂ ⇄ formyl-methanofuran steps in an archaeon, on the formyl-H₄MPT → formate step
(as the Fhc complex) in a bacterium that has `fae` and no Mcr.

**Modules, direction and key genes.** A module in `synergy_completeness.tsv` is a
list of slots: `requires` genes, plus `requires_any` choices written `a+b|c`
("a and b, or c"). Two rules can rule a module out — status `absent`, the reason
under `forbids_violated`, "ruled out" in the figures:

- **direction** — the three methanogenesis modules need Mcr `methanogenic`,
  `reverse_methanogenesis_aom` needs `reverse`;
- **key genes** — what defines the module (`nod` for nitrite-dependent methane
  oxidation, a substrate methyltransferase for methylotrophic methanogenesis, `hdrD`
  for acetoclastic, a monooxygenase for aerobic methane oxidation …). Without them,
  the genes a module shares with its neighbours do not make it "partial".

The figure and report scripts are shared, byte-identical, with the sister pipelines;
only `workflow/scripts/_domain.py` (palette, labels, file names) and
`workflow/scripts/_cycle_model.py` (the cycle diagram) are methane-specific, plus the
rules in `apply_rules.py` and `compute_complex_completeness.py` described above and
the ruled-out message in `make_gap_analysis.py`.

## Pathways & targets

Seven process modules / 83 markers (see [`Info-methane.md`](Info-methane.md) for
the per-gene atlas, KO anchors, measured gate calibration and trap rationale):

1. **Mcr / Mtr core** — mcrABG(CD), mcrA_anme, mtrA–H
2. **CO₂ reduction (H₄MPT C1 module)** — fwdABCD, ftr, mch, mtd, hmd, mer, fdhAB
3. **Acetoclastic** — ackA, pta, acs, cdhABCDE
4. **Methylotrophic methanogenesis** — mtaABC, mtmBC, mtbBC, mttBC, mtbA, mtsAB
5. **Hdr / hydrogenases** — hdrABC, hdrDE, mvhAGD, frhABG
6. **Aerobic methane oxidation** — pmoABC, mmoXYZBCD, mxaFI, xoxF, nod
7. **Formaldehyde oxidation / C1 assimilation** — fae, mtdB, gfa, frmAB, fdh, hxlAB,
   sgaA, hprA, gckA, mtkAB, mcl

## Known limits

Found or confirmed by the validation campaign (`validation/REPORT.md`):

- **Truncated McrA in fragmented MAGs is not called.** An McrA gene cut at a contig
  end scores under the K00399 threshold, so the genome gets no `mcrA` call and no
  direction (4 of the 67 genomes of Mcr-carrying lineages in the GTDB-500 set). Its
  `mcrB` / `mcrG` are reported normally since the amendment of 2026-10-05 (before it
  they were mislabelled as alkyl-CoM reductase subunits; `validation/CHANGELOG.md`).
  Calling the fragment itself needs a "direction unresolved" state and is left for a
  separate cycle.
- **Non-euryarchaeal McrA** can fall under the KOfam threshold: the full-length McrA
  of one *Ca.* Methanomethylicus genome of the GTDB-500 set scores 578 (threshold
  775.5) and its `mcrB` / `mcrG` are then labelled as alkyl-CoM reductase subunits,
  while the *Ca.* M. mesodigestus V2 MAG of the MAG study is called. Separate cycle.
- **A MAG that lacks the gene cannot be called**: 6 GTDB genomes and 2 of the 12
  study MAGs have no Mcr subunit A in the assembly.
- **ANME-3 is not cleanly separable** from methylotrophic Methanosarcinaceae on McrA
  sequence; the model is set for precision, and two of five left-out ANME-3 species
  fell under its threshold. The direction call is about the clade, not about what the
  cell is doing on the day.
- **Divergent Cu-monooxygenase paralogues are not called**: pxmA (gammaproteobacteria,
  *Methylocystis*) and verrucomicrobial pmoA3; a genome carrying only such a copy has
  `pmoA` `disqualified`.
- **Not tested out-of-genus** (every genus is a seed or training genus): ANME-2d,
  NC10, alpha-proteobacterial pMMO, acetoclastic methanogenesis.
- **Residual gene-level errors** (47 of 3,173 cells): `acs` against its paralogues,
  `fdh` / `fdhA` adjudication and selenocysteine formate dehydrogenases whose halves
  do not add up to the threshold, `frhB` (F420-binding subunits of other complexes),
  `hdrD`, the broad `sgaA` / `hprA`, corrinoid proteins `mttC` / `mtbC`.
- **Genomic potential, not physiology.** *Methanosarcina acetivorans* scores complete
  for hydrogenotrophic methanogenesis (it carries `frh` but does not grow on H₂ + CO₂).
- **Recoded genes** (selenocysteine, pyrrolysine). From genome input they are called
  when one half passes on its own or the two halves can be joined (above). From a
  proteome they depend on the annotation.
- **Ubiquitous genes** (`ackA`, `pta`, `acs`, `fdh`, `frmAB`, `hdrA`) light up in most
  genomes; read them through the modules, which need Mcr or another key gene.

## Validation

| what | where |
|---|---|
| panel roster, split, rationale (frozen 2026-10-04) | `validation/panel.tsv`, `validation/PANEL_PLAN.md` |
| ground truths (KEGG; curated function, with evidence per cell) | `validation/ground_truth.tsv`, `curated_function_gt.tsv`, `curated_cells.tsv`, `phenotype_gt.tsv` |
| hardening log (training genomes only) | `validation/HARDENING.md` |
| scoring, seed-leak and independence audits, gate | `validation/score_mcycle.py`, `detect_seed_leakage.py`, `trap_independence.py`, `test_regression.py` |
| comparator benchmark (pre-registered) | `validation/benchmark/` |
| GTDB-500 concordance and clade checks | `comparators/gtdb500_m/` |
| MAG study | `validation/metagenomes/` |
| every change after a freeze | `validation/CHANGELOG.md` |

```bash
make ref-panel validate-panel     # fetch the 49 genomes, check each is what the roster says
make regression                   # run the panel, score, 8-row gate (14 checks)
make regression-score             # re-score existing results_ref/
MCYCLE_SCOPE=train python validation/score_mcycle.py    # one side of the split
GT_FILE=ground_truth.tsv python validation/score_mcycle.py   # KEGG contrast
```

## Build / test

```bash
python workflow/scripts/build_hmm_db.py     # KOfam KO profiles (+ custom HMMs / Pfam fallbacks if any)
python workflow/scripts/build_blast_db.py    # curated UniProt seeds → DIAMOND DBs
python validation/check_smoke.py             # re-check existing mcycle_results/ against the expectations
```

**KOfam is pinned in the repository.** The 90 KOfam profiles and their thresholds are
built from `resources/kofam_pinned/`, the release of 2026-05-24 the tool was validated
on. genome.jp serves KOfam from a rolling URL and keeps no old releases — the release
of 2026-09-29 has a different threshold for 85 of the 90 KOs — so a database built
from a fresh download is not the validated one. `build_hmm_db.py --upstream` downloads
the current release (~1.5 GB) anyway; re-validate (`make regression`) before trusting it.

**What a clone does not contain.** Pipeline results, the downloaded genomes, and the
third-party tools and databases of the comparator benchmark (METABOLIC, DRAM, MCycDB,
GTDB proteomes). `make ref-panel` downloads the 49 reference genomes into `../ref_panel`,
next to the clone. The scripts under `comparators/`, the study configs
`config/config_{ref,gtdb500,mags}.yaml` and the clade-model builders
(`harvest_clade_refs.py`, `build_clade_hmms.py`) are the record of how the validation
was run: they carry paths of the machine it ran on and need editing to be re-run
elsewhere. Their outputs — the tables under `validation/` and `comparators/`, and the
clade HMMs under `targets/` — are committed.
