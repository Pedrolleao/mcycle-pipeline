# mcycle-pipeline

Maps MAGs / isolate proteomes to their participation in the **methane cycle**:
methanogenesis (CO₂-reducing, acetoclastic, methylotrophic), anaerobic methane
oxidation by reverse methanogenesis, aerobic methane oxidation, and the fate of the
C1 units that follow. Third sister of `Nitrogen_Cycle/ncycle-pipeline` and
`Sulfur_Cycle/scycle-pipeline`, on the same engine: detection is **KO-primary**
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

What these numbers do not cover is in *Known limits* below — in particular McrA genes
truncated at contig ends in fragmented MAGs.

## Run

```bash
cd mcycle-pipeline
python run.py --input <dir-of-.faa-or-.fna> --cores 8
# First run creates the `cycle-pipeline` conda env from envs/mcycle.yaml if it does
# not exist yet (shared with the nitrogen and sulfur sister tools). To use another
# env with the same dependencies:
#   MCYCLE_ENV=<env-name> python run.py --input <dir> --skip-db-setup
```

`run.py` auto-detects protein (`.faa`) vs nucleotide (`.fna`, → Prodigal) input,
builds the databases on first run (and again whenever `config/targets.yaml` is newer
than them), then dispatches Snakemake.

```bash
make smoke          # fetch the smoke panel, run it, check the expected calls (~1 min)
```

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
     `disqualified` when the genome has an McrA homologue that fails the McrA call
     (alkyl-coenzyme M reductases of alkane oxidizers); `mmoY/Z/B/C/D` need `mmoX`.
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
   cross-sample `multisample_matrix.tsv`, figures, and an interactive `report.html`.

## Outputs

All paths are under `paths.results_dir` (`results/` by default). Figures are written
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
| `multisample_matrix.tsv` | genomes × (targets, complexes, modules) |
| `multisample_heatmap.svg/.png` | overview dot grid; genomes ordered by gene-content similarity |
| `figures/pathway_<pathway>.svg/.png` | one dot grid per pathway |
| `figures/complexes.svg/.png`, `figures/synergies.svg/.png` | complex / process-module completeness |
| `figures/mcycle_maps.svg/.png` | every genome's methane-cycle map side by side (up to 48 genomes) |
| `report.html` | self-contained interactive report (no network needed): the gene grid and the complex / module grid with hover evidence, row search / ordering, and a per-genome panel with the cycle map, locus maps and the full calls table. Light and dark themes. Every figure in it (gene grid, complex / module grid, cycle map, each locus map) has a **Save PNG (300 dpi)** button: it downloads that figure as currently shown — row filter and order, hidden pathways, selected genome, light or dark theme — with its title and legend, rendered at 300 dpi (a grid too large for a browser canvas is saved at the highest resolution that fits, and says so). The page follows the group's *Simple Terminal* design system (`design/Simple`): JetBrains Mono, hairline `[ bracketed ]` frames, its dark palette or its Light variant according to the system theme, with a LIGHT / DARK selector in the top-right corner to pin either. The font is inlined from `workflow/scripts/fonts/` (SIL OFL 1.1, licence alongside), so the report looks the same offline and the PNG export uses it too; pathway colours stay the validated palette of the static figures. |

**Reading the glyphs** (same in every figure and in `report.html`): solid disc =
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

Seven process modules / 83 markers (see [`../Info-methane.md`](../Info-methane.md) for
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

- **Truncated McrA in fragmented MAGs.** An McrA gene cut at a contig end scores under
  the K00399 threshold; the genome then gets no `mcrA` call, and — worse — its `mcrB` /
  `mcrG` are `disqualified` as alkyl-CoM reductase subunits, because the rule that
  recognizes alkane oxidizers sees "an McrA homologue that failed". 4 of the 67 genomes of Mcr-carrying
  lineages in the GTDB-500 set. Not fixed: the tool was frozen for the
  validation. First item of `ROADMAP.md`.
- **Non-euryarchaeal McrA** can fall under the KOfam threshold: the full-length McrA
  of one *Ca.* Methanomethylicus genome of the GTDB-500 set scores 578 (threshold
  775.5), while the *Ca.* M. mesodigestus V2 MAG of the MAG study is called.
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
python validation/check_smoke.py             # re-check existing results/ against the expectations
```

The KOfam cache (`resources/.cache/`) is symlinked to the ncycle-pipeline download
(~1.5 GB profiles.tar.gz, pinned release 2026-05-24, verified by SHA-256 at build
time) to avoid a re-fetch; delete the links to make the build download its own copy.
