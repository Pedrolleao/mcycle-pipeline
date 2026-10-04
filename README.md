# mcycle-pipeline

Maps MAGs / isolate proteomes to their participation in the **methane cycle**:
methanogenesis (CO₂-reducing, acetoclastic, methylotrophic), anaerobic methane
oxidation by reverse methanogenesis, aerobic methane oxidation, and the fate of the
C1 units that follow. Third sister of `Nitrogen_Cycle/ncycle-pipeline` and
`Sulfur_Cycle/scycle-pipeline`, on the same engine: detection is **KO-primary**
(KOfam HMMs + adaptive per-KO thresholds) with DIAMOND-BLAST gating for the homology
traps. Covers **83 targets** (across 7 process modules), **17 obligatory complexes**,
and **11 process-completeness synergies**.

**Status: reference layer built, runs end-to-end, smoke-tested.** The engine is the
sulfur pipeline's; the methane biology is authored fresh (`config/targets.yaml`,
[`../Info-methane.md`](../Info-methane.md)). On the 18-genome smoke panel — run from
genome sequence, so every genome has gene coordinates and locus maps — all 168
expectations hold (`make smoke`): each of five methanogens is complete for its own
methanogenesis type(s), the ANME-2d genome is called *reverse*, four aerobic
methanotrophs are called by the right monooxygenase, a methylotroph without one is
not, and three ammonia oxidizers have their amoABC **disqualified** as pmoABC.
**This is not an accuracy estimate** — the panel is small and the two BLAST gates were
calibrated on it. The validation campaign the sister tools went through (reference
panel, ground truth, hold-out, comparators) has not been run; see
[`ROADMAP.md`](ROADMAP.md).

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
2. **HMM scan** — `hmmscan` against `resources/hmm/mcycle_targets.hmm`, the **KOfam
   profile HMM for every KO in `config/targets.yaml`** (90 profiles). Thresholds
   (KOfam `ko_list`, per-target `ko_tc` overrides) live in
   `resources/hmm/tc_cutoffs.tsv`. The custom-HMM tier of the sister tools is wired
   in but no methane clade HMM has been trained yet.
3. **BLAST gating** — `diamond blastp` against **curated UniProt seeds**
   (`resources/blast_db/`), per-target `blast_identity_min`, tagged `>{target_id}||{acc}`.
4. **Calls** — `apply_rules.py` integrates evidence per target (signature precedence
   custom HMM > KO > Pfam). Three things are methane-specific:
   - **`pmoA`** needs a hit ≥ 70 % to a methanotroph pmoA seed, otherwise it is
     `disqualified` — KOfam alone calls the amoA of *Nitrosococcus* a methane
     monooxygenase. **`pmoB` / `pmoC`** share their KOs with amoB / amoC and follow
     the pmoA call.
   - **`mcrA_anme`** is McrA evaluated against ANME-clade seeds (≥ 80 %). If it
     passes, the genome's Mcr direction is `reverse`, otherwise `methanogenic`; the
     call is tagged on the mcrA / mcrB / mcrG `evidence_source`
     (e.g. `ko|mcr_reverse`).
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

- **Marine ANME are called methanogenic.** Only ANME-2d (*Ca.* Methanoperedens) has
  a seed; ANME-1 / -2a / -2c / -3 do not.
- **Genomic potential, not physiology.** *Methanosarcina acetivorans* scores complete
  for hydrogenotrophic methanogenesis (it carries `frh` but does not grow on H₂ + CO₂).
- **Recoded genes** (selenocysteine, pyrrolysine). From genome input they are called
  when one half passes on its own or the two halves can be joined (above); a gene
  whose halves do not add up to the threshold is still missed. From a proteome they
  depend on the annotation — the NCBI proteome of *M. luminyensis* lacks `mtmB`,
  `mtbB` and `mttB`, which are found from its genome sequence.
- **Gene calls move a few scores.** Prodigal's start codon can differ from the
  annotation's; `mtaA` of *M. luminyensis* lost 0.6 bit that way and needed a
  threshold override (`ko_tc`). Expect the same at other knife-edge thresholds.
- **Untested clades**: NC10 (`nod`, pmoA at the KO threshold), USC pmoA, alkane
  oxidizers with alkyl-CoM reductase beyond the one sequence checked.
- **Ubiquitous genes** (`ackA`, `pta`, `acs`, `fdh`, `frmAB`, `hdrA`) light up in most
  genomes; read them through the modules, which need Mcr or another key gene.

## Build / test

```bash
python workflow/scripts/build_hmm_db.py     # KOfam KO profiles (+ custom HMMs / Pfam fallbacks if any)
python workflow/scripts/build_blast_db.py    # curated UniProt seeds → DIAMOND DBs
python validation/check_smoke.py             # re-check existing results/ against the expectations
```

The KOfam cache (`resources/.cache/`) is symlinked to the ncycle-pipeline download
(~1.5 GB profiles.tar.gz, pinned release 2026-05-24, verified by SHA-256 at build
time) to avoid a re-fetch; delete the links to make the build download its own copy.
