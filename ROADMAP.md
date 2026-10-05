# mcycle-pipeline — roadmap

Where the methane tool stands against the sister-tool contract
(`Unify_Tools/shared/sister-tool-contract.md`). The validation campaign was run on
2026-10-04; the phase-by-phase record is [`WORKPLAN.md`](WORKPLAN.md), the results
[`validation/REPORT.md`](validation/REPORT.md).

| Dimension | State |
|---|---|
| Launcher (`run.py`, shared env, no legacy modes) | done; `MCYCLE_CONFIG` selects a per-study config |
| KOfam cache, SHA-pinned (same release as N and S) | done — cache is a symlink to the nitrogen download |
| DB staleness guard (`_stale()` vs targets.yaml) | done |
| Custom-HMM provenance (`targets/<id>/manifest.yaml`) | done — 9 clade HMMs (pmoA x 4, mcrA_anme x 5) from GTDB-typed genomes |
| Smoke panel + expectations (`make smoke`) | done — 18 genomes, 168 expectations |
| Dual ground truth (KEGG + curated overlay) | done — 3,034 / 3,173 cells |
| Held-out split + seed-leak detection | done — 27 training / 22 hold-out, 0 leaked cells |
| Trap-independence audit | done — 46 independent trap positives, precision 1.000 |
| Bootstrap-CI regression gate + CI workflow | done — 14 / 14 checks; `.github/workflows/regression.yml` (template until there is a remote) |
| Domain-DB comparator | done — MCycDB 2021 |
| Benchmark vs KofamScan / METABOLIC / DRAM / MCycDB | done — pre-registered, 8 / 8 contrasts |
| GTDB concordance (500 genomes) | done — four tools incl. METABOLIC; `comparators/gtdb500_m/GTDB500_REPORT.md` |
| Orthogonal check of the clade calls | done — genome lineage + gene-tree placement |
| MAG realism study | done — 12 MAGs from nucleotide input |

## Done after the campaign

- **2026-10-05, amendment (commit `c7bf634`).** mcrB / mcrG of a MAG whose McrA is a
  gene fragment are no longer labelled as alkyl-CoM reductase subunits. Minimal by
  decision: 8 calls change in the GTDB-500 set, none on the reference panel
  (`validation/CHANGELOG.md`, `validation/REPORT.md` section 6).

## Next cycle — new call semantics, own validation

Decided 2026-10-04: these two are NOT amendments. Each changes what the tool can
report, so each gets a design, a frozen test set chosen before it is run, and a
re-run of the gate and of the GTDB check.

1. **Call the truncated McrA.** A fragment at a contig end (4 GTDB-500 genomes;
   1.62-1.86 bits per aligned position against 1.39 for a full-length McrA at the
   threshold) should be an `mcrA` call tagged `partial`. The catch: the ANME clade
   thresholds are calibrated on full-length proteins, so a truncated ANME McrA would
   fail them and default to `methanogenic` — a wrong-direction call, of which there
   are none today. It needs a third state, "direction unresolved", through
   `resolve_mcr_direction`, the synergies, the cycle map and the report, and
   length-aware clade thresholds if the direction is to be called at all.
2. **Non-euryarchaeal McrA** (Verstraetearchaeota and relatives) that score under the
   KOfam threshold at full length (*Ca.* Methanomethylicus GCA_024464205: 578 against
   775.5) and are therefore indistinguishable, by K00399, from an alkyl-CoM reductase
   (642). Needs a canonical-McrA clade model built like the ANME ones
   (`harvest_clade_refs.py`, `build_clade_hmms.py`), with the alkyl-CoM reductases as
   negatives.

## Open, smaller

1. **Ground-truth review of the hold-out candidates** listed in the changelog
   (`mtmB` / `mtbB` of *M. burtonii*, `ftr` of *M. thermophilum*), with independent
   evidence only.
2. **pxmA / pmoA3**: decide whether to report them as their own target.
3. **Residual KO-level confusions** (`acs`, `fdh` / `fdhA` without FdhB, `frhB`, `hdrD`).
4. Archive `resources/blast_db/` and the clade-reference genome list with a release
   (`../mcycle_blast_db_2026-10-04.tar.gz`, SHA-256 beside it).
