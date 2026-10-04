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
| GTDB concordance (500 genomes) | see `comparators/gtdb500_m/` |
| Orthogonal check of the clade calls | done — genome lineage + gene-tree placement |
| MAG realism study | done — 12 MAGs from nucleotide input |

## Open, in order

1. **Truncated McrA.** Tell a gene fragment from an alkyl-CoM reductase before
   disqualifying `mcrB` / `mcrG`: a fragment scores ~1.8 bits per aligned profile
   position on K00399, a full-length alkyl-CoM reductase or non-canonical McrA ~1.0-1.1
   (measured on the GTDB-500 set). Call the fragment `mcrA` with a `partial` tag. This
   changes the tool after the hold-out was scored: re-run `make regression`, report
   before and after, and log it in `validation/CHANGELOG.md`.
2. **Non-euryarchaeal McrA** (Verstraetearchaeota and relatives) under the KOfam
   threshold: needs a canonical-McrA clade model, built like the ANME ones.
3. **Ground-truth review of the hold-out candidates** listed in the changelog
   (`mtmB` / `mtbB` of *M. burtonii*, `ftr` of *M. thermophilum*), with independent
   evidence only.
4. **pxmA / pmoA3**: decide whether to report them as their own target.
5. **Residual KO-level confusions** (`acs`, `fdh` / `fdhA` without FdhB, `frhB`, `hdrD`).
6. Archive `resources/blast_db/` and the clade-reference genome list with a release
   (`../mcycle_blast_db_2026-10-04.tar.gz`, SHA-256 beside it).
