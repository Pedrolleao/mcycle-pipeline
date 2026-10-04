# mcycle-pipeline — roadmap

Where the methane tool stands against the sister-tool contract
(`Unify_Tools/shared/sister-tool-contract.md`) that the nitrogen and sulfur pipelines
satisfy. Done = the operational half; open = the whole validation half.

| Dimension | State |
|---|---|
| Launcher (`run.py`, shared env, no legacy modes) | done |
| KOfam cache, SHA-pinned (same release as N and S) | done — cache is a symlink to the nitrogen download |
| DB staleness guard (`_stale()` vs targets.yaml) | done |
| Custom-HMM provenance (`targets/<id>/manifest.yaml`) | schema in place, no HMM trained |
| Smoke panel + expectations (`make smoke`) | done — 18 genomes (genome sequence), 168 expectations |
| Dual ground truth (KEGG + curated overlay) | open |
| Held-out split + seed-leak detection | open |
| Trap-independence audit | open |
| Bootstrap-CI regression gate + CI workflow | open |
| GTDB concordance (500 genomes) | open |
| MAG realism study | open |
| Domain-DB comparator | open — MCycDB is the counterpart of NCycDB / SCycDB |
| Benchmark vs KofamScan / METABOLIC / DRAM | open |

## Next steps, in order

The executable plan, phase by phase, is [`WORKPLAN.md`](WORKPLAN.md); the list below is
its outline.

1. **Close the seed gaps** listed in `../Info-methane.md` (hardening backlog): marine
   ANME McrA, and decoy genomes for `nod`, NC10 pmoA, `mmoX`, `xoxF`, alkyl-CoM
   reductase. Train clade HMMs (pmoA vs amoA; ANME McrA) where a BLAST identity gate
   is too brittle — both current gates sit in gaps only 6–10 identity points from the
   nearest decoy.
2. **Reference panel** (~40 genomes: methanogen orders incl. Methanocellales,
   Methanomicrobiales, Methanonatronarchaeia, Bathy-/Verstraetearchaeota-type Mcr;
   ANME; methanotroph clades; methylotrophs; decoys) with a train / hold-out split.
   Port `build_ground_truth.py`, `score_*.py`, `logo_cv.py`, `trap_independence.py`,
   `test_regression.py` from scycle-pipeline (they are domain-neutral bar file names).
3. **Comparators and scale**: KofamScan, METABOLIC, DRAM, MCycDB; GTDB-500
   concordance; a MAG study on nucleotide input.
4. Only then quote accuracy figures. Until then the honest claim is the smoke panel.
