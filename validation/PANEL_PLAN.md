# PANEL_PLAN — reference panel of mcycle-pipeline

**Roster and train / hold-out split FROZEN on 2026-10-04**, before any ground truth was
built and before the pipeline was run on any of the added genomes. The machine-readable
roster is `validation/panel.tsv`; this file explains it. Later changes are allowed only
as dated additions in `validation/CHANGELOG.md`; the `split` of a genome never changes.

- **49 genomes: 27 training, 22 hold-out.** All are genome sequence (`.fna`), fetched
  by assembly accession with `validation/fetch_ref_panel.sh` into `../ref_panel/`.
- Identity QC (`validation/validate_panel.py --online`, 2026-10-04): 49 / 49 pass —
  every file's deflines name the expected organism, and every file holds exactly the
  total sequence length of the NCBI assembly in the roster. Table: `panel_qc.tsv`.
- Every accession was resolved through the NCBI Datasets API, every KEGG organism code
  through `rest.kegg.jp/get/gn:<code>` (the entry names the same assembly), every GTDB
  genus through the r232 taxonomy tables, and every DOI through the Crossref API, on
  the day of writing. Nothing in the roster is recalled from memory.

## How the split was made

A genome is **training** if any of these holds, otherwise it may be hold-out:

1. it is one of the 18 smoke-panel genomes (thresholds were tuned while looking at them);
2. its genus supplies a BLAST seed (`config/targets.yaml`, `blast_refs_uniprot`):
   *Methanothermobacter, Methanosarcina, Methanocaldococcus, Methanococcus,
   Methanopyrus, Methanomassiliicoccus, Methanoperedens, Methylococcus,
   Methylotuvimicrobium, Methylomonas, Methylocystis, Methylocapsa, Methylacidiphilum,
   Methylomirabilis, Methylosinus*;
3. a sequence of it was measured during the build (*Ca.* Ethanoperedens: its ethyl-CoM
   reductase was scored against K00399 and the McrA seeds);
4. it is earmarked as a seed or calibration source for M3 (ANME1_G60, the butane
   oxidizers, *P. putida*, *Ca.* Methylomirabilis oxyfera).

"Genus" means the NCBI genus **and** the GTDB r232 genus; a hold-out genome shares
neither with any training genome or seed. For uncultured lineages only the GTDB genus
exists (ANME1_G37 is `g__ANME1a`, ANME1_G60 is `g__G60ANME1`). The check was run
mechanically: no hold-out `taxa` token matches a seed organism in
`resources/blast_db/*.fasta`.

**Binding on M3 and later:** no seed, no HMM training sequence and no threshold
calibration may come from a hold-out genus (either naming). A sequence set pulled from
UniProt or NCBI for an HMM must be filtered against the hold-out genera before use.

## Training genomes (27)

| genome | assembly | contigs | GTDB genus | KEGG | role | published physiology | reference |
|---|---|---|---|---|---|---|---|
| Mmaripaludis_S2 | GCF_000011585.1 | 1 | Methanococcus | mmp | hydrogenotrophic methanogen (Methanococcales) | CH4 from H2/CO2 and formate | Hendrickson 2004, doi:10.1128/jb.186.20.6956-6969.2004 |
| Mthermautotrophicus_DeltaH | GCF_000008645.1 | 1 | Methanothermobacter | mth | hydrogenotrophic methanogen (Methanobacteriales) | CH4 from H2/CO2 | Smith 1997, doi:10.1128/jb.179.22.7135-7155.1997 |
| Macetivorans_C2A | GCF_000007345.1 | 1 | Methanosarcina | mac | acetoclastic + methylotrophic methanogen (Methanosarcinales) | CH4 from acetate, methanol, methylamines; no growth on H2/CO2 | Galagan 2002, doi:10.1101/gr.223902 |
| Msoehngenii_GP6 | GCF_000204415.1 | 2 | Methanothrix | mcj | obligate acetoclastic methanogen | CH4 from acetate only | Barber 2011, doi:10.1128/jb.05031-11 |
| Mluminyensis_B10 | GCF_000308215.1 | 26 | Methanomassiliicoccus | tbd | H2-dependent methylotrophic methanogen (Methanomassiliicoccales) | CH4 from methanol + H2 | Dridi 2012, doi:10.1099/ijs.0.033712-0 |
| Mnitroreducens_ANME2d | GCA_000685155.1 | 10 | Methanoperedens | - | anaerobic methanotroph (ANME-2d) | nitrate-dependent anaerobic methane oxidation; Mcr in reverse | Haroon 2013, doi:10.1038/nature12375 |
| Mcapsulatus_Bath | GCF_000008325.1 | 1 | Methylococcus | mca | aerobic methanotroph (gamma, type X); pmoA + mmoX seed strain | pMMO + sMMO; RuMP | Ward 2004, doi:10.1371/journal.pbio.0020303 |
| Mfumariolicum_SolV | GCA_000953475.1 | 1 | Methylacidiphilum | tbd | aerobic methanotroph (Verrucomicrobia) | pMMO; lanthanide methanol dehydrogenase only; no sMMO | Pol 2007, doi:10.1038/nature06222 |
| Mtrichosporium_OB3b | GCF_002752655.1 | 4 | Methylosinus | mtw | aerobic methanotroph (alpha, type II); mmoX seed species | pMMO + sMMO; serine cycle | Stein 2010, doi:10.1128/jb.01144-10 |
| Msilvestris_BL2 | GCF_000021745.1 | 1 | Methylocapsa | msl | facultative aerobic methanotroph (alpha), sMMO only | sMMO, no pMMO | Chen 2010, doi:10.1128/jb.00506-10 |
| Mextorquens_AM1 | GCF_000022685.1 | 5 | Methylobacterium | mea | methylotroph, not a methanotroph | methanol and methylamine; no methane oxidation | Vuilleumier 2009, doi:10.1371/journal.pone.0005584 |
| Neuropaea_ATCC19718 | GCF_000009145.1 | 1 | Nitrosomonas | neu | decoy: beta-AOB (amoCAB) | ammonia oxidizer; no methane-dependent growth | Chain 2003, doi:10.1128/jb.185.21.6496.2003 |
| Noceani_ATCC19707 | GCF_000012805.1 | 2 | Nitrosococcus | noc | decoy: gamma-AOB (amoA closest to pmoA) | ammonia oxidizer; no methane-dependent growth | Klotz 2006, doi:10.1128/aem.00463-06 |
| Nmaritimus_SCM1 | GCF_000018465.1 | 1 | Nitrosopumilus | nmr | decoy: AOA (archaeal amoA) | ammonia oxidizer | Walker 2010, doi:10.1073/pnas.0913533107 |
| Afulgidus_DSM4304 | GCF_000008665.1 | 1 | Archaeoglobus | afu | decoy: archaeal sulfate reducer with H4MPT C1 genes + CODH/ACS, no Mcr | sulfate reducer; not a methanogen | Klenk 1997, doi:10.1038/37052 |
| Dvulgaris_Hildenborough | GCF_000195755.1 | 2 | Nitratidesulfovibrio | dvu | decoy: sulfate reducer (Hdr-like, CODH) | sulfate reducer | Heidelberg 2004, doi:10.1038/nbt959 |
| Ecoli_K12_MG1655 | GCF_000005845.2 | 1 | G047199095 | eco | decoy: ack / pta / acs + frmAB only | no C1 metabolism beyond formaldehyde detoxification | Blattner 1997, doi:10.1126/science.277.5331.1453 |
| Spneumoniae_ref | GCF_001457635.1 | 1 | Streptococcus | - | negative (strain NCTC7465; the citation is the species genome, strain TIGR4) | nothing methane-related | Tettelin 2001, doi:10.1126/science.1061217 |
| Mkandleri_AV19 | GCF_000007185.1 | 1 | Methanopyrus | mka | hydrogenotrophic methanogen (Methanopyrales); mcrA seed strain | CH4 from H2/CO2 only | Kurr 1991, doi:10.1007/bf00262992 |
| Mthermoacetophila_PT | GCF_000014945.1 | 1 | Methanothrix_B | mtp | obligate acetoclastic methanogen, second Methanothrix | CH4 from acetate only | Kamagata 1991, doi:10.1099/00207713-41-2-191 |
| ANME1_G60 | GCA_003194435.1 | 38 | G60ANME1 | - | anaerobic methanotroph, ANME-1 (60 C consortium); intended seed source for marine-ANME McrA in M3 | anaerobic methane oxidation with a sulfate-reducing partner; Mcr in reverse | Krukenberg 2018, doi:10.1111/1462-2920.14077 |
| Ethanoperedens_E50 | GCA_905171685.1 | 1 | EX4572-44 | - | decoy: anaerobic ethane oxidizer (ethyl-CoM reductase); its EcrA was scored in the smoke phase | ethane oxidation; not methane | Hahn 2020, doi:10.1128/mbio.00600-20 |
| Methylocystis_SC2 | GCF_000304315.1 | 1 | Methylocystis | msc | aerobic methanotroph (alpha), pMMO only; pmoA seed genus | two pMMO isozymes; no sMMO | Baani 2008, doi:10.1073/pnas.0702643105 |
| Moxyfera | GCF_000091165.1 | 1 | Methylomirabilis | mox | nitrite-dependent methanotroph (NC10); pmoA seed genus | pMMO + NO dismutase (intra-aerobic) | Ettwig 2010, doi:10.1038/nature08883 |
| Tbutanivorans_NBRC103042 | GCF_001591165.1 | 156 | Thauera | - | decoy: butane oxidizer with a soluble di-iron butane monooxygenase (closest relative of sMMO); type strain | butane oxidation; sBMO is not a methane monooxygenase | Cooley 2009, doi:10.1099/mic.0.028175-0 |
| Nocardioides_CF8 | GCF_000389985.1 | 1 | Nocardioides | - | decoy: butane oxidizer with a copper membrane monooxygenase (pBMO) | butane oxidation; pBMO is not a methane monooxygenase | Sayavedra-Soto 2011, doi:10.1111/j.1758-2229.2010.00239.x |
| Pputida_KT2440 | GCF_000007565.2 | 1 | Aquipseudomonas | ppu | decoy: PQQ alcohol dehydrogenases PedE (Ca) and PedH (lanthanide) | not a methylotroph; PedH is an ethanol dehydrogenase | Wehrmann 2017, doi:10.1128/mbio.00570-17 |

## Hold-out genomes (22)

| genome | assembly | contigs | GTDB genus | KEGG | role | published physiology | reference |
|---|---|---|---|---|---|---|---|
| Mbsmithii_ATCC35061 | GCF_000016525.1 | 1 | Methanocatella | msi | hydrogenotrophic methanogen (Methanobacteriales, human gut) | CH4 from H2/CO2 and formate | Samuel 2007, doi:10.1073/pnas.0704189104 |
| Msstadtmanae_DSM3091 | GCF_000012545.1 | 1 | Methanosphaera | mst | H2 + methanol specialist (Methanobacteriales) | CH4 only from methanol + H2; cannot reduce CO2 to CH4 | Fricke 2006, doi:10.1128/jb.188.2.642-658.2006 |
| Mhungatei_JF1 | GCF_000013445.1 | 1 | Methanospirillum | mhu | hydrogenotrophic methanogen (Methanomicrobiales) | CH4 from H2/CO2 and formate | Gunsalus 2016, doi:10.1186/s40793-015-0124-8 |
| Mcpaludicola_SANAE | GCF_000011005.1 | 1 | Methanocella | mpd | hydrogenotrophic methanogen (Methanocellales) | CH4 from H2/CO2 and formate | Sakai 2008, doi:10.1099/ijs.0.65571-0 |
| Mburtonii_DSM6242 | GCF_000013725.1 | 1 | Methanococcoides | mbu | obligate methylotrophic methanogen (Methanosarcinales) | CH4 from methylamines and methanol; not from H2/CO2 or acetate | Allen 2009, doi:10.1038/ismej.2009.45 |
| Malvi_Mx1201 | GCF_000300255.2 | 1 | Methanomethylophilus | max | H2-dependent methylotrophic methanogen (Methanomassiliicoccales, second family) | CH4 from methanol / methylamines + H2 | Borrel 2012, doi:10.1128/jb.01867-12 |
| Mnthermophilum_AMET1 | GCF_002153915.1 | 8 | Methanonatronarchaeum | - | methyl-reducing methanogen (Methanonatronarchaeia) | CH4 from C1-methylated compounds with formate or H2 | Sorokin 2017, doi:10.1038/nmicrobiol.2017.81 |
| ANME1_G37 | GCA_003194425.1 | 5 | ANME1a | - | anaerobic methanotroph, ANME-1 (37 C consortium); GTDB genus differs from ANME1_G60 | anaerobic methane oxidation with a sulfate-reducing partner; Mcr in reverse | Krukenberg 2018, doi:10.1111/1462-2920.14077 |
| Sbutanivorans_BOX1 | GCA_001766825.1 | 16 | Syntropharchaeum | - | decoy: anaerobic butane oxidizer (alkyl-CoM reductases) | butane oxidation; not methane | Laso-Pérez 2016, doi:10.1038/nature20152 |
| Malbum_BG8 | GCF_000214275.2 | 2 | Methylomicrobium | - | aerobic methanotroph (gamma, type I) | pMMO; no sMMO | Kits 2013, doi:10.1128/genomea.00170-13 |
| Mcmarinum_S8 | GCF_003584645.1 | 1 | Methylocaldum | mmai | aerobic methanotroph (gamma, Methylococcaceae, thermotolerant) | obligate methanotroph; MMO complement to be read from the paper in M2 | Takeuchi 2014, doi:10.1099/ijs.0.063503-0 |
| Mfstellata_AR4 | GCF_000385335.1 | 1 | Methyloferula | - | aerobic methanotroph (alpha), sMMO only | sMMO, no pMMO | Vorobev 2011, doi:10.1099/ijs.0.028118-0 |
| Methylacidimicrobium_AP8 | GCF_903064525.1 | 1 | Methylacidimicrobium | meap | aerobic methanotroph (Verrucomicrobia, second genus) | methane and H2 oxidizer; pMMO | Picone 2021, doi:10.3389/fmicb.2021.637762 |
| Mbflagellatus_KT | GCF_000013705.1 | 1 | Methylobacillus | mfa | obligate methylotroph, not a methanotroph | methanol and methylamine; no methane oxidation | Chistoserdova 2007, doi:10.1128/jb.00045-07 |
| Nsmultiformis_ATCC25196 | GCF_000196355.1 | 4 | Nitrosospira | nmu | decoy: beta-AOB, second genus | ammonia oxidizer | Norton 2008, doi:10.1128/aem.02722-07 |
| Ngterrae_TAO100 | GCF_002356115.1 | 2 | Nitrosoglobus | ntt | decoy: gamma-AOB, second genus | ammonia oxidizer | Hayatsu 2017, doi:10.1038/ismej.2016.191 |
| Ninopinata_ENR4 | GCF_001458695.1 | 1 | Nitrospira_D | nio | decoy: comammox Nitrospira (amoA clade A) | complete ammonia oxidizer | Daims 2015, doi:10.1038/nature16461 |
| Mchubuense_NBB4 | GCF_000266905.1 | 3 | Mycobacterium | mcb | decoy: hydrocarbon oxidizer with a copper membrane monooxygenase (HMO) and soluble di-iron monooxygenases | ethene / alkane oxidation; not methane | Coleman 2011, doi:10.1038/ismej.2011.98 |
| Goxydans_621H | GCF_000011685.1 | 6 | Gluconobacter | gox | decoy: PQQ-dependent membrane alcohol dehydrogenases | incomplete oxidation of alcohols and sugars; not a methylotroph | Prust 2005, doi:10.1038/nbt1062 |
| Awoodii_DSM1030 | GCF_000247605.1 | 1 | Acetobacterium | awo | decoy: acetogen (Wood-Ljungdahl pathway, CODH/ACS, methyltransferases) | acetate from H2/CO2 and methyl groups; no methane | Poehlein 2012, doi:10.1371/journal.pone.0033439 |
| Dautotrophicum_HRM2 | GCF_000020365.1 | 2 | Desulforapulum | dat | decoy: completely oxidizing sulfate reducer (Wood-Ljungdahl pathway, Hdr-like) | sulfate reducer | Strittmatter 2009, doi:10.1111/j.1462-2920.2008.01825.x |
| Fplacidus_DSM10642 | GCF_000025505.1 | 1 | Ferroglobus | fpl | decoy: Archaeoglobi without Mcr, second genus | Fe(II) oxidizer / nitrate reducer; not a methanogen | Hafenbradl 1996, doi:10.1007/s002030050388 |

`KEGG`: organism code whose genome entry names this assembly; `-` = not in KEGG (its
ground truth is curated by hand in M2); `tbd` = smoke genome whose code is resolved in M2.
The physiology column is what the cited paper reports, and is the source for
`phenotype_gt.tsv` in M2 — it is not a prediction of what the tool will call.

## Guild coverage

| guild | training | hold-out |
|---|---|---|
| hydrogenotrophic methanogens | Methanococcales, Methanobacteriales, Methanopyrales | Methanobacteriales (*Methanobrevibacter*), Methanomicrobiales, Methanocellales |
| acetoclastic methanogens | *Methanosarcina*, *Methanothrix* x 2 | — (both acetoclastic genera are training) |
| methylotrophic methanogens | *Methanosarcina* | *Methanococcoides* (obligate) |
| H2-dependent methyl reducers | *Methanomassiliicoccus* | *Methanosphaera*, *Methanomethylophilus*, *Methanonatronarchaeum* |
| ANME | ANME-2d, ANME-1 (G60) | ANME-1 (G37, other GTDB genus) |
| alkyl-CoM reductase decoys | *Ca.* Ethanoperedens | *Ca.* Syntrophoarchaeum |
| methanotrophs, gamma | *Methylococcus* (pMMO + sMMO) | *Methylomicrobium* (pMMO), *Methylocaldum* |
| methanotrophs, alpha | *Methylosinus* (both), *Methylocella* (sMMO), *Methylocystis* (pMMO) | *Methyloferula* (sMMO only) |
| methanotrophs, Verrucomicrobia | *Methylacidiphilum* | *Methylacidimicrobium* |
| methanotrophs, NC10 | *Ca.* Methylomirabilis | — |
| methylotrophs without MMO | *Methylorubrum* | *Methylobacillus* |
| ammonia-oxidizer decoys | beta-AOB, gamma-AOB, AOA | beta-AOB, gamma-AOB (second genera), comammox |
| Cu-monooxygenase hydrocarbon decoys | *Nocardioides* sp. CF8 | *Mycolicibacterium chubuense* NBB4 |
| di-iron monooxygenase decoys | *Thauera butanivorans* | *M. chubuense* NBB4 |
| PQQ alcohol-dehydrogenase decoys | *Pseudomonas putida* | *Gluconobacter oxydans* |
| Wood-Ljungdahl / Hdr decoys | *Archaeoglobus*, *Desulfovibrio* | *Acetobacterium*, *Desulforapulum*, *Ferroglobus* |
| general negatives | *E. coli*, *S. pneumoniae* | (the decoys above) |

## What the hold-out cannot test, and why

- **ANME-2d.** The only seeded genus is *Methanoperedens*; the second genus-level
  lineage in GTDB (`g__Methanoperedens_A`) has no closed or near-closed genome — the
  best assemblies are MAGs of 180-270 scaffolds. It goes to the MAG study (M8).
- **Marine ANME-2a / -2c / -3.** Best assemblies are MAGs of 58-230 contigs; they go
  to M8. ANME-1 is in the panel because two good genomes of different GTDB genera
  exist (G37ANME1, 5 scaffolds in the file; G60ANME1, 38).
- **Alpha-proteobacterial pMMO.** Every alpha methanotroph genus with pMMO
  (*Methylosinus, Methylocystis, Methylocapsa*) is a seed genus; the hold-out alpha
  methanotroph is sMMO-only. Alpha pmoA is tested out-of-genus only through the
  leave-one-clade-out cross-validation of M3.
- **NC10** (`nod`, NC10 pmoA). One genus with genomes, and it is a seed genus.
- **Acetoclastic methanogenesis.** *Methanosarcina* and *Methanothrix* are the only
  acetoclastic genera and both are training.
- **Methanopyrales, Methanococcales.** *Methanopyrus* is a single seed genus;
  Methanococcales is represented in training only.

## Deviations from the work plan

- 49 genomes instead of ~45 (22 hold-out instead of ~20): each trap needed a decoy on
  both sides of the split.
- The ANME-2d second lineage and marine ANME-2 / -3 are moved to M8 for lack of a
  high-quality genome (the plan allowed this for marine ANME; it is extended to ANME-2d).
- Five hold-out genomes are not in KEGG (ANME1_G37, Sbutanivorans_BOX1,
  Mnthermophilum_AMET1, Malbum_BG8, Mfstellata_AR4), and six training genomes are not
  either (the ANME-2d and ANME-1 genomes, *Ca.* Ethanoperedens, *Thauera butanivorans*,
  *Nocardioides* sp. CF8, *S. pneumoniae* NCTC7465). M2 has to curate their ground
  truth from the literature and the NCBI annotation.
