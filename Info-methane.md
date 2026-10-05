# Info-methane.md — the methane-cycle atlas behind `mcycle-pipeline`

Narrative rationale for every target, complex, and synergy in
`config/targets.yaml`. Companion to `Info-nitrogen.md` and
`Info-sulfur.md`; the detection machinery is identical (KO-primary KOfam HMMs +
adaptive thresholds, an optional custom-HMM tier, DIAMOND-BLAST gating).
All KO numbers were verified against the cached KOfam `ko_list` (release 2026-05-24),
all Pfam ids against the InterPro API, and all BLAST seed accessions against the
UniProt REST API (2026-10-04).

The methane cycle has two halves that meet at CH₄. **Methanogenic archaea** make it —
from CO₂ + H₂ / formate, from acetate, or from methylated compounds — through one
terminal enzyme, methyl-coenzyme M reductase (Mcr). **Methanotrophs** consume it:
aerobic bacteria with a methane monooxygenase (CH₄ → CH₃OH → HCHO → CO₂ / biomass),
and anaerobic archaea (ANME) that run the methanogenesis pathway backwards. The
pipeline resolves a genome's role across **7 process modules / 83 marker genes**.

---

## 1. Mcr / Mtr core — methanogenesis and reverse methanogenesis

Common to every methanogen and every ANME archaeon.

| gene | KO | role |
|---|---|---|
| mcrA / mcrB / mcrG | K00399 / K00401 / K00402 | methyl-coenzyme M reductase α / β / γ (methyl-CoM + CoB ⇄ CH₄ + CoM-S-S-CoB) |
| mcrC / mcrD | K03421 / K03422 | Mcr operon proteins (activation / assembly; accessory) |
| mcrA_anme | K00399 + seeds | the same McrA, called against **ANME-clade** seeds — the direction marker (trap 1) |
| mtrA–H | K00577–K00584 | Na⁺-pumping methyl-H₄MPT:CoM methyltransferase (mtrF / mtrG reported, not required) |

## 2. CO₂ reduction — the H₄MPT C1 module (CO₂ → methyl-H₄MPT)

| gene | KO | role |
|---|---|---|
| fwdA / B / C / D | K00200 / K00201 / K00202 / K00203 | formylmethanofuran dehydrogenase (CO₂ → formyl-MFR) — **A/B/C shared with bacterial FhcABC** |
| ftr | K00672 | formyl-MFR:H₄MPT formyltransferase — **shared with bacterial FhcD** |
| mch | K01499 | methenyl-H₄MPT cyclohydrolase — **shared with methylotrophs** |
| mtd / hmd | K00319 / K13942 | methylene-H₄MPT dehydrogenase (F₄₂₀) / H₂-forming alternative |
| mer | K00320 | methylene-H₄MPT reductase |
| fdhA / fdhB | K22516 / K00125 | F₄₂₀-dependent formate dehydrogenase (formate as electron donor) |

fwdF/G/H are not targets: their KOs (K00204, K00205, K11260) are generic
"4Fe-4S ferredoxin" profiles.

## 3. Acetoclastic methanogenesis (acetate → CH₄ + CO₂)

| gene | KO | role |
|---|---|---|
| ackA / pta | K00925 / K00625, K13788 | acetate activation in *Methanosarcina* — ⚠️ ubiquitous in bacteria |
| acs | K01895 | acetate activation in *Methanothrix* — ⚠️ ubiquitous |
| cdhA / B / C / D / E | K00192 / K00195 / K00193 / K00194 / K00197 | CODH/ACS complex (cleaves acetyl-CoA; **synthesizes** it in autotrophs, acetogens, sulfate reducers) |

`pta` carries a threshold override (`ko_tc: 350`): the genuine Pta of
*M. acetivorans* scores 378 on K00625, under the KOfam threshold of 391.4.

## 4. Methylotrophic methanogenesis

| gene | KO | role |
|---|---|---|
| mtaB / mtaC / mtaA | K04480 / K14081 / K14080 | methanol:corrinoid methyltransferase, corrinoid protein, methylcobamide:CoM MT |
| mtmB / mtmC | K16176 / K16177 | monomethylamine MT + corrinoid protein |
| mtbB / mtbC | K16178 / K16179 | dimethylamine MT + corrinoid protein |
| mttB / mttC | K14083 / K14084 | trimethylamine MT + corrinoid protein |
| mtbA | K14082 | methylcobamide:CoM MT (methylamine-specific) |
| mtsA / mtsB | K16954 / K16955 | methylthiol:CoM MT system (DMS, methanethiol) |

⚠️ MtmB / MtbB / MttB contain **pyrrolysine** (an in-frame UAG codon) — see
"Recoded genes" below. `mtaA` carries a threshold override (`ko_tc: 375`, KOfam
382.7): from genome sequence the *M. luminyensis* MtaA starts 8 residues late and
scores 382.1.

## 5. Energy conservation — heterodisulfide reductases and hydrogenases

| gene | KO | role |
|---|---|---|
| hdrA / B / C | K03388, K22480 / K03389, K22481 / K03390, K22482 | soluble, electron-bifurcating heterodisulfide reductase — ⚠️ HdrA-like proteins are widespread outside methanogens |
| hdrD / hdrE | K08264 / K08265 | membrane-bound heterodisulfide reductase (cytochrome-containing Methanosarcinales) |
| mvhA / G / D | K14126 / K14128 / K14127 | F₄₂₀-non-reducing hydrogenase (electron donor to HdrABC) |
| frhA / B / G | K00440 / K00441 / K00443 | F₄₂₀-reducing hydrogenase |

## 6. Aerobic methane oxidation (CH₄ → CH₃OH → HCHO)

| gene | KO | role |
|---|---|---|
| pmoA / B / C | K10944 (+K28504) / K10945 / K10946 | particulate (copper) methane monooxygenase — ⚠️ trap 2 |
| mmoX / Y / Z | K16157 / K16158 / K16159 | soluble (di-iron) methane monooxygenase hydroxylase |
| mmoB / C / D | K16160 / K16161 / K16162 | sMMO regulatory protein / reductase / MmoD |
| mxaF / mxaI | K14028 / K14029 | calcium-dependent methanol dehydrogenase |
| xoxF | K23995 | lanthanide-dependent methanol dehydrogenase |
| nod | K27148 | nitric oxide dismutase (NC10 "intra-aerobic" methane oxidation) — unvalidated |

## 7. Formaldehyde oxidation and C1 assimilation

| gene | KO | role |
|---|---|---|
| fae | K10713, K13812 | formaldehyde-activating enzyme (HCHO + H₄MPT → methylene-H₄MPT) |
| mtdB | K10714 | NAD(P)-dependent methylene-H₄MPT dehydrogenase (bacterial) |
| gfa / frmA / frmB | K03396 / K00121 / K01070 | glutathione-dependent formaldehyde oxidation |
| fdh | K00123, K00122 | NAD⁺ / quinone-linked formate dehydrogenase — ⚠️ widespread |
| hxlA / hxlB | K08093, K13831, K13812 / K08094, K13831 | RuMP cycle: hexulose-6-P synthase / isomerase |
| sgaA, hprA, gckA | K00830, K00018, K11529 | serine cycle (⚠️ broad KOs) |
| mtkA / mtkB, mcl | K14067 / K08692, K08691 | malate thiokinase, malyl-CoA lyase (serine cycle) |

The bacterial H₄MPT route also needs `mch` and the Fhc complex — the `mch`, `ftr`
and `fwdA/B/C` targets of module 2.

---

## Trap 1 — Mcr: methanogenesis vs reverse methanogenesis (ANME)

Mcr is the **same enzyme** in methanogens and in anaerobic methanotrophic archaea,
which oxidize methane by running the whole pathway backwards. Gene content does not
separate them — ANME-2 carry the complete methanogenesis pathway. This is the
structure of the sulfur tool's `dsrAB` trap, but here the sequence does carry signal:
ANME McrA form their own clades.

`mcrA_anme` is the same protein evaluated a second time, against **five ANME clade
HMMs** (ANME-1, -2a/b, -2c, -2d, -3; `targets/mcrA_anme/`), each with
its own threshold. Until the hardening of 2026-10-04 it was a BLAST gate at ≥ 80 %
identity to a single ANME-2d seed — the measurements below are from that stage and
show why the gate was brittle. `apply_rules.resolve_mcr_direction` then tags
`mcrA`/`mcrB`/`mcrG` as `…|mcr_reverse` when the gate passes, `…|mcr_methanogenic`
otherwise; the methanogenesis modules and `reverse_methanogenesis_aom` each need
their own direction and are ruled out under the other. Measured on 2026-10-04:

| comparison | identity to the ANME-2d seed (A0A284VQI9) |
|---|---|
| *Ca.* Methanoperedens, second genome (the panel's ANME-2d) | 96.1 % |
| *Methanosarcina barkeri* / *acetivorans* McrA | 70.1 % / 68.5 % |
| *Methanothrix soehngenii* McrA | 63.1 % |
| other methanogen orders | 55–61 % |

**Clade models (2026-10-04).** UniProt has no function-verified, full-length McrA for
the marine ANME clades, so the training sequences were taken from genomes instead:
37 McrA sequences from genomes of the ANME families of GTDB r232 (positives) and 80 from methanogen
and alkyl-CoM-reductase genomes (negatives), one model per clade, thresholds set by
leave-one-genus-out (`validation/HARDENING.md`). ANME-1, -2a/b, -2c and -2d separate
from methanogen McrA; **ANME-3 does not separate cleanly** from methylotrophic
Methanosarcinaceae (*Methanomethylovorans*, *Methanolobus*), and its model is set for
precision. The call is about the clade, not about what the cell is doing: trace
methane oxidation by methanogens and methane production by ANME have both been
reported.

**Alkyl-CoM reductases** (short-chain-alkane oxidizers) are divergent McrA
homologues. KOfam already excludes the one tested: the ethyl-CoM reductase of
*Ca.* Ethanoperedens (A0A7R9R780) scores 643 on K00399 (threshold 775.5; canonical
McrA ≥ 1000, ANME-2d ≈ 915) and is 42–45 % identical to every McrA seed, under the
50 % floor of the `mcrA` BLAST corroboration.

## Trap 2 — pmoABC vs amoABC

Particulate methane monooxygenase and ammonia monooxygenase are one enzyme family —
the mirror image of the nitrogen tool's `amoA` trap. KOfam does not separate them:

| protein | K10944 "methane MO A" (TC 353.0) | K28504 "ammonia MO A" (TC 273.2) |
|---|---|---|
| *M. capsulatus* pmoA | 418 | 276 |
| *M. fumariolicum* pmoA1 / pmoA2 | 405 / 397 | 274 / 273 |
| *M. trichosporium* pmoA | 399 | 249 |
| ***Nitrosococcus oceani* amoA** (γ-AOB) | **384** | 260 |
| *Nitrosomonas europaea* amoA (β-AOB) | 323 | 482 |
| *Nitrosopumilus maritimus* amoA (AOA) | — | 280 |

A γ-proteobacterial ammonia oxidizer passes the *methane* monooxygenase KO. The clade
call is made on subunit A by **four clade HMMs** (γ, α, verrucomicrobial, NC10;
`targets/pmoA/`), trained on 76 PmoA sequences from methanotroph genomes, with the
AmoA of β-AOB, γ-AOB and AOA and a butane monooxygenase as negatives; every left-out
methanotroph genus still passes its clade model (54 / 54 γ, 41 / 41 α). Before
2026-10-04 the call was a BLAST gate at ≥ 70 % identity to a seed. Measured gap then: the best
decoy is *N. oceani* amoA at 63.5 %; true pmoA reach 78–99 % to a seed of their own
clade with their own strain left out. `pmoB`/`pmoC` share K10945 / K10946 verbatim
with amoB / amoC, so `apply_rules.gate_pmo_subunits` lets them follow the pmoA call:
without a confirmed pmoA they are `disqualified`.

Not covered: the divergent `pxmA` and verrucomicrobial `pmoA3` paralogues (function
unresolved; pmoA3 scores under the K10944 threshold), the atmospheric-methane clades
(USCγ, most of USCα), and the hydrocarbon monooxygenases of the same family. NC10
pmoA sits at the threshold (358 vs 353) — a knife edge with no test genome yet.

## Trap 3 — the shared H₄MPT C1 module

`fwdA/B/C`, `ftr` and `mch` are the same KOs in methanogens (CO₂ → methyl, reductive)
and in methylotrophic bacteria, where FhcABCD + Mch oxidize formaldehyde to formate.
The genes are called once; direction is assigned above them:

- synergy layer — `hydrogenotrophic_methanogenesis` needs Mcr (key gene);
  `h4mpt_formaldehyde_oxidation` needs the bacterial `mtdB` (key gene);
- cycle map — in a genome with `fae` and no Mcr the shared genes light the
  formyl-H₄MPT → formate step and `mch` points upwards.

*Archaeoglobus fulgidus* (sulfate reducer; fwd, ftr, mch, mtd, mer and CODH/ACS but
no Mcr) is the panel's decoy for this trap.

## Recoded genes — selenocysteine and pyrrolysine

Methanogens read two stop codons as amino acids: UGA as selenocysteine (formate
dehydrogenase, hydrogenases, Fwd, Hdr in *Methanococcus*) and UAG as pyrrolysine (the
methylamine methyltransferases). A gene caller stops at them, so from genome sequence
such a gene arrives as two neighbouring ORFs. Measured on the smoke panel (2026-10-04):

| gene | halves (aa) | scores | KOfam threshold | outcome |
|---|---|---|---|---|
| *M. maripaludis* fdhA, copy 1 | 133 + 528 | 221 + 782 | 989.5 | joined |
| *M. maripaludis* fdhA, copy 2 | 132 + 539 | 216 + 837 | 989.5 | joined |
| *M. luminyensis* / *M. acetivorans* mtmB, mtbB, mttB | — | one half passes alone | — | called without joining |

`apply_rules.join_inframe_stops` scores two ORFs as one gene only when they are
neighbours on one strand, the second continues **in the reading frame of the first**
across at most 150 nt, they hit **consecutive, non-overlapping** parts of the same KO
profile, and their scores together reach the threshold. The conditions matter:
adjacent ORFs hitting one profile are common in the panel (tandem paralogues hit the
same range twice; genuinely two-gene or frameshifted arrangements are out of frame or
overlap), and only four pairs in 18 genomes pass — the two fdhA above, one of them
again as the generic `fdh`, and an interrupted mtaA-like gene in *M. acetivorans*
that changes no call. Every such call is tagged `joined:<other ORF>`: sequence alone
cannot tell recoding from a nonsense mutation.

Switching the panel from annotated proteomes to genome sequence changed 5 of 1,494
gene calls: `mtmB` / `mtbB` / `mttB` appeared in *M. luminyensis* (its NCBI proteome
lacks them, although it grows on methylamines), and `acs` appeared in
*M. thermautotrophicus* and *M. trichosporium* (full-length Prodigal genes well above
threshold; why the proteomes lacked them was not traced).

### Other shared signatures
- **ack / pta / acs, cdh, hdrABC** occur far outside methanogens — they count only
  inside modules that also require Mcr.
- **acetoclastic vs autotrophic**: a hydrogenotroph with CODH/ACS + acs (acetate
  assimilation) has every acetoclastic gene except one. The module's key gene is
  `hdrD`: membrane-bound HdrDE belongs to the cytochrome-containing Methanosarcinales,
  the only acetoclasts.
- **hxlA / hxlB** run in reverse in archaea (pentose-phosphate synthesis), so every
  methanogen carries them; `rump_assimilation` therefore also needs a methanol
  dehydrogenase.
- **mmoX** shares the di-iron monooxygenase fold with propane / butane / phenol /
  toluene monooxygenases; K16157 is specific and two reviewed seeds corroborate.

---

## Obligatory complexes & process synergies

**Complexes (17):** mcr (mcrABG), mtr (mtrABCDE+H), fwd (fwdABCD), hdr_abc, hdr_de,
mvh_hydrogenase, frh_hydrogenase, f420_formate_dh, codh_acs (cdhABCDE),
methanol / monomethylamine / dimethylamine / trimethylamine methyltransferase, pmmo,
smmo_hydroxylase, mxa_methanol_dh, malate_thiokinase.

**Synergies (11).** Slots are `requires` genes plus `requires_any` choices; `key`
entries define the module — without them it is ruled out, so shared genes alone never
make a module look half-present.

| module | slots | key |
|---|---|---|
| hydrogenotrophic_methanogenesis | fwdAB, ftr, mch, mer, mtrAH, mcrABG; mtd \| hmd; frhABG \| fdhAB | mcrA |
| acetoclastic_methanogenesis | cdhCDE, mtrA, mcrABG, hdrD; ackA+pta \| acs | mcrA, hdrD |
| methylotrophic_methanogenesis | mcrABG; mtaABC \| mtbA+mtmBC \| mtbA+mtbBC \| mtbA+mttBC \| mtsAB | mcrA + a substrate MT |
| reverse_methanogenesis_aom | mcrABG, mtrAH, mch, ftr, fwdAB — direction `reverse` | mcrA |
| heterodisulfide_recycling | mcrA; hdrABC \| hdrDE | mcrA |
| aerobic_methane_oxidation | pmoABC \| mmoXYZ; mxaFI \| xoxF | pmoA \| mmoX |
| nitrite_dependent_methane_oxidation | pmoABC, nod | nod |
| methanol_oxidation | mxaFI \| xoxF | — |
| h4mpt_formaldehyde_oxidation | fae, mtdB, mch, ftr, fwdA | mtdB |
| rump_assimilation | hxlAB; mxaFI \| xoxF | mxaF \| xoxF |
| serine_cycle | sgaA, hprA, mtkAB, mcl; mxaFI \| xoxF | mxaF \| xoxF |

---

## Provenance of the BLAST seeds

| target | accessions | what they are |
|---|---|---|
| mcrA | P11558, P58815, P07962, Q58256, P07961, Q49605 | reviewed McrA: *Methanothermobacter marburgensis* I + II, *Methanosarcina barkeri*, *Methanocaldococcus jannaschii*, *Methanococcus vannielii*, *Methanopyrus kandleri* |
| mcrA | R9TAZ9, A0A284VQI9 | *Methanomassiliicoccus intestinalis*; *Ca.* Methanoperedens (ANME-2d) |
| mcrA_anme | A0A284VQI9 | *Ca.* Methanoperedens nitratireducens McrA, 567 aa |
| pmoA | Q607G3 (reviewed), G4SZ63, G0A0H7 | γ: *Methylococcus capsulatus* Bath, *Methylotuvimicrobium alcaliphilum* 20Z, *Methylomonas methanica* MC09 |
| pmoA | A0A3G8M5C3, A0A1I3Z8P3 | α: *Methylocystis rosea*, *Methylocapsa palsarum* |
| pmoA | A9QPE3, A9QPE6 | Verrucomicrobia: *Methylacidiphilum infernorum* V4 pmoA1, pmoA2 |
| pmoA | A0A2T4U1K1 | NC10: *Ca.* Methylomirabilis limnetica |
| mmoX | P22869, P27353 | reviewed sMMO α: *M. capsulatus* Bath, *Methylosinus trichosporium* |

Seeds were chosen so that the panel's test genomes are hit through a *different*
strain or species wherever one exists (*M. trichosporium* OB3b through *Methylocystis*,
*M. fumariolicum* through *M. infernorum*, the ANME-2d genome through a second
*Methanoperedens*). Two exceptions are in-sample: *M. capsulatus* Bath (pmoA, mmoX)
and *M. trichosporium* (mmoX P27353). The seeds are now the nearest-reference
annotation for `pmoA` and `mcrA_anme` (decided by the clade HMMs) and the
corroboration for `mcrA` and `mmoX`. Empty
`blast_refs_uniprot` marks targets called on their KO alone; seeds are curated,
**never guessed**.

## After the validation campaign (2026-10-04)

The campaign of `WORKPLAN.md` was run; results are in
`validation/REPORT.md`. What changed in the reference layer, and why,
is logged gene by gene in `validation/HARDENING.md`:

- clade HMMs for `pmoA` and `mcrA_anme` (above);
- `mcrB` / `mcrG` follow McrA when the genome carries an alkyl-CoM reductase (a full
  homologue under the threshold; a gene fragment of a canonical McrA does not count);
  `mmoY/Z/B/C/D` follow `mmoX`;
- `fdhA` vs `fdh` resolved through the F420-binding FdhB;
- `mvhD` needs to cover half of the protein (HdrA–MvhD fusions pass K14127);
- threshold overrides, each with the measured gap: `hmd`, `fwdD`, `hdrB`, `hdrC`, `fdh`.

Decoy margins measured on training genomes: sBMO of *Thauera butanivorans* scores 884
on K16157 (threshold 1168.7, true MmoX 1178–1201); PedH of *Pseudomonas putida* 548
on K23995 (859.1); the ethyl-CoM reductase of *Ca.* Ethanoperedens 642 on K00399
(775.5); NC10 `nod` 1351 against 609 for the nearest qNor-type protein.

Open after the campaign (ordered in `ROADMAP.md`):

- **Truncated McrA in fragmented MAGs** is not called (no `mcrA`, no direction). Since
  the amendment of 2026-10-05 its `mcrB` / `mcrG` are reported normally; before, the
  alkyl-CoM rule mislabelled them. Calling the fragment is a separate cycle.
- **Non-euryarchaeal McrA** can fall under the K00399 threshold at full length and is
  then indistinguishable, by that profile, from an alkyl-CoM reductase. Separate cycle.
- **ANME-3 vs methylotrophic Methanosarcinaceae** is at the limit of what McrA
  sequence resolves.
- **pxmA and verrucomicrobial pmoA3** are recognized as outside the pmoA clades and
  not called; whether to report them as a target is undecided.
- **Known over-call**: `hydrogenotrophic_methanogenesis` is complete for
  *M. acetivorans*, which carries but does not use Frh. Genomic potential, not
  physiology.
- **Recoded genes the joining rule still misses**: halves that do not add up to the
  threshold (the selenocysteine formate dehydrogenases of *E. coli* and
  *D. vulgaris* against the generic `fdh` profile), and splits into more than two ORFs.
- **Not yet targets**: coenzyme M / F₄₃₀ biosynthesis, Fpo / Vht / Ech / Rnf energy
  converters, the multiheme cytochromes of ANME electron export, methylamine
  oxidation in bacteria, NAD-dependent methanol dehydrogenase.
