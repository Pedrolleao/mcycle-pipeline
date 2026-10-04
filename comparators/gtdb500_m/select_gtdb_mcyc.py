#!/usr/bin/env python3
"""
select_gtdb_mcyc.py — the reproducible 500-genome GTDB slice for the METHANE
cross-tool CONCORDANCE study (no ground truth: it measures where tools agree and
disagree at scale).

Mirrors the nitrogen and sulfur selections. The 380-genome cross-phylum BACKBONE is
reused verbatim from ncycle's gtdb500 selection, so the three sister studies share
it; only the 120 methane-ENRICHED species representatives differ: methanogens of
every class, the ANME families, alkane-oxidizing archaea with alkyl-CoM reductase,
the aerobic methanotroph lineages, and the ammonia oxidizers that populate the
pmoA / amoA trap (5 random species representatives per clade, seed 1234).

Input : bac120_taxonomy.tsv + ar53_taxonomy.tsv (GTDB r232; ncycle's gtdb_pilot)
        + ncycle's gtdb500/selection.tsv (the shared backbone rows)
Output: selection.tsv (accession, ncbi_acc, stratum, clade, gtdb_taxonomy)
"""
from __future__ import annotations
import csv, random, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
NCYCLE = Path("/home/dmin/Grants/Nitrogen_Cycle/ncycle-pipeline/comparators")
BAC = NCYCLE / "gtdb_pilot" / "bac120_taxonomy.tsv"
AR = NCYCLE / "gtdb_pilot" / "ar53_taxonomy.tsv"
N_SELECTION = NCYCLE / "gtdb500" / "selection.tsv"   # reuse its 380 backbone rows
SEED = 1234
DEF_OUT = HERE / "selection.tsv"
DEF_PER_CLADE = 5
DEF_ENRICHED_CAP = 120        # cap enriched at this many; backbone fills to TARGET_TOTAL
DEF_TARGET_TOTAL = 500

# clade label -> GTDB lineage substring(s); a genome matches when any is in its lineage
ENRICH = {
    # methanogens
    "Methanobacteria": ["c__Methanobacteria"], "Methanococci": ["c__Methanococci"],
    "Methanomicrobia": ["c__Methanomicrobia"], "Methanocellia": ["c__Methanocellia"],
    "Methanosarcinaceae": ["f__Methanosarcinaceae"], "Methanotrichaceae": ["f__Methanotrichaceae"],
    "Methanomassiliicoccales": ["o__Methanomassiliicoccales"],
    "Methanofastidiosales": ["o__Methanofastidiosales"],
    "Methanomethylicaceae": ["f__Methanomethylicaceae"],          # non-euryarchaeal Mcr
    # anaerobic methanotrophs
    "ANME-1": ["f__Methanospirareceae"], "ANME-2ab": ["f__Methanocomedenaceae"],
    "ANME-2c": ["f__Methanogasteraceae"], "ANME-2d": ["f__Methanoperedenaceae"],
    "ANME-3": ["g__Methanovorans"],
    # alkane oxidizers with alkyl-CoM reductase
    "alkane_AcrA": ["f__EX4572-44", "f__Syntropharchaeaceae", "g__Alkanophaga",
                    "c__Methanoliparia"],
    # aerobic methanotrophs
    "MOB_Methylomonadaceae": ["f__Methylomonadaceae"], "MOB_Methylococcaceae": ["f__Methylococcaceae"],
    "MOB_alpha": ["g__Methylocystis", "g__Methylosinus", "g__Methylocapsa"],
    "MOB_Verrucomicrobia": ["f__Methylacidiphilaceae"], "MOB_NC10": ["f__Methylomirabilaceae"],
    # ammonia oxidizers (pmoA / amoA trap)
    "AOB_beta": ["f__Nitrosomonadaceae"], "AOB_gamma": ["f__Nitrosococcaceae"],
    "AOA": ["f__Nitrosopumilaceae", "f__Nitrososphaeraceae"],
    "Nitrospira_comammox": ["g__Nitrospira_D"],
}


def ncbi_acc(gtdb_acc: str) -> str:
    return gtdb_acc.split("_", 1)[1] if gtdb_acc[:3] in ("RS_", "GB_") else gtdb_acc


def load(path: Path) -> list[tuple[str, str]]:
    out = []
    with open(path) as fh:
        for line in fh:
            acc, lin = line.rstrip("\n").split("\t")
            out.append((acc, lin))
    return out


def species_reps(rows):
    by_sp = defaultdict(list)
    for acc, lin in rows:
        by_sp[lin.split(";")[-1]].append((acc, lin))
    reps = {}
    for sp, cands in by_sp.items():
        cands.sort(key=lambda x: (0 if x[0].startswith("RS_") else 1, x[0]))
        reps[sp] = cands[0]
    return reps


def load_backbone(path: Path) -> list[dict]:
    """Reuse ncycle's 380 backbone rows verbatim (shared cross-phylum set)."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("stratum") == "backbone":
                out.append(r)
    return out


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEF_OUT)
    ap.add_argument("--per-clade", type=int, default=DEF_PER_CLADE)
    ap.add_argument("--enriched-cap", type=int, default=DEF_ENRICHED_CAP)
    ap.add_argument("--target-total", type=int, default=DEF_TARGET_TOTAL)
    args = ap.parse_args()

    rng = random.Random(SEED)
    reps = species_reps(load(BAC) + load(AR))
    print(f"[select] {len(reps)} species-reps", file=sys.stderr)

    chosen: dict[str, dict] = {}
    # ---- enriched stratum (methane clades) ----
    for clade, pat in ENRICH.items():
        if len(chosen) >= args.enriched_cap:
            break
        hits = sorted((acc, lin) for sp, (acc, lin) in reps.items()
                      if any(p + ";" in lin + ";" for p in pat))
        rng.shuffle(hits)
        for acc, lin in hits[:args.per_clade]:
            if acc not in chosen and len(chosen) < args.enriched_cap:
                chosen[acc] = {"accession": acc, "ncbi_acc": ncbi_acc(acc),
                               "stratum": "enriched", "clade": clade, "gtdb_taxonomy": lin}
    n_enriched = len(chosen)
    print(f"[select] enriched: {n_enriched} methane genomes across {len(ENRICH)} clades", file=sys.stderr)

    # ---- backbone stratum: REUSE ncycle's 380 backbone rows verbatim ----
    backbone = load_backbone(N_SELECTION)
    added = 0
    for r in backbone:
        if len(chosen) >= args.target_total:
            break
        if r["accession"] not in chosen:
            chosen[r["accession"]] = r
            added += 1
    print(f"[select] backbone: {added} reused from ncycle (shared cross-phylum set)", file=sys.stderr)

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "ncbi_acc", "stratum", "clade", "gtdb_taxonomy"],
                           delimiter="\t")
        w.writeheader()
        for rec in chosen.values():
            w.writerow(rec)
    print(f"[select] wrote {args.out} ({len(chosen)} genomes: {n_enriched} enriched + {added} backbone)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
