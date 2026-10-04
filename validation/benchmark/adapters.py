#!/usr/bin/env python3
"""
adapters.py — load the ground truth and every tool's calls into one vocabulary
(the 83 mcycle targets), exactly as fixed in prereg.md (2026-10-04).

Each loader returns {(genome, target): bool}. Comparator loaders read a normalized
TSV written by comparators/build_{metabolic,dram,mcycdb}_tsv.py and return {} when
the file is missing, so the benchmark runs with whatever is available.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation"))
import score_mcycle as S          # noqa: E402  (ground truth, results dir, TRAP, HOLDOUT)
import compare_kofam as K         # noqa: E402  (raw KofamScan baseline)

BENCH = ROOT / "validation" / "benchmark"
TRAP = S.TRAP
HOLDOUT = S.HOLDOUT

# Step resolution: step -> diagnostic markers (present when any is present).
STEP_DIAG: dict[str, list[str]] = {
    "mcr": ["mcrA"],
    "mtr": ["mtrA"],
    "co2_to_methyl": ["fwdA"],
    "acetoclastic": ["cdhC"],
    "methylotrophic": ["mtaB", "mtmB", "mtbB", "mttB"],
    "heterodisulfide": ["hdrA", "hdrD"],
    "reverse_methanogenesis": ["mcrA_anme"],
    "pmmo": ["pmoA"],
    "smmo": ["mmoX"],
    "methanol_oxidation": ["mxaF", "xoxF"],
    "nod": ["nod"],
    "formaldehyde_h4mpt": ["mtdB"],
    "rump": ["hxlA"],
    "serine_cycle": ["mtkA"],
}


def _targets() -> list[dict]:
    return yaml.safe_load(open(ROOT / "config" / "targets.yaml"))["targets"]


def all_target_ids() -> list[str]:
    return [t["id"] for t in _targets()]


def target_kos() -> dict[str, list[str]]:
    return {t["id"]: list(t.get("ko") or []) for t in _targets()}


def load_truth() -> tuple[dict[tuple[str, str], str], set[str]]:
    truth = S.load_truth()
    return truth, {g for g, _ in truth}


def _with_anme(pred: dict[tuple[str, str], bool], genomes: set[str]) -> dict:
    """Tools without a direction call: mcrA_anme = their mcrA call (prereg)."""
    for g in genomes:
        pred[(g, "mcrA_anme")] = pred.get((g, "mcrA"), False)
    return pred


def _read(path: Path | None) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ── mcycle and raw KofamScan (no extra input) ────────────────────────────────

def load_mcycle(genomes: set[str]) -> dict[tuple[str, str], bool]:
    calls = S.load_calls(genomes)
    return {(g, t): calls[g].get(t) in S.PRESENT for g in calls for t in all_target_ids()}


def load_kofam(genomes: set[str]) -> dict[tuple[str, str], bool]:
    # shared KOs count for every target that lists them, so mcrA_anme = K00399 already
    return K.raw_kofam(sorted(genomes), _targets(), K.stock_thresholds())


# ── METABOLIC v4.0 ───────────────────────────────────────────────────────────
# metabolic.tsv: genome, kind, id.  kind = ko   -> a KO of KEGG_identifier_result
#                                   kind = gene -> a worksheet-1 row reported Present
# Targets with a dedicated gene-named row are decided by that row alone.
METABOLIC_GENE_ROWS = {"pmoA", "pmoB", "pmoC", "mcrA", "mcrB", "mcrC", "mmoB", "mmoD",
                       "mxaF", "fae", "frmA", "cdhD", "cdhE"}


def load_metabolic(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    rows = _read(path)
    if not rows:
        return {}
    kos: dict[str, set[str]] = {}
    genes: dict[str, set[str]] = {}
    for r in rows:
        (kos if r["kind"] == "ko" else genes).setdefault(r["genome"], set()).add(r["id"])
    pred = {}
    for g in genomes:
        for tid, tk in target_kos().items():
            if tid in METABOLIC_GENE_ROWS:
                pred[(g, tid)] = tid in genes.get(g, set())
            else:
                pred[(g, tid)] = any(k in kos.get(g, set()) for k in tk)
    return _with_anme(pred, genomes)


# ── DRAM v1.4.6 ──────────────────────────────────────────────────────────────
# dram.tsv: genome, function, present — the distillate (product.tsv) functions of
# the category "Methanogenesis and methanotrophy". A present function is expanded to
# the targets below; every other target is predicted absent.
DRAM_FUNCTION_MAP: dict[str, list[str]] = {
    "key functional gene": ["mcrA", "mcrB", "mcrG"],
    "acetate => methane, pt 1": ["acs"],
    "acetate => methane, pt 2": ["ackA"],
    "acetate => methane, pt 3": ["pta"],
    "methanol => methane": ["mtaB"],
    "trimethylamine => dimethylamine": ["mttB"],
    "dimethylamine => monomethylamine": ["mtbB"],
    "monomethylamine => ammonia": ["mtmB"],
    "putative but not defining co2 => methane": ["fwdA"],
    "methane => methanol, with oxygen (pmo)": ["pmoA", "pmoB", "pmoC"],
    "methane => methanol, with oxygen (mmo)": ["mmoX", "mmoY", "mmoZ", "mmoB", "mmoC", "mmoD"],
}


def load_dram(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    rows = _read(path)
    if not rows:
        return {}
    hit: dict[str, set[str]] = {}
    for r in rows:
        targets = DRAM_FUNCTION_MAP.get(r["function"].strip().lower())
        if targets and str(r["present"]).lower() in ("1", "true", "yes", "present"):
            hit.setdefault(r["genome"], set()).update(targets)
    pred = {(g, tid): tid in hit.get(g, set()) for g in genomes for tid in all_target_ids()}
    return _with_anme(pred, genomes)


# ── MCycDB 2021 ──────────────────────────────────────────────────────────────
# mcycdb.tsv: genome, family, count. Family -> target of the same name, plus the
# renames below. amoA / amoB / amoC never count as pmo. No family: fwdA, fwdD, sgaA.
MCYC_RENAME: dict[str, list[str]] = {
    "hdrA1": ["hdrA"], "hdrA2": ["hdrA"], "hdrB1": ["hdrB"], "hdrB2": ["hdrB"],
    "hdrC1": ["hdrC"], "hdrC2": ["hdrC"],
    "xoxF1": ["xoxF"], "xoxF2": ["xoxF"], "xoxF4": ["xoxF"], "xoxF5": ["xoxF"],
    "fmdB": ["fwdB"], "fmdC": ["fwdC"], "fae-hps": ["fae", "hxlA"], "fghA": ["frmB"],
    "gck": ["gckA"], "fdoG": ["fdh"], "fdhF": ["fdh"],
}


def load_mcycdb(path: Path | None, genomes: set[str]) -> dict[tuple[str, str], bool]:
    rows = _read(path)
    if not rows:
        return {}
    ids = set(all_target_ids())
    hit: dict[str, set[str]] = {}
    for r in rows:
        if float(r["count"] or 0) <= 0:
            continue
        fam = r["family"]
        for tid in MCYC_RENAME.get(fam, [fam] if fam in ids else []):
            hit.setdefault(r["genome"], set()).add(tid)
    pred = {(g, tid): tid in hit.get(g, set()) for g in genomes for tid in ids}
    return _with_anme(pred, genomes)


# Targets each comparator can express at all (coverage sensitivity, prereg).
def covered_targets() -> dict[str, set[str]]:
    ids = set(all_target_ids())
    with_ko = {t for t, k in target_kos().items() if k}
    dram = {t for ts in DRAM_FUNCTION_MAP.values() for t in ts}
    mcyc = (ids - {"fwdA", "fwdD", "sgaA", "mcrA_anme"})
    fam_file = ROOT / "comparators" / "MCyc" / "id2gene.map"
    if fam_file.exists():
        fams = {l.split("\t")[1] for l in open(fam_file) if "\t" in l}
        mcyc = {t for t in ids if t in fams} | {t for f, ts in MCYC_RENAME.items()
                                               if f in fams for t in ts}
    return {"kofam": with_ko - {"mcrA_anme"}, "metabolic": with_ko - {"mcrA_anme"},
            "dram": dram, "mcycdb": mcyc}


LOADERS = {
    "mcycle": (load_mcycle, False),
    "kofam": (load_kofam, False),
    "metabolic": (load_metabolic, True),
    "dram": (load_dram, True),
    "mcycdb": (load_mcycdb, True),
}
