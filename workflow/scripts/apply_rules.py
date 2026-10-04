#!/usr/bin/env python3
"""
apply_rules.py — combine hmmscan + DIAMOND-blastp evidence, emit one row per
target into calls/mcycle_calls.tsv.

Same engine as the nitrogen / sulfur sister pipelines:
  - evidence per protein: hmmscan domtblout (KOfam KO, custom HMM, Pfam) and
    DIAMOND hits against curated UniProt references tagged for a target
    (header convention: >{target_id}||{acc})
  - status set: confirmed | domain-only | narrow-no-IPR | disqualified | absent

Methane-specific rules:
  - join_inframe_stops:    (before the evaluation, nucleotide input only) a gene
                           that Prodigal split at an in-frame stop codon —
                           selenocysteine / pyrrolysine recoding — is scored as one
  - resolve_mcr_direction: methanogenesis vs reverse methanogenesis (ANME)
  - gate_pmo_subunits:     pmoB / pmoC follow the clade-gated pmoA call
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

KO_RE = re.compile(r"^K\d{4,}$")


# ───────────────────────────── parsers ───────────────────────────────────────

def parse_hmmscan_domtbl(path: Path, tc_cutoffs: dict[str, float],
                         qcov_cutoffs: dict[str, float] | None = None,
                         below: dict[str, dict] | None = None,
                        ) -> dict[str, dict]:
    """Return {protein_id: {pfam: {pfam_id: best_evalue},
                             custom: {target_id: best_evalue}}}.

    Pfam profiles are identified by target_acc starting with "PF" (e.g.
    PF13435.10). Custom HMMs (built via build_custom_hmms.py with
    `hmmbuild --name target_id`) carry the target_id in target_name and
    typically have no ACC field (`-`). A target may have several clade models,
    named `target_id__clade` (build_clade_hmms.py), each with its own threshold;
    a protein passing any of them is a custom hit of the target, and
    out[protein]["clade"][target_id] names the best-scoring clade.

    `qcov_cutoffs` is a per-target {target_id: min_query_coverage} dict.
    When set, domains with env-coverage of the query protein below the
    threshold are dropped — prevents small-domain HMMs (e.g. 104-aa hcnA)
    from false-positive-hitting the matching domain of a large multidomain
    protein (~950 aa Fe-S-containing oxidoreductase).

    `below`, when given, is filled with the KOfam hits that miss their
    threshold: {KO: {protein_id: [score, evalue, hmm_from, hmm_to]}} (the
    profile range over all domains) — the input of join_inframe_stops.
    """
    qcov_cutoffs = qcov_cutoffs or {}
    out: dict[str, dict] = defaultdict(
        lambda: {"pfam": {}, "custom": {}, "ko": {}, "clade": {}})
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            target_name = parts[0]
            target_acc = parts[1]                 # "PF#####.NN" or "-"
            qlen = int(parts[5])
            query_name = parts[3]
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            env_from = int(parts[19])
            env_to = int(parts[20])
            if target_acc.startswith("PF"):
                key = target_acc.split(".")[0]
                bucket = "pfam"
            elif KO_RE.match(target_name) or KO_RE.match(target_acc):
                # KOfam profile — name (or acc) is the KO number, e.g. K02588.
                key = target_name if KO_RE.match(target_name) else target_acc
                bucket = "ko"
            else:
                key = target_name
                bucket = "custom"
            tc = tc_cutoffs.get(key)
            # Apply TC cutoff if available, else trust the e-value filter that
            # hmmscan already applied via -E.
            if tc is not None and full_score < tc:
                if below is not None and bucket == "ko":
                    rec = below.setdefault(key, {}).setdefault(
                        query_name, [full_score, full_evalue, 10 ** 9, 0])
                    rec[2] = min(rec[2], int(parts[15]))
                    rec[3] = max(rec[3], int(parts[16]))
                continue
            # Apply per-target query-coverage filter (defends against
            # small-domain HMMs over-calling on multidomain proteins).
            min_qcov = qcov_cutoffs.get(key)
            if min_qcov is not None and qlen > 0:
                qcov = (env_to - env_from + 1) / qlen
                if qcov < min_qcov:
                    continue
            clade = None
            if bucket == "custom" and "__" in key:
                # clade model `<target>__<clade>`: counts for the target, the best
                # clade is remembered for the evidence column
                key, clade = key.split("__", 1)
            entry = out[query_name][bucket]
            prev = entry.get(key)
            if prev is None or full_evalue < prev:
                entry[key] = full_evalue
                if clade:
                    out[query_name]["clade"][key] = clade
    return out


def parse_blast_tsv(path: Path) -> dict[str, dict[str, list[tuple[str, float, float]]]]:
    """Return {protein_id: {target_id: [(uniprot_acc, pident, evalue), …]}}.
    Subject header convention: target_id||uniprot_acc."""
    out: dict[str, dict] = defaultdict(lambda: defaultdict(list))
    if not path.exists() or path.stat().st_size == 0:
        return out
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qid = parts[0]
            sid = parts[1]
            if "||" not in sid:
                continue
            target_id, acc = sid.split("||", 1)
            pident = float(parts[2])
            evalue = float(parts[10])
            out[qid][target_id].append((acc.split(" ", 1)[0], pident, evalue))
    return out


def load_tc_cutoffs(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2]:
                try:
                    out[parts[0]] = float(parts[2])
                except ValueError:
                    pass
    return out


# ───────────────────────────── per-target evaluation ─────────────────────────

def evaluate_target(target: dict,
                    pfam_hits: dict[str, dict[str, float]],
                    custom_hits: dict[str, dict[str, float]],
                    ko_hits: dict[str, dict[str, float]],
                    blast_hits: dict[str, dict],
                    pident_default: float, qcov_default: float,
                    custom_clade: dict[str, dict[str, str]] | None = None,
                    ) -> tuple[str, dict] | None:
    """Return (status, evidence_dict) for the BEST protein matching this target,
    or None if nothing matches.

    Status set:
      confirmed     — sequence signature met (Pfam OR custom HMM) AND BLAST hit
                       (or signature-only when no fallback is configured)
      domain-only   — signature met, no BLAST hit (only meaningful when fallback set)
      narrow-no-IPR — BLAST hit, signature not met
      disqualified  — requires_blast_for_confirmation set, signature met, BLAST missing

    `clade_hmm_required` targets (pmoA, mcrA_anme) are decided by their clade HMMs
    alone: a protein passing a clade model is `confirmed`; a protein with only the
    family signature (KO) is `disqualified` — the family is there, the clade is
    not. The BLAST seeds are then reported as the nearest reference, not required.
    """
    tid = target["id"]
    is_custom = bool(target.get("custom_hmm"))
    req_pfams = set(target.get("pfam") or [])
    target_kos = set(target.get("ko") or [])
    logic = (target.get("pfam_logic") or "any").lower()   # any | all
    fallback = bool(target.get("blast_fallback"))
    requires_blast_for_confirmation = bool(
        target.get("requires_blast_for_confirmation"))
    pident_min = float(target.get("blast_identity_min")
                       or pident_default)
    clade_gate = bool(target.get("clade_hmm_required"))
    custom_clade = custom_clade or {}

    # Collect per-protein evidence summaries.
    candidates = []
    proteins = (set(pfam_hits.keys()) | set(custom_hits.keys())
                | set(ko_hits.keys()) | set(blast_hits.keys()))
    for prot in proteins:
        p_hits = pfam_hits.get(prot, {})
        c_hits = custom_hits.get(prot, {})
        k_hits = ko_hits.get(prot, {})
        b_hits_all = blast_hits.get(prot, {}).get(tid, [])
        b_hits = [(a, pi, ev) for (a, pi, ev) in b_hits_all
                  if pi >= pident_min]
        # Sequence-signature check. Precedence: custom HMM > KOfam KO > Pfam.
        # The custom-HMM tier is preferred when a clade model has been built
        # (hardening phase); until then KO is the primary signature, so a
        # `custom_hmm: true` target with no built model still resolves via KO.
        sig_source = None
        sig_evalue = None
        present_pfams: set[str] = set()
        present_kos: set[str] = set()
        if is_custom and tid in c_hits:
            sig_source = "custom-hmm"
            sig_evalue = c_hits[tid]
        if sig_source is None and target_kos:
            present_kos = target_kos & set(k_hits.keys())
            if present_kos:
                sig_source = "ko"
                sig_evalue = min(k_hits[k] for k in present_kos)
        if sig_source is None and req_pfams:
            present_pfams = req_pfams & set(p_hits.keys())
            ok = ((logic == "all" and present_pfams == req_pfams)
                  or (logic == "any" and bool(present_pfams)))
            if ok:
                sig_source = "pfam"
                sig_evalue = min(p_hits[p] for p in present_pfams)
        sig_ok = sig_source is not None
        blast_ok = bool(b_hits)
        if not (sig_ok or blast_ok):
            continue
        # Status.
        if clade_gate:
            if not sig_ok:
                continue
            status = "confirmed" if sig_source == "custom-hmm" else "disqualified"
        elif requires_blast_for_confirmation and sig_ok and not blast_ok:
            status = "disqualified"
        elif sig_ok and (fallback and blast_ok or not fallback):
            status = "confirmed"
        elif sig_ok and fallback and not blast_ok:
            status = "domain-only"
        elif not sig_ok and blast_ok:
            # BLAST-only hit with NO sequence signature. For any target that has a
            # primary signature tier (custom HMM, KOfam KO, or Pfam), a BLAST-only hit
            # means that signature was rejected → cross-reactivity to a generic
            # homologue (often a UniProt entry auto-fetched by a shared gene name),
            # not a real call (S6, 2026-05-25 hardening). Only true Tier-3 targets
            # (no custom HMM, no KO, no Pfam — e.g. otr) legitimately call on BLAST
            # alone, where narrow-no-IPR remains a meaningful "BLAST says yes" status.
            if is_custom or target_kos or req_pfams:
                continue
            status = "narrow-no-IPR"
        else:
            continue
        evidence_source = sig_source if sig_ok else "blast"
        best_b = min(b_hits, key=lambda x: x[2]) if b_hits else None
        candidates.append({
            "protein_id": prot,
            "status": status,
            "evidence_source": evidence_source,
            "pfam_hits": (sorted(present_kos) if sig_source == "ko"
                          else sorted(present_pfams) if sig_source == "pfam"
                          else [f"{tid}__{custom_clade[prot][tid]}"
                                if tid in custom_clade.get(prot, {}) else tid]
                          if sig_source == "custom-hmm"
                          else []),
            "pfam_evalue": sig_evalue,
            "blast_acc": best_b[0] if best_b else "",
            "blast_pident": best_b[1] if best_b else "",
            "blast_evalue": best_b[2] if best_b else "",
        })
    if not candidates:
        return None
    # Prefer confirmed > domain-only > narrow-no-IPR > disqualified.
    rank = {"confirmed": 0, "domain-only": 1,
            "narrow-no-IPR": 2, "disqualified": 3}
    # Tie-breakers after (status, e-value): higher BLAST identity, then protein
    # id in natural order. Paralog copies often tie on the first two, and
    # `proteins` is a set, so without these the reported copy changed from run
    # to run. Status is unaffected either way.
    ranked = sorted(candidates,
                    key=lambda c: (rank[c["status"]], _evalue_key(c),
                                   -float(c["blast_pident"] or 0),
                                   _natural_key(c["protein_id"])))
    best = ranked[0]
    # Other proteins that reach the same status for this target (paralog
    # copies). Reported, not scored.
    best["other_copies"] = sorted((c["protein_id"] for c in ranked[1:]
                                   if c["status"] == best["status"]),
                                  key=_natural_key)
    return best["status"], best


# ───────────────────────────── read-mode shortcut ────────────────────────────

def evaluate_target_from_reads(target: dict,
                               read_presence: dict[str, dict]) -> tuple[str, dict] | None:
    """For Mode B, presence is per-target boolean from parse_coverage.py."""
    tid = target["id"]
    pr = read_presence.get(tid)
    if pr is None:
        return None
    if not pr.get("present"):
        return None
    status = "confirmed"   # if breadth + depth + identity met, treat as confirmed
    return status, {
        "protein_id": "",
        "status": status,
        "pfam_hits": [],
        "pfam_evalue": "",
        "blast_acc": pr.get("source_acc", ""),
        "blast_pident": pr.get("identity", ""),
        "blast_evalue": "",
        "breadth_pct": pr.get("breadth_pct"),
        "depth_median": pr.get("depth_median"),
    }


def load_read_presence(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not path or not Path(path).exists():
        return out
    import csv
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            tid = row["target_id"]
            out[tid] = {
                "present": row["present"] in ("True", "true", "1"),
                "source_acc": row.get("source_acc", ""),
                "identity": float(row["identity_pct"] or 0.0),
                "breadth_pct": float(row["breadth_pct"] or 0.0),
                "depth_median": float(row["depth_median"] or 0.0),
            }
    return out


def _evalue_key(c: dict) -> float:
    """E-value used to rank candidates: the signature (HMM) e-value, else the
    BLAST one. An e-value of exactly 0.0 is the BEST possible score, so it must
    not be treated as missing (the old `a or b or 1e9` chain ranked it last)."""
    for v in (c["pfam_evalue"], c["blast_evalue"]):
        if v not in (None, ""):
            return float(v)
    return 1e9


def _natural_key(s: str) -> list:
    """Sort key that orders `contig_2` before `contig_10`."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def parse_prodigal_strands(faa_path: Path | None) -> dict[str, tuple[str, int, int, str]]:
    """{protein_id: (contig, start, end, strand)} from Prodigal .faa headers, for
    the gene-coordinate columns of mcycle_calls.tsv. {} for a pre-called proteome."""
    out: dict[str, tuple[str, int, int, str]] = {}
    if not faa_path or not Path(faa_path).exists():
        return out
    with open(faa_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = [p.strip() for p in line[1:].split("#")]
            if len(parts) < 4:
                continue
            protid = parts[0].split()[0]
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            strand = "-" if parts[3] == "-1" else "+"
            out[protid] = (protid.rsplit("_", 1)[0], start, end, strand)
    return out


# ───────────────────────────── methane-specific rules ────────────────────────

MAX_STOP_GAP_NT = 150      # between the two halves of a gene split at an in-frame stop
MAX_HMM_OVERLAP = 10       # profile positions the two halves may share


def join_inframe_stops(below: dict[str, dict], tc_cutoffs: dict[str, float],
                       loci: dict[str, tuple[str, int, int, str]]
                       ) -> dict[tuple[str, str], tuple[float, str]]:
    """Score as ONE gene the two halves of a gene that Prodigal split at an
    in-frame stop codon.

    Methanogens read UGA as selenocysteine (formate dehydrogenase, hydrogenases,
    Fwd, Hdr) and UAG as pyrrolysine (the methylamine methyltransferases). A gene
    caller stops there, so the gene arrives as two neighbouring ORFs, each
    scoring under the KOfam threshold. Two ORFs are joined only if they look like
    exactly that:
      - neighbours on one contig and strand (consecutive Prodigal gene numbers);
      - the gap between them is 0-150 nt and a multiple of 3, i.e. the second ORF
        continues in the reading frame of the first (a frameshifted or genuinely
        two-gene arrangement is out of frame or overlapping);
      - both hit the same KO profile below its threshold, the upstream half on the
        N-terminal part of the profile and the downstream half on what follows
        (ranges in order, sharing at most 10 positions) — tandem paralogues hit
        the same range twice and are not joined;
      - the two scores together reach the threshold.
    Returns {(protein_id, KO): (evalue, partner_id)} for the higher-scoring half,
    which then stands for the gene. Sequence alone cannot tell recoding from a
    nonsense mutation, so every such call is tagged `joined:<partner>` in the
    evidence_source. Needs gene coordinates: {} for a pre-called proteome."""
    out: dict[tuple[str, str], tuple[float, str]] = {}
    if not loci:
        return out
    for ko, hits in below.items():
        tc = tc_cutoffs.get(ko)
        if tc is None:
            continue
        for prot, (score, evalue, h_from, h_to) in hits.items():
            contig, _, num = prot.rpartition("_")
            if not num.isdigit():
                continue
            nxt = f"{contig}_{int(num) + 1}"
            if nxt not in hits or prot not in loci or nxt not in loci:
                continue
            (_, s1, e1, strand1), (_, s2, e2, strand2) = loci[prot], loci[nxt]
            gap = s2 - e1 - 1
            if strand1 != strand2 or not 0 <= gap <= MAX_STOP_GAP_NT or gap % 3:
                continue
            # `prot` lies at the lower coordinates: upstream on +, downstream on -.
            up, down = (prot, nxt) if strand1 == "+" else (nxt, prot)
            if hits[up][3] - hits[down][2] > MAX_HMM_OVERLAP or hits[up][2] >= hits[down][2]:
                continue
            if score + hits[nxt][0] < tc:
                continue
            best, other = (prot, nxt) if score >= hits[nxt][0] else (nxt, prot)
            out[(best, ko)] = (hits[best][1], other)
    return out


PRESENT_STATUSES = ("confirmed", "domain-only", "narrow-no-IPR")


def _present(rows_by_id: dict[str, dict], tid: str) -> bool:
    r = rows_by_id.get(tid)
    return r is not None and r["status"] in PRESENT_STATUSES


def resolve_mcr_direction(rows_by_id: dict[str, dict]) -> str | None:
    """Call the metabolic DIRECTION of an Mcr-carrying genome.

    Methyl-coenzyme M reductase is the SAME enzyme in methanogens (methyl-CoM →
    CH4) and in anaerobic methanotrophic archaea (ANME; CH4 → methyl-CoM, the
    whole methanogenesis pathway run backwards). Gene content does not separate
    the two — ANME-2 carry the complete pathway — but ANME McrA form their own
    clades, so the call is made on the McrA sequence itself:
      reverse      — the `mcrA_anme` target is confirmed, i.e. the genome's McrA
                     is >= the gate identity to an ANME-clade seed.
      methanogenic — McrA present, not an ANME clade.
    Returns 'methanogenic' | 'reverse', or None if mcrA is absent.

    Limits (see Info-methane.md): only ANME-2d (Ca. Methanoperedens) is seeded,
    so marine ANME-1 / -2a / -2c / -3 are still called methanogenic; and the
    call is about the clade, not the physiology of the cell on the day."""
    if not _present(rows_by_id, "mcrA"):
        return None
    return "reverse" if _present(rows_by_id, "mcrA_anme") else "methanogenic"


def gate_pmo_subunits(rows_by_id: dict[str, dict]) -> list[str]:
    """pmoB / pmoC share their KOs (K10945 / K10946) verbatim with amoB / amoC of
    ammonia oxidizers, so on their own they say "copper membrane monooxygenase",
    not "methane". The clade is decided on subunit A (pmoA, BLAST-gated to
    methanotroph seeds): unless pmoA is confirmed in the same genome, pmoB / pmoC
    are demoted to `disqualified`. Returns the ids that were demoted."""
    if (rows_by_id.get("pmoA") or {}).get("status") == "confirmed":
        return []
    demoted = []
    for tid in ("pmoB", "pmoC"):
        r = rows_by_id.get(tid)
        if r and r["status"] in PRESENT_STATUSES:
            r["status"] = "disqualified"
            r["evidence_source"] = f"{r['evidence_source'] or 'ko'}|no_pmoA"
            demoted.append(tid)
    return demoted


def _demote(rows_by_id: dict[str, dict], tids: tuple[str, ...], tag: str) -> list[str]:
    demoted = []
    for tid in tids:
        r = rows_by_id.get(tid)
        if r and r["status"] in PRESENT_STATUSES:
            r["status"] = "disqualified"
            r["evidence_source"] = f"{r['evidence_source'] or 'ko'}|{tag}"
            demoted.append(tid)
    return demoted


ACR_MIN_FRACTION = 0.5     # of the K00399 threshold: "an McrA homologue is there"


def gate_mcr_subunits(rows_by_id: dict[str, dict], below: dict[str, dict],
                      tc_cutoffs: dict[str, float]) -> list[str]:
    """Alkyl-coenzyme M reductases of alkane-oxidizing archaea are Mcr homologues.
    Their alpha subunit fails K00399 (Ca. Ethanoperedens EcrA: 642, threshold 775.5)
    but the beta subunit passes K00401 (586, threshold 504.3). mcrB / mcrG are
    therefore demoted to `disqualified` when the genome has no McrA call AND
    carries an McrA homologue under the threshold (>= half of it) — the signature
    of an alkyl-CoM reductase. A genome with mcrB / mcrG and no McrA homologue at
    all (a MAG that lost the contig) keeps its calls."""
    if _present(rows_by_id, "mcrA"):
        return []
    tc = tc_cutoffs.get("K00399")
    homologue = tc is not None and any(
        rec[0] >= ACR_MIN_FRACTION * tc for rec in below.get("K00399", {}).values())
    return _demote(rows_by_id, ("mcrB", "mcrG"), "no_mcrA_acr_like") if homologue else []


def gate_mmo_subunits(rows_by_id: dict[str, dict]) -> list[str]:
    """The reductase, regulatory and assembly components of soluble methane
    monooxygenase (and the beta / gamma chains of its hydroxylase) have close
    relatives in every other soluble di-iron monooxygenase and among unrelated
    oxidoreductases (a lone FNR-type reductase of Methylocystis sp. SC2 scores 410
    on K16161, threshold 407.7). What makes them sMMO is the hydroxylase alpha
    chain: without mmoX in the genome, mmoY / mmoZ / mmoB / mmoC / mmoD are demoted
    to `disqualified`."""
    if _present(rows_by_id, "mmoX"):
        return []
    return _demote(rows_by_id, ("mmoY", "mmoZ", "mmoB", "mmoC", "mmoD"), "no_mmoX")


FDHA_PARTNER_FRACTION = 0.70   # of the K22516 threshold, when fdhB is in the genome


def resolve_fdh(rows_by_id: dict[str, dict], below: dict[str, dict],
                tc_cutoffs: dict[str, float], ko_hits: dict[str, dict],
                joined: dict) -> list[str]:
    """Separate the F420-dependent formate dehydrogenase of methanogens (fdhA,
    K22516) from the NAD / quinone-linked one (fdh, K00123).

    The two alpha subunits are one family and the profiles overlap: archaeal FdhA
    score 0.72-0.95 of the K22516 threshold while passing K00123, bacterial FdhA
    score 0.76-0.86 on K22516. The F420-binding beta subunit FdhB (K00125) exists
    only in the F420-dependent enzyme. So, in a genome WITH fdhB:
      - an FdhA-family protein scoring >= 0.70 of the K22516 threshold is fdhA
        (evidence `ko|fdhB_partner` when it is under the threshold itself);
      - that protein is not also the generic `fdh`: if every protein behind the
        `fdh` call is such an FdhA, `fdh` is demoted to `disqualified`.
    Without fdhB nothing changes. Returns a list of notes for the log."""
    if not _present(rows_by_id, "fdhB"):
        return []
    tc = tc_cutoffs.get("K22516")
    if tc is None:
        return []
    notes = []
    family = {p for p, rec in below.get("K22516", {}).items()
              if rec[0] >= FDHA_PARTNER_FRACTION * tc}
    family |= {p for p, kos in ko_hits.items() if "K22516" in kos}
    family |= {partner for (p, ko), (_e, partner) in joined.items() if ko == "K22516"}
    r = rows_by_id.get("fdhA")
    if r is not None and r["status"] not in PRESENT_STATUSES and family:
        best = max((p for p in family if p in below.get("K22516", {})),
                   key=lambda p: below["K22516"][p][0], default=None)
        if best is not None:
            r.update(status="confirmed", evidence_source="ko|fdhB_partner",
                     protein_id=best, pfam_hits="K22516",
                     best_pfam_evalue=below["K22516"][best][1],
                     other_copies=";".join(sorted((family - {best}) & set(below["K22516"]),
                                                  key=_natural_key)))
            notes.append(f"fdhA called on {best} (K22516 "
                         f"{below['K22516'][best][0]:.0f} of {tc:.0f}, fdhB present)")
    g = rows_by_id.get("fdh")
    if g is not None and g["status"] in PRESENT_STATUSES and _present(rows_by_id, "fdhA"):
        behind = {g["protein_id"], *filter(None, g.get("other_copies", "").split(";"))}
        if behind <= family:
            g["status"] = "disqualified"
            g["evidence_source"] = f"{g['evidence_source'] or 'ko'}|f420_fdhA"
            notes.append("fdh disqualified (its protein is the F420-dependent FdhA)")
    return notes


# ───────────────────────────── main ──────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", choices=["protein", "read"], required=True)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--hmm", type=Path)
    ap.add_argument("--blast-unstable", type=Path)
    ap.add_argument("--blast-gated", type=Path)
    ap.add_argument("--read-presence", type=Path,
                    help="TSV from parse_coverage.py (Mode B)")
    ap.add_argument("--tc-cutoffs", type=Path,
                    default=Path("resources/hmm/tc_cutoffs.tsv"))
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pident-default", type=float, default=30.0)
    ap.add_argument("--qcov-default", type=float, default=50.0)
    ap.add_argument("--gene-coords", type=Path,
                    help="proteome .faa; if it carries Prodigal coordinate headers "
                         "(nucleotide/MAG input), the contig / start / end / strand "
                         "of every called gene are reported. Empty for pre-called "
                         "proteomes.")
    args = ap.parse_args()

    with open(args.targets) as fh:
        targets = yaml.safe_load(fh)["targets"]

    if args.mode == "protein":
        tc = load_tc_cutoffs(args.tc_cutoffs)
        # Per-target hmm_min_qcov override from targets.yaml — defends small
        # domain HMMs (e.g. 104-aa hcnA) from hitting domains of large
        # multidomain proteins.
        qcov_cutoffs = {key: float(t["hmm_min_qcov"])
                        for t in targets if t.get("hmm_min_qcov") is not None
                        for key in [t["id"], *(t.get("ko") or [])]}
        below: dict[str, dict] = {}
        hmm_parsed = (parse_hmmscan_domtbl(args.hmm, tc, qcov_cutoffs, below)
                      if args.hmm else {})
        # Genes split at an in-frame stop (selenocysteine / pyrrolysine) count once.
        loci = parse_prodigal_strands(args.gene_coords)
        joined = join_inframe_stops(below, tc, loci)
        for (prot, ko), (evalue, partner) in joined.items():
            hmm_parsed[prot]["ko"][ko] = evalue
            print(f"[apply_rules] {args.sample}: {ko} joined across an in-frame stop: "
                  f"{prot} + {partner}", file=sys.stderr)
        blast_hits = {}
        for p in (args.blast_unstable, args.blast_gated):
            if p and p.exists():
                for prot, td in parse_blast_tsv(p).items():
                    blast_hits.setdefault(prot, {}).update(td)
        # Split per signature tier: pfam_hits[prot]={pfam_id: evalue};
        # custom_hits[prot]={tid: evalue}; ko_hits[prot]={KO: evalue}.
        pfam_hits_flat = {p: d["pfam"] for p, d in hmm_parsed.items()}
        custom_hits_flat = {p: d["custom"] for p, d in hmm_parsed.items()}
        ko_hits_flat = {p: d["ko"] for p, d in hmm_parsed.items()}
        clade_flat = {p: d["clade"] for p, d in hmm_parsed.items() if d.get("clade")}

        rows = []
        for t in targets:
            res = evaluate_target(t, pfam_hits_flat, custom_hits_flat,
                                  ko_hits_flat, blast_hits,
                                  args.pident_default, args.qcov_default,
                                  custom_clade=clade_flat)
            if res is None:
                rows.append({
                    "target_id": t["id"], "name": t["name"],
                    "category": t["category"], "complex": t.get("complex", ""),
                    "status": "absent", "evidence_source": "",
                    "protein_id": "", "pfam_hits": "",
                    "best_pfam_evalue": "", "blast_acc": "",
                    "blast_pident": "", "blast_evalue": "",
                    "other_copies": "",
                })
                continue
            status, ev = res
            rows.append({
                "target_id": t["id"], "name": t["name"],
                "category": t["category"], "complex": t.get("complex", ""),
                "status": status,
                "evidence_source": ev.get("evidence_source", ""),
                "protein_id": ev["protein_id"],
                "pfam_hits": ",".join(ev["pfam_hits"]),
                "best_pfam_evalue": ev["pfam_evalue"],
                "blast_acc": ev["blast_acc"],
                "blast_pident": ev["blast_pident"],
                "blast_evalue": ev["blast_evalue"],
                "other_copies": ";".join(ev.get("other_copies", [])),
            })

        by_id = {r["target_id"]: r for r in rows}

        # Mark the calls that rest on two ORFs joined across an in-frame stop.
        for r in rows:
            for ko in r["pfam_hits"].split(","):
                if (r["protein_id"], ko) in joined:
                    r["evidence_source"] += f"|joined:{joined[(r['protein_id'], ko)][1]}"
                    break

        # pmoB / pmoC follow the clade-gated pmoA call (shared KOs with amoB / amoC).
        demoted = gate_pmo_subunits(by_id)
        if demoted:
            print(f"[apply_rules] {args.sample}: {', '.join(demoted)} disqualified "
                  "(no confirmed pmoA — amoB/amoC-type copper monooxygenase)",
                  file=sys.stderr)

        # Subunits that mean something only next to the subunit that defines the
        # enzyme (training-panel hardening, 2026-10-04).
        for what, ids in (("alkyl-CoM reductase subunits", gate_mcr_subunits(by_id, below, tc)),
                          ("no mmoX — not a soluble methane monooxygenase",
                           gate_mmo_subunits(by_id))):
            if ids:
                print(f"[apply_rules] {args.sample}: {', '.join(ids)} disqualified ({what})",
                      file=sys.stderr)
        for note in resolve_fdh(by_id, below, tc, ko_hits_flat, joined):
            print(f"[apply_rules] {args.sample}: {note}", file=sys.stderr)

        # Mcr direction call: methanogenesis vs reverse methanogenesis (ANME).
        # Annotates the mcrA / mcrB / mcrG evidence_source (e.g. "ko|mcr_reverse");
        # the methanogenesis vs reverse_methanogenesis_aom synergies and the cycle
        # map read the tag back.
        direction = resolve_mcr_direction(by_id)
        if direction:
            for tid in ("mcrA", "mcrB", "mcrG"):
                r = by_id.get(tid)
                if r and r["status"] in PRESENT_STATUSES:
                    r["evidence_source"] = f"{r['evidence_source'] or 'blast'}|mcr_{direction}"
            print(f"[apply_rules] {args.sample}: Mcr direction = {direction}",
                  file=sys.stderr)

        # Gene coordinates of each called protein (nucleotide/MAG input only; the
        # columns stay empty for a pre-called proteome), then the other copies.
        # Appended last so existing readers of the first 12 columns are unaffected.
        for r in rows:
            contig, start, end, strand = loci.get(r["protein_id"], ("", "", "", ""))
            copies = r.pop("other_copies", "")
            r.update(contig=contig, start=start, end=end, strand=strand,
                     other_copies=copies)

    elif args.mode == "read":
        rp = load_read_presence(args.read_presence)
        rows = []
        for t in targets:
            res = evaluate_target_from_reads(t, rp)
            if res is None:
                rows.append({
                    "target_id": t["id"], "name": t["name"],
                    "category": t["category"], "complex": t.get("complex", ""),
                    "status": "absent", "evidence_source": "",
                    "protein_id": "", "pfam_hits": "",
                    "best_pfam_evalue": "", "blast_acc": "",
                    "blast_pident": "", "blast_evalue": "",
                    "breadth_pct": rp.get(t["id"], {}).get("breadth_pct", ""),
                    "depth_median": rp.get(t["id"], {}).get("depth_median", ""),
                })
                continue
            _, ev = res
            rows.append({
                "target_id": t["id"], "name": t["name"],
                "category": t["category"], "complex": t.get("complex", ""),
                "status": ev["status"],
                "evidence_source": "read-coverage",
                "protein_id": "", "pfam_hits": "",
                "best_pfam_evalue": "",
                "blast_acc": ev["blast_acc"],
                "blast_pident": ev["blast_pident"], "blast_evalue": "",
                "breadth_pct": ev.get("breadth_pct", ""),
                "depth_median": ev.get("depth_median", ""),
            })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    n_status = defaultdict(int)
    for r in rows:
        n_status[r["status"]] += 1
    print(f"[apply_rules] {args.sample}: " +
          ", ".join(f"{k}={v}" for k, v in sorted(n_status.items())),
          file=sys.stderr)


if __name__ == "__main__":
    main()
