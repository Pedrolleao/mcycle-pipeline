#!/usr/bin/env python3
"""
compute_complex_completeness.py — given calls/mcycle_calls.tsv and the
complexes/synergies blocks of targets.yaml, emit:

  calls/complex_completeness.tsv
      complex_id  completeness  status  members_present  members_missing  application

  calls/synergy_completeness.tsv
      synergy_id  completeness  status  requires_present  requires_missing
      forbids_satisfied  forbids_violated  n_present  n_total  benefit

A module is a list of slots. Every gene under `requires:` is one slot. Every
entry under `requires_any:` is one slot too, satisfied by any ONE of its
alternatives; an alternative is a gene, or a list of genes that must all be
present (e.g. `[[ackA, pta], acs]` = "ack + pta, or acs"). A satisfied slot is
written to requires_present as the alternative that satisfied it ("ackA+pta");
an unsatisfied one to requires_missing as the whole choice ("ackA+pta|acs").

The methanogenesis modules and the reverse-methanogenesis (ANME) module need the
same Mcr, so they are told apart by the Mcr direction that apply_rules.py tags on
the mcrA / mcrB / mcrG evidence_source (e.g. `ko|mcr_reverse`): when a direction
is called, it counts as one extra slot, and a module that needs the OPPOSITE
direction is ruled out (status absent, the direction listed under
forbids_violated). Without Mcr there is no call and the modules are scored on
gene presence alone.

A module may also name `key:` entries — what defines it (e.g. `nod` for
nitrite-dependent methane oxidation; Mcr plus a substrate methyltransferase
for methylotrophic methanogenesis). Every entry must be present; an entry that
is a list is satisfied by any one of its genes. If an entry is not met the
module is ruled out the same way (status absent, "no key gene (a|b)" under
forbids_violated), so that shared genes alone (Mcr, pMMO, ...) never make a
module look half-present.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import yaml

PRESENT_STATUSES = {"confirmed", "domain-only", "narrow-no-IPR"}

# Module → the Mcr direction it needs. (Kept here rather than in targets.yaml:
# it mirrors resolve_mcr_direction in apply_rules.py, not a marker definition.)
MCR_DIRECTION_OF = {
    "hydrogenotrophic_methanogenesis": "methanogenic",
    "acetoclastic_methanogenesis":     "methanogenic",
    "methylotrophic_methanogenesis":   "methanogenic",
    "reverse_methanogenesis_aom":      "reverse",
}


def mcr_direction(path: Path) -> str | None:
    """'methanogenic' | 'reverse' when apply_rules.py made a direction call
    for this genome, else None (no Mcr)."""
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["target_id"] in ("mcrA", "mcrB", "mcrG"):
                m = re.search(r"mcr_(methanogenic|reverse)\b",
                              row.get("evidence_source", ""))
                if m:
                    return m.group(1)
    return None


def _alt_genes(alt) -> list[str]:
    return [alt] if isinstance(alt, str) else list(alt)


def score_slots(synergy: dict, calls: dict[str, str]) -> tuple[list[str], list[str]]:
    """(present, missing) slot labels of one module — see the module docstring."""
    def has(gene: str) -> bool:
        return calls.get(gene) in PRESENT_STATUSES
    present = [g for g in synergy.get("requires") or [] if has(g)]
    missing = [g for g in synergy.get("requires") or [] if not has(g)]
    for slot in synergy.get("requires_any") or []:
        hit = next((a for a in slot if all(has(g) for g in _alt_genes(a))), None)
        if hit is not None:
            present.append("+".join(_alt_genes(hit)))
        else:
            missing.append("|".join("+".join(_alt_genes(a)) for a in slot))
    return present, missing


def load_calls(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["target_id"]] = row["status"]
    return out


def status_for(completeness: float) -> str:
    if completeness >= 0.999:
        return "complete"
    if completeness >= 0.5:
        return "partial"
    return "absent"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path)
    ap.add_argument("--complexes-out", required=True, type=Path)
    ap.add_argument("--synergies-out", required=True, type=Path)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.targets))
    complexes = cfg.get("complexes", {})
    synergies = cfg.get("synergies", [])
    calls = load_calls(args.calls)
    direction = mcr_direction(args.calls)

    # Complexes.
    args.complexes_out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.complexes_out, "w") as fh:
        fh.write("complex_id\tcompleteness\tstatus\tmembers_present\t"
                 "members_missing\tn_present\tn_total\tapplication\n")
        for cid, body in complexes.items():
            members = body["members"]
            present = [m for m in members
                       if calls.get(m) in PRESENT_STATUSES]
            missing = [m for m in members if m not in present]
            comp = len(present) / len(members) if members else 0.0
            fh.write("\t".join([
                cid, f"{comp:.3f}", status_for(comp),
                ",".join(present), ",".join(missing),
                str(len(present)), str(len(members)),
                body.get("application", ""),
            ]) + "\n")

    # Synergies.
    with open(args.synergies_out, "w") as fh:
        fh.write("synergy_id\tcompleteness\tstatus\trequires_present\t"
                 "requires_missing\tforbids_satisfied\tforbids_violated\t"
                 "n_present\tn_total\tbenefit\n")
        for s in synergies:
            present, missing = score_slots(s, calls)
            # Direction rule (Mcr modules only, and only when a direction is called).
            satisfied, violated = [], []
            want = MCR_DIRECTION_OF.get(s["name"])
            if want and direction:
                (satisfied if direction == want else violated).append(
                    f"mcr={direction}")
            # Key-gene rule: none of the defining genes present → ruled out.
            unmet = ["|".join(_alt_genes(k)) for k in s.get("key") or []
                     if not any(calls.get(g) in PRESENT_STATUSES
                                for g in _alt_genes(k))]
            if unmet:
                violated.append(f"no key gene ({'; '.join(unmet)})")
            n_present = len(present) + len(satisfied)
            n_total = len(present) + len(missing) + len(satisfied) + len(violated)
            comp = n_present / n_total if n_total else 0.0
            status = "absent" if violated else status_for(comp)
            fh.write("\t".join([
                s["name"], f"{comp:.3f}", status,
                ",".join(present), ",".join(missing),
                ",".join(satisfied), ",".join(violated),
                str(n_present), str(n_total),
                s.get("benefit", ""),
            ]) + "\n")


if __name__ == "__main__":
    main()
