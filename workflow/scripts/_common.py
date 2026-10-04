"""Shared helpers for mcycle-pipeline scripts."""

from __future__ import annotations

from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "config.yaml"
TARGETS = ROOT / "config" / "targets.yaml"


def load_config() -> dict:
    with open(CONFIG) as fh:
        return yaml.safe_load(fh)


def load_targets() -> dict:
    with open(TARGETS) as fh:
        return yaml.safe_load(fh)


def unique_pfam_ids(targets: list[dict]) -> list[str]:
    seen, out = set(), []
    for t in targets:
        for pf in t.get("pfam") or []:
            if pf not in seen:
                seen.add(pf)
                out.append(pf)
    return out


def blast_fallback_targets(targets: list[dict]) -> list[dict]:
    return [t for t in targets if t.get("blast_fallback")]


def blast_gated_targets(targets: list[dict]) -> list[dict]:
    """Targets routed to the blast_gated_refs DIAMOND DB.

    Two categories collapse here:
    1. No Pfam at all — BLAST is the only signature besides the KO.
    2. Has Pfam but `requires_blast_for_confirmation: true` (the clade gates:
       pmoA, mcrA_anme). Their BLAST refs are the diagnostic anchors that
       decide the clade; they must be in the DIAMOND DB or the gate silently
       disqualifies every hit.
    """
    return [t for t in targets
            if not (t.get("pfam") or [])
            or t.get("requires_blast_for_confirmation")]
