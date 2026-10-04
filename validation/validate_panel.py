#!/usr/bin/env python3
"""
validate_panel.py — identity QC gate for the reference panel.

Ported from ncycle-pipeline/validation/validate_panel.py, where a proteome filed
under the wrong organism invalidated three ground-truth corrections. The methane
panel is genome sequence (.fna), so the checks are made on the contigs:

  1. defline check (offline) — NCBI genomic FASTA deflines name the organism.
     A genome passes when >= 95 % of its contigs carry one of the roster's `taxa`
     tokens. Deflines that name no organism at all (bare accessions) are agnostic
     and are left to check 2.
  2. assembly check (--online) — the NCBI Datasets report for the roster accession
     must name an organism matching the same tokens, and the file must hold exactly
     the assembly's total sequence length.

Exit status: 0 = all pass, 1 = any failure, 2 = --strict and a genome could not
be confirmed by either check.

    python validation/validate_panel.py                 # offline, ../ref_panel
    python validation/validate_panel.py --online        # + NCBI assembly check
    python validation/validate_panel.py --online --write validation/panel_qc.tsv
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PANEL = (ROOT / ".." / "ref_panel").resolve()
DEFAULT_ROSTER = ROOT / "validation" / "panel.tsv"
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"

CONTIG_MATCH_FLOOR = 0.95   # >= 95 % of organism-naming deflines must match


def read_roster(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                   delimiter="\t"))


def scan_fasta(path: Path) -> tuple[list[str], int]:
    """Deflines (without '>') and the total sequence length of a FASTA file."""
    deflines, total = [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                deflines.append(line[1:].rstrip())
            else:
                total += len(line.strip())
    return deflines, total


def has_token(text: str, tokens: list[str]) -> bool:
    low = text.lower()
    return any(re.search(r"(?<![a-z0-9])" + re.escape(t) + r"(?![a-z])", low)
               for t in tokens)


def ncbi_reports(accessions: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for i in range(0, len(accessions), 50):
        url = API + ",".join(accessions[i:i + 50]) + "/dataset_report?page_size=100"
        with urllib.request.urlopen(url, timeout=90) as r:
            for rep in json.load(r).get("reports", []):
                out[rep["accession"]] = rep
        time.sleep(0.4)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Identity QC of the reference panel")
    ap.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    ap.add_argument("--roster", type=Path, default=DEFAULT_ROSTER)
    ap.add_argument("--online", action="store_true",
                    help="also check each file against its NCBI assembly report")
    ap.add_argument("--strict", action="store_true",
                    help="fail when a genome is confirmed by neither check")
    ap.add_argument("--write", type=Path, help="write the QC table to this TSV")
    args = ap.parse_args()

    roster = read_roster(args.roster)
    reports = ncbi_reports([r["source"] for r in roster]) if args.online else {}

    rows, n_fail, n_unknown = [], 0, 0
    for r in roster:
        name, acc = r["name"], r["source"]
        tokens = [t for t in r["taxa"].split(";") if t]
        f = args.panel / f"{name}.{r['input']}"
        row = {"name": name, "accession": acc, "split": r["split"], "n_seqs": "",
               "total_bp": "", "defline_match": "", "ncbi_organism": "",
               "ncbi_bp": "", "status": "", "reason": ""}
        if not f.is_file() or f.stat().st_size == 0:
            row.update(status="fail", reason="file missing or empty")
            rows.append(row); n_fail += 1
            continue

        deflines, total = scan_fasta(f)
        # a defline names an organism when it has text after the sequence accession
        naming = [d for d in deflines if len(d.split(None, 1)) > 1]
        hits = sum(has_token(d, tokens) for d in naming)
        row.update(n_seqs=len(deflines), total_bp=total,
                   defline_match=f"{hits}/{len(naming)}")
        problems, confirmed = [], False
        if naming:
            if hits / len(naming) >= CONTIG_MATCH_FLOOR:
                confirmed = True
            else:
                wrong = next(d for d in naming if not has_token(d, tokens))
                problems.append(f"deflines name another organism, e.g. '{wrong[:70]}'")

        if args.online:
            rep = reports.get(acc)
            if rep is None:
                problems.append("accession not found at NCBI")
            else:
                org = rep["organism"]["organism_name"]
                iso = rep["organism"].get("infraspecific_names", {})
                label = " ".join([org, iso.get("strain", ""), iso.get("isolate", "")])
                ncbi_bp = int(rep["assembly_stats"]["total_sequence_length"])
                row.update(ncbi_organism=org, ncbi_bp=ncbi_bp)
                if not has_token(label, tokens):
                    problems.append(f"NCBI names the assembly '{org}'")
                elif ncbi_bp != total:
                    problems.append(f"length {total} != assembly length {ncbi_bp}")
                else:
                    confirmed = True
                status = rep["assembly_info"].get("assembly_status", "")
                if status not in ("current", ""):
                    problems.append(f"assembly status is '{status}'")

        if problems:
            row.update(status="fail", reason="; ".join(problems)); n_fail += 1
        elif confirmed:
            row.update(status="pass", reason="")
        else:
            row.update(status="unknown",
                       reason="no organism in deflines; run with --online")
            n_unknown += 1
        rows.append(row)

    mark = {"pass": "  PASS", "fail": "x FAIL", "unknown": "? WARN"}
    for row in rows:
        print(f"{mark[row['status']]}  {row['name']:<28s} {row['accession']:<16s} "
              f"{row['n_seqs']!s:>4} seqs {row['total_bp']!s:>9} bp  "
              f"deflines {row['defline_match']:<9s} {row['reason']}")
    extra = sorted(p.name for p in args.panel.glob("*.f[an]a")
                   if p.stem not in {r["name"] for r in roster})
    if extra:
        print(f"\nnot in the roster: {', '.join(extra)}")
    n_pass = sum(r["status"] == "pass" for r in rows)
    print(f"\n[validate_panel] {n_pass} pass / {n_fail} fail / {n_unknown} unknown "
          f"of {len(rows)} ({'online' if args.online else 'offline'})")

    if args.write:
        with open(args.write, "w") as fh:
            fh.write(f"# identity QC of the reference panel, "
                     f"{datetime.date.today().isoformat()} "
                     f"({'defline + NCBI assembly check' if args.online else 'defline check'})\n")
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    if n_fail or extra:
        return 1
    if args.strict and n_unknown:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
