"""The methane-cycle diagram: compounds, reaction steps, and how a genome's
calls light the steps up.

Single source of truth for the static cycle maps (make_cycle_map.py) and
the interactive one in report.html (make_html_report.py ships the geometry
computed here as JSON, so the browser only draws polylines).

A step is `complete` when every gene of at least one route is present
(confirmed or domain-only — the same rule as complex completeness), `partial`
when some route gene is present but no route is whole, otherwise `absent`.
This is a drawing-level summary; the authoritative per-process calls remain
complex_completeness.tsv and synergy_completeness.tsv.

The diagram is a loop around a shared ladder of H₄MPT-bound C1 intermediates.
Methanogens run it downwards (CO₂ → CH₄, right-hand side); aerobic
methanotrophs and methylotrophs climb the other side (CH₄ → CH₃OH → HCHO →
formate → CO₂). Two things depend on the genome:

* the Mcr direction called by apply_rules.py — in a `reverse` (ANME) genome
  the whole methanogenesis chain is drawn pointing from CH₄ to CO₂;
* which organism type owns the shared fwdABC / ftr / mch genes — in a
  bacterial methylotroph (fae present, no Mcr) they are the Fhc complex and
  Mch of formaldehyde oxidation, so they light the formyl-H₄MPT → formate
  step instead of the CO₂ ⇄ formyl-methanofuran steps, and mch points upwards.
"""

from __future__ import annotations

import math
import re

# Canvas (data units, y up).
XLIM = (-3.0, 181.0)
YLIM = (-3.0, 106.0)

# id: (label, x, y, width, height)
NODES = {
    "CO2":   ("CO₂",                92.0, 99.0, 10.0, 7.0),
    "FMFR":  ("formyl-MFR",        128.0, 92.0, 19.0, 7.0),
    "FH4":   ("formyl-H₄MPT",       92.0, 78.0, 22.0, 7.0),
    "MEN":   ("methenyl-H₄MPT",     92.0, 60.0, 25.0, 7.0),
    "MLN":   ("methylene-H₄MPT",    92.0, 42.0, 26.0, 7.0),
    "MH4":   ("methyl-H₄MPT",      130.0, 34.0, 22.0, 7.0),
    "MCOM":  ("methyl-CoM",        130.0, 14.0, 18.0, 7.0),
    "CH4":   ("CH₄",                70.0,  8.0, 10.0, 7.0),
    "MEOH":  ("CH₃OH",              40.0, 26.0, 13.0, 7.0),
    "HCHO":  ("HCHO",               40.0, 48.0, 13.0, 7.0),
    "BIO":   ("biomass",             9.0, 48.0, 16.0, 7.0),
    "FOR":   ("formate",            56.0, 90.0, 14.0, 7.0),
    "AC":    ("acetate",           166.0, 92.0, 13.0, 7.0),
    "ACCOA": ("acetyl-CoA",        166.0, 74.0, 18.0, 7.0),
    "HDS":   ("CoM-S-S-CoB",       164.0, 50.0, 24.0, 7.0),
    "THIOL": ("CoM-SH + CoB-SH",   164.0, 30.0, 28.0, 7.0),
    "MA":    ("methylamines · DMS", 164.0, 12.0, 30.0, 7.0),
}

# Short display names.
DISPLAY: dict[str, str] = {"mcrA_anme": "mcrA(ANME)"}

_K = "mcr_core"
_R = "co2_reduction"
_A = "acetoclastic"
_M = "methylotrophic_methanogenesis"
_E = "energy_conservation"
_O = "aerobic_methane_oxidation"
_F = "formaldehyde_c1"

# arrows: (from, to, bow, (dx0, dy0), (dx1, dy1)) — bow is the offset of the
#   curve's control point to the LEFT of the direction of travel; the two
#   offsets shift the start / end anchor off the node centre (parallel arrows).
# routes: alternative gene sets, any one of which performs the step.
# lines:  how the genes are written next to the arrow (may include accessory
#   genes that no route requires, e.g. mcrC/D, mvhAGD, frhABG).
# label:  (x, y, ha) of the first label line; further lines stack downwards.
# reversible: part of the methanogenesis chain — drawn the other way round in
#   a reverse-Mcr (ANME) genome.
# hide_in / only_in: genome modes (see genome_context) in which the step is
#   switched off / the only modes in which it is evaluated.
STEPS = [
    {"id": "fwd", "name": "CO₂ ⇄ formyl-methanofuran (formylmethanofuran dehydrogenase)",
     "pathway": _R, "reversible": True, "hide_in": {"oxidative"},
     "arrows": [("CO2", "FMFR", 0.0, (0, 0), (0, 0))],
     "routes": [["fwdA", "fwdB", "fwdC"]],
     "lines": [["fwdA", "fwdB"], ["fwdC", "fwdD"]],
     "label": (139.5, 93.8, "left")},
    {"id": "ftr", "name": "Formyl-methanofuran ⇄ formyl-H₄MPT (formyltransferase)",
     "pathway": _R, "reversible": True, "hide_in": {"oxidative"},
     "arrows": [("FMFR", "FH4", 0.0, (0, 0), (0, 0))],
     "routes": [["ftr"]], "lines": [["ftr"]],
     "label": (114.0, 81.5, "left")},
    {"id": "fhc", "name": "Formyl-H₄MPT → formate (Fhc complex: FhcABC + FhcD, the "
                          "bacterial homologues of fwdABC + ftr)",
     "pathway": _F, "only_in": {"oxidative"},
     "arrows": [("FH4", "FOR", 0.0, (0, 0), (0, 0))],
     "routes": [["ftr", "fwdA", "fwdB", "fwdC"]],
     "lines": [["ftr", "fwdA"], ["fwdB", "fwdC"]],
     "label": (71.0, 79.0, "right")},
    {"id": "fdh", "name": "Formate → CO₂ (formate dehydrogenase)", "pathway": _F,
     "arrows": [("FOR", "CO2", 0.0, (0, 0), (0, 0))],
     "routes": [["fdh"], ["fdhA", "fdhB"]],
     "lines": [["fdh", "fdhA", "fdhB"]],
     "label": (68.0, 99.8, "center")},
    {"id": "mch", "name": "Formyl-H₄MPT ⇄ methenyl-H₄MPT (cyclohydrolase)",
     "pathway": _R, "reversible": True, "up_in": {"oxidative"},
     "arrows": [("FH4", "MEN", 0.0, (0, 0), (0, 0))],
     "routes": [["mch"]], "lines": [["mch"]],
     "label": (95.0, 69.0, "left")},
    {"id": "mtd", "name": "Methenyl-H₄MPT ⇄ methylene-H₄MPT (F₄₂₀-dependent Mtd, or H₂-forming Hmd)",
     "pathway": _R, "reversible": True, "hide_in": {"oxidative"},
     "arrows": [("MEN", "MLN", 0.0, (2.4, 0), (2.4, 0))],
     "routes": [["mtd"], ["hmd"]], "lines": [["mtd", "hmd"]],
     "label": (97.5, 51.0, "left")},
    {"id": "mtdB", "name": "Methylene-H₄MPT → methenyl-H₄MPT (NAD(P)-dependent MtdB)",
     "pathway": _F,
     "arrows": [("MLN", "MEN", 0.0, (-2.4, 0), (-2.4, 0))],
     "routes": [["mtdB"]], "lines": [["mtdB"]],
     "label": (86.5, 51.0, "right")},
    {"id": "mer", "name": "Methylene-H₄MPT ⇄ methyl-H₄MPT (Mer; F₄₂₀H₂ from Frh)",
     "pathway": _R, "reversible": True,
     "arrows": [("MLN", "MH4", 0.0, (0, 0), (0, 0))],
     "routes": [["mer"]], "lines": [["mer"], ["frhA", "frhB", "frhG"]],
     "label": (110.0, 47.5, "left")},
    {"id": "mtr", "name": "Methyl-H₄MPT ⇄ methyl-CoM (Mtr, Na⁺-pumping)",
     "pathway": _K, "reversible": True,
     "arrows": [("MH4", "MCOM", 0.0, (0, 0), (0, 0))],
     "routes": [["mtrA", "mtrB", "mtrC", "mtrD", "mtrE", "mtrH"]],   # mtrF / mtrG: see targets.yaml
     "lines": [["mtrA", "mtrB", "mtrC", "mtrD"], ["mtrE", "mtrF", "mtrG", "mtrH"]],
     "label": (126.5, 27.6, "right")},
    {"id": "mcr", "name": "Methyl-CoM ⇄ CH₄ (methyl-coenzyme M reductase)",
     "pathway": _K, "reversible": True,
     "arrows": [("MCOM", "CH4", 0.0, (0, -1.6), (0, 0))],
     "routes": [["mcrA", "mcrB", "mcrG"]],
     "lines": [["mcrA", "mcrB", "mcrG"], ["mcrC", "mcrD", "mcrA_anme"]],
     "label": (104.0, 6.3, "center")},
    {"id": "hdr", "name": "Heterodisulfide → CoM-SH + CoB-SH (heterodisulfide reductase)",
     "pathway": _E,
     "arrows": [("HDS", "THIOL", 0.0, (0, 0), (0, 0))],
     "routes": [["hdrA", "hdrB", "hdrC"], ["hdrD", "hdrE"]],
     "lines": [["hdrA", "hdrB", "hdrC"], ["mvhA", "mvhG", "mvhD"], ["hdrD", "hdrE"]],
     "label": (160.5, 43.6, "right")},
    {"id": "ack", "name": "Acetate → acetyl-CoA (Ack + Pta, or Acs)", "pathway": _A,
     "arrows": [("AC", "ACCOA", 0.0, (0, 0), (0, 0))],
     "routes": [["ackA", "pta"], ["acs"]],
     "lines": [["ackA", "pta"], ["acs"]],
     "label": (162.5, 84.8, "right")},
    {"id": "cdh", "name": "Acetyl-CoA ⇄ methyl-H₄MPT + CO₂ (CODH/ACS complex)",
     "pathway": _A,
     "arrows": [("ACCOA", "MH4", 0.0, (0, 0), (0, 0))],
     "routes": [["cdhA", "cdhB", "cdhC", "cdhD", "cdhE"]],
     "lines": [["cdhA", "cdhB", "cdhC"], ["cdhD", "cdhE"]],
     "label": (143.0, 60.0, "right")},
    {"id": "mta", "name": "Methanol → methyl-CoM (MtaABC)", "pathway": _M,
     "arrows": [("MEOH", "MCOM", 0.0, (0, 0), (0, 2.0))],
     "routes": [["mtaA", "mtaB", "mtaC"]],
     "lines": [["mtaA", "mtaB", "mtaC"]],
     "label": (78.0, 25.6, "center")},
    {"id": "mtm", "name": "Methylamines → methyl-CoM (MtmBC / MtbBC / MttBC + MtbA)",
     "pathway": _M,
     "arrows": [("MA", "MCOM", 0.0, (0, 1.7), (0, 1.7))],
     "routes": [["mtbA", "mtmB", "mtmC"], ["mtbA", "mtbB", "mtbC"],
                ["mtbA", "mttB", "mttC"]],
     "lines": [["mtbA", "mtmB", "mtmC"], ["mtbB", "mtbC", "mttB", "mttC"]],
     "label": (179.0, 23.2, "right")},
    {"id": "mts", "name": "Methyl sulfides → methyl-CoM (MtsAB)", "pathway": _M,
     "arrows": [("MA", "MCOM", 0.0, (0, -1.7), (0, -1.7))],
     "routes": [["mtsA", "mtsB"]], "lines": [["mtsA", "mtsB"]],
     "label": (179.0, 5.0, "right")},
    {"id": "mmo", "name": "CH₄ → CH₃OH (particulate or soluble methane monooxygenase)",
     "pathway": _O,
     "arrows": [("CH4", "MEOH", 0.0, (0, 0), (0, 0))],
     "routes": [["pmoA", "pmoB", "pmoC"], ["mmoX", "mmoY", "mmoZ"]],
     "lines": [["pmoA", "pmoB", "pmoC"], ["mmoX", "mmoY", "mmoZ"],
               ["mmoB", "mmoC", "nod"]],
     "label": (50.0, 10.6, "right")},
    {"id": "mdh", "name": "CH₃OH → HCHO (MxaFI or XoxF methanol dehydrogenase)",
     "pathway": _O,
     "arrows": [("MEOH", "HCHO", 0.0, (0, 0), (0, 0))],
     "routes": [["mxaF", "mxaI"], ["xoxF"]],
     "lines": [["mxaF", "mxaI"], ["xoxF"]],
     "label": (43.5, 38.8, "left")},
    {"id": "fae", "name": "HCHO → methylene-H₄MPT (formaldehyde-activating enzyme)",
     "pathway": _F, "hide_in": {"methanogenic", "reverse"},
     "arrows": [("HCHO", "MLN", 0.0, (0, 0), (0, 0))],
     "routes": [["fae"]], "lines": [["fae"]],
     "label": (61.0, 49.8, "center")},
    {"id": "gsh", "name": "HCHO → formate (glutathione-dependent: Gfa, FrmA, FrmB)",
     "pathway": _F,
     "arrows": [("HCHO", "FOR", 0.0, (3.0, 0), (0, 0))],
     "routes": [["frmA", "frmB"]], "lines": [["gfa"], ["frmA", "frmB"]],
     "label": (53.5, 69.5, "left")},
    {"id": "rump", "name": "HCHO → biomass (ribulose-monophosphate cycle)",
     "pathway": _F, "hide_in": {"methanogenic", "reverse"},
     "arrows": [("HCHO", "BIO", 0.0, (0, 1.7), (0, 1.7))],
     "routes": [["hxlA", "hxlB"]], "lines": [["hxlA", "hxlB"]],
     "label": (25.2, 54.2, "center")},
    {"id": "serine", "name": "C1 → biomass (serine cycle)", "pathway": _F,
     "hide_in": {"methanogenic", "reverse"},
     "arrows": [("HCHO", "BIO", 0.0, (0, -1.7), (0, -1.7))],
     "routes": [["sgaA", "hprA", "mtkA", "mtkB", "mcl"]],
     "lines": [["sgaA", "hprA", "gckA"], ["mtkA", "mtkB", "mcl"]],
     "label": (25.2, 42.0, "center")},
]

# Figure sizing used by make_cycle_map.py (inches / points).
FIG = {"single_w": 11.4, "grid_w": 8.4, "single_fs": 9.0, "grid_fs": 6.6,
       "grid_cols": 2}

LINE_STEP = 3.6          # distance between stacked label lines (data units)
HEAD_LEN, HEAD_W = 2.4, 1.25
NODE_PAD = 1.3           # gap between an arrow end and the node box


# ───────────────────────────── evaluation ────────────────────────────────────

def _slot_genes(slot) -> list[str]:
    return [slot] if isinstance(slot, str) else list(slot)


def _token(slot, codes: dict[str, int], muted: bool = False) -> dict:
    """The gene to print for a slot: the best-supported alternative (first one
    when none is present). `muted` prints it as absent (the step is switched
    off in this genome's mode and the gene is shown on another step)."""
    genes = _slot_genes(slot)
    rank = {2: 0, 1: 1, -1: 2, 0: 3}
    best = min(genes, key=lambda g: (rank[codes.get(g, 0)], genes.index(g)))
    return {"gene": best, "label": DISPLAY.get(best, best),
            "code": 0 if muted else codes.get(best, 0)}


_PRESENT = ("confirmed", "domain-only", "narrow-no-IPR")


def genome_context(calls: dict[str, dict]) -> dict:
    """Two facts the gene list alone does not show.

    mcr_direction — Mcr works in either direction and apply_rules.py calls
        which (ANME-clade McrA ⇒ reverse, otherwise methanogenic), tagging the
        mcrA / mcrB / mcrG evidence_source, e.g. `ko|mcr_reverse`. Read it back.
    mode — who owns the shared H₄MPT genes: `methanogenic` / `reverse` follow
        the Mcr direction; `oxidative` is a bacterial methylotroph (no Mcr, fae
        present), where fwdABC + ftr are the Fhc complex."""
    ctx: dict = {}
    for tid in ("mcrA", "mcrB", "mcrG"):
        m = re.search(r"mcr_(methanogenic|reverse)",
                      (calls.get(tid) or {}).get("evidence_source", ""))
        if m:
            ctx["mcr_direction"] = m.group(1)
            ctx["mode"] = m.group(1)
            return ctx
    if (calls.get("fae") or {}).get("status") in _PRESENT:
        ctx["mode"] = "oxidative"
    return ctx


def context_note(ctx: dict) -> str:
    """One-line description of the context, for subtitles ('' if none)."""
    d = (ctx or {}).get("mcr_direction")
    if d == "reverse":
        return "Mcr direction: reverse (ANME-type methane oxidation)"
    if d == "methanogenic":
        return "Mcr direction: methanogenic"
    if (ctx or {}).get("mode") == "oxidative":
        return "H₄MPT C1 module: oxidative (Fae / Fhc, no Mcr)"
    return ""


def evaluate(codes: dict[str, int], ctx: dict | None = None) -> dict[str, dict]:
    """{step_id: {state, reversed, lines:[[token,…],…]}} for one genome."""
    mode = (ctx or {}).get("mode")
    out: dict[str, dict] = {}
    for st in STEPS:
        def present(slot) -> bool:
            return any(codes.get(g, 0) in (1, 2) for g in _slot_genes(slot))
        off = (mode in st.get("hide_in", ())
               or ("only_in" in st and mode not in st["only_in"]))
        route_full = any(all(present(s) for s in r) for r in st["routes"])
        any_gene = any(present(s) for r in st["routes"] for s in r)
        state = ("absent" if off else "complete" if route_full
                 else "partial" if any_gene else "absent")
        out[st["id"]] = {
            "state": state,
            "reversed": bool((st.get("reversible") and mode == "reverse")
                             or mode in st.get("up_in", ())),
            "lines": [[_token(s, codes, off) for s in line] for line in st["lines"]],
        }
    return out


# ───────────────────────────── geometry ──────────────────────────────────────

def _inside(pt, node, pad) -> bool:
    _, x, y, w, h = NODES[node]
    return abs(pt[0] - x) <= w / 2 + pad and abs(pt[1] - y) <= h / 2 + pad


def arrow_geometry(frm, to, bow, off0, off1, n=48) -> dict:
    """Polyline (trimmed at both node boxes) + arrowhead triangle."""
    p0 = (NODES[frm][1] + off0[0], NODES[frm][2] + off0[1])
    p2 = (NODES[to][1] + off1[0], NODES[to][2] + off1[1])
    dx, dy = p2[0] - p0[0], p2[1] - p0[1]
    d = math.hypot(dx, dy) or 1.0
    c = ((p0[0] + p2[0]) / 2 - bow * dy / d, (p0[1] + p2[1]) / 2 + bow * dx / d)
    pts = []
    for i in range(n + 1):
        t = i / n
        a, b, e = (1 - t) ** 2, 2 * (1 - t) * t, t ** 2
        pts.append((a * p0[0] + b * c[0] + e * p2[0],
                    a * p0[1] + b * c[1] + e * p2[1]))
    pts = [p for p in pts
           if not _inside(p, frm, NODE_PAD) and not _inside(p, to, NODE_PAD)]
    if len(pts) < 2:
        return {"line": [], "head": []}
    tip = pts[-1]
    # Direction at the tip, taken a few samples back for stability.
    ref = pts[max(0, len(pts) - 4)]
    ux, uy = tip[0] - ref[0], tip[1] - ref[1]
    u = math.hypot(ux, uy) or 1.0
    ux, uy = ux / u, uy / u
    base = (tip[0] - ux * HEAD_LEN, tip[1] - uy * HEAD_LEN)
    head = [tip,
            (base[0] - uy * HEAD_W, base[1] + ux * HEAD_W),
            (base[0] + uy * HEAD_W, base[1] - ux * HEAD_W)]
    # End the shaft at the arrowhead base so a thick line never pokes through.
    shaft = [p for p in pts
             if (p[0] - tip[0]) * ux + (p[1] - tip[1]) * uy <= -HEAD_LEN * 0.8]
    shaft.append(base)
    return {"line": shaft, "head": head}


def layout() -> dict:
    """Everything a renderer needs, JSON-serialisable."""
    r = lambda p: [round(p[0], 2), round(p[1], 2)]          # noqa: E731
    steps = []
    for st in STEPS:
        def geom(specs):
            out = []
            for a in specs:
                g = arrow_geometry(*a)
                out.append({"line": [r(p) for p in g["line"]],
                            "head": [r(p) for p in g["head"]]})
            return out
        # A reversible step also carries the same arrows drawn the other way.
        rev = [(to, frm, -bow, o1, o0) for frm, to, bow, o0, o1 in st["arrows"]] \
            if st.get("reversible") or st.get("up_in") else []
        steps.append({"id": st["id"], "name": st["name"],
                      "pathway": st["pathway"], "arrows": geom(st["arrows"]),
                      "arrows_rev": geom(rev), "label": list(st["label"])})
    return {
        "xlim": list(XLIM), "ylim": list(YLIM), "line_step": LINE_STEP,
        "wide": (XLIM[1] - XLIM[0]) > 130,
        "nodes": [{"id": k, "label": v[0], "x": v[1], "y": v[2],
                   "w": v[3], "h": v[4]} for k, v in NODES.items()],
        "steps": steps,
    }
