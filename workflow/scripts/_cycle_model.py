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

The diagram is drawn as a metro map: each pathway is a line in its colour,
each compound a station, each gene a stop on the line of the step it performs.
The shared ladder of H₄MPT-bound C1 intermediates is the straight trunk from
CO₂ (left) to CH₄ (right); methanogens travel it rightwards, aerobic
methanotrophs and methylotrophs return along the bottom (CH₄ → CH₃OH → HCHO →
formate → CO₂); acetate, methylamines and the heterodisulfide cycle feed in
from the top. Two things depend on the genome:

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
XLIM = (0.0, 244.0)
YLIM = (0.0, 118.0)

# Stations.  id: (label, x, y, (label dx, label dy, ha, va))
# Lines meet at the station centre; the label offset keeps the name clear of
# the lines that arrive there.
NODES = {
    "CO2":   ("CO₂",                16.0,  60.0, (0.0,   6.0, "center", "bottom")),
    "FMFR":  ("formyl-MFR",         46.0,  60.0, (0.0,  -6.0, "center", "top")),
    "FH4":   ("formyl-H₄MPT",       76.0,  60.0, (4.5, -11.0, "left",   "top")),
    "MEN":   ("methenyl-H₄MPT",    106.0,  60.0, (0.0,  -6.0, "center", "top")),
    "MLN":   ("methylene-H₄MPT",   136.0,  60.0, (4.0, -11.0, "left",   "top")),
    "MH4":   ("methyl-H₄MPT",      166.0,  60.0, (0.0,  -6.0, "center", "top")),
    "MCOM":  ("methyl-CoM",        196.0,  60.0, (-4.0, -11.0, "right",  "top")),
    "CH4":   ("CH₄",               226.0,  60.0, (6.5,   0.0, "left",   "center")),
    "MEOH":  ("CH₃OH",             196.0,  30.0, (-4.5, -5.5, "right",  "top")),
    "HCHO":  ("HCHO",              136.0,  30.0, (-4.5, -5.5, "right",  "top")),
    "BIO":   ("biomass",           106.0,  10.0, (-6.0,  0.0, "right",  "center")),
    "FOR":   ("formate",            46.0,  30.0, (0.0,  -6.0, "center", "top")),
    "AC":    ("acetate",           166.0, 104.0, (-6.0,  0.0, "right",  "center")),
    "ACCOA": ("acetyl-CoA",        166.0,  84.0, (-6.0,  0.0, "right",  "center")),
    "HDS":   ("CoM-S-S-CoB",       226.0, 108.0, (0.0,   6.0, "center", "bottom")),
    "THIOL": ("CoM-SH + CoB-SH",   226.0,  84.0, (0.0,  -4.5, "center", "top")),
    "MA":    ("methylamines · DMS", 196.0, 94.0, (0.0,   6.0, "center", "bottom")),
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

# path:   waypoints of the line from the first station to the last, octilinear
#   (horizontal / vertical / 45°); the stations' centres are the two ends.
# offset: shift of the whole path to the left of the direction of travel, for
#   parallel strands between the same two stations (0 = on the centre line).
# side:   where the stop labels go — +1 left of the direction of travel, -1
#   right of it, 0 alternating (dense steps); a list gives one value per
#   straight stretch of the path.
# routes: alternative gene sets, any one of which performs the step.
# lines:  the genes drawn as stops, in order (may include accessory genes that
#   no route requires, e.g. mcrC/D, mvhAGD, frhABG).
# reversible: part of the methanogenesis chain — travelled the other way round
#   in a reverse-Mcr (ANME) genome (the chevron moves to the other end).
# hide_in / only_in: genome modes (see genome_context) in which the step is
#   switched off / the only modes in which it is evaluated.
STEPS = [
    {"id": "fwd", "name": "CO₂ ⇄ formyl-methanofuran (formylmethanofuran dehydrogenase)",
     "pathway": _R, "reversible": True, "hide_in": {"oxidative"},
     "path": ["CO2", "FMFR"], "side": 1,
     "routes": [["fwdA", "fwdB", "fwdC"]],
     "lines": [["fwdA", "fwdB"], ["fwdC", "fwdD"]]},
    {"id": "ftr", "name": "Formyl-methanofuran ⇄ formyl-H₄MPT (formyltransferase)",
     "pathway": _R, "reversible": True, "hide_in": {"oxidative"},
     "path": ["FMFR", "FH4"], "side": 1,
     "routes": [["ftr"]], "lines": [["ftr"]]},
    {"id": "fhc", "name": "Formyl-H₄MPT → formate (Fhc complex: FhcABC + FhcD, the "
                          "bacterial homologues of fwdABC + ftr)",
     "pathway": _F, "only_in": {"oxidative"},
     "path": ["FH4", "FOR"], "side": 1,
     "routes": [["ftr", "fwdA", "fwdB", "fwdC"]],
     "lines": [["ftr", "fwdA"], ["fwdB", "fwdC"]]},
    {"id": "fdh", "name": "Formate → CO₂ (formate dehydrogenase)", "pathway": _F,
     "path": ["FOR", "CO2"], "side": -1,
     "routes": [["fdh"], ["fdhA", "fdhB"]],
     "lines": [["fdh", "fdhA", "fdhB"]]},
    {"id": "mch", "name": "Formyl-H₄MPT ⇄ methenyl-H₄MPT (cyclohydrolase)",
     "pathway": _R, "reversible": True, "up_in": {"oxidative"},
     "path": ["FH4", "MEN"], "side": 1,
     "routes": [["mch"]], "lines": [["mch"]]},
    {"id": "mtd", "name": "Methenyl-H₄MPT ⇄ methylene-H₄MPT (F₄₂₀-dependent Mtd, or H₂-forming Hmd)",
     "pathway": _R, "reversible": True, "hide_in": {"oxidative"},
     "path": ["MEN", "MLN"], "side": 1,
     "routes": [["mtd"], ["hmd"]], "lines": [["mtd", "hmd"]]},
    {"id": "mtdB", "name": "Methylene-H₄MPT → methenyl-H₄MPT (NAD(P)-dependent MtdB)",
     "pathway": _F,
     "path": ["MLN", "MEN"], "offset": -2.6, "side": -1,
     "routes": [["mtdB"]], "lines": [["mtdB"]]},
    {"id": "mer", "name": "Methylene-H₄MPT ⇄ methyl-H₄MPT (Mer; F₄₂₀H₂ from Frh)",
     "pathway": _R, "reversible": True,
     "path": ["MLN", "MH4"], "side": 1,
     "routes": [["mer"]], "lines": [["mer"], ["frhA", "frhB", "frhG"]]},
    {"id": "mtr", "name": "Methyl-H₄MPT ⇄ methyl-CoM (Mtr, Na⁺-pumping)",
     "pathway": _K, "reversible": True,
     "path": ["MH4", "MCOM"], "side": 1,
     "routes": [["mtrA", "mtrB", "mtrC", "mtrD", "mtrE", "mtrH"]],   # mtrF / mtrG: see targets.yaml
     "lines": [["mtrA", "mtrB", "mtrC", "mtrD"], ["mtrE", "mtrF", "mtrG", "mtrH"]]},
    {"id": "mcr", "name": "Methyl-CoM ⇄ CH₄ (methyl-coenzyme M reductase)",
     "pathway": _K, "reversible": True,
     "path": ["MCOM", "CH4"], "side": 1,
     "routes": [["mcrA", "mcrB", "mcrG"]],
     "lines": [["mcrA", "mcrB", "mcrG"], ["mcrC", "mcrD", "mcrA_anme"]]},
    {"id": "hdr", "name": "Heterodisulfide → CoM-SH + CoB-SH (heterodisulfide reductase)",
     "pathway": _E,
     "path": ["HDS", "THIOL"], "side": 0,
     "routes": [["hdrA", "hdrB", "hdrC"], ["hdrD", "hdrE"]],
     "lines": [["hdrA", "hdrB", "hdrC"], ["mvhA", "mvhG", "mvhD"], ["hdrD", "hdrE"]]},
    {"id": "ack", "name": "Acetate → acetyl-CoA (Ack + Pta, or Acs)", "pathway": _A,
     "path": ["AC", "ACCOA"], "side": 1,
     "routes": [["ackA", "pta"], ["acs"]],
     "lines": [["ackA", "pta"], ["acs"]]},
    {"id": "cdh", "name": "Acetyl-CoA ⇄ methyl-H₄MPT + CO₂ (CODH/ACS complex)",
     "pathway": _A,
     "path": ["ACCOA", "MH4"], "side": 1,
     "routes": [["cdhA", "cdhB", "cdhC", "cdhD", "cdhE"]],
     "lines": [["cdhA", "cdhB", "cdhC"], ["cdhD", "cdhE"]]},
    {"id": "mta", "name": "Methanol → methyl-CoM (MtaABC)", "pathway": _M,
     "path": ["MEOH", "MCOM"], "side": -1,
     "routes": [["mtaA", "mtaB", "mtaC"]],
     "lines": [["mtaA", "mtaB", "mtaC"]]},
    {"id": "mtm", "name": "Methylamines → methyl-CoM (MtmBC / MtbBC / MttBC + MtbA)",
     "pathway": _M,
     "path": ["MA", "MCOM"], "offset": 1.6, "side": 1,
     "routes": [["mtbA", "mtmB", "mtmC"], ["mtbA", "mtbB", "mtbC"],
                ["mtbA", "mttB", "mttC"]],
     "lines": [["mtbA", "mtmB", "mtmC"], ["mtbB", "mtbC", "mttB", "mttC"]]},
    {"id": "mts", "name": "Methyl sulfides → methyl-CoM (MtsAB)", "pathway": _M,
     "path": ["MA", "MCOM"], "offset": -1.6, "side": -1,
     "routes": [["mtsA", "mtsB"]], "lines": [["mtsA", "mtsB"]]},
    {"id": "mmo", "name": "CH₄ → CH₃OH (particulate or soluble methane monooxygenase)",
     "pathway": _O,
     "path": ["CH4", (226.0, 42.0), (214.0, 30.0), "MEOH"], "side": [-1, -1, 1],
     "routes": [["pmoA", "pmoB", "pmoC"], ["mmoX", "mmoY", "mmoZ"]],
     "lines": [["pmoA", "pmoB", "pmoC"], ["mmoX", "mmoY", "mmoZ"],
               ["mmoB", "mmoC", "mmoD", "nod"]]},
    {"id": "mdh", "name": "CH₃OH → HCHO (MxaFI or XoxF methanol dehydrogenase)",
     "pathway": _O,
     "path": ["MEOH", "HCHO"], "side": -1,
     "routes": [["mxaF", "mxaI"], ["xoxF"]],
     "lines": [["mxaF", "mxaI"], ["xoxF"]]},
    {"id": "fae", "name": "HCHO → methylene-H₄MPT (formaldehyde-activating enzyme)",
     "pathway": _F, "hide_in": {"methanogenic", "reverse"},
     "path": ["HCHO", "MLN"], "side": -1,
     "routes": [["fae"]], "lines": [["fae"]]},
    {"id": "gsh", "name": "HCHO → formate (glutathione-dependent: Gfa, FrmA, FrmB)",
     "pathway": _F,
     "path": ["HCHO", "FOR"], "side": -1,
     "routes": [["frmA", "frmB"]], "lines": [["gfa"], ["frmA", "frmB"]]},
    {"id": "rump", "name": "HCHO → biomass (ribulose-monophosphate cycle)",
     "pathway": _F, "hide_in": {"methanogenic", "reverse"},
     "path": ["HCHO", (136.0, 22.0), (124.0, 10.0), "BIO"], "offset": -1.6, "side": -1,
     "routes": [["hxlA", "hxlB"]], "lines": [["hxlA", "hxlB"]]},
    {"id": "serine", "name": "C1 → biomass (serine cycle)", "pathway": _F,
     "hide_in": {"methanogenic", "reverse"},
     "path": ["HCHO", (136.0, 22.0), (124.0, 10.0), "BIO"], "offset": 1.6, "side": 1,
     "routes": [["sgaA", "hprA", "mtkA", "mtkB", "mcl"]],
     "lines": [["sgaA", "hprA", "gckA"], ["mtkA", "mtkB", "mcl"]]},
]

# Figure sizing used by make_cycle_map.py (inches / points).
FIG = {"single_w": 11.4, "grid_w": 8.4, "single_fs": 9.0, "grid_fs": 6.6,
       "grid_cols": 2}

# Metro geometry (data units).
LINE_W = 1.5             # stroke width of a line
STATION_R = 2.5          # station disc radius (ring drawn on it)
STOP_R = 1.05            # stop disc radius
STOP_MARGIN = (3.8, 5.2) # free track after the first / before the last station
CORNER_MARGIN = 3.0      # no stop this close to a bend
STOP_SPACING_MAX = 9.0   # stops on a long stretch are grouped, not spread out
LABEL_GAP = 2.3          # stop centre → label
CHEVRON = 1.9            # direction chevron size


# ───────────────────────────── evaluation ────────────────────────────────────

def _slot_genes(slot) -> list[str]:
    return [slot] if isinstance(slot, str) else list(slot)


def _token(slot, codes: dict[str, int]) -> dict:
    """The gene to print for a slot: the best-supported alternative (first one
    when none is present). The code is always the gene's real call — also on a
    step that is switched off in this genome's mode (see `off` in evaluate()),
    so the map never shows a called gene as absent."""
    genes = _slot_genes(slot)
    rank = {2: 0, 1: 1, -1: 2, 0: 3}
    best = min(genes, key=lambda g: (rank[codes.get(g, 0)], genes.index(g)))
    return {"gene": best, "label": DISPLAY.get(best, best),
            "code": codes.get(best, 0)}


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
    """{step_id: {state, reversed, off, lines:[[token,…],…]}} for one genome.

    `off` is '' for a step evaluated normally, otherwise the reason the step is
    not drawn as active in this genome's mode (hide_in / only_in). An off step
    has state `absent` whatever its genes, but its gene tokens keep their real
    calls — a present gene is never printed as absent; the renderers mark the
    step and state the reason instead."""
    mode = (ctx or {}).get("mode")
    out: dict[str, dict] = {}
    for st in STEPS:
        def present(slot) -> bool:
            return any(codes.get(g, 0) in (1, 2) for g in _slot_genes(slot))
        off = ""
        if mode in st.get("hide_in", ()):
            off = f"not drawn in a {mode} genome"
        elif "only_in" in st and mode not in st["only_in"]:
            off = "drawn only in " + " / ".join(f"{m} genomes" for m in sorted(st["only_in"]))
        route_full = any(all(present(s) for s in r) for r in st["routes"])
        any_gene = any(present(s) for r in st["routes"] for s in r)
        state = ("absent" if off else "complete" if route_full
                 else "partial" if any_gene else "absent")
        out[st["id"]] = {
            "state": state,
            "off": off,
            "reversed": bool((st.get("reversible") and mode == "reverse")
                             or mode in st.get("up_in", ())),
            "lines": [[_token(s, codes) for s in line] for line in st["lines"]],
        }
    return out


# ───────────────────────────── geometry ──────────────────────────────────────

def _pt(p):
    """A waypoint: a station id or an (x, y) pair."""
    return (NODES[p][1], NODES[p][2]) if isinstance(p, str) else tuple(p)


def _offset_polyline(pts, d):
    """Shift a polyline by d to the left of its direction of travel (y up);
    corners are joined by the intersection of the shifted segments."""
    if not d or len(pts) < 2:
        return list(pts)
    segs = []
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        L = math.hypot(x1 - x0, y1 - y0) or 1.0
        nx, ny = -(y1 - y0) / L, (x1 - x0) / L
        segs.append(((x0 + nx * d, y0 + ny * d), (x1 + nx * d, y1 + ny * d)))
    out = [segs[0][0]]
    for (a0, a1), (b0, b1) in zip(segs, segs[1:]):
        ax, ay = a1[0] - a0[0], a1[1] - a0[1]
        bx, by = b1[0] - b0[0], b1[1] - b0[1]
        den = ax * by - ay * bx
        if abs(den) < 1e-9:
            out.append(a1)
            continue
        t = ((b0[0] - a0[0]) * by - (b0[1] - a0[1]) * bx) / den
        out.append((a0[0] + ax * t, a0[1] + ay * t))
    out.append(segs[-1][1])
    return out


def _along(pts, dist):
    """Point, unit tangent and stretch index at a distance along the polyline."""
    left = dist
    for k, ((x0, y0), (x1, y1)) in enumerate(zip(pts, pts[1:])):
        L = math.hypot(x1 - x0, y1 - y0)
        if left <= L or (x1, y1) == pts[-1]:
            t = min(max(left / L, 0.0), 1.0) if L else 0.0
            return (x0 + (x1 - x0) * t, y0 + (y1 - y0) * t), \
                (((x1 - x0) / L, (y1 - y0) / L) if L else (1.0, 0.0)), k
        left -= L
    return pts[-1], (1.0, 0.0), len(pts) - 2


def _length(pts) -> float:
    return sum(math.hypot(x1 - x0, y1 - y0) for (x0, y0), (x1, y1) in zip(pts, pts[1:]))


def _stop_distances(pts, n: int) -> list[float]:
    """Distances along the polyline of n evenly spread stops, keeping clear of
    the two stations and of every bend; a long stretch groups its stops around
    its middle instead of spreading them out."""
    m0, m1 = STOP_MARGIN
    usable, at = [], 0.0
    nseg = len(pts) - 1
    for k, ((x0, y0), (x1, y1)) in enumerate(zip(pts, pts[1:])):
        L = math.hypot(x1 - x0, y1 - y0)
        a = at + (m0 if k == 0 else CORNER_MARGIN)
        b = at + L - (m1 if k == nseg - 1 else CORNER_MARGIN)
        if b > a:
            usable.append((a, b))
        at += L
    U = sum(b - a for a, b in usable) or 1.0
    span = min(U, n * STOP_SPACING_MAX)
    lead = (U - span) / 2
    out = []
    for i in range(n):
        u = lead + span * (i + 0.5) / n
        for a, b in usable:
            if u <= b - a + 1e-9:
                out.append(a + u)
                break
            u -= b - a
        else:
            out.append(usable[-1][1] if usable else L / 2)
    return out


def _orient(tangent) -> str:
    tx, ty = tangent
    if abs(ty) > 0.9:
        return "v"
    if abs(tx) > 0.9:
        return "h"
    return "d"


def step_geometry(st) -> dict:
    """Polyline, stops (one per gene token, in `lines` order) and the chevron
    positions of one step. Stop labels: on a vertical or diagonal stretch the
    text is horizontal beside the stop; on a horizontal stretch it is written
    at 45° so dense stops do not collide (the metro-map convention)."""
    pts = _offset_polyline([_pt(p) for p in st["path"]], st.get("offset", 0.0))
    L = _length(pts)
    tokens = [g for line in st["lines"] for g in line]
    n = len(tokens)
    stops = []
    sides = st.get("side", 1)
    for i, d in enumerate(_stop_distances(pts, n)):
        (x, y), (tx, ty), k = _along(pts, d)
        side = sides[k] if isinstance(sides, (list, tuple)) else sides
        side = side or (1 if i % 2 == 0 else -1)
        nx, ny = -ty * side, tx * side          # left normal × side
        o = _orient((tx, ty))
        if o == "h":
            # 45° label: starts above the stop reading up-right, or ends below
            # it reading up to the stop.
            up = 1 if ny > 0 else -1
            lx, ly = x + 0.55 * up, y + LABEL_GAP * up
            stops.append({"x": x, "y": y, "lx": lx, "ly": ly, "rot": 45,
                          "ha": "left" if up > 0 else "right", "va": "center"})
        else:
            lx, ly = x + nx * LABEL_GAP, y + ny * LABEL_GAP
            ha = "left" if nx > 0.3 else "right" if nx < -0.3 else "center"
            va = "center" if o == "v" else ("bottom" if ny > 0 else "top")
            stops.append({"x": x, "y": y, "lx": lx, "ly": ly, "rot": 0,
                          "ha": ha, "va": va})

    # Chevron tips: just before each end station, pointing along the travel.
    def chevron(at_end: bool):
        d = L - STATION_R - 1.0 if at_end else STATION_R + 1.0
        (x, y), (tx, ty), _ = _along(pts, d)
        if not at_end:
            tx, ty = -tx, -ty
        return {"x": x, "y": y, "tx": tx, "ty": ty}
    return {"path": pts, "stops": stops, "chevron": chevron(True),
            "chevron_rev": chevron(False)}


def layout() -> dict:
    """Everything a renderer needs, JSON-serialisable."""
    r = lambda v: round(v, 2)                                 # noqa: E731
    steps = []
    for st in STEPS:
        g = step_geometry(st)
        steps.append({
            "id": st["id"], "name": st["name"], "pathway": st["pathway"],
            "path": [[r(x), r(y)] for x, y in g["path"]],
            "stops": [{k: (r(v) if isinstance(v, float) else v) for k, v in s.items()}
                      for s in g["stops"]],
            "chevron": {k: r(v) for k, v in g["chevron"].items()},
            "chevron_rev": {k: r(v) for k, v in g["chevron_rev"].items()},
        })
    return {
        "xlim": list(XLIM), "ylim": list(YLIM),
        "wide": (XLIM[1] - XLIM[0]) > 130,
        "geom": {"line_w": LINE_W, "station_r": STATION_R, "stop_r": STOP_R,
                 "chevron": CHEVRON},
        "nodes": [{"id": k, "label": v[0], "x": v[1], "y": v[2],
                   "ldx": v[3][0], "ldy": v[3][1], "ha": v[3][2], "va": v[3][3]}
                  for k, v in NODES.items()],
        "steps": steps,
    }
