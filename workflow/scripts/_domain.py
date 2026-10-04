"""Everything specific to the METHANE cycle that the figures and report need:
pathway order, palette, labels and output file names. The nitrogen and sulfur
sister pipelines have their own _domain.py and _cycle_model.py; every other
figure / report script is identical between the three.

Palette: the pathway hues are the first seven slots of the same validated
palette the sister pipelines use, assigned in CAT_ORDER (the order the pathway
blocks sit next to each other; adjacent CVD ΔE >= 8, normal-vision ΔE >= 15,
both modes). Three light-mode hues are below 3:1 on white, so colour never
carries identity alone: every coloured mark has a text label or a block header.
"""

CYCLE_LETTER = "CH₄"
CYCLE_NAME = "methane"
CALLS_TSV = "mcycle_calls.tsv"
LOCI_TSV = "mcycle_loci.tsv"

CAT_ORDER = [
    "mcr_core", "co2_reduction", "acetoclastic",
    "methylotrophic_methanogenesis", "energy_conservation",
    "aerobic_methane_oxidation", "formaldehyde_c1",
]

PATHWAY_COLOR = {
    "mcr_core":                      "#2a78d6",   # blue
    "co2_reduction":                 "#eb6834",   # orange
    "acetoclastic":                  "#1baf7a",   # aqua
    "methylotrophic_methanogenesis": "#eda100",   # yellow
    "energy_conservation":           "#e87ba4",   # magenta
    "aerobic_methane_oxidation":     "#008300",   # green
    "formaldehyde_c1":               "#4a3aa7",   # violet
}
# Same hues stepped for a dark surface (used by report.html only).
PATHWAY_COLOR_DARK = {
    "mcr_core":                      "#3987e5",
    "co2_reduction":                 "#d95926",
    "acetoclastic":                  "#199e70",
    "methylotrophic_methanogenesis": "#c98500",
    "energy_conservation":           "#d55181",
    "aerobic_methane_oxidation":     "#008300",
    "formaldehyde_c1":               "#9085e9",
}
PATHWAY_LABEL = {
    "mcr_core":                      "Mcr / Mtr core (and reverse, ANME)",
    "co2_reduction":                 "CO₂ reduction (H₄MPT C1 module)",
    "acetoclastic":                  "Acetoclastic methanogenesis",
    "methylotrophic_methanogenesis": "Methylotrophic methanogenesis",
    "energy_conservation":           "Hdr / hydrogenases",
    "aerobic_methane_oxidation":     "Aerobic methane oxidation",
    "formaldehyde_c1":               "Formaldehyde oxidation / C1 uptake",
}
# Shorter form for the angled block headers of the overview grid.
PATHWAY_SHORT = dict(
    PATHWAY_LABEL,
    mcr_core="Mcr / Mtr core",
    co2_reduction="CO₂ reduction",
    acetoclastic="Acetoclastic",
    methylotrophic_methanogenesis="Methylotrophic",
    aerobic_methane_oxidation="Aerobic CH₄ oxidation",
    formaldehyde_c1="C1 oxidation / uptake",
)

# Tidy-ups applied to complex / module ids when they are shown as labels.
PRETTY_REPLACE = [("hdr abc", "HdrABC"), ("hdr de", "HdrDE"), ("mcr", "Mcr"),
                  ("mtr", "Mtr"), ("fwd", "Fwd"), ("mvh", "Mvh"), ("frh", "Frh"),
                  ("f420", "F₄₂₀"), ("codh acs", "CODH/ACS"), ("pmmo", "pMMO"),
                  ("smmo", "sMMO"), ("mxa ", "Mxa "), (" dh", " dehydrogenase"),
                  ("aom", "(AOM)"), ("h4mpt", "H₄MPT"), ("rump", "RuMP")]

# Legend text for a module that is ruled out by an exclusion rule.
RULED_OUT_LABEL = "ruled out — key gene absent, or Mcr runs the other way"
