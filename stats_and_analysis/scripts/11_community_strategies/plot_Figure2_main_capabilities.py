#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict
from html import escape

SOURCE = Path(
    "results/ecological_strategies/community_weighted_GRESP/"
    "community_capabilities_20_samples_long.tsv"
)

OUT = Path(
    "results/functional_annotation/figure2_functional_reorganisation/"
    "Figure2_main_capabilities.svg"
)

# Each panel has a biologically motivated selection.
PANELS = [
    ("BS", "A  |  Bulk soil",
     ["NAT", "BESE", "BARE", "DEGR"], [
        ("Osmotic", "Capability_Compatible_solute_uptake",
         "Compatible-solute uptake"),
        ("Carbon_P", "Capability_PolyP_turnover",
         "Polyphosphate turnover"),
        ("Iron", "Iron_complex_uptake",
         "Iron-complex uptake"),
        ("Colonisation", "Direct_plant_interaction",
         "ACC deaminase potential"),
        ("Redox", "Capability_Reduced_sulfur_oxidation",
         "Reduced-sulfur oxidation"),
        ("Redox", "Capability_DNRA",
         "Nitrite to ammonium (NirBD/NrfAH)"),
        ("Redox", "Capability_NO_reduction",
         "NO reduction"),
        ("Redox", "Capability_N2O_reduction",
         "N2O reduction"),
        ("Nitrogen", "Capability_Nitrogen_fixation",
         "Nitrogen fixation"),
        ("Nitrogen", "Capability_Assimilatory_nitrate_reduction",
         "Assimilatory nitrate reduction"),
        ("Carbon_P", "Capability_Carbon_fixation",
         "Carbon fixation"),
    ]),

    ("RZ", "B  |  Rhizosphere",
     ["NAT", "BESE", "BARE"], [
        ("Osmotic", "Capability_Trehalose",
         "Trehalose potential"),
        ("Carbon_P", "Capability_PolyP_turnover",
         "Polyphosphate turnover"),
        ("Nitrogen", "Capability_Ammonium_uptake",
         "Ammonium uptake"),
        ("Iron", "Iron_complex_uptake",
         "Iron-complex uptake"),
        ("Redox", "Capability_LowO2_respiration",
         "Low-O2 respiration"),
        ("Redox", "Capability_Nitrate_respiration",
         "Nitrate respiration"),
    ])
]

# Load unfiltered abundance-weighted capability prevalence.
values = defaultdict(list)
seen = set()

with SOURCE.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["Filter"] != "Unfiltered":
            continue

        key = (
            r["Habitat"], r["Site"], r["Treatment"],
            r["GRESP"], r["Capability"]
        )

        if key in seen:
            raise ValueError(f"Duplicate observation: {key}")

        seen.add(key)

        group_key = (
            r["Habitat"], r["Treatment"],
            r["GRESP"], r["Capability"]
        )

        values[group_key].append(
            float(r["Weighted_prevalence_pct"])
        )

# Continuous colour map (shared across both panels).
def heat_color(value):
    t = max(0, min(1, value / 100))

    stops = [
        (0.0, (248, 247, 241)),
        (0.5, (149, 190, 211)),
        (1.0, (35, 64, 130))
    ]

    for i in range(2):
        a, c1 = stops[i]
        b, c2 = stops[i + 1]

        if t <= b:
            u = (t - a) / (b - a)
            rgb = [
                round(x + u * (y - x))
                for x, y in zip(c1, c2)
            ]
            return "#" + "".join(f"{v:02x}" for v in rgb)

def rect(x, y, w, h, color):
    svg.append(
        f'<rect x="{x}" y="{y}" width="{w}" '
        f'height="{h}" fill="{color}"/>'
    )

def text(x, y, label, size=13, weight="normal",
         color="#202631", anchor="start"):
    svg.append(
        f'<text x="{x}" y="{y}" '
        f'font-family="Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" '
        f'fill="{color}" text-anchor="{anchor}">'
        f'{escape(str(label))}</text>'
    )

SVG_W = 1090
SVG_H = 805
X0 = 545
CELL_W = 120
CELL_H = 29

svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" '
    f'width="{SVG_W}" height="{SVG_H}" '
    f'viewBox="0 0 {SVG_W} {SVG_H}">',
    f'<rect width="{SVG_W}" height="{SVG_H}" fill="white"/>'
]

text(35, 38, "Microbial ecological capabilities across restoration environments",
     22, "bold")

text(35, 64,
     "Mean MAG-abundance-weighted capability prevalence (%)",
     14, color="#526071")

# Shared legend
LEGEND_X = 620
LEGEND_Y = 83

for i in range(301):
    rect(LEGEND_X + i, LEGEND_Y, 1, 13,
         heat_color(i / 3))

for tick in [0, 25, 50, 75, 100]:
    x = LEGEND_X + tick * 3
    text(x, 111, str(tick), 11, anchor="middle")

group_colors = {
    "Osmotic": "#047A99",
    "Redox": "#DE013A",
    "Nitrogen": "#F77D2E",
    "Carbon_P": "#032349",
    "Iron": "#8E5EA2",
    "Colonisation": "#2E7D32"
}

panel_y = {"BS": 149, "RZ": 542}

for habitat, title, treatments, capabilities in PANELS:

    start_y = panel_y[habitat]

    text(35, start_y - 37, title, 18, "bold")

    for j, treatment in enumerate(treatments):
        x = X0 + j * CELL_W + CELL_W / 2
        label = "DEGRADED" if treatment == "DEGR" else treatment
        text(x, start_y - 12, label, 13, "bold",
             anchor="middle")

    for i, (group, capability, label) in enumerate(capabilities):

        y = start_y + i * CELL_H

        # Small category marker
        rect(38, y + 4, 5, CELL_H - 8, group_colors[group])

        text(56, y + 20, label, 13)

        for j, treatment in enumerate(treatments):
            key = (habitat, treatment, group, capability)

            if key not in values:
                raise ValueError(f"Missing capability: {key}")

            vals = values[key]
            mean = sum(vals) / len(vals)

            x = X0 + j * CELL_W

            rect(x, y, CELL_W, CELL_H, heat_color(mean))

            fg = "white" if mean >= 67 else "#202631"
            text(x + CELL_W / 2, y + 20,
                 f"{mean:.1f}", 12,
                 color=fg, anchor="middle")

# Explanatory notes
text(35, 758,
     "Cells show site means; BARE bulk soil n=2, all other groups n=3.",
     12, color="#526071")

text(35, 778,
     "Nitrite-to-ammonium row represents NirBD/NrfAH potential, not a confirmed complete DNRA pathway.",
     12, color="#526071")

svg.append("</svg>")

OUT.write_text("\n".join(svg))

print("SUCCESS")
print("Bulk-soil capabilities:", len(PANELS[0][3]))
print("Rhizosphere capabilities:", len(PANELS[1][3]))
print("Output:", OUT)
