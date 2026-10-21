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
    "Figure2B_refined_BS_RZ.svg"
)

PANELS = [
    ("BS", "i. Bulk soil", ["NAT", "BESE", "BARE", "DEGR"], [
        ("Osmotic", "Capability_Compatible_solute_uptake", "Compatible-solute uptake"),
        ("Carbon_P", "Capability_PolyP_turnover", "Polyphosphate turnover"),
        ("Iron", "Iron_complex_uptake", "Iron-complex uptake"),
        ("Colonisation", "Direct_plant_interaction", "ACC deaminase"),
        ("Redox", "Capability_Reduced_sulfur_oxidation", "Reduced-sulfur oxidation"),
        ("Redox", "Capability_DNRA", "Nitrite to ammonium"),
        ("Redox", "Capability_NO_reduction", "NO reduction"),
        ("Redox", "Capability_N2O_reduction", "N2O reduction"),
        ("Nitrogen", "Capability_Nitrogen_fixation", "Nitrogen fixation"),
        ("Nitrogen", "Capability_Assimilatory_nitrate_reduction",
         "Assimilatory nitrate reduction"),
        ("Carbon_P", "Capability_Carbon_fixation", "Carbon fixation"),
    ]),
    ("RZ", "ii. Rhizosphere", ["NAT", "BESE", "BARE"], [
        ("Osmotic", "Capability_Trehalose", "Trehalose potential"),
        ("Carbon_P", "Capability_PolyP_turnover", "Polyphosphate turnover"),
        ("Nitrogen", "Capability_Ammonium_uptake", "Ammonium uptake"),
        ("Iron", "Iron_complex_uptake", "Iron-complex uptake"),
        ("Redox", "Capability_LowO2_respiration", "Low-O2 respiration"),
        ("Redox", "Capability_Nitrate_respiration", "Nitrate respiration"),
    ])
]

values = defaultdict(list)

with SOURCE.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["Filter"] == "Unfiltered":
            key = (
                r["Habitat"], r["Treatment"],
                r["GRESP"], r["Capability"]
            )
            values[key].append(float(r["Weighted_prevalence_pct"]))

def colour(v):
    t = min(1, max(0, v / 100))
    stops = [
        (0.0, (248, 247, 241)),
        (0.5, (149, 190, 211)),
        (1.0, (35, 64, 130))
    ]
    for (a, c1), (b, c2) in zip(stops[:-1], stops[1:]):
        if t <= b:
            q = (t-a)/(b-a)
            rgb = [
                round(x + q*(y-x))
                for x, y in zip(c1, c2)
            ]
            return "#" + "".join(f"{x:02x}" for x in rgb)
    return "#234082"

svg = [
    '<svg xmlns="http://www.w3.org/2000/svg" '
    'width="1260" height="610" viewBox="0 0 1260 610">',
    '<rect width="1260" height="610" fill="white"/>'
]

def text(x, y, value, size=13, bold=False, fill="#243040", anchor="start"):
    svg.append(
        f'<text x="{x}" y="{y}" font-family="Arial" '
        f'font-size="{size}" font-weight="{"bold" if bold else "normal"}" '
        f'fill="{fill}" text-anchor="{anchor}">{escape(str(value))}</text>'
    )

def rect(x, y, w, h, fill):
    svg.append(
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}"/>'
    )

group_colours = {
    "Osmotic": "#047A99",
    "Redox": "#DE013A",
    "Nitrogen": "#F77D2E",
    "Carbon_P": "#032349",
    "Iron": "#8E5EA2",
    "Colonisation": "#2E7D32"
}

text(30, 39, "Metabolic capabilities across restoration environments",
     22, True)
text(30, 65,
     "Mean abundance-weighted prevalence among recovered MAGs (%)",
     14, fill="#526071")

# Continuous colour legend
for i in range(301):
    rect(820+i, 77, 1, 16, colour(i/3))

for tick in [0,25,50,75,100]:
    text(820+tick*3, 112, str(tick), 11, anchor="middle")

positions = {
    "BS": (30, 335),
    "RZ": (680, 990)
}

for habitat, title, treatments, capabilities in PANELS:
    label_x, heat_x = positions[habitat]
    cell_w = 73
    cell_h = 29
    y0 = 207

    text(label_x, 158, title, 18, True)

    for j, treatment in enumerate(treatments):
        name = "DEG" if treatment == "DEGR" else treatment
        text(
            heat_x+j*cell_w+cell_w/2,
            190, name, 13, True, anchor="middle"
        )

    for i, (group, capability, label) in enumerate(capabilities):
        y = y0+i*cell_h

        rect(label_x, y+5, 5, 19, group_colours[group])
        text(label_x+13, y+20, label, 12)

        for j, treatment in enumerate(treatments):
            key = (habitat, treatment, group, capability)

            if key not in values:
                raise ValueError(f"Missing: {key}")

            v = values[key]
            expected = 2 if habitat == "BS" and treatment == "BARE" else 3

            if len(v) != expected:
                raise ValueError(f"Unexpected number of sites: {key}, n={len(v)}")

            mean = sum(v)/len(v)
            x = heat_x+j*cell_w

            rect(x, y, cell_w, cell_h, colour(mean))
            text(
                x+cell_w/2, y+20, f"{mean:.1f}",
                12, fill="white" if mean >= 67 else "#243040",
                anchor="middle"
            )

text(30, 555,
     "Descriptive site means: n = 3 pooled-DNA metagenomes per treatment × habitat, except BARE BS (n = 2).",
     12, fill="#526071")

text(30, 574,
     "Nitrite-to-ammonium potential uses NirBD/NrfAH markers; it does not establish a complete DNRA pathway.",
     12, fill="#526071")

svg.append("</svg>")
OUT.write_text("\n".join(svg))

print("SUCCESS")
print("Bulk soil:", len(PANELS[0][3]), "capabilities")
print("Rhizosphere:", len(PANELS[1][3]), "capabilities")
print("Output:", OUT)
