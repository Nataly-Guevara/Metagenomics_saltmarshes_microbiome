#!/usr/bin/env python3
import csv
from pathlib import Path
from html import escape

source = Path(
    "results/ecological_strategies/community_weighted_GRESP/"
    "community_GRESP_20_samples.tsv"
)
outdir = Path(
    "results/functional_annotation/figure2_functional_reorganisation"
)
outdir.mkdir(parents=True, exist_ok=True)

strategies = [
    ("Osmotic", "Osmotic"),
    ("Redox", "Redox"),
    ("Nitrogen", "Nitrogen"),
    ("Carbon_P", "Carbon + P"),
    ("Iron", "Iron"),
    ("Colonisation", "Colonisation")
]

with source.open() as f:
    records = [
        r for r in csv.DictReader(f, delimiter="\t")
        if r["Filter"] == "Unfiltered"
    ]

assert len(records) == 20

lookup = {}
for r in records:
    key = (r["Habitat"], r["Site"], r["Treatment"])
    if key in lookup:
        raise ValueError(f"Duplicate sample: {key}")
    lookup[key] = r

ranges = {}
for col, label in strategies:
    vals = [float(r[col]) for r in records]
    ranges[col] = (min(vals), max(vals))

colors = {
    "NAT": "#248A48",
    "BESE": "#F15A52",
    "BARE": "#D39A35",
    "DEGR": "#794A9B"
}

svg = [
    '<svg xmlns="http://www.w3.org/2000/svg" '
    'width="1200" height="1040" viewBox="0 0 1200 1040">',
    '<rect width="1200" height="1040" fill="white"/>'
]

def txt(x, y, value, size=14, weight="normal",
        fill="#222222", anchor="start"):
    svg.append(
        f'<text x="{x}" y="{y}" font-family="Arial" '
        f'font-size="{size}" font-weight="{weight}" '
        f'fill="{fill}" text-anchor="{anchor}">'
        f'{escape(str(value))}</text>'
    )

def rect(x, y, w, h, fill):
    svg.append(
        f'<rect x="{x}" y="{y}" width="{w}" '
        f'height="{h}" fill="{fill}"/>'
    )

def heat_color(value, minimum, maximum):
    t = ((value - minimum) / (maximum - minimum)
         if maximum > minimum else 0.5)
    low = (238, 247, 250)
    high = (18, 101, 137)
    rgb = [
        round(a + t * (b - a))
        for a, b in zip(low, high)
    ]
    return "#" + "".join(f"{v:02X}" for v in rgb), t

txt(45, 39, "Genome-resolved ecological strategies",
    23, "bold")
txt(45, 64,
    "Community abundance-weighted GRESP scores across saltmarsh environments",
    14, fill="#555555")

x0 = 335
cw = 135
row_h = 33

txt(50, 109, "Site", 14, "bold")
txt(165, 109, "Treatment", 14, "bold")

for i, (col, label) in enumerate(strategies):
    txt(x0 + i * cw + cw / 2, 109,
        label, 13, "bold", anchor="middle")

y = 140
count = 0

for habitat, title in [
    ("BS", "BULK SOIL"),
    ("RZ", "RHIZOSPHERE")
]:
    txt(45, y, title, 16, "bold")
    y += 24

    for site in ["ACP", "BASS", "TRUG"]:
        for treatment in ["NAT", "BESE", "BARE", "DEGR"]:
            key = (habitat, site, treatment)
            if key not in lookup:
                continue

            r = lookup[key]

            txt(50, y + 21, site, 13)
            rect(165, y + 5, 12, 19, colors[treatment])
            txt(185, y + 21, treatment, 13)

            for i, (col, label) in enumerate(strategies):
                value = float(r[col])
                minimum, maximum = ranges[col]
                color, t = heat_color(value, minimum, maximum)
                x = x0 + i * cw

                rect(x, y, cw - 3, row_h - 3, color)
                txt(
                    x + (cw - 3) / 2,
                    y + 21,
                    f"{value:.3f}",
                    12,
                    fill="white" if t > 0.65 else "#222222",
                    anchor="middle"
                )

            y += row_h
            count += 1

        y += 7

    y += 32

txt(45, y + 10,
    "Colour scaled independently within each strategy; numbers show original 0–1 scores.",
    12, fill="#555555")
txt(45, y + 29,
    "TRUG BARE bulk soil was not sequenced. Degraded rhizosphere is not applicable.",
    12, fill="#555555")

svg.append("</svg>")

output = outdir / "Figure2_GRESP_heatmap.svg"
output.write_text("\n".join(svg))

print("SUCCESS")
print("Samples plotted:", count)
print("Strategies:", len(strategies))
print("Output:", output)
