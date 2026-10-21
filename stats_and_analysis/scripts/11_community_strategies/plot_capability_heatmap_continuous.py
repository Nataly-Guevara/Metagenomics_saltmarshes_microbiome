#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict
from html import escape

SOURCE = Path(
    "results/ecological_strategies/community_weighted_GRESP/"
    "community_capabilities_20_samples_long.tsv"
)

OUTDIR = Path(
    "results/functional_annotation/figure2_functional_reorganisation"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# 1. Display labels and row order
# -----------------------------

GROUPS = [
    ("Osmotic", [
        ("Capability_Trehalose", "Trehalose synthesis / use"),
        ("Capability_Compatible_solute_uptake", "Betaine / proline uptake"),
        ("Capability_Ectoine_biosynthesis", "Ectoine synthesis"),
        ("Capability_Mechanosensitive_protection", "Mechanosensitive protection"),
        ("Capability_Ion_homeostasis", "Na+/H+ and K+ homeostasis"),
    ]),
    ("Redox", [
        ("Capability_Aerobic_respiration", "Aerobic respiration"),
        ("Capability_LowO2_respiration", "Low-O2 respiration"),
        ("Capability_Nitrate_respiration", "Nitrate respiration"),
        ("Capability_Nitrite_to_NO", "Nitrite to NO"),
        ("Capability_DNRA", "Nitrite to ammonium (NirBD/NrfAH)"),
        ("Capability_NO_reduction", "NO reduction"),
        ("Capability_N2O_reduction", "N2O reduction"),
        ("Capability_Dissimilatory_sulfate_reduction_core", "DsrAB sulfate-reduction core*"),
        ("Capability_Reduced_sulfur_oxidation", "Reduced-sulfur oxidation"),
        ("Capability_Sulfite_oxidation", "Sulfite oxidation"),
    ]),
    ("Nitrogen", [
        ("Capability_Nitrogen_fixation", "Nitrogen fixation"),
        ("Capability_Ammonium_uptake", "Ammonium uptake"),
        ("Capability_Assimilatory_nitrate_reduction", "Assimilatory nitrate reduction"),
        ("Capability_Urea_utilisation", "Urea utilisation"),
        ("Capability_DNRA", "Nitrate to ammonium pathway"),
        ("Capability_Complete_denitrification", "Complete denitrification"),
        ("Capability_Nitrification", "Nitrification"),
        ("Capability_Anammox", "Anammox"),
    ]),
    ("Carbon_P", [
        ("Capability_Carbon_fixation", "Carbon fixation"),
        ("Capability_C1_utilisation", "C1 utilisation"),
        ("Capability_High_affinity_phosphate_uptake", "High-affinity phosphate uptake"),
        ("Capability_Organic_P_mobilisation", "Organic P mobilisation"),
        ("Capability_Phosphonate_acquisition", "Phosphonate acquisition"),
        ("Capability_PolyP_turnover", "Polyphosphate turnover"),
    ]),
    ("Iron", [
        ("Ferrous_iron_uptake", "Ferrous iron uptake"),
        ("Direct_ferric_iron_uptake", "Direct ferric iron uptake"),
        ("Siderophore_biosynthesis", "Siderophore biosynthesis"),
        ("Iron_complex_uptake", "Iron-complex uptake"),
    ]),
    ("Colonisation", [
        ("Chemotaxis", "Chemotaxis"),
        ("Motility", "Motility"),
        ("Attachment_biofilm_matrix", "Attachment / biofilm matrix"),
        ("Direct_plant_interaction", "Direct plant interaction"),
    ]),
]

COLUMN_ORDER = [
    ("BS", "NAT", "NATURAL\nBS"),
    ("BS", "BESE", "BESE\nBS"),
    ("BS", "BARE", "BARE\nBS"),
    ("BS", "DEGR", "DEGRADED\nBS"),
    ("RZ", "NAT", "NATURAL\nRZ"),
    ("RZ", "BESE", "BESE\nRZ"),
    ("RZ", "BARE", "BARE\nRZ"),
]

GROUP_LABELS = {
    "Osmotic": "Osmotic and\nsalinity adaptation",
    "Redox": "Redox and\nO2 / sulfur cycling",
    "Nitrogen": "Nitrogen\ncycling",
    "Carbon_P": "Carbon and\nphosphorus",
    "Iron": "Iron\nacquisition",
    "Colonisation": "Colonisation /\nplant interaction",
}

GROUP_COLORS = {
    "Osmotic": "#7B2CBF",
    "Redox": "#D62828",
    "Nitrogen": "#008B8B",
    "Carbon_P": "#1565C0",
    "Iron": "#E67E22",
    "Colonisation": "#2E7D32",
}

# -----------------------------
# 2. Read data and compute means
# -----------------------------

rows = []
with SOURCE.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        if r["Filter"] == "Unfiltered":
            rows.append(r)

means = defaultdict(list)

for r in rows:
    key = (r["Habitat"], r["Treatment"], r["GRESP"], r["Capability"])
    means[key].append(float(r["Weighted_prevalence_pct"]))

mean_value = {}
for key, vals in means.items():
    mean_value[key] = sum(vals) / len(vals)

# Also save matrix as TSV
matrix_out = OUTDIR / "Supplementary_capability_heatmap_matrix.tsv"
with matrix_out.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["GRESP", "Capability"] + [c[2].replace("\n", " ") for c in COLUMN_ORDER])
    for group, items in GROUPS:
        for capability, label in items:
            row = [group, label]
            for habitat, treatment, _ in COLUMN_ORDER:
                row.append(round(mean_value.get((habitat, treatment, group, capability), 0.0), 3))
            writer.writerow(row)

# -----------------------------
# 3. SVG settings
# -----------------------------

cell_w = 92
cell_h = 24
header_h = 88
left_group_w = 300
left_label_w = 360
top_title_h = 60
legend_h = 24
legend_w = 260
gap_between_blocks = 10  # small gap only between BS and RZ
bottom_margin = 55
right_margin = 25

n_rows = sum(len(items) for _, items in GROUPS)
heat_x = left_group_w + left_label_w
heat_y = top_title_h + header_h

heat_w = cell_w * 4 + gap_between_blocks + cell_w * 3
heat_h = cell_h * n_rows

svg_w = heat_x + heat_w + right_margin
svg_h = heat_y + heat_h + bottom_margin

svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
    f'<rect width="{svg_w}" height="{svg_h}" fill="white"/>'
]

def text(x, y, s, size=14, fill="#222", weight="normal", anchor="start", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    svg.append(
        f'<text x="{x}" y="{y}" font-family="Arial" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{transform}>{escape(str(s))}</text>'
    )

def rect(x, y, w, h, fill, stroke="none", sw=0):
    svg.append(
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'
    )

def line(x1, y1, x2, y2, color="#cccccc", width=1):
    svg.append(
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>'
    )

def heat_color(value):
    # 0 to 100 scale
    t = max(0.0, min(1.0, value / 100.0))
    low = (245, 241, 220)    # pale cream
    high = (54, 48, 146)     # deep purple-blue
    rgb = [round(a + t * (b - a)) for a, b in zip(low, high)]
    return "#" + "".join(f"{v:02X}" for v in rgb)

# -----------------------------
# 4. Title and legend
# -----------------------------

text(28, 34, "Functional annotation patterns of ecological capabilities", size=24, weight="bold")
text(28, 54, "Values are mean abundance-weighted prevalence in recovered MAGs (%)", size=13, fill="#555555")

legend_x = 28
legend_y = 72
text(legend_x, legend_y - 8, "Abundance-weighted capability prevalence (%) ", size=12, weight="bold")

for i in range(legend_w):
    value = 100 * i / legend_w
    rect(legend_x + i, legend_y, 1, legend_h, heat_color(value))

for tick_val in [0, 20, 40, 60, 80, 100]:
    tx = legend_x + legend_w * tick_val / 100
    line(tx, legend_y + legend_h, tx, legend_y + legend_h + 5, "#666666", 1)
    text(tx, legend_y + legend_h + 18, str(tick_val), size=11, anchor="middle", rotate=30)

# -----------------------------
# 5. Column headers
# -----------------------------

# BS block
x_positions = []
x = heat_x
for idx, (habitat, treatment, label) in enumerate(COLUMN_ORDER):
    x_positions.append(x)
    if idx == 3:
        x += cell_w + gap_between_blocks
    else:
        x += cell_w

for x, (_, _, label) in zip(x_positions, COLUMN_ORDER):
    parts = label.split("\n")
    text(x + cell_w / 2, heat_y - 28, parts[0], size=12, weight="bold", anchor="middle")
    text(x + cell_w / 2, heat_y - 10, parts[1], size=12, anchor="middle")

text(heat_x + cell_w * 2, heat_y - 54, "Bulk soil", size=14, weight="bold", anchor="middle")
text(heat_x + cell_w * 4 + gap_between_blocks + cell_w * 1.5, heat_y - 54, "Rhizosphere", size=14, weight="bold", anchor="middle")

# -----------------------------
# 6. Heatmap body
# -----------------------------

current_row = 0

for group, items in GROUPS:
    start_row = current_row

    for capability, label in items:
        y = heat_y + current_row * cell_h

        # row label
        text(left_group_w + 8, y + 17, label, size=12)

        # cells
        for (x, (habitat, treatment, _)) in zip(x_positions, COLUMN_ORDER):
            val = mean_value.get((habitat, treatment, group, capability), None)

            # missing combination shouldn't happen, but if so leave blank grey
            if val is None:
                fill = "#eeeeee"
                display = ""
            else:
                fill = heat_color(val)
                display = f"{val:.1f}"

            rect(x, y, cell_w, cell_h, fill)

            # choose text colour
            txt_fill = "white" if (val is not None and val >= 58) else "#222222"
            if display:
                text(x + cell_w / 2, y + 17, display, size=11, fill=txt_fill, anchor="middle")

        current_row += 1

    end_row = current_row

    # group label
    group_mid_y = heat_y + ((start_row + end_row) / 2) * cell_h
    text(left_group_w - 10, group_mid_y, GROUP_LABELS[group], size=12,
         fill=GROUP_COLORS[group], weight="bold", anchor="end")

    # horizontal separator between groups
    if group != GROUPS[-1][0]:
        sep_y = heat_y + current_row * cell_h
        line(left_group_w, sep_y, heat_x + heat_w, sep_y, "#dddddd", 1.5)

# outline around BS and RZ blocks
rect(heat_x, heat_y, cell_w * 4, heat_h, fill="none", stroke="#cccccc", sw=1)
rect(heat_x + cell_w * 4 + gap_between_blocks, heat_y, cell_w * 3, heat_h, fill="none", stroke="#cccccc", sw=1)

# -----------------------------
# 7. Footnote
# -----------------------------

foot_y = heat_y + heat_h + 28
text(28, foot_y, "Means are descriptive site averages (most groups n=3 sites; BARE bulk soil n=2 because TRUG BARE BS is missing).", size=11, fill="#555555")
text(28, foot_y + 16, "* DsrAB sulfate-reduction core indicates marker presence; reductive vs oxidative directionality remains unresolved and should be interpreted cautiously.", size=11, fill="#555555")

svg.append("</svg>")

out_svg = OUTDIR / "Supplementary_capability_heatmap_continuous.svg"
out_svg.write_text("\n".join(svg))

print("SUCCESS")
print("Rows:", n_rows)
print("Columns:", len(COLUMN_ORDER))
print("SVG:", out_svg)
print("Matrix:", matrix_out)
