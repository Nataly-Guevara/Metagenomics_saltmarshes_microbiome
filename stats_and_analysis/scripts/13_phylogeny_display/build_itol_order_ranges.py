#!/usr/bin/env python3

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

TREE_FILE = Path(
    "results/phylogeny_bacteria_352/"
    "bacteria_352_bac120.tree"
)

TAXONOMY_FILE = Path(
    "results/phylogeny_bacteria_352/"
    "tree_taxonomy_display_352.tsv"
)

OUTPUT_FILE = Path(
    "results/phylogeny_bacteria_352/"
    "itol_order_coloured_clades_352.txt"
)

SUMMARY_FILE = Path(
    "results/phylogeny_bacteria_352/"
    "order_coloured_clade_summary_352.tsv"
)

# Do not highlight isolated or very small groups.
MIN_CLADE_SIZE = 5

# Order, display colour and publication-style letter.
ORDER_STYLES = {
    "Solirubrobacterales": ("a", "#1f77b4"),
    "Longimicrobiales": ("b", "#ff7f0e"),
    "Acidimicrobiales": ("c", "#2ca02c"),
    "Nitriliruptorales": ("d", "#d62728"),
    "Spongiisociales": ("e", "#9467bd"),
    "Promineifilales": ("f", "#8c564b"),
    "Multivoradales": ("g", "#e377c2"),
    "Flavobacteriales": ("h", "#7f7f7f"),
    "Rhodothermales": ("i", "#bcbd22"),
    "Palauibacterales": ("j", "#17becf"),
    "Cyanobacteriales": ("k", "#393b79"),
    "Xanthomonadales": ("l", "#637939"),
    "Desulfobacterales": ("m", "#8c6d31"),
    "Rhizobiales": ("n", "#843c39"),
    "Rhodobacterales": ("o", "#7b4173"),
    "Pseudomonadales": ("p", "#3182bd"),
}


def hex_to_rgba(hex_colour, alpha=0.20):
    """Convert #RRGGBB to an iTOL-compatible RGBA colour."""
    value = hex_colour.lstrip("#")

    red = int(value[0:2], 16)
    green = int(value[2:4], 16)
    blue = int(value[4:6], 16)

    return f"rgba({red},{green},{blue},{alpha})"


if not TREE_FILE.exists():
    sys.exit(f"ERROR: tree file not found:\n{TREE_FILE}")

if not TAXONOMY_FILE.exists():
    sys.exit(f"ERROR: taxonomy file not found:\n{TAXONOMY_FILE}")

# ---------------------------------------------------------
# Read MAG-to-order assignments
# ---------------------------------------------------------

mag_to_order = {}

with TAXONOMY_FILE.open("r", newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    required_columns = {"Final_MAG", "Order_display"}

    if reader.fieldnames is None:
        sys.exit("ERROR: taxonomy table has no header.")

    missing_columns = required_columns - set(reader.fieldnames)

    if missing_columns:
        sys.exit(
            "ERROR: taxonomy table is missing columns: "
            + ", ".join(sorted(missing_columns))
        )

    for row in reader:
        mag = row["Final_MAG"].strip()
        order = row["Order_display"].strip()

        if mag:
            mag_to_order[mag] = order

# ---------------------------------------------------------
# Extract terminal tips in displayed Newick order
# ---------------------------------------------------------

tree_text = TREE_FILE.read_text().strip()

tree_tips = re.findall(
    r"(?<=[(,])([^():,;\s]+)(?=:)",
    tree_text,
)

if len(tree_tips) != 352:
    sys.exit(
        f"ERROR: expected 352 tree tips, but extracted "
        f"{len(tree_tips)}."
    )

if len(set(tree_tips)) != len(tree_tips):
    sys.exit("ERROR: duplicate tree-tip identifiers detected.")

missing_taxonomy = [
    tip
    for tip in tree_tips
    if tip not in mag_to_order
]

if missing_taxonomy:
    sys.exit(
        f"ERROR: {len(missing_taxonomy)} tree tips are missing "
        "from the taxonomy table."
    )

# ---------------------------------------------------------
# Identify contiguous displayed clades for selected orders
# ---------------------------------------------------------

runs_by_order = defaultdict(list)

current_order = None
current_tips = []


def save_current_run():
    if (
        current_order in ORDER_STYLES
        and len(current_tips) >= MIN_CLADE_SIZE
    ):
        runs_by_order[current_order].append(
            list(current_tips)
        )


for tip in tree_tips:
    order = mag_to_order[tip]

    if order == current_order:
        current_tips.append(tip)
    else:
        save_current_run()
        current_order = order
        current_tips = [tip]

save_current_run()

# ---------------------------------------------------------
# Prepare legend information
# ---------------------------------------------------------

displayed_orders = [
    order
    for order in ORDER_STYLES
    if runs_by_order.get(order)
]

legend_shapes = []
legend_colours = []
legend_labels = []

for order in displayed_orders:
    letter, solid_colour = ORDER_STYLES[order]

    legend_shapes.append("1")
    legend_colours.append(solid_colour)
    legend_labels.append(f"{letter}. {order}")

# ---------------------------------------------------------
# Write iTOL coloured-clade range dataset
# ---------------------------------------------------------

with OUTPUT_FILE.open("w") as outfile:
    outfile.write("DATASET_RANGE\n")
    outfile.write("SEPARATOR TAB\n")
    outfile.write("DATASET_LABEL\tMajor order clades\n")
    outfile.write("COLOR\t#808080\n")

    # Coloured sector rather than bracket.
    outfile.write("RANGE_TYPE\tbox\n")

    # Highlight the phylogenetic clade itself.
    outfile.write("RANGE_COVER\tclade\n")

    # Do not extend through the external rings or datasets.
    outfile.write("COVER_DATASETS\t0\n")

    # Keep MAG labels readable rather than colouring over them.
    outfile.write("COVER_LABELS\t0\n")

    outfile.write("FIT_LABELS\t0\n")
    outfile.write("UNROOTED_SMOOTH\tsmooth\n")

    # Display publication-style letter inside each clade.
    outfile.write("SHOW_LABELS\t1\n")
    outfile.write("LABEL_POSITION\tcenter-center\n")
    outfile.write("LABELS_VERTICAL\t0\n")
    outfile.write("ROTATE_ALIGNED_LABELS\t1\n")
    outfile.write("STRAIGHT_LABELS\t0\n")
    outfile.write("LABEL_ROTATION\t0\n")
    outfile.write("LABEL_SHIFT_X\t0\n")
    outfile.write("LABEL_SHIFT_Y\t0\n")
    outfile.write("LABEL_AUTO_COLOR\t0\n")
    outfile.write("LABEL_SIZE_FACTOR\t1.4\n")
    outfile.write("LABEL_OUTLINE_WIDTH\t2\n")
    outfile.write("LABEL_OUTLINE_COLOR\t#ffffff\n")

    # Small gap between adjacent coloured sectors.
    outfile.write("VERTICAL_SHRINK\t0.04\n")

    if displayed_orders:
        outfile.write("LEGEND_TITLE\tOrder\n")
        outfile.write("LEGEND_SCALE\t1\n")
        outfile.write(
            "LEGEND_SHAPES\t"
            + "\t".join(legend_shapes)
            + "\n"
        )
        outfile.write(
            "LEGEND_COLORS\t"
            + "\t".join(legend_colours)
            + "\n"
        )
        outfile.write(
            "LEGEND_LABELS\t"
            + "\t".join(legend_labels)
            + "\n"
        )

    outfile.write("DATA\n")

    for order in displayed_orders:
        letter, solid_colour = ORDER_STYLES[order]
        transparent_fill = hex_to_rgba(
            solid_colour,
            alpha=0.20,
        )

        order_runs = runs_by_order[order]

        for tips in order_runs:
            start_tip = tips[0]
            end_tip = tips[-1]

            fields = [
                start_tip,
                end_tip,
                transparent_fill,
                "",
                solid_colour,
                "solid",
                "1",
                letter,
                "#000000",
                "1.2",
                "bold",
            ]

            outfile.write("\t".join(fields) + "\n")

# ---------------------------------------------------------
# Write validation summary
# ---------------------------------------------------------

with SUMMARY_FILE.open("w", newline="") as outfile:
    fieldnames = [
        "Letter",
        "Order",
        "Clade_number",
        "MAG_count",
        "Start_MAG",
        "End_MAG",
        "Colour",
    ]

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()

    for order in displayed_orders:
        letter, colour = ORDER_STYLES[order]

        for clade_number, tips in enumerate(
            runs_by_order[order],
            start=1,
        ):
            writer.writerow(
                {
                    "Letter": letter,
                    "Order": order,
                    "Clade_number": clade_number,
                    "MAG_count": len(tips),
                    "Start_MAG": tips[0],
                    "End_MAG": tips[-1],
                    "Colour": colour,
                }
            )

total_clades = sum(
    len(runs_by_order[order])
    for order in displayed_orders
)

highlighted_mags = sum(
    len(tips)
    for order in displayed_orders
    for tips in runs_by_order[order]
)

print("Finished successfully.")
print(f"Tree tips processed: {len(tree_tips)}")
print(f"Orders highlighted: {len(displayed_orders)}")
print(f"Coloured clades created: {total_clades}")
print(f"MAGs inside coloured clades: {highlighted_mags}")
print(f"Minimum clade size: {MIN_CLADE_SIZE}")
print(f"iTOL dataset: {OUTPUT_FILE}")
print(f"Validation summary: {SUMMARY_FILE}")
