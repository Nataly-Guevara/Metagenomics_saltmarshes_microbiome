#!/usr/bin/env python3

"""
Build coordinated iTOL ecological-strategy symbol rings.

CURRENT STRATEGIES
------------------
1. Osmotic adaptation
2. Redox flexibility
3. Nitrogen acquisition and transformation
4. Carbon + phosphorus resource potential

Later we will add:
5. Iron acquisition
6. Colonisation and plant-growth-promotion potential

The script uses the FINAL curated per-MAG strategy scores already
generated in the ecological-strategy analyses.

Visual logic
------------
- iTOL DATASET_SYMBOL
- circle symbol for every MAG
- all strategy scores are on a common 0-1 scale
- symbol size increases linearly with strategy breadth
- colour intensity increases with strategy breadth
- each strategy has its own base colour
- score 0 is represented by very pale grey
- identical score values have identical symbol sizes across strategies

The previous osmotic iTOL dataset used approximately:
    size = 1 + 9 * score

We retain this rule for all ecological strategies.
"""

import csv
import re
from pathlib import Path


# ============================================================
# PATHS
# ============================================================

BASE = Path(
    "/path/to/metagenomics_project"
)

TREE = (
    BASE
    / "results/phylogeny_bacteria_352/"
      "bacteria_352_bac120.tree"
)

OUTPUT_DIR = (
    BASE
    / "results/ecological_strategies/"
      "itol_strategy_rings"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

COMBINED_MATRIX = (
    OUTPUT_DIR
    / "ecological_strategy_scores_352_MAGs.tsv"
)


# ============================================================
# FINAL STRATEGY DEFINITIONS
# ============================================================

STRATEGIES = {

    "osmotic": {
        "label":
            "Osmotic adaptation",

        "table":
            BASE
            / "results/ecological_strategies/"
              "osmotic_adaptation/tables/"
              "osmotic_strategy_final_373_MAGs.tsv",

        "score_column":
            "Osmotic_strategy_score",

        "color":
            "#047A99",

        "output":
            OUTPUT_DIR
            / "itol_ring_1_osmotic_adaptation.txt",
    },


    "redox": {
        "label":
            "Redox flexibility",

        "table":
            BASE
            / "results/ecological_strategies/"
              "redox_flexibility/tables/"
              "redox_strategy_final_v2_373_MAGs.tsv",

        "score_column":
            "Redox_strategy_score",

        "color":
            "#DE013A",

        "output":
            OUTPUT_DIR
            / "itol_ring_2_redox_flexibility.txt",
    },


    "nitrogen": {
        "label":
            "Nitrogen acquisition & transformation",

        "table":
            BASE
            / "results/ecological_strategies/"
              "nitrogen_acquisition_transformation/tables/"
              "nitrogen_strategy_final_373_MAGs.tsv",

        "score_column":
            "Nitrogen_strategy_score",

        "color":
            "#F77D2E",

        "output":
            OUTPUT_DIR
            / "itol_ring_3_nitrogen_strategy.txt",
    },


    "carbon_phosphorus": {
        "label":
            "Carbon + phosphorus",

        "table":
            BASE
            / "results/ecological_strategies/"
              "carbon_phosphorus_resource/tables/"
              "carbon_phosphorus_strategy_final_373_MAGs.tsv",

        "score_column":
            "Carbon_Phosphorus_strategy_score",

        "color":
            "#032349",

        "output":
            OUTPUT_DIR
            / "itol_ring_4_carbon_phosphorus.txt",
    },
}


# ============================================================
# VISUAL SETTINGS
# ============================================================

PALE_COLOR = "#F2F2F2"

MIN_SYMBOL_SIZE = 1.0
MAX_SYMBOL_SIZE = 10.0

ITOL_MAXIMUM_SIZE = 28
SYMBOL_SPACING = 5

# iTOL symbol code:
# 2 = circle
SYMBOL_CODE = 2


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def hex_to_rgb(hex_color):
    """
    Convert #RRGGBB to integer RGB tuple.
    """

    value = hex_color.lstrip("#")

    return tuple(
        int(
            value[i:i + 2],
            16
        )
        for i in (
            0,
            2,
            4,
        )
    )


def rgb_to_hex(rgb):
    """
    Convert RGB tuple to #RRGGBB.
    """

    return (
        "#"
        + "".join(
            f"{max(0, min(255, int(round(v)))):02X}"
            for v in rgb
        )
    )


def interpolate_color(
    low_color,
    high_color,
    score
):
    """
    Linear colour interpolation between pale grey at 0
    and the strategy-specific colour at 1.
    """

    score = max(
        0.0,
        min(
            1.0,
            float(score)
        )
    )

    low = hex_to_rgb(
        low_color
    )

    high = hex_to_rgb(
        high_color
    )

    rgb = tuple(
        low[i]
        +
        score
        * (
            high[i]
            - low[i]
        )
        for i in range(3)
    )

    return rgb_to_hex(
        rgb
    )


def score_to_size(score):
    """
    Common symbol-size scaling across all six strategy axes.

    score = 0   -> size 1
    score = 0.5 -> size 5.5
    score = 1   -> size 10
    """

    score = max(
        0.0,
        min(
            1.0,
            float(score)
        )
    )

    return (
        MIN_SYMBOL_SIZE
        +
        (
            MAX_SYMBOL_SIZE
            - MIN_SYMBOL_SIZE
        )
        * score
    )


def read_strategy_table(
    table_path,
    score_column
):
    """
    Read Final_MAG and one strategy-score column.
    """

    if not table_path.exists():

        raise FileNotFoundError(
            f"Strategy table not found: {table_path}"
        )

    scores = {}

    with table_path.open(
        "r",
        newline=""
    ) as infile:

        reader = csv.DictReader(
            infile,
            delimiter="\t"
        )

        if reader.fieldnames is None:

            raise ValueError(
                f"No header in {table_path}"
            )

        required = {
            "Final_MAG",
            score_column,
        }

        missing = (
            required
            - set(
                reader.fieldnames
            )
        )

        if missing:

            raise ValueError(
                f"Missing columns in {table_path}: "
                + ", ".join(
                    sorted(
                        missing
                    )
                )
            )

        for row in reader:

            mag = row[
                "Final_MAG"
            ].strip()

            value = row[
                score_column
            ].strip()

            if not mag:
                continue

            if value == "":
                continue

            scores[
                mag
            ] = float(
                value
            )

    return scores


# ============================================================
# EXTRACT 352 BACTERIAL TREE TIPS
# ============================================================

if not TREE.exists():

    raise FileNotFoundError(
        f"Tree not found: {TREE}"
    )


newick = (
    TREE
    .read_text()
    .strip()
)


matches = re.findall(
    r"(?<=[(,])([^():,]+)(?=:)",
    newick
)


tree_tips = []

seen_tree_tips = set()


for item in matches:

    tip = (
        item
        .strip()
        .strip("'")
        .strip('"')
    )

    if (
        tip
        and
        tip not in seen_tree_tips
    ):

        tree_tips.append(
            tip
        )

        seen_tree_tips.add(
            tip
        )


if len(
    tree_tips
) != 352:

    print(
        "WARNING: expected 352 bacterial tree tips, "
        f"detected {len(tree_tips)}."
    )


# ============================================================
# READ ALL FOUR FINAL STRATEGY TABLES
# ============================================================

strategy_scores = {}


for key, info in (
    STRATEGIES.items()
):

    scores = read_strategy_table(
        info[
            "table"
        ],
        info[
            "score_column"
        ]
    )

    strategy_scores[
        key
    ] = scores


# ============================================================
# QC: CHECK TREE-TIP MATCHING
# ============================================================

print()
print(
    "=" * 72
)

print(
    "iTOL ECOLOGICAL STRATEGY RING PREPARATION"
)

print(
    "=" * 72
)

print()

print(
    f"Bacterial tree tips: "
    f"{len(tree_tips)}"
)

print()


for key, info in (
    STRATEGIES.items()
):

    scores = strategy_scores[
        key
    ]

    matched = [
        tip
        for tip in tree_tips
        if tip in scores
    ]

    missing = [
        tip
        for tip in tree_tips
        if tip not in scores
    ]

    extras = [
        mag
        for mag in scores
        if mag not in seen_tree_tips
    ]


    print(
        f"{info['label']}:"
    )

    print(
        f"  MAGs in final table: "
        f"{len(scores)}"
    )

    print(
        f"  Tree tips matched:   "
        f"{len(matched)}"
    )

    print(
        f"  Tree tips missing:   "
        f"{len(missing)}"
    )

    print(
        f"  Non-tree MAGs:       "
        f"{len(extras)}"
    )


    if missing:

        print(
            "  WARNING - missing tree tips:"
        )

        for tip in missing[
            :10
        ]:

            print(
                f"    {tip}"
            )

        if len(
            missing
        ) > 10:

            print(
                f"    ... plus "
                f"{len(missing) - 10} more"
            )

    print()


# ============================================================
# REQUIRE COMPLETE MATCHING BEFORE WRITING
# ============================================================

for key, info in (
    STRATEGIES.items()
):

    scores = strategy_scores[
        key
    ]

    missing = [
        tip
        for tip in tree_tips
        if tip not in scores
    ]

    if missing:

        raise ValueError(
            f"{info['label']} is missing "
            f"{len(missing)} bacterial tree MAGs. "
            "iTOL files were not written."
        )


# ============================================================
# WRITE COORDINATED iTOL SYMBOL RINGS
# ============================================================

for key, info in (
    STRATEGIES.items()
):

    label = info[
        "label"
    ]

    base_color = info[
        "color"
    ]

    output = info[
        "output"
    ]

    scores = strategy_scores[
        key
    ]


    mid_color = interpolate_color(
        PALE_COLOR,
        base_color,
        0.5
    )


    with output.open(
        "w"
    ) as outfile:

        outfile.write(
            "DATASET_SYMBOL\n"
        )

        outfile.write(
            "SEPARATOR TAB\n"
        )

        outfile.write(
            f"DATASET_LABEL\t"
            f"{label}\n"
        )

        outfile.write(
            f"COLOR\t"
            f"{base_color}\n"
        )

        outfile.write(
            f"MAXIMUM_SIZE\t"
            f"{ITOL_MAXIMUM_SIZE}\n"
        )

        outfile.write(
            f"SYMBOL_SPACING\t"
            f"{SYMBOL_SPACING}\n"
        )

        outfile.write(
            "SHOW_LABELS\t1\n"
        )

        outfile.write(
            "LABEL_SIZE_FACTOR\t1\n"
        )

        outfile.write(
            "LABEL_ROTATION\t0\n"
        )

        outfile.write(
            "LABEL_ALIGN_TO_TREE\t1\n"
        )

        outfile.write(
            f"LEGEND_TITLE\t"
            f"{label} score\n"
        )

        outfile.write(
            "LEGEND_SHAPES\t"
            "2\t2\t2\n"
        )

        outfile.write(
            "LEGEND_COLORS\t"
            f"{PALE_COLOR}\t"
            f"{mid_color}\t"
            f"{base_color}\n"
        )

        outfile.write(
            "LEGEND_LABELS\t"
            "0.0\t0.5\t1.0\n"
        )

        outfile.write(
            "LEGEND_SHAPE_SCALES\t"
            "0.10\t0.55\t1.00\n"
        )

        outfile.write(
            "LEGEND_SHAPE_INVERT\t"
            "0\t0\t0\n"
        )

        outfile.write(
            "\nDATA\n"
        )


        for mag in tree_tips:

            score = scores[
                mag
            ]

            size = score_to_size(
                score
            )

            color = interpolate_color(
                PALE_COLOR,
                base_color,
                score
            )

            label_text = (
                f"{label}: "
                f"{score:.3f}"
            )

            outfile.write(
                f"{mag}\t"
                f"{SYMBOL_CODE}\t"
                f"{size:.4f}\t"
                f"{color}\t"
                f"1\t"
                f"-1\t"
                f"{label_text}\n"
            )


# ============================================================
# WRITE COMBINED 352-MAG STRATEGY MATRIX
# ============================================================

matrix_fields = [
    "Final_MAG",
    "Osmotic_adaptation",
    "Redox_flexibility",
    "Nitrogen_acquisition_transformation",
    "Carbon_phosphorus_resource",
]


with COMBINED_MATRIX.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=matrix_fields,
        delimiter="\t",
    )

    writer.writeheader()


    for mag in tree_tips:

        writer.writerow({
            "Final_MAG":
                mag,

            "Osmotic_adaptation":
                f"{strategy_scores['osmotic'][mag]:.4f}",

            "Redox_flexibility":
                f"{strategy_scores['redox'][mag]:.4f}",

            "Nitrogen_acquisition_transformation":
                f"{strategy_scores['nitrogen'][mag]:.4f}",

            "Carbon_phosphorus_resource":
                f"{strategy_scores['carbon_phosphorus'][mag]:.4f}",
        })


# ============================================================
# PRINT FINAL SUMMARY
# ============================================================

print(
    "=" * 72
)

print(
    "iTOL STRATEGY RINGS CREATED SUCCESSFULLY"
)

print(
    "=" * 72
)

print()


for number, (
    key,
    info
) in enumerate(
    STRATEGIES.items(),
    start=1
):

    scores = [
        strategy_scores[
            key
        ][
            mag
        ]
        for mag in tree_tips
    ]

    mean_score = (
        sum(
            scores
        )
        / len(
            scores
        )
    )

    minimum = min(
        scores
    )

    maximum = max(
        scores
    )

    print(
        f"Ring {number}: "
        f"{info['label']}"
    )

    print(
        f"  Colour: "
        f"{info['color']}"
    )

    print(
        f"  Mean score: "
        f"{mean_score:.3f}"
    )

    print(
        f"  Range: "
        f"{minimum:.3f} - "
        f"{maximum:.3f}"
    )

    print(
        f"  File: "
        f"{info['output']}"
    )

    print()


print(
    f"Combined score matrix:"
)

print(
    f"  {COMBINED_MATRIX}"
)

print()

print(
    "All four rings use:"
)

print(
    "  - the same 352 bacterial MAGs"
)

print(
    "  - a common 0-1 score scale"
)

print(
    "  - identical score-to-symbol-size scaling"
)

print(
    "  - strategy-specific colour gradients"
)

print()

print(
    "Future rings:"
)

print(
    "  Ring 5: Iron acquisition"
)

print(
    "  Ring 6: Colonisation / PGP potential"
)

print(
    "=" * 72
)
