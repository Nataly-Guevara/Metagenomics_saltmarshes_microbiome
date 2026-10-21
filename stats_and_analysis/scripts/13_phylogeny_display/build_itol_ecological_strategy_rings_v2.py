#!/usr/bin/env python3

"""
Build coordinated iTOL ecological-strategy rings - V2.

VISUAL GRAMMAR
==============

Broad ecological dimension -> symbol shape

Environmental tolerance:
    circle
    - Osmotic adaptation
    - Redox flexibility

Resource acquisition:
    square
    - Nitrogen acquisition and transformation
    - Carbon + phosphorus resource potential
    - Iron acquisition [future]

Host interaction:
    triangle
    - Colonisation / plant-growth-promotion potential [future]


Individual ecological strategy -> fixed colour

Strategy breadth score -> symbol size

The strategy score ranges from 0 to 1.

IMPORTANT
---------
Colour does NOT encode score.
Colour identifies the ecological strategy.

Symbol size is the only quantitative visual encoding.

All strategy rings use exactly the same score-to-size scaling.
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
      "itol_strategy_rings_v2"
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
# iTOL SYMBOL CODES
# ============================================================

# iTOL:
# 1 = square
# 2 = circle
# 4 = right triangle

SQUARE = 1
CIRCLE = 2
TRIANGLE = 4


# ============================================================
# FINAL STRATEGY DEFINITIONS
# ============================================================

STRATEGIES = {

    "osmotic": {

        "label":
            "Osmotic adaptation",

        "dimension":
            "Environmental tolerance",

        "symbol":
            CIRCLE,

        "color":
            "#047A99",

        "table":
            BASE
            / "results/ecological_strategies/"
              "osmotic_adaptation/tables/"
              "osmotic_strategy_final_373_MAGs.tsv",

        "score_column":
            "Osmotic_strategy_score",

        "output":
            OUTPUT_DIR
            / "itol_ring_1_osmotic_adaptation.txt",
    },


    "redox": {

        "label":
            "Redox flexibility",

        "dimension":
            "Environmental tolerance",

        "symbol":
            CIRCLE,

        "color":
            "#DE013A",

        "table":
            BASE
            / "results/ecological_strategies/"
              "redox_flexibility/tables/"
              "redox_strategy_final_v2_373_MAGs.tsv",

        "score_column":
            "Redox_strategy_score",

        "output":
            OUTPUT_DIR
            / "itol_ring_2_redox_flexibility.txt",
    },


    "nitrogen": {

        "label":
            "Nitrogen acquisition & transformation",

        "dimension":
            "Resource acquisition",

        "symbol":
            SQUARE,

        "color":
            "#F77D2E",

        "table":
            BASE
            / "results/ecological_strategies/"
              "nitrogen_acquisition_transformation/tables/"
              "nitrogen_strategy_final_373_MAGs.tsv",

        "score_column":
            "Nitrogen_strategy_score",

        "output":
            OUTPUT_DIR
            / "itol_ring_3_nitrogen_strategy.txt",
    },


    "carbon_phosphorus": {

        "label":
            "Carbon + phosphorus",

        "dimension":
            "Resource acquisition",

        "symbol":
            SQUARE,

        "color":
            "#032349",

        "table":
            BASE
            / "results/ecological_strategies/"
              "carbon_phosphorus_resource/tables/"
              "carbon_phosphorus_strategy_final_373_MAGs.tsv",

        "score_column":
            "Carbon_Phosphorus_strategy_score",

        "output":
            OUTPUT_DIR
            / "itol_ring_4_carbon_phosphorus.txt",
    },
}


# ============================================================
# FUTURE VISUAL ASSIGNMENTS
# ============================================================
#
# Iron acquisition:
#
#     dimension = Resource acquisition
#     symbol = SQUARE
#     suggested colour = #8E5EA2
#
#
# Colonisation / PGP:
#
#     dimension = Host interaction
#     symbol = TRIANGLE
#     suggested colour = #3A923A
#
# ============================================================


# ============================================================
# VISUAL SETTINGS
# ============================================================

MIN_SYMBOL_SIZE = 1.0
MAX_SYMBOL_SIZE = 10.0

ITOL_MAXIMUM_SIZE = 28

SYMBOL_SPACING = 5


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def score_to_size(score):

    """
    Convert strategy breadth score (0-1)
    into a common iTOL symbol size.

    0.0 -> 1.0
    0.2 -> 2.8
    0.4 -> 4.6
    0.5 -> 5.5
    0.6 -> 6.4
    0.8 -> 8.2
    1.0 -> 10.0
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
    Read Final_MAG and final strategy score.
    """

    if not table_path.exists():

        raise FileNotFoundError(
            f"Strategy table not found: "
            f"{table_path}"
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
                f"No header found in "
                f"{table_path}"
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
                f"Missing columns in "
                f"{table_path}: "
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


            if (
                not mag
                or value == ""
            ):
                continue


            scores[
                mag
            ] = float(
                value
            )


    return scores


# ============================================================
# READ BACTERIAL TREE
# ============================================================

if not TREE.exists():

    raise FileNotFoundError(
        f"Bacterial tree not found: "
        f"{TREE}"
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
        "WARNING: expected 352 bacterial "
        f"tree tips, found {len(tree_tips)}."
    )


# ============================================================
# READ FINAL STRATEGY TABLES
# ============================================================

strategy_scores = {}


for key, info in (
    STRATEGIES.items()
):

    strategy_scores[
        key
    ] = read_strategy_table(
        info[
            "table"
        ],
        info[
            "score_column"
        ]
    )


# ============================================================
# QC TREE-TIP MATCHING
# ============================================================

print()
print(
    "=" * 72
)

print(
    "iTOL ECOLOGICAL STRATEGY RINGS V2"
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
        f"{info['label']}"
    )

    print(
        f"  Dimension:          "
        f"{info['dimension']}"
    )

    print(
        f"  Symbol code:        "
        f"{info['symbol']}"
    )

    print(
        f"  Colour:             "
        f"{info['color']}"
    )

    print(
        f"  MAGs in table:      "
        f"{len(scores)}"
    )

    print(
        f"  Tree tips matched:  "
        f"{len(matched)}"
    )

    print(
        f"  Tree tips missing:  "
        f"{len(missing)}"
    )

    print(
        f"  Non-tree MAGs:      "
        f"{len(extras)}"
    )

    print()


# ============================================================
# REQUIRE FULL MATCH BEFORE WRITING
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
            f"{len(missing)} tree MAGs. "
            "No iTOL datasets written."
        )


# ============================================================
# WRITE iTOL SYMBOL RINGS
# ============================================================

for key, info in (
    STRATEGIES.items()
):

    label = info[
        "label"
    ]

    dimension = info[
        "dimension"
    ]

    symbol = info[
        "symbol"
    ]

    color = info[
        "color"
    ]

    output = info[
        "output"
    ]

    scores = strategy_scores[
        key
    ]


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
            f"{color}\n"
        )

        outfile.write(
            f"MAXIMUM_SIZE\t"
            f"{ITOL_MAXIMUM_SIZE}\n"
        )

        outfile.write(
            f"SYMBOL_SPACING\t"
            f"{SYMBOL_SPACING}\n"
        )


        # ----------------------------------------------
        # We don't need MAG labels beside every symbol.
        # The tree itself contains MAG/taxon information.
        # ----------------------------------------------

        outfile.write(
            "SHOW_LABELS\t0\n"
        )


        # ----------------------------------------------
        # Legend
        # ----------------------------------------------

        outfile.write(
            f"LEGEND_TITLE\t"
            f"{label}\n"
        )


        outfile.write(
            f"LEGEND_SHAPES\t"
            f"{symbol}\t"
            f"{symbol}\t"
            f"{symbol}\n"
        )


        outfile.write(
            f"LEGEND_COLORS\t"
            f"{color}\t"
            f"{color}\t"
            f"{color}\n"
        )


        outfile.write(
            "LEGEND_LABELS\t"
            "Low (0.0)\t"
            "Intermediate (0.5)\t"
            "High (1.0)\n"
        )


        outfile.write(
            "LEGEND_SHAPE_SCALES\t"
            "0.10\t"
            "0.55\t"
            "1.00\n"
        )


        outfile.write(
            "LEGEND_SHAPE_INVERT\t"
            "0\t0\t0\n"
        )


        outfile.write(
            "\nDATA\n"
        )


        # ----------------------------------------------
        # One symbol for every bacterial MAG
        # ----------------------------------------------

        for mag in tree_tips:

            score = scores[
                mag
            ]


            size = score_to_size(
                score
            )


            label_text = (
                f"{label}: "
                f"{score:.3f}"
            )


            outfile.write(
                f"{mag}\t"
                f"{symbol}\t"
                f"{size:.4f}\t"
                f"{color}\t"
                f"1\t"
                f"-1\t"
                f"{label_text}\n"
            )


# ============================================================
# WRITE COMBINED STRATEGY MATRIX
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
# PRINT SCORE SUMMARY
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
        f"  Ecological dimension: "
        f"{info['dimension']}"
    )

    print(
        f"  Symbol: "
        f"{info['symbol']}"
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
    "Visual grammar:"
)

print()

print(
    "  CIRCLE = Environmental tolerance"
)

print(
    "    Osmotic adaptation"
)

print(
    "    Redox flexibility"
)

print()

print(
    "  SQUARE = Resource acquisition"
)

print(
    "    Nitrogen acquisition"
)

print(
    "    Carbon + phosphorus"
)

print(
    "    Iron acquisition [future]"
)

print()

print(
    "  TRIANGLE = Host interaction"
)

print(
    "    Colonisation / PGP [future]"
)

print()

print(
    "  Colour = strategy identity"
)

print(
    "  Symbol size = strategy breadth score (0-1)"
)

print(
    "  Symbol size scaling is identical across all rings."
)

print()

print(
    f"Combined score matrix:"
)

print(
    f"  {COMBINED_MATRIX}"
)

print(
    "=" * 72
)
