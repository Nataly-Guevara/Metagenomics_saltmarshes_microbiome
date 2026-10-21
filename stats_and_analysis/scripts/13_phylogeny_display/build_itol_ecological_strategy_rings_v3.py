#!/usr/bin/env python3

"""
Build coordinated iTOL ecological-strategy rings - V3.

PURPOSE
-------
Generate a cleaner Figure 1 phylogenetic display of ecological
strategy breadth across the 352 bacterial MAGs.

BIOLOGICAL SCORES ARE NOT CHANGED.

Only the graphical representation is changed.

VISUAL GRAMMAR
==============

Broad ecological dimension -> shape

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
    - Colonisation / PGP [future]


Individual strategy -> fixed colour

Strategy breadth -> symbol size


THRESHOLDED DISPLAY
===================

score = 0
    no symbol
    = none of the predefined capabilities detected

0 < score <= 0.25
    tiny symbol
    = limited capability breadth

0.25 < score <= 0.50
    medium symbol
    = intermediate capability breadth

score > 0.50
    large symbol
    = broad capability repertoire

Within the high category, symbol size still increases slightly
toward score = 1.

The exact continuous score is retained in the iTOL tooltip.
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
      "itol_strategy_rings_v3"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

COMBINED_MATRIX = (
    OUTPUT_DIR
    / "ecological_strategy_scores_352_MAGs.tsv"
)

LEGEND_FILE = (
    OUTPUT_DIR
    / "figure1_strategy_legend.txt"
)


# ============================================================
# iTOL SYMBOL CODES
# ============================================================

# iTOL symbol codes:
#
# 1 = square
# 2 = circle
# 4 = right triangle

SQUARE = 1
CIRCLE = 2
TRIANGLE = 4


# ============================================================
# STRATEGY DEFINITIONS
# ============================================================

STRATEGIES = {

    "osmotic": {

        "label":
            "Osmotic adaptation",

        "short_label":
            "Osmotic",

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

        "short_label":
            "Redox",

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

        "short_label":
            "Nitrogen",

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

        "short_label":
            "Carbon + P",

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
# FUTURE ASSIGNMENTS
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
# DISPLAY THRESHOLDS
# ============================================================

LOW_THRESHOLD = 0.25
INTERMEDIATE_THRESHOLD = 0.50

# Symbol sizes.
#
# These are deliberately smaller than V2 because osmotic
# adaptation was visually dominating the circular tree.

LOW_SIZE = 1.5
INTERMEDIATE_SIZE = 4.0

HIGH_MIN_SIZE = 6.0
HIGH_MAX_SIZE = 8.0

ITOL_MAXIMUM_SIZE = 22

SYMBOL_SPACING = 5


# ============================================================
# DISPLAY FUNCTIONS
# ============================================================

def breadth_class(score):

    """
    Convert continuous strategy score to an interpretive class.
    """

    score = float(
        score
    )

    if score == 0:

        return "Absent"

    elif score <= LOW_THRESHOLD:

        return "Low"

    elif score <= INTERMEDIATE_THRESHOLD:

        return "Intermediate"

    else:

        return "High"


def score_to_size(score):

    """
    Thresholded visual size.

    score = 0
        returns None -> no symbol written

    0 < score <= 0.25
        tiny symbol

    0.25 < score <= 0.50
        medium symbol

    > 0.50
        large symbol, increasing toward score = 1
    """

    score = float(
        score
    )

    if score <= 0:

        return None


    if score <= LOW_THRESHOLD:

        return LOW_SIZE


    if score <= INTERMEDIATE_THRESHOLD:

        return INTERMEDIATE_SIZE


    # ----------------------------------------------
    # High breadth:
    #
    # 0.50 -> approximately HIGH_MIN_SIZE
    # 1.00 -> HIGH_MAX_SIZE
    # ----------------------------------------------

    high_fraction = (
        score
        - INTERMEDIATE_THRESHOLD
    ) / (
        1.0
        - INTERMEDIATE_THRESHOLD
    )

    return (
        HIGH_MIN_SIZE
        +
        high_fraction
        * (
            HIGH_MAX_SIZE
            - HIGH_MIN_SIZE
        )
    )


# ============================================================
# READ FINAL STRATEGY TABLE
# ============================================================

def read_strategy_table(
    table_path,
    score_column
):

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


            score = float(
                value
            )


            if (
                score < 0
                or score > 1
            ):

                raise ValueError(
                    f"Score outside 0-1 range: "
                    f"{mag} | {score}"
                )


            scores[
                mag
            ] = score


    return scores


# ============================================================
# READ 352-TIP BACTERIAL TREE
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
# READ STRATEGY SCORES
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
# TREE MATCHING QC
# ============================================================

print()
print(
    "=" * 72
)

print(
    "iTOL ECOLOGICAL STRATEGY RINGS V3"
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
# REQUIRE COMPLETE TREE MATCH
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
# WRITE iTOL RINGS
# ============================================================

for key, info in (
    STRATEGIES.items()
):

    label = info[
        "label"
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


        # No labels beside each symbol.
        outfile.write(
            "SHOW_LABELS\t0\n"
        )


        # ----------------------------------------------------
        # IMPORTANT:
        #
        # No individual legend is added here.
        #
        # The final publication figure will use ONE shared
        # legend explaining:
        #
        # shape = ecological dimension
        # colour = strategy identity
        # size = strategy breadth
        # ----------------------------------------------------


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


            # ----------------------------------------------
            # TRUE ABSENCE:
            #
            # no symbol is written.
            # ----------------------------------------------

            if size is None:

                continue


            category = breadth_class(
                score
            )


            tooltip = (
                f"{label}: "
                f"{score:.3f} "
                f"({category} breadth)"
            )


            outfile.write(
                f"{mag}\t"
                f"{symbol}\t"
                f"{size:.4f}\t"
                f"{color}\t"
                f"1\t"
                f"-1\t"
                f"{tooltip}\n"
            )


# ============================================================
# WRITE COMBINED MATRIX
# ============================================================

matrix_fields = [
    "Final_MAG",

    "Osmotic_adaptation",
    "Osmotic_class",

    "Redox_flexibility",
    "Redox_class",

    "Nitrogen_acquisition_transformation",
    "Nitrogen_class",

    "Carbon_phosphorus_resource",
    "Carbon_phosphorus_class",
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

        osmotic = strategy_scores[
            "osmotic"
        ][
            mag
        ]

        redox = strategy_scores[
            "redox"
        ][
            mag
        ]

        nitrogen = strategy_scores[
            "nitrogen"
        ][
            mag
        ]

        carbon_p = strategy_scores[
            "carbon_phosphorus"
        ][
            mag
        ]


        writer.writerow({

            "Final_MAG":
                mag,

            "Osmotic_adaptation":
                f"{osmotic:.4f}",

            "Osmotic_class":
                breadth_class(
                    osmotic
                ),

            "Redox_flexibility":
                f"{redox:.4f}",

            "Redox_class":
                breadth_class(
                    redox
                ),

            "Nitrogen_acquisition_transformation":
                f"{nitrogen:.4f}",

            "Nitrogen_class":
                breadth_class(
                    nitrogen
                ),

            "Carbon_phosphorus_resource":
                f"{carbon_p:.4f}",

            "Carbon_phosphorus_class":
                breadth_class(
                    carbon_p
                ),
        })


# ============================================================
# WRITE SHARED FIGURE-LEGEND TEXT
# ============================================================

with LEGEND_FILE.open(
    "w"
) as outfile:

    outfile.write(
        "FIGURE 1 ECOLOGICAL STRATEGY LEGEND\n"
    )

    outfile.write(
        "===================================\n\n"
    )

    outfile.write(
        "SHAPE = broad ecological dimension\n\n"
    )

    outfile.write(
        "Circle   Environmental tolerance\n"
    )

    outfile.write(
        "Square   Resource acquisition\n"
    )

    outfile.write(
        "Triangle Host interaction\n\n"
    )


    outfile.write(
        "COLOUR = ecological strategy\n\n"
    )

    for key, info in (
        STRATEGIES.items()
    ):

        outfile.write(
            f"{info['color']}  "
            f"{info['label']}\n"
        )


    outfile.write(
        "\nFuture:\n"
    )

    outfile.write(
        "#8E5EA2  Iron acquisition\n"
    )

    outfile.write(
        "#3A923A  Colonisation / PGP\n\n"
    )


    outfile.write(
        "SYMBOL SIZE = standardized ecological capability breadth\n\n"
    )

    outfile.write(
        "No symbol   score = 0; no predefined capability detected\n"
    )

    outfile.write(
        "Tiny        0 < score <= 0.25; low breadth\n"
    )

    outfile.write(
        "Medium      0.25 < score <= 0.50; intermediate breadth\n"
    )

    outfile.write(
        "Large       score > 0.50; high breadth\n\n"
    )

    outfile.write(
        "Exact continuous scores are retained in the underlying "
        "data and iTOL tooltips.\n"
    )


# ============================================================
# PRINT FINAL QC AND DISTRIBUTIONS
# ============================================================

print(
    "=" * 72
)

print(
    "iTOL STRATEGY RINGS V3 CREATED SUCCESSFULLY"
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


    class_counts = {
        "Absent": 0,
        "Low": 0,
        "Intermediate": 0,
        "High": 0,
    }


    for score in scores:

        class_counts[
            breadth_class(
                score
            )
        ] += 1


    visible_symbols = (
        len(
            scores
        )
        - class_counts[
            "Absent"
        ]
    )


    print(
        f"Ring {number}: "
        f"{info['label']}"
    )

    print(
        f"  Dimension: "
        f"{info['dimension']}"
    )

    print(
        f"  Colour: "
        f"{info['color']}"
    )

    print(
        f"  Absent / no symbol: "
        f"{class_counts['Absent']}"
    )

    print(
        f"  Low / tiny:         "
        f"{class_counts['Low']}"
    )

    print(
        f"  Intermediate:       "
        f"{class_counts['Intermediate']}"
    )

    print(
        f"  High / large:       "
        f"{class_counts['High']}"
    )

    print(
        f"  Visible symbols:    "
        f"{visible_symbols}/"
        f"{len(scores)}"
    )

    print(
        f"  File: "
        f"{info['output']}"
    )

    print()


print(
    "Thresholds:"
)

print(
    "  Absent:       score = 0"
)

print(
    "  Low:          0 < score <= 0.25"
)

print(
    "  Intermediate: 0.25 < score <= 0.50"
)

print(
    "  High:         score > 0.50"
)

print()


print(
    "IMPORTANT:"
)

print(
    "  No symbol means true absence of all predefined "
    "capabilities for that strategy."
)

print(
    "  Low non-zero scores remain visible as tiny symbols."
)

print(
    "  Biological strategy scores have not been changed."
)

print()


print(
    f"Combined matrix:"
)

print(
    f"  {COMBINED_MATRIX}"
)

print(
    f"Shared legend:"
)

print(
    f"  {LEGEND_FILE}"
)

print(
    "=" * 72
)
