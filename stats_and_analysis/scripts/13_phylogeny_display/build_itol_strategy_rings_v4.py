#!/usr/bin/env python3

"""
Build coordinated iTOL ecological-strategy rings - V4.

V4 changes ONLY the graphical symbol sizes.

BIOLOGICAL SCORES AND THRESHOLDS ARE IDENTICAL TO V3.

Visual categories
-----------------
score = 0
    no symbol

0 < score <= 0.25
    low breadth, tiny symbol

0.25 < score <= 0.50
    intermediate breadth, medium symbol

score > 0.50
    high breadth, larger symbol

Within the high category, symbol size increases gently with
the exact continuous strategy score.

Visual grammar
--------------
Circle   = Environmental tolerance
Square   = Resource acquisition
Triangle = Host interaction

Colour = individual ecological strategy
Size   = standardized ecological capability breadth
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
      "itol_strategy_rings_v4"
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
# DISPLAY THRESHOLDS
# ============================================================

LOW_THRESHOLD = 0.25
INTERMEDIATE_THRESHOLD = 0.50


# ============================================================
# SMALLER SYMBOL SIZES FOR V4
# ============================================================

LOW_SIZE = 0.8

INTERMEDIATE_SIZE = 2.2

HIGH_MIN_SIZE = 3.5

HIGH_MAX_SIZE = 5.0


# Smaller overall iTOL maximum than V3
ITOL_MAXIMUM_SIZE = 16

SYMBOL_SPACING = 5


# ============================================================
# DISPLAY FUNCTIONS
# ============================================================

def breadth_class(score):

    score = float(score)

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
    Convert continuous score to V4 plotting size.

    Biological thresholds are unchanged from V3.

    score = 0
        no symbol

    0 < score <= 0.25
        size = 0.8

    0.25 < score <= 0.50
        size = 2.2

    score > 0.50
        size ranges continuously from 3.5 to 5.0
    """

    score = float(score)


    if score <= 0:

        return None


    if score <= LOW_THRESHOLD:

        return LOW_SIZE


    if score <= INTERMEDIATE_THRESHOLD:

        return INTERMEDIATE_SIZE


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
# READ STRATEGY TABLE
# ============================================================

def read_strategy_table(
    table_path,
    score_column
):

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
                f"No header found in {table_path}"
            )


        required = {
            "Final_MAG",
            score_column,
        }


        missing = (
            required
            - set(reader.fieldnames)
        )


        if missing:

            raise ValueError(
                f"Missing columns in {table_path}: "
                + ", ".join(sorted(missing))
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


            score = float(value)


            if (
                score < 0
                or score > 1
            ):

                raise ValueError(
                    f"Score outside 0-1 range: "
                    f"{mag} | {score}"
                )


            scores[mag] = score


    return scores


# ============================================================
# READ 352-TIP BACTERIAL TREE
# ============================================================

if not TREE.exists():

    raise FileNotFoundError(
        f"Bacterial tree not found: {TREE}"
    )


newick = TREE.read_text().strip()


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
        and tip not in seen_tree_tips
    ):

        tree_tips.append(tip)

        seen_tree_tips.add(tip)


if len(tree_tips) != 352:

    print(
        "WARNING: expected 352 bacterial "
        f"tree tips, found {len(tree_tips)}."
    )


# ============================================================
# READ STRATEGY SCORES
# ============================================================

strategy_scores = {}


for key, info in STRATEGIES.items():

    strategy_scores[key] = read_strategy_table(
        info["table"],
        info["score_column"]
    )


# ============================================================
# TREE MATCHING QC
# ============================================================

print()

print("=" * 72)

print(
    "iTOL ECOLOGICAL STRATEGY RINGS V4"
)

print("=" * 72)

print()

print(
    f"Bacterial tree tips: {len(tree_tips)}"
)

print()


for key, info in STRATEGIES.items():

    scores = strategy_scores[key]


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

for key, info in STRATEGIES.items():

    scores = strategy_scores[key]


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

for key, info in STRATEGIES.items():

    label = info["label"]

    symbol = info["symbol"]

    color = info["color"]

    output = info["output"]

    scores = strategy_scores[key]


    with output.open("w") as outfile:


        outfile.write(
            "DATASET_SYMBOL\n"
        )

        outfile.write(
            "SEPARATOR TAB\n"
        )

        outfile.write(
            f"DATASET_LABEL\t{label}\n"
        )

        outfile.write(
            f"COLOR\t{color}\n"
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
            "SHOW_LABELS\t0\n"
        )


        outfile.write(
            "\nDATA\n"
        )


        for mag in tree_tips:

            score = scores[mag]

            size = score_to_size(
                score
            )


            # True absence:
            # no symbol written

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

        osmotic = (
            strategy_scores[
                "osmotic"
            ][mag]
        )

        redox = (
            strategy_scores[
                "redox"
            ][mag]
        )

        nitrogen = (
            strategy_scores[
                "nitrogen"
            ][mag]
        )

        carbon_p = (
            strategy_scores[
                "carbon_phosphorus"
            ][mag]
        )


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
# WRITE SHARED LEGEND TEXT
# ============================================================

with LEGEND_FILE.open("w") as outfile:

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


    for key, info in STRATEGIES.items():

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
        "SYMBOL SIZE = standardized ecological "
        "capability breadth\n\n"
    )

    outfile.write(
        "No symbol   score = 0\n"
    )

    outfile.write(
        "Tiny        0 < score <= 0.25\n"
    )

    outfile.write(
        "Medium      0.25 < score <= 0.50\n"
    )

    outfile.write(
        "Large       score > 0.50\n"
    )


# ============================================================
# FINAL QC
# ============================================================

print(
    "=" * 72
)

print(
    "iTOL STRATEGY RINGS V4 CREATED SUCCESSFULLY"
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
        strategy_scores[key][mag]
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
            breadth_class(score)
        ] += 1


    visible_symbols = (
        len(scores)
        - class_counts["Absent"]
    )


    print(
        f"Ring {number}: "
        f"{info['label']}"
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
    "V4 plotting sizes:"
)

print(
    f"  Low:          {LOW_SIZE}"
)

print(
    f"  Intermediate: {INTERMEDIATE_SIZE}"
)

print(
    f"  High:         "
    f"{HIGH_MIN_SIZE}-{HIGH_MAX_SIZE}"
)

print()

print(
    "Biological thresholds are unchanged from V3."
)

print(
    "Only the graphical symbol sizes were reduced."
)

print()

print(
    f"Combined matrix:"
)

print(
    f"  {COMBINED_MATRIX}"
)

print()

print(
    f"Shared legend:"
)

print(
    f"  {LEGEND_FILE}"
)

print(
    "=" * 72
)
