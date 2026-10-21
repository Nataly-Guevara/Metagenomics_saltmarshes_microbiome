#!/usr/bin/env python3

"""
Create a preliminary iTOL square-symbol dataset for redox flexibility.

Each MAG is represented by one filled square.

Visual encoding
---------------
Shape:
    Square: iTOL symbol 1

Size:
    Proportional to the preliminary genome-level redox score.

Colour:
    Score 0.0 = #F2F2F2
    Score 0.5 = #8FB7D9
    Score 1.0 = #1F4E79

Position:
    -2, placing redox flexibility beside the osmotic-adaptation
    circle previously assigned to external position -1.

Redox capabilities
------------------
1. Nitrate respiration
2. Nitrite and/or nitric-oxide reduction
3. DNRA
4. Nitrous-oxide reduction
5. Dissimilatory sulfate reduction
6. Sulfide oxidation
7. Sulfite oxidation
8. Low-oxygen respiration

Input
-----
redox_capability_presence_preliminary.tsv

Outputs
-------
1. itol_redox_flexibility_square_symbols_preliminary.txt
2. itol_redox_flexibility_square_symbols_audit.tsv

Important
---------
- This visualization is preliminary.
- No individual-MAG completeness correction is applied.
- The script validates the 0-8 count and 0-1 score directly from
  the eight binary capability calls.
- Existing osmotic and redox audit files are not modified.
- Only the Python standard library is used.
"""

import csv
import math
import os
import sys
from collections import Counter


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = "/path/to/metagenomics_project"

INPUT_FILE = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/redox_flexibility/tables/"
    "redox_capability_presence_preliminary.tsv",
)

ITOL_OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/redox_flexibility/itol",
)

TABLE_OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/redox_flexibility/tables",
)

ITOL_OUTPUT_FILE = os.path.join(
    ITOL_OUTPUT_DIR,
    "itol_redox_flexibility_square_symbols_preliminary.txt",
)

AUDIT_OUTPUT_FILE = os.path.join(
    TABLE_OUTPUT_DIR,
    "itol_redox_flexibility_square_symbols_audit.tsv",
)


# ---------------------------------------------------------------------
# Redox capability columns
# ---------------------------------------------------------------------

CAPABILITY_COLUMNS = [
    "Redox_Nitrate_respiration_preliminary",
    "Redox_Nitrite_NO_reduction_preliminary",
    "Redox_DNRA_preliminary",
    "Redox_Nitrous_oxide_reduction_preliminary",
    "Redox_Sulfate_reduction_preliminary",
    "Redox_Sulfide_oxidation_preliminary",
    "Redox_Sulfite_oxidation_preliminary",
    "Redox_Low_oxygen_respiration_strict",
]

CAPABILITY_LABELS = [
    "Nitrate respiration",
    "Nitrite and nitric-oxide reduction",
    "DNRA",
    "Nitrous-oxide reduction",
    "Dissimilatory sulfate reduction",
    "Sulfide oxidation",
    "Sulfite oxidation",
    "Low-oxygen respiration",
]


# ---------------------------------------------------------------------
# Input score columns
# ---------------------------------------------------------------------

ID_COLUMN = "Final_MAG"

COUNT_COLUMN = (
    "Redox_preliminary_capability_count_0_8"
)

SCORE_COLUMN = (
    "Redox_preliminary_strategy_score_0_1"
)

COMPLETENESS_COLUMN = "Completeness"


# ---------------------------------------------------------------------
# iTOL settings
# ---------------------------------------------------------------------

SYMBOL = 1
FILL = 1
EXTERNAL_POSITION = -2

MIN_DATA_SIZE = 1.0
MAX_DATA_SIZE = 10.0

ITOL_MAXIMUM_SIZE = 28

DATASET_LABEL = "Redox flexibility preliminary"
DATASET_COLOR = "#1F4E79"

COLOR_MIN = "#F2F2F2"
COLOR_MID = "#8FB7D9"
COLOR_MAX = "#1F4E79"


# ---------------------------------------------------------------------
# Parsing and validation
# ---------------------------------------------------------------------

def clean_text(value):
    """Return a stripped string."""
    if value is None:
        return ""

    return str(value).strip()


def parse_float(value, column, mag):
    """Parse a finite numerical value."""
    cleaned = clean_text(value)

    try:
        number = float(cleaned)
    except ValueError:
        sys.exit(
            "ERROR: Non-numeric value encountered.\n"
            f"MAG: {mag}\n"
            f"Column: {column}\n"
            f"Value: {value}"
        )

    if not math.isfinite(number):
        sys.exit(
            "ERROR: Non-finite value encountered.\n"
            f"MAG: {mag}\n"
            f"Column: {column}\n"
            f"Value: {value}"
        )

    return number


def parse_binary(value, column, mag):
    """Parse a binary field containing zero or one."""
    number = parse_float(
        value,
        column,
        mag,
    )

    if number == 0:
        return 0

    if number == 1:
        return 1

    sys.exit(
        "ERROR: Capability field is not binary.\n"
        f"MAG: {mag}\n"
        f"Column: {column}\n"
        f"Value: {value}"
    )


def parse_integer(value, column, mag):
    """Parse an integer-valued field."""
    number = parse_float(
        value,
        column,
        mag,
    )

    integer = int(round(number))

    if number != integer:
        sys.exit(
            "ERROR: Expected an integer value.\n"
            f"MAG: {mag}\n"
            f"Column: {column}\n"
            f"Value: {value}"
        )

    return integer


def validate_required_columns(fieldnames):
    """Confirm all required input columns are present."""
    if fieldnames is None:
        sys.exit(
            "ERROR: Input table does not contain a header."
        )

    required = {
        ID_COLUMN,
        COUNT_COLUMN,
        SCORE_COLUMN,
        COMPLETENESS_COLUMN,
    }

    required.update(
        CAPABILITY_COLUMNS
    )

    missing = sorted(
        required.difference(fieldnames)
    )

    if missing:
        sys.exit(
            "ERROR: Required columns are missing:\n"
            + "\n".join(missing)
        )


# ---------------------------------------------------------------------
# Colour and size conversion
# ---------------------------------------------------------------------

def hex_to_rgb(hex_colour):
    """Convert a hexadecimal colour to an RGB tuple."""
    cleaned = hex_colour.lstrip("#")

    if len(cleaned) != 6:
        raise ValueError(
            f"Invalid hexadecimal colour: {hex_colour}"
        )

    return tuple(
        int(cleaned[index:index + 2], 16)
        for index in (0, 2, 4)
    )


def rgb_to_hex(rgb):
    """Convert an RGB tuple to uppercase hexadecimal colour."""
    bounded = tuple(
        max(
            0,
            min(
                255,
                int(round(value)),
            ),
        )
        for value in rgb
    )

    return "#{:02X}{:02X}{:02X}".format(
        bounded[0],
        bounded[1],
        bounded[2],
    )


def interpolate_rgb(start_rgb, end_rgb, fraction):
    """Linearly interpolate between two RGB colours."""
    fraction = max(
        0.0,
        min(
            1.0,
            fraction,
        ),
    )

    return tuple(
        start + (end - start) * fraction
        for start, end in zip(
            start_rgb,
            end_rgb,
        )
    )


def score_to_colour(score):
    """
    Convert a score from 0 to 1 into the requested three-point
    colour gradient.
    """
    min_rgb = hex_to_rgb(
        COLOR_MIN
    )

    mid_rgb = hex_to_rgb(
        COLOR_MID
    )

    max_rgb = hex_to_rgb(
        COLOR_MAX
    )

    if score <= 0.5:
        local_fraction = score / 0.5

        interpolated = interpolate_rgb(
            min_rgb,
            mid_rgb,
            local_fraction,
        )

    else:
        local_fraction = (
            score - 0.5
        ) / 0.5

        interpolated = interpolate_rgb(
            mid_rgb,
            max_rgb,
            local_fraction,
        )

    return rgb_to_hex(
        interpolated
    )


def score_to_size(score):
    """
    Convert score 0-1 to a positive relative symbol size.

    Score 0 receives size 1 so zero-score MAGs remain visible.
    Score 1 receives size 10.
    """
    return (
        MIN_DATA_SIZE
        + score
        * (
            MAX_DATA_SIZE
            - MIN_DATA_SIZE
        )
    )


# ---------------------------------------------------------------------
# Read, recalculate and validate scores
# ---------------------------------------------------------------------

def read_redox_scores():
    """
    Recalculate the eight-capability count and proportional score
    for every MAG, then derive the square properties.
    """
    try:
        handle = open(
            INPUT_FILE,
            "r",
            encoding="utf-8-sig",
            newline="",
        )
    except OSError as exc:
        sys.exit(
            "ERROR: Could not open input file:\n"
            f"{INPUT_FILE}\n"
            f"{exc}"
        )

    output_rows = []
    seen_mags = set()

    with handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        validate_required_columns(
            reader.fieldnames
        )

        for row in reader:
            mag = clean_text(
                row.get(ID_COLUMN)
            )

            if not mag:
                continue

            if mag in seen_mags:
                sys.exit(
                    "ERROR: Duplicate MAG identifier found:\n"
                    f"{mag}"
                )

            seen_mags.add(
                mag
            )

            capability_values = []

            for column in CAPABILITY_COLUMNS:
                capability_values.append(
                    parse_binary(
                        row[column],
                        column,
                        mag,
                    )
                )

            recalculated_count = sum(
                capability_values
            )

            recalculated_score = (
                recalculated_count
                / len(CAPABILITY_COLUMNS)
            )

            existing_count = parse_integer(
                row[COUNT_COLUMN],
                COUNT_COLUMN,
                mag,
            )

            existing_score = parse_float(
                row[SCORE_COLUMN],
                SCORE_COLUMN,
                mag,
            )

            count_matches = (
                recalculated_count == existing_count
            )

            score_matches = math.isclose(
                recalculated_score,
                existing_score,
                rel_tol=0.0,
                abs_tol=0.00005,
            )

            if not count_matches or not score_matches:
                sys.exit(
                    "ERROR: Recalculated redox score does not match "
                    "the preliminary presence matrix.\n"
                    f"MAG: {mag}\n"
                    f"Recalculated count: {recalculated_count}\n"
                    f"Existing count: {existing_count}\n"
                    f"Recalculated score: {recalculated_score:.4f}\n"
                    f"Existing score: {existing_score:.4f}"
                )

            completeness = parse_float(
                row[COMPLETENESS_COLUMN],
                COMPLETENESS_COLUMN,
                mag,
            )

            square_size = score_to_size(
                recalculated_score
            )

            square_colour = score_to_colour(
                recalculated_score
            )

            output_row = {
                "Final_MAG": mag,
                "Redox_capability_count_0_8": (
                    recalculated_count
                ),
                "Redox_strategy_score_0_1": (
                    recalculated_score
                ),
                "Square_relative_size": square_size,
                "Square_colour": square_colour,
                "Completeness": completeness,
                "Count_matches_preliminary_matrix": (
                    "Yes" if count_matches else "No"
                ),
                "Score_matches_preliminary_matrix": (
                    "Yes" if score_matches else "No"
                ),
                "Individual_MAG_completeness_correction": (
                    "Not_applied"
                ),
            }

            for column, value in zip(
                CAPABILITY_COLUMNS,
                capability_values,
            ):
                output_row[column] = value

            output_rows.append(
                output_row
            )

    if not output_rows:
        sys.exit(
            "ERROR: No MAG score rows were recovered."
        )

    return output_rows


# ---------------------------------------------------------------------
# Write iTOL square-symbol dataset
# ---------------------------------------------------------------------

def write_itol_symbol_file(rows):
    """
    Write the iTOL DATASET_SYMBOL file.

    Data columns:
        ID
        symbol
        size
        colour
        fill
        position
        label
    """
    with open(
        ITOL_OUTPUT_FILE,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        handle.write("DATASET_SYMBOL\n")
        handle.write("SEPARATOR TAB\n")

        handle.write(
            f"DATASET_LABEL\t{DATASET_LABEL}\n"
        )

        handle.write(
            f"COLOR\t{DATASET_COLOR}\n"
        )

        handle.write(
            f"MAXIMUM_SIZE\t{ITOL_MAXIMUM_SIZE}\n"
        )

        handle.write("SYMBOL_SPACING\t5\n")
        handle.write("SHOW_LABELS\t1\n")
        handle.write("LABEL_SIZE_FACTOR\t1\n")
        handle.write("LABEL_ROTATION\t0\n")
        handle.write("LABEL_ALIGN_TO_TREE\t1\n")

        handle.write(
            "LEGEND_TITLE\t"
            "Redox flexibility score preliminary\n"
        )

        handle.write(
            "LEGEND_SHAPES\t1\t1\t1\n"
        )

        handle.write(
            f"LEGEND_COLORS\t"
            f"{COLOR_MIN}\t{COLOR_MID}\t{COLOR_MAX}\n"
        )

        handle.write(
            "LEGEND_LABELS\t0.0\t0.5\t1.0\n"
        )

        handle.write(
            "LEGEND_SHAPE_SCALES\t0.10\t0.55\t1.00\n"
        )

        handle.write(
            "LEGEND_SHAPE_INVERT\t0\t0\t0\n"
        )

        handle.write("\n")
        handle.write("DATA\n")

        for row in rows:
            mag = row["Final_MAG"]

            score = row[
                "Redox_strategy_score_0_1"
            ]

            square_size = row[
                "Square_relative_size"
            ]

            square_colour = row[
                "Square_colour"
            ]

            label = (
                "Redox flexibility preliminary: "
                f"{score:.3f}"
            )

            handle.write(
                f"{mag}\t"
                f"{SYMBOL}\t"
                f"{square_size:.4f}\t"
                f"{square_colour}\t"
                f"{FILL}\t"
                f"{EXTERNAL_POSITION}\t"
                f"{label}\n"
            )


# ---------------------------------------------------------------------
# Write audit table
# ---------------------------------------------------------------------

def write_audit(rows):
    """Write the score-to-symbol conversion audit."""
    columns = (
        [
            "Final_MAG",
        ]
        + CAPABILITY_COLUMNS
        + [
            "Redox_capability_count_0_8",
            "Redox_strategy_score_0_1",
            "Square_relative_size",
            "Square_colour",
            "Completeness",
            "Count_matches_preliminary_matrix",
            "Score_matches_preliminary_matrix",
            "Individual_MAG_completeness_correction",
        ]
    )

    with open(
        AUDIT_OUTPUT_FILE,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            output_row = dict(
                row
            )

            output_row[
                "Redox_strategy_score_0_1"
            ] = (
                f"{row['Redox_strategy_score_0_1']:.4f}"
            )

            output_row[
                "Square_relative_size"
            ] = (
                f"{row['Square_relative_size']:.4f}"
            )

            output_row[
                "Completeness"
            ] = (
                f"{row['Completeness']:.4f}"
            )

            writer.writerow(
                output_row
            )


# ---------------------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------------------

def print_summary(rows):
    """Print output paths and score-to-square mapping."""
    score_counts = Counter(
        row["Redox_capability_count_0_8"]
        for row in rows
    )

    score_examples = {}

    for row in rows:
        count = row[
            "Redox_capability_count_0_8"
        ]

        if count not in score_examples:
            score_examples[count] = {
                "score": row[
                    "Redox_strategy_score_0_1"
                ],
                "size": row[
                    "Square_relative_size"
                ],
                "colour": row[
                    "Square_colour"
                ],
            }

    print(
        "Preliminary redox square-symbol export completed."
    )
    print()
    print(f"MAGs exported: {len(rows)}")
    print("iTOL symbol: 1 (square)")
    print("External symbol position: -2")
    print(
        "Individual MAG completeness correction applied: No"
    )
    print()
    print(
        f"iTOL square dataset: {ITOL_OUTPUT_FILE}"
    )
    print(
        f"Square conversion audit: {AUDIT_OUTPUT_FILE}"
    )
    print()
    print("Score-to-square mapping")
    print("-----------------------")

    for count in range(0, 9):
        example = score_examples.get(
            count
        )

        if example is None:
            continue

        print(
            f"{count}/8 "
            f"(score {example['score']:.3f}): "
            f"{score_counts[count]} MAGs, "
            f"relative size {example['size']:.2f}, "
            f"colour {example['colour']}"
        )

    print()
    print(
        "This is a preliminary visualization pending final review "
        "of the redox capability rules."
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    os.makedirs(
        ITOL_OUTPUT_DIR,
        exist_ok=True,
    )

    os.makedirs(
        TABLE_OUTPUT_DIR,
        exist_ok=True,
    )

    rows = read_redox_scores()

    write_itol_symbol_file(
        rows
    )

    write_audit(
        rows
    )

    print_summary(
        rows
    )


if __name__ == "__main__":
    main()

