#!/usr/bin/env python3

"""
Create an iTOL circle-symbol dataset for osmotic adaptation.

Each MAG is represented by one filled circle.

Visual encoding
---------------
Shape:
    Circle (iTOL symbol 2)

Size:
    Proportional to the genome-level osmotic strategy score.

Colour:
    Score 0.0 = #F2F2F2
    Score 0.5 = #8FB7D9
    Score 1.0 = #1F4E79

Position:
    -1, meaning the first external symbol column outside the tree.

Input
-----
osmotic_ricci_capability_matrix.tsv

Output
------
itol_osmotic_adaptation_circle_symbols.txt

Audit
-----
itol_osmotic_adaptation_circle_symbols_audit.tsv

Important
---------
- No individual-MAG completeness correction is applied.
- Existing heatmap and capability files are not modified.
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
    "results/ecological_strategies/osmotic_adaptation/tables/"
    "osmotic_ricci_capability_matrix.tsv",
)

ITOL_OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/itol",
)

TABLE_OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/tables",
)

ITOL_OUTPUT_FILE = os.path.join(
    ITOL_OUTPUT_DIR,
    "itol_osmotic_adaptation_circle_symbols.txt",
)

AUDIT_OUTPUT_FILE = os.path.join(
    TABLE_OUTPUT_DIR,
    "itol_osmotic_adaptation_circle_symbols_audit.tsv",
)


# ---------------------------------------------------------------------
# Input columns
# ---------------------------------------------------------------------

ID_COLUMN = "Final_MAG"
SCORE_COLUMN = "Osmotic_strategy_score_0_1"
COUNT_COLUMN = "Osmotic_capability_count_0_5"
COMPLETENESS_COLUMN = "Completeness"


# ---------------------------------------------------------------------
# iTOL symbol settings
# ---------------------------------------------------------------------

SYMBOL = 2
FILL = 1
EXTERNAL_POSITION = -1

# Relative sizes written in each data row.
# iTOL scales the largest value to MAXIMUM_SIZE below.
MIN_DATA_SIZE = 1.0
MAX_DATA_SIZE = 10.0

# Largest displayed symbol size in iTOL.
ITOL_MAXIMUM_SIZE = 28

DATASET_LABEL = "Osmotic adaptation"
DATASET_COLOR = "#1F4E79"


# ---------------------------------------------------------------------
# Three-point colour gradient
# ---------------------------------------------------------------------

COLOR_MIN = "#F2F2F2"
COLOR_MID = "#8FB7D9"
COLOR_MAX = "#1F4E79"


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


def hex_to_rgb(hex_colour):
    """Convert hexadecimal colour to an RGB tuple."""
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
        max(0, min(255, int(round(value))))
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
        min(1.0, fraction),
    )

    return tuple(
        start + (end - start) * fraction
        for start, end in zip(start_rgb, end_rgb)
    )


def score_to_colour(score):
    """
    Convert a score from 0 to 1 into the requested three-point
    colour gradient.

    0.0 to 0.5:
        COLOR_MIN -> COLOR_MID

    0.5 to 1.0:
        COLOR_MID -> COLOR_MAX
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
    Convert score 0-1 into a positive relative symbol size.

    Score 0.0 receives size 1.
    Score 1.0 receives size 10.

    A nonzero minimum ensures that MAGs with score zero remain visible.
    """
    return (
        MIN_DATA_SIZE
        + score * (
            MAX_DATA_SIZE - MIN_DATA_SIZE
        )
    )


def validate_required_columns(fieldnames):
    """Confirm required columns are available."""
    if fieldnames is None:
        sys.exit(
            "ERROR: Input table does not contain a header."
        )

    required = {
        ID_COLUMN,
        SCORE_COLUMN,
        COUNT_COLUMN,
        COMPLETENESS_COLUMN,
    }

    missing = sorted(
        required.difference(fieldnames)
    )

    if missing:
        sys.exit(
            "ERROR: Required columns are missing:\n"
            + "\n".join(missing)
        )


def read_scores():
    """
    Read osmotic scores and calculate circle properties.

    Returns one row per MAG.
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

            score = parse_float(
                row[SCORE_COLUMN],
                SCORE_COLUMN,
                mag,
            )

            if score < 0 or score > 1:
                sys.exit(
                    "ERROR: Osmotic strategy score is outside 0-1.\n"
                    f"MAG: {mag}\n"
                    f"Score: {score}"
                )

            capability_count = parse_integer(
                row[COUNT_COLUMN],
                COUNT_COLUMN,
                mag,
            )

            if capability_count < 0 or capability_count > 5:
                sys.exit(
                    "ERROR: Capability count is outside 0-5.\n"
                    f"MAG: {mag}\n"
                    f"Count: {capability_count}"
                )

            expected_score = (
                capability_count / 5.0
            )

            if not math.isclose(
                score,
                expected_score,
                rel_tol=0.0,
                abs_tol=0.00005,
            ):
                sys.exit(
                    "ERROR: Score does not match capability count.\n"
                    f"MAG: {mag}\n"
                    f"Count: {capability_count}\n"
                    f"Expected score: {expected_score:.4f}\n"
                    f"Observed score: {score:.4f}"
                )

            completeness = parse_float(
                row[COMPLETENESS_COLUMN],
                COMPLETENESS_COLUMN,
                mag,
            )

            symbol_size = score_to_size(
                score
            )

            symbol_colour = score_to_colour(
                score
            )

            output_rows.append(
                {
                    "Final_MAG": mag,
                    "Osmotic_capability_count_0_5": (
                        capability_count
                    ),
                    "Osmotic_strategy_score_0_1": (
                        score
                    ),
                    "Circle_relative_size": (
                        symbol_size
                    ),
                    "Circle_colour": (
                        symbol_colour
                    ),
                    "Completeness": (
                        completeness
                    ),
                    "Individual_MAG_completeness_correction": (
                        "Not_applied"
                    ),
                }
            )

    if not output_rows:
        sys.exit(
            "ERROR: No MAG score rows were recovered."
        )

    return output_rows


def write_itol_symbol_file(rows):
    """
    Write the iTOL DATASET_SYMBOL file.

    Required fields per MAG:
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
            "LEGEND_TITLE\tOsmotic adaptation score\n"
        )
        handle.write(
            "LEGEND_SHAPES\t2\t2\t2\n"
        )
        handle.write(
            f"LEGEND_COLORS\t{COLOR_MIN}\t{COLOR_MID}\t{COLOR_MAX}\n"
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
                "Osmotic_strategy_score_0_1"
            ]

            symbol_size = row[
                "Circle_relative_size"
            ]

            symbol_colour = row[
                "Circle_colour"
            ]

            label = (
                "Osmotic adaptation: "
                f"{score:.1f}"
            )

            handle.write(
                f"{mag}\t"
                f"{SYMBOL}\t"
                f"{symbol_size:.4f}\t"
                f"{symbol_colour}\t"
                f"{FILL}\t"
                f"{EXTERNAL_POSITION}\t"
                f"{label}\n"
            )


def write_audit(rows):
    """Write the score-to-symbol conversion audit."""
    columns = [
        "Final_MAG",
        "Osmotic_capability_count_0_5",
        "Osmotic_strategy_score_0_1",
        "Circle_relative_size",
        "Circle_colour",
        "Completeness",
        "Individual_MAG_completeness_correction",
    ]

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
        )

        writer.writeheader()

        for row in rows:
            output_row = dict(
                row
            )

            output_row[
                "Osmotic_strategy_score_0_1"
            ] = (
                f"{row['Osmotic_strategy_score_0_1']:.4f}"
            )

            output_row[
                "Circle_relative_size"
            ] = (
                f"{row['Circle_relative_size']:.4f}"
            )

            output_row[
                "Completeness"
            ] = (
                f"{row['Completeness']:.4f}"
            )

            writer.writerow(
                output_row
            )


def print_summary(rows):
    """Print output locations and score-level symbol properties."""
    score_counts = Counter(
        row["Osmotic_capability_count_0_5"]
        for row in rows
    )

    score_examples = {}

    for row in rows:
        count = row[
            "Osmotic_capability_count_0_5"
        ]

        if count not in score_examples:
            score_examples[count] = {
                "score": row[
                    "Osmotic_strategy_score_0_1"
                ],
                "size": row[
                    "Circle_relative_size"
                ],
                "colour": row[
                    "Circle_colour"
                ],
            }

    print("Osmotic circle-symbol export completed.")
    print()
    print(f"MAGs exported: {len(rows)}")
    print("iTOL symbol: 2 (circle)")
    print("External symbol position: -1")
    print(
        "Individual MAG completeness correction applied: No"
    )
    print()
    print(
        f"iTOL circle dataset: {ITOL_OUTPUT_FILE}"
    )
    print(
        f"Circle conversion audit: {AUDIT_OUTPUT_FILE}"
    )
    print()
    print("Score-to-circle mapping")
    print("-----------------------")

    for count in range(0, 6):
        example = score_examples.get(
            count
        )

        if example is None:
            continue

        print(
            f"{count}/5 "
            f"(score {example['score']:.1f}): "
            f"{score_counts[count]} MAGs, "
            f"relative size {example['size']:.2f}, "
            f"colour {example['colour']}"
        )


def main():
    os.makedirs(
        ITOL_OUTPUT_DIR,
        exist_ok=True,
    )

    os.makedirs(
        TABLE_OUTPUT_DIR,
        exist_ok=True,
    )

    rows = read_scores()

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

