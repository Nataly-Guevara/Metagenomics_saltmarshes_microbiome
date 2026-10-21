#!/usr/bin/env python3

"""
Create Ricci-style genome-level osmotic capability matrices and
iTOL datasets.

Method
------
For each MAG, five biologically curated capabilities are represented
as binary variables:

    1. Trehalose biosynthesis
    2. Compatible-solute uptake
    3. Ectoine biosynthesis
    4. Mechanosensitive-channel protection
    5. Sodium homeostasis

For MAG i and capability j:

    M_ij = 1 if the accepted biological rule is satisfied
           0 otherwise

The genome-level osmotic strategy score is:

    Osmotic score = sum(M_ij) / 5

This follows the module-based strategy used by Ricci et al.:

- accepted marker genes are grouped into biological capabilities;
- each capability is converted to binary presence/absence;
- capabilities receive equal weight;
- individual MAG scores are not completeness-corrected.

Completeness correction is reserved for later group-level prevalence
summaries and is not applied to the iTOL tree rings.

Inputs
------
osmotic_capability_scores_final_with_taxonomy.tsv

Outputs
-------
1. osmotic_ricci_capability_matrix.tsv

   MAG-level binary capability matrix, raw count, proportional score,
   completeness and selected taxonomy.

2. itol_osmotic_strategy_score_heatmap.txt

   Primary iTOL heatmap ring containing one continuous score per MAG.

3. itol_osmotic_subcapabilities_binary.txt

   Optional iTOL binary dataset displaying the five individual
   osmotic sub-capabilities.

4. osmotic_itol_export_audit.tsv

   Audit confirming score consistency and iTOL identifier uniqueness.

Existing files are not overwritten except files with these exact
new output names when this script is rerun.

Only the Python standard library is used.
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
    "osmotic_capability_scores_final_with_taxonomy.tsv",
)

TABLE_OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/tables",
)

ITOL_OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/itol",
)

RICCI_MATRIX_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "osmotic_ricci_capability_matrix.tsv",
)

ITOL_SCORE_OUTPUT = os.path.join(
    ITOL_OUTPUT_DIR,
    "itol_osmotic_strategy_score_heatmap.txt",
)

ITOL_BINARY_OUTPUT = os.path.join(
    ITOL_OUTPUT_DIR,
    "itol_osmotic_subcapabilities_binary.txt",
)

AUDIT_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "osmotic_itol_export_audit.tsv",
)


# ---------------------------------------------------------------------
# Final Ricci-style capability definitions
# ---------------------------------------------------------------------

CAPABILITIES = [
    {
        "number": 1,
        "source_column": "Osmotic_Trehalose_biosynthesis",
        "output_column": "Trehalose_biosynthesis",
        "itol_label": "Trehalose biosynthesis",
    },
    {
        "number": 2,
        "source_column": "Osmotic_Compatible_solute_uptake",
        "output_column": "Compatible_solute_uptake",
        "itol_label": "Compatible-solute uptake",
    },
    {
        "number": 3,
        "source_column": "Osmotic_Ectoine_biosynthesis",
        "output_column": "Ectoine_biosynthesis",
        "itol_label": "Ectoine biosynthesis",
    },
    {
        "number": 4,
        "source_column": "Osmotic_Mechanosensitive_channel",
        "output_column": "Mechanosensitive_channel_protection",
        "itol_label": "Mechanosensitive channels",
    },
    {
        "number": 5,
        "source_column": "Osmotic_Sodium_homeostasis",
        "output_column": "Sodium_homeostasis",
        "itol_label": "Sodium homeostasis",
    },
]


# ---------------------------------------------------------------------
# Metadata retained in the Ricci-style matrix
# ---------------------------------------------------------------------

METADATA_COLUMNS = [
    "Completeness",
    "Contamination",
    "Quality_class",
    "Site",
    "Treatment",
    "Compartment",
    "GTDB_domain",
    "GTDB_phylum",
    "GTDB_class",
    "GTDB_order",
    "GTDB_family",
    "GTDB_genus",
    "GTDB_species",
]


# ---------------------------------------------------------------------
# Existing score columns used for validation
# ---------------------------------------------------------------------

EXISTING_COUNT_COLUMN = "Osmotic_capability_count_0_5"
EXISTING_SCORE_COLUMN = "Osmotic_capability_score_0_1"


# ---------------------------------------------------------------------
# Basic parsing and validation
# ---------------------------------------------------------------------

def clean_text(value):
    """Return a stripped string, treating missing values as empty."""
    if value is None:
        return ""

    return str(value).strip()


def parse_float(value, column, mag):
    """Parse a finite numeric value."""
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
    """Parse a binary value represented as zero or one."""
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
    """Confirm that all required input columns exist."""
    if fieldnames is None:
        sys.exit(
            "ERROR: Input file does not contain a header."
        )

    required = {
        "Final_MAG",
        EXISTING_COUNT_COLUMN,
        EXISTING_SCORE_COLUMN,
    }

    required.update(
        capability["source_column"]
        for capability in CAPABILITIES
    )

    required.update(
        METADATA_COLUMNS
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
# Read and verify the MAG-level matrix
# ---------------------------------------------------------------------

def read_and_calculate_matrix():
    """
    Read the final taxonomy-enriched matrix.

    Recalculate the Ricci-style score directly from the five binary
    capabilities and verify that it agrees with the previously
    generated score columns.
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
    audit_rows = []
    seen_mags = set()

    with handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        validate_required_columns(
            reader.fieldnames
        )

        for input_row in reader:
            mag = clean_text(
                input_row.get("Final_MAG")
            )

            if not mag:
                continue

            if mag in seen_mags:
                sys.exit(
                    "ERROR: Duplicate Final_MAG identifier found:\n"
                    f"{mag}"
                )

            seen_mags.add(mag)

            capability_values = {}

            for capability in CAPABILITIES:
                source_column = capability["source_column"]
                output_column = capability["output_column"]

                capability_values[output_column] = parse_binary(
                    input_row[source_column],
                    source_column,
                    mag,
                )

            recalculated_count = sum(
                capability_values.values()
            )

            recalculated_score = (
                recalculated_count / len(CAPABILITIES)
            )

            existing_count = parse_integer(
                input_row[EXISTING_COUNT_COLUMN],
                EXISTING_COUNT_COLUMN,
                mag,
            )

            existing_score = parse_float(
                input_row[EXISTING_SCORE_COLUMN],
                EXISTING_SCORE_COLUMN,
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
                    "ERROR: Recalculated Ricci-style score does not "
                    "match the existing final score.\n"
                    f"MAG: {mag}\n"
                    f"Recalculated count: {recalculated_count}\n"
                    f"Existing count: {existing_count}\n"
                    f"Recalculated score: {recalculated_score:.4f}\n"
                    f"Existing score: {existing_score:.4f}"
                )

            output_row = {
                "Final_MAG": mag,
            }

            for column in METADATA_COLUMNS:
                output_row[column] = clean_text(
                    input_row.get(column)
                )

            for capability in CAPABILITIES:
                output_column = capability["output_column"]

                output_row[output_column] = (
                    capability_values[output_column]
                )

            output_row[
                "Osmotic_capability_count_0_5"
            ] = recalculated_count

            output_row[
                "Osmotic_strategy_score_0_1"
            ] = f"{recalculated_score:.4f}"

            output_row[
                "Genome_level_completeness_correction"
            ] = "Not_applied"

            output_rows.append(
                output_row
            )

            audit_rows.append(
                {
                    "Final_MAG": mag,
                    "Unique_tree_identifier": "Yes",
                    "Capability_count_recalculated": (
                        recalculated_count
                    ),
                    "Capability_count_existing": (
                        existing_count
                    ),
                    "Count_matches": (
                        "Yes" if count_matches else "No"
                    ),
                    "Strategy_score_recalculated": (
                        f"{recalculated_score:.4f}"
                    ),
                    "Strategy_score_existing": (
                        f"{existing_score:.4f}"
                    ),
                    "Score_matches": (
                        "Yes" if score_matches else "No"
                    ),
                    "Completeness": clean_text(
                        input_row.get("Completeness")
                    ),
                    "Individual_MAG_correction_applied": "No",
                }
            )

    if not output_rows:
        sys.exit(
            "ERROR: No MAG rows were recovered from the input."
        )

    return output_rows, audit_rows


# ---------------------------------------------------------------------
# Write the Ricci-style matrix
# ---------------------------------------------------------------------

def write_ricci_matrix(rows):
    """Write the final binary capability matrix."""
    output_columns = (
        ["Final_MAG"]
        + METADATA_COLUMNS
        + [
            capability["output_column"]
            for capability in CAPABILITIES
        ]
        + [
            "Osmotic_capability_count_0_5",
            "Osmotic_strategy_score_0_1",
            "Genome_level_completeness_correction",
        ]
    )

    with open(
        RICCI_MATRIX_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_columns,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


# ---------------------------------------------------------------------
# Write the primary iTOL strategy ring
# ---------------------------------------------------------------------

def write_itol_score_dataset(rows):
    """
    Write one continuous heatmap ring showing the osmotic strategy
    score from zero to one.

    Tree identifiers must exactly match the Final_MAG leaf labels
    in the phylogenomic tree.
    """
    with open(
        ITOL_SCORE_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        handle.write("DATASET_HEATMAP\n")
        handle.write("SEPARATOR TAB\n")
        handle.write("DATASET_LABEL\tOsmotic adaptation\n")
        handle.write("COLOR\t#3B6FB6\n")
        handle.write("FIELD_LABELS\tOsmotic adaptation\n")
        handle.write("FIELD_COLORS\t#3B6FB6\n")
        handle.write("COLOR_MIN\t#F2F2F2\n")
        handle.write("COLOR_MAX\t#1F4E79\n")
        handle.write("COLOR_MID\t#8FB7D9\n")
        handle.write("USE_MID_COLOR\t1\n")
        handle.write("MID_VALUE\t0.5\n")
        handle.write("MINIMUM\t0\n")
        handle.write("MAXIMUM\t1\n")
        handle.write("MARGIN\t1\n")
        handle.write("STRIP_WIDTH\t50\n")
        handle.write("SHOW_INTERNAL\t0\n")
        handle.write("SHOW_LABELS\t0\n")
        handle.write("\n")
        handle.write("DATA\n")

        for row in rows:
            mag = row["Final_MAG"]
            score = row["Osmotic_strategy_score_0_1"]

            handle.write(
                f"{mag}\t{score}\n"
            )


# ---------------------------------------------------------------------
# Write optional binary sub-capability dataset
# ---------------------------------------------------------------------

def write_itol_binary_dataset(rows):
    """
    Write an optional five-field iTOL binary dataset.

    This is useful for checking which component capabilities generate
    the overall osmotic score. It is not required as one of the six
    final ecological-strategy rings.
    """
    field_labels = [
        capability["itol_label"]
        for capability in CAPABILITIES
    ]

    field_colors = [
        "#355C7D",
        "#4C78A8",
        "#72A0C1",
        "#9BBFD4",
        "#C3D9E6",
    ]

    field_shapes = [
        "1"
        for _ in CAPABILITIES
    ]

    with open(
        ITOL_BINARY_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        handle.write("DATASET_BINARY\n")
        handle.write("SEPARATOR TAB\n")
        handle.write(
            "DATASET_LABEL\tOsmotic sub-capabilities\n"
        )
        handle.write("COLOR\t#3B6FB6\n")
        handle.write(
            "FIELD_LABELS\t"
            + "\t".join(field_labels)
            + "\n"
        )
        handle.write(
            "FIELD_COLORS\t"
            + "\t".join(field_colors)
            + "\n"
        )
        handle.write(
            "FIELD_SHAPES\t"
            + "\t".join(field_shapes)
            + "\n"
        )
        handle.write("MARGIN\t1\n")
        handle.write("SHOW_LABELS\t0\n")
        handle.write("\n")
        handle.write("DATA\n")

        for row in rows:
            values = [
                str(
                    row[capability["output_column"]]
                )
                for capability in CAPABILITIES
            ]

            handle.write(
                row["Final_MAG"]
                + "\t"
                + "\t".join(values)
                + "\n"
            )


# ---------------------------------------------------------------------
# Write audit output
# ---------------------------------------------------------------------

def write_audit(audit_rows):
    """Write score and identifier validation results."""
    columns = [
        "Final_MAG",
        "Unique_tree_identifier",
        "Capability_count_recalculated",
        "Capability_count_existing",
        "Count_matches",
        "Strategy_score_recalculated",
        "Strategy_score_existing",
        "Score_matches",
        "Completeness",
        "Individual_MAG_correction_applied",
    ]

    with open(
        AUDIT_OUTPUT,
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

        for row in audit_rows:
            writer.writerow(row)


# ---------------------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------------------

def print_summary(rows):
    """Print output paths and strategy-score distributions."""
    score_distribution = Counter(
        int(row["Osmotic_capability_count_0_5"])
        for row in rows
    )

    capability_prevalence = {}

    for capability in CAPABILITIES:
        column = capability["output_column"]

        positive = sum(
            int(row[column])
            for row in rows
        )

        capability_prevalence[
            capability["itol_label"]
        ] = positive

    total_mags = len(rows)

    print("Ricci-style osmotic iTOL export completed.")
    print()
    print(f"MAGs exported: {total_mags}")
    print("Binary capabilities per MAG: 5")
    print(
        "Individual MAG completeness correction applied: No"
    )
    print()
    print(
        f"Ricci capability matrix: {RICCI_MATRIX_OUTPUT}"
    )
    print(
        f"iTOL osmotic strategy ring: {ITOL_SCORE_OUTPUT}"
    )
    print(
        "iTOL sub-capability dataset: "
        f"{ITOL_BINARY_OUTPUT}"
    )
    print(
        f"Export audit: {AUDIT_OUTPUT}"
    )
    print()
    print("Capability prevalence")
    print("---------------------")

    for capability in CAPABILITIES:
        label = capability["itol_label"]
        positive = capability_prevalence[label]

        prevalence = (
            positive / total_mags * 100.0
        )

        print(
            f"{label:<38}"
            f"{positive:>5} MAGs "
            f"({prevalence:>6.2f}%)"
        )

    print()
    print("Osmotic strategy-score distribution")
    print("-----------------------------------")

    for count in range(0, 6):
        mag_count = score_distribution.get(
            count,
            0,
        )

        prevalence = (
            mag_count / total_mags * 100.0
        )

        print(
            f"{count}/5 "
            f"(score {count / 5.0:.1f}): "
            f"{mag_count} MAGs "
            f"({prevalence:.2f}%)"
        )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    os.makedirs(
        TABLE_OUTPUT_DIR,
        exist_ok=True,
    )

    os.makedirs(
        ITOL_OUTPUT_DIR,
        exist_ok=True,
    )

    rows, audit_rows = read_and_calculate_matrix()

    write_ricci_matrix(
        rows
    )

    write_itol_score_dataset(
        rows
    )

    write_itol_binary_dataset(
        rows
    )

    write_audit(
        audit_rows
    )

    print_summary(
        rows
    )


if __name__ == "__main__":
    main()
