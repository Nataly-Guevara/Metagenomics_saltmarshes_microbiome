#!/usr/bin/env python3

"""
Summarise final osmotic-adaptation capabilities and scores by
environmental and taxonomic groups.

Input
-----
osmotic_capability_scores_final_with_taxonomy.tsv

Outputs
-------
1. osmotic_capability_prevalence_by_group.tsv

   Long-format observed and Ricci-style completeness-corrected
   prevalence for each of the five final binary capabilities.

2. osmotic_score_state_prevalence_by_group.tsv

   Long-format prevalence of each exact score state:
       0/5, 1/5, ..., 5/5
       0/4, 1/4, ..., 4/4

   These corrected score-state prevalences are exploratory because
   incomplete genome recovery can move a MAG between score states.

3. osmotic_score_summary_by_group.tsv

   Uncorrected descriptive statistics for:
       - raw 0-5 capability count;
       - proportional 0-1 capability score;
       - channel-excluded 0-4 count;
       - channel-excluded proportional score.

4. osmotic_group_summary_audit.tsv

   Group sizes, completeness ranges, missing-value counts and flags
   for groups where corrected prevalence should be interpreted
   cautiously.

Methodological distinction
--------------------------
Genome-level capability scores are not corrected for completeness.

For a binary capability j in group g:

    observed prevalence (%) =
        positive MAGs / total MAGs * 100

    corrected prevalence (%) =
        observed prevalence / mean completeness (%) * 100

This follows the group-level correction described from the
Ricci et al. supplementary spreadsheet.

The corrected value is not capped at 100 because values above 100
are useful diagnostics showing that the correction has exceeded
the biologically possible prevalence range.

Only the Python standard library is used.
"""

import csv
import math
import os
import statistics
import sys
from collections import defaultdict


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = "/path/to/metagenomics_project"

INPUT_FILE = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/tables/"
    "osmotic_capability_scores_final_with_taxonomy.tsv",
)

OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/tables",
)

CAPABILITY_PREVALENCE_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_capability_prevalence_by_group.tsv",
)

SCORE_STATE_PREVALENCE_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_score_state_prevalence_by_group.tsv",
)

SCORE_SUMMARY_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_score_summary_by_group.tsv",
)

GROUP_AUDIT_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_group_summary_audit.tsv",
)


# ---------------------------------------------------------------------
# Grouping variables
# ---------------------------------------------------------------------

GROUP_COLUMNS = [
    "Site",
    "Treatment",
    "Compartment",
    "Quality_class",
    "GTDB_phylum",
    "GTDB_order",
]


# ---------------------------------------------------------------------
# Final binary capabilities
# ---------------------------------------------------------------------

CAPABILITIES = [
    {
        "column": "Osmotic_Trehalose_biosynthesis",
        "label": "Trehalose biosynthesis",
    },
    {
        "column": "Osmotic_Compatible_solute_uptake",
        "label": "Compatible-solute uptake",
    },
    {
        "column": "Osmotic_Ectoine_biosynthesis",
        "label": "Ectoine biosynthesis",
    },
    {
        "column": "Osmotic_Mechanosensitive_channel",
        "label": "Mechanosensitive-channel protection",
    },
    {
        "column": "Osmotic_Sodium_homeostasis",
        "label": "Sodium homeostasis",
    },
]


# ---------------------------------------------------------------------
# Score definitions
# ---------------------------------------------------------------------

SCORE_VARIABLES = [
    {
        "column": "Osmotic_capability_count_0_5",
        "label": "Primary capability count",
        "scale": "0-5",
        "possible_values": list(range(0, 6)),
        "integer": True,
    },
    {
        "column": "Osmotic_capability_score_0_1",
        "label": "Primary proportional score",
        "scale": "0-1",
        "possible_values": None,
        "integer": False,
    },
    {
        "column": (
            "Osmotic_capability_count_without_channels_0_4"
        ),
        "label": "Capability count without channels",
        "scale": "0-4",
        "possible_values": list(range(0, 5)),
        "integer": True,
    },
    {
        "column": (
            "Osmotic_capability_score_without_channels_0_1"
        ),
        "label": "Proportional score without channels",
        "scale": "0-1",
        "possible_values": None,
        "integer": False,
    },
]


SCORE_STATE_DEFINITIONS = [
    {
        "column": "Osmotic_capability_count_0_5",
        "score_name": "Primary capability count",
        "denominator": 5,
        "states": list(range(0, 6)),
    },
    {
        "column": (
            "Osmotic_capability_count_without_channels_0_4"
        ),
        "score_name": "Capability count without channels",
        "denominator": 4,
        "states": list(range(0, 5)),
    },
]


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def clean_text(value):
    """Convert a value to a stripped string."""
    if value is None:
        return ""

    return str(value).strip()


def validate_columns(fieldnames, required_columns):
    """Stop if one or more required columns are absent."""
    if fieldnames is None:
        sys.exit(
            "ERROR: Input table does not contain a header."
        )

    missing = sorted(
        set(required_columns).difference(fieldnames)
    )

    if missing:
        sys.exit(
            "ERROR: Required columns are missing from the input:\n"
            + "\n".join(missing)
        )


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
            "ERROR: Non-finite numerical value encountered.\n"
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
        "ERROR: Non-binary capability value encountered.\n"
        f"MAG: {mag}\n"
        f"Column: {column}\n"
        f"Value: {value}"
    )


def parse_integer(value, column, mag):
    """Parse an integer-valued score."""
    number = parse_float(
        value,
        column,
        mag,
    )

    rounded = int(round(number))

    if number != rounded:
        sys.exit(
            "ERROR: Non-integer score encountered.\n"
            f"MAG: {mag}\n"
            f"Column: {column}\n"
            f"Value: {value}"
        )

    return rounded


def percentile(values, proportion):
    """
    Calculate a linearly interpolated percentile.

    proportion must range from 0 to 1.
    """
    if not values:
        return float("nan")

    sorted_values = sorted(values)

    if len(sorted_values) == 1:
        return sorted_values[0]

    position = proportion * (
        len(sorted_values) - 1
    )

    lower_index = int(
        math.floor(position)
    )

    upper_index = int(
        math.ceil(position)
    )

    if lower_index == upper_index:
        return sorted_values[lower_index]

    fraction = position - lower_index

    lower_value = sorted_values[lower_index]
    upper_value = sorted_values[upper_index]

    return (
        lower_value
        + fraction * (upper_value - lower_value)
    )


def safe_percentage(numerator, denominator):
    """Calculate a percentage, returning NaN for zero denominator."""
    if denominator == 0:
        return float("nan")

    return (
        numerator / denominator * 100.0
    )


def corrected_prevalence(
    observed_prevalence_percent,
    mean_completeness_percent,
):
    """
    Apply the Ricci-style group-level correction:

        corrected prevalence =
            observed prevalence / mean completeness * 100
    """
    if mean_completeness_percent <= 0:
        return float("nan")

    return (
        observed_prevalence_percent
        / mean_completeness_percent
        * 100.0
    )


def format_number(value, decimals=4):
    """Format a number while preserving missing values."""
    if value is None:
        return ""

    if isinstance(value, float) and math.isnan(value):
        return ""

    return f"{value:.{decimals}f}"


def group_size_flag(group_size):
    """Flag groups with potentially unstable estimates."""
    if group_size < 5:
        return "Very_small_N_lt5"

    if group_size < 10:
        return "Small_N_lt10"

    return "Adequate_N_ge10"


def correction_flag(corrected_value):
    """Flag corrected prevalence above the biological maximum."""
    if isinstance(corrected_value, float):
        if math.isnan(corrected_value):
            return "Correction_not_calculated"

    if corrected_value > 100:
        return "Corrected_prevalence_above_100"

    return "Corrected_prevalence_0_to_100"


# ---------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------

def read_input():
    """
    Read and validate the taxonomy-enriched MAG score matrix.
    """
    required_columns = [
        "Final_MAG",
        "Completeness",
    ]

    required_columns.extend(
        GROUP_COLUMNS
    )

    required_columns.extend(
        capability["column"]
        for capability in CAPABILITIES
    )

    required_columns.extend(
        score["column"]
        for score in SCORE_VARIABLES
    )

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

    rows = []
    seen_mags = set()

    with handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        validate_columns(
            reader.fieldnames,
            required_columns,
        )

        for input_row in reader:
            mag = clean_text(
                input_row.get("Final_MAG")
            )

            if not mag:
                continue

            if mag in seen_mags:
                sys.exit(
                    "ERROR: Duplicate MAG in input matrix:\n"
                    f"{mag}"
                )

            seen_mags.add(mag)

            completeness = parse_float(
                input_row["Completeness"],
                "Completeness",
                mag,
            )

            if completeness <= 0 or completeness > 100:
                sys.exit(
                    "ERROR: Completeness must be greater than zero "
                    "and no greater than 100.\n"
                    f"MAG: {mag}\n"
                    f"Completeness: {completeness}"
                )

            parsed_row = {
                "Final_MAG": mag,
                "Completeness": completeness,
            }

            for group_column in GROUP_COLUMNS:
                group_value = clean_text(
                    input_row.get(group_column)
                )

                if not group_value:
                    group_value = "Unclassified_or_missing"

                parsed_row[group_column] = group_value

            for capability in CAPABILITIES:
                column = capability["column"]

                parsed_row[column] = parse_binary(
                    input_row[column],
                    column,
                    mag,
                )

            for score in SCORE_VARIABLES:
                column = score["column"]

                if score["integer"]:
                    parsed_row[column] = parse_integer(
                        input_row[column],
                        column,
                        mag,
                    )
                else:
                    parsed_row[column] = parse_float(
                        input_row[column],
                        column,
                        mag,
                    )

            rows.append(parsed_row)

    if not rows:
        sys.exit(
            "ERROR: No MAG rows were recovered from the input."
        )

    return rows


# ---------------------------------------------------------------------
# Group construction
# ---------------------------------------------------------------------

def build_groups(rows):
    """
    Build:

        group_column
            -> group_value
                -> list of MAG rows
    """
    groups = {
        group_column: defaultdict(list)
        for group_column in GROUP_COLUMNS
    }

    for row in rows:
        for group_column in GROUP_COLUMNS:
            group_value = row[group_column]

            groups[group_column][group_value].append(
                row
            )

    return groups


# ---------------------------------------------------------------------
# Capability prevalence
# ---------------------------------------------------------------------

def calculate_capability_prevalence(groups):
    """
    Calculate observed and corrected prevalence for each capability
    in every group.
    """
    output_rows = []

    for group_column in GROUP_COLUMNS:
        group_values = groups[group_column]

        for group_value in sorted(group_values):
            group_rows = group_values[group_value]

            total_mags = len(group_rows)

            completeness_values = [
                row["Completeness"]
                for row in group_rows
            ]

            mean_completeness = statistics.mean(
                completeness_values
            )

            for capability in CAPABILITIES:
                capability_column = capability["column"]

                positive_mags = sum(
                    row[capability_column]
                    for row in group_rows
                )

                observed = safe_percentage(
                    positive_mags,
                    total_mags,
                )

                corrected = corrected_prevalence(
                    observed,
                    mean_completeness,
                )

                output_rows.append(
                    {
                        "Grouping_variable": group_column,
                        "Group": group_value,
                        "Group_MAGs": total_mags,
                        "Group_size_flag": group_size_flag(
                            total_mags
                        ),
                        "Capability": capability["label"],
                        "Capability_column": capability_column,
                        "Positive_MAGs": positive_mags,
                        "Negative_MAGs": (
                            total_mags - positive_mags
                        ),
                        "Observed_prevalence_percent": (
                            observed
                        ),
                        "Mean_completeness_percent": (
                            mean_completeness
                        ),
                        "Completeness_fraction": (
                            mean_completeness / 100.0
                        ),
                        "Corrected_prevalence_percent": (
                            corrected
                        ),
                        "Correction_flag": correction_flag(
                            corrected
                        ),
                        "Correction_method": (
                            "Observed_prevalence_percent / "
                            "Mean_completeness_percent * 100"
                        ),
                    }
                )

    return output_rows


def write_capability_prevalence(output_rows):
    """Write long-format capability prevalence."""
    columns = [
        "Grouping_variable",
        "Group",
        "Group_MAGs",
        "Group_size_flag",
        "Capability",
        "Capability_column",
        "Positive_MAGs",
        "Negative_MAGs",
        "Observed_prevalence_percent",
        "Mean_completeness_percent",
        "Completeness_fraction",
        "Corrected_prevalence_percent",
        "Correction_flag",
        "Correction_method",
    ]

    with open(
        CAPABILITY_PREVALENCE_OUTPUT,
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

        for row in output_rows:
            formatted_row = dict(row)

            for column in [
                "Observed_prevalence_percent",
                "Mean_completeness_percent",
                "Corrected_prevalence_percent",
            ]:
                formatted_row[column] = format_number(
                    row[column],
                    decimals=4,
                )

            formatted_row[
                "Completeness_fraction"
            ] = format_number(
                row["Completeness_fraction"],
                decimals=6,
            )

            writer.writerow(formatted_row)


# ---------------------------------------------------------------------
# Score-state prevalence
# ---------------------------------------------------------------------

def calculate_score_state_prevalence(groups):
    """
    Calculate prevalence of exact count states such as 0/5 or 3/4.

    Completeness-corrected score-state prevalences are labelled
    exploratory because missing capabilities may shift MAGs between
    states rather than simply removing a binary marker.
    """
    output_rows = []

    for group_column in GROUP_COLUMNS:
        group_values = groups[group_column]

        for group_value in sorted(group_values):
            group_rows = group_values[group_value]

            total_mags = len(group_rows)

            mean_completeness = statistics.mean(
                row["Completeness"]
                for row in group_rows
            )

            for definition in SCORE_STATE_DEFINITIONS:
                score_column = definition["column"]
                score_name = definition["score_name"]
                denominator = definition["denominator"]

                for state in definition["states"]:
                    state_mags = sum(
                        1
                        for row in group_rows
                        if row[score_column] == state
                    )

                    observed = safe_percentage(
                        state_mags,
                        total_mags,
                    )

                    corrected = corrected_prevalence(
                        observed,
                        mean_completeness,
                    )

                    output_rows.append(
                        {
                            "Grouping_variable": group_column,
                            "Group": group_value,
                            "Group_MAGs": total_mags,
                            "Group_size_flag": group_size_flag(
                                total_mags
                            ),
                            "Score_variable": score_name,
                            "Score_column": score_column,
                            "Score_state_count": state,
                            "Score_denominator": denominator,
                            "Score_state_label": (
                                f"{state}/{denominator}"
                            ),
                            "MAGs_in_score_state": state_mags,
                            "Observed_prevalence_percent": (
                                observed
                            ),
                            "Mean_completeness_percent": (
                                mean_completeness
                            ),
                            "Corrected_prevalence_percent": (
                                corrected
                            ),
                            "Correction_flag": correction_flag(
                                corrected
                            ),
                            "Interpretation": (
                                "Exploratory_corrected_score_state_"
                                "prevalence"
                            ),
                        }
                    )

    return output_rows


def write_score_state_prevalence(output_rows):
    """Write long-format exact score-state prevalence."""
    columns = [
        "Grouping_variable",
        "Group",
        "Group_MAGs",
        "Group_size_flag",
        "Score_variable",
        "Score_column",
        "Score_state_count",
        "Score_denominator",
        "Score_state_label",
        "MAGs_in_score_state",
        "Observed_prevalence_percent",
        "Mean_completeness_percent",
        "Corrected_prevalence_percent",
        "Correction_flag",
        "Interpretation",
    ]

    with open(
        SCORE_STATE_PREVALENCE_OUTPUT,
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

        for row in output_rows:
            formatted_row = dict(row)

            for column in [
                "Observed_prevalence_percent",
                "Mean_completeness_percent",
                "Corrected_prevalence_percent",
            ]:
                formatted_row[column] = format_number(
                    row[column],
                    decimals=4,
                )

            writer.writerow(formatted_row)


# ---------------------------------------------------------------------
# Score summaries
# ---------------------------------------------------------------------

def calculate_score_summaries(groups):
    """
    Calculate uncorrected descriptive summaries of genome-level
    scores within each group.
    """
    output_rows = []

    for group_column in GROUP_COLUMNS:
        group_values = groups[group_column]

        for group_value in sorted(group_values):
            group_rows = group_values[group_value]

            total_mags = len(group_rows)

            completeness_values = [
                row["Completeness"]
                for row in group_rows
            ]

            mean_completeness = statistics.mean(
                completeness_values
            )

            for score_definition in SCORE_VARIABLES:
                score_column = score_definition["column"]

                values = [
                    float(row[score_column])
                    for row in group_rows
                ]

                mean_value = statistics.mean(
                    values
                )

                median_value = statistics.median(
                    values
                )

                standard_deviation = (
                    statistics.stdev(values)
                    if len(values) > 1
                    else 0.0
                )

                standard_error = (
                    standard_deviation
                    / math.sqrt(len(values))
                    if values
                    else float("nan")
                )

                output_rows.append(
                    {
                        "Grouping_variable": group_column,
                        "Group": group_value,
                        "Group_MAGs": total_mags,
                        "Group_size_flag": group_size_flag(
                            total_mags
                        ),
                        "Score_variable": (
                            score_definition["label"]
                        ),
                        "Score_column": score_column,
                        "Score_scale": (
                            score_definition["scale"]
                        ),
                        "Mean": mean_value,
                        "Median": median_value,
                        "Standard_deviation": (
                            standard_deviation
                        ),
                        "Standard_error": standard_error,
                        "Minimum": min(values),
                        "First_quartile": percentile(
                            values,
                            0.25,
                        ),
                        "Third_quartile": percentile(
                            values,
                            0.75,
                        ),
                        "Maximum": max(values),
                        "Mean_completeness_percent": (
                            mean_completeness
                        ),
                        "Completeness_correction_applied": (
                            "No"
                        ),
                        "Reason_not_corrected": (
                            "Genome-level composite scores are "
                            "reported directly; completeness "
                            "correction is applied to group-level "
                            "binary prevalence."
                        ),
                    }
                )

    return output_rows


def write_score_summaries(output_rows):
    """Write long-format descriptive score summaries."""
    columns = [
        "Grouping_variable",
        "Group",
        "Group_MAGs",
        "Group_size_flag",
        "Score_variable",
        "Score_column",
        "Score_scale",
        "Mean",
        "Median",
        "Standard_deviation",
        "Standard_error",
        "Minimum",
        "First_quartile",
        "Third_quartile",
        "Maximum",
        "Mean_completeness_percent",
        "Completeness_correction_applied",
        "Reason_not_corrected",
    ]

    with open(
        SCORE_SUMMARY_OUTPUT,
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

        for row in output_rows:
            formatted_row = dict(row)

            for column in [
                "Mean",
                "Median",
                "Standard_deviation",
                "Standard_error",
                "Minimum",
                "First_quartile",
                "Third_quartile",
                "Maximum",
                "Mean_completeness_percent",
            ]:
                formatted_row[column] = format_number(
                    row[column],
                    decimals=6,
                )

            writer.writerow(formatted_row)


# ---------------------------------------------------------------------
# Group audit
# ---------------------------------------------------------------------

def calculate_group_audit(groups):
    """
    Summarise group sizes and completeness characteristics.
    """
    output_rows = []

    for group_column in GROUP_COLUMNS:
        group_values = groups[group_column]

        for group_value in sorted(group_values):
            group_rows = group_values[group_value]

            completeness_values = [
                row["Completeness"]
                for row in group_rows
            ]

            total_mags = len(group_rows)

            output_rows.append(
                {
                    "Grouping_variable": group_column,
                    "Group": group_value,
                    "Group_MAGs": total_mags,
                    "Group_size_flag": group_size_flag(
                        total_mags
                    ),
                    "Mean_completeness_percent": (
                        statistics.mean(
                            completeness_values
                        )
                    ),
                    "Median_completeness_percent": (
                        statistics.median(
                            completeness_values
                        )
                    ),
                    "Minimum_completeness_percent": min(
                        completeness_values
                    ),
                    "Maximum_completeness_percent": max(
                        completeness_values
                    ),
                    "Completeness_standard_deviation": (
                        statistics.stdev(
                            completeness_values
                        )
                        if total_mags > 1
                        else 0.0
                    ),
                    "MAGs_below_50_percent_complete": sum(
                        1
                        for value in completeness_values
                        if value < 50
                    ),
                    "MAGs_50_to_below_70_percent_complete": sum(
                        1
                        for value in completeness_values
                        if 50 <= value < 70
                    ),
                    "MAGs_70_to_below_90_percent_complete": sum(
                        1
                        for value in completeness_values
                        if 70 <= value < 90
                    ),
                    "MAGs_at_least_90_percent_complete": sum(
                        1
                        for value in completeness_values
                        if value >= 90
                    ),
                }
            )

    return output_rows


def write_group_audit(output_rows):
    """Write group-size and completeness audit."""
    columns = [
        "Grouping_variable",
        "Group",
        "Group_MAGs",
        "Group_size_flag",
        "Mean_completeness_percent",
        "Median_completeness_percent",
        "Minimum_completeness_percent",
        "Maximum_completeness_percent",
        "Completeness_standard_deviation",
        "MAGs_below_50_percent_complete",
        "MAGs_50_to_below_70_percent_complete",
        "MAGs_70_to_below_90_percent_complete",
        "MAGs_at_least_90_percent_complete",
    ]

    with open(
        GROUP_AUDIT_OUTPUT,
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

        for row in output_rows:
            formatted_row = dict(row)

            for column in [
                "Mean_completeness_percent",
                "Median_completeness_percent",
                "Minimum_completeness_percent",
                "Maximum_completeness_percent",
                "Completeness_standard_deviation",
            ]:
                formatted_row[column] = format_number(
                    row[column],
                    decimals=4,
                )

            writer.writerow(formatted_row)


# ---------------------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------------------

def print_summary(
    rows,
    groups,
    capability_rows,
    score_state_rows,
):
    """
    Print output paths and key diagnostics.
    """
    number_of_groups = sum(
        len(groups[group_column])
        for group_column in GROUP_COLUMNS
    )

    very_small_groups = sum(
        1
        for group_column in GROUP_COLUMNS
        for group_rows in groups[group_column].values()
        if len(group_rows) < 5
    )

    corrected_above_100 = sum(
        1
        for row in capability_rows
        if row["Corrected_prevalence_percent"] > 100
    )

    corrected_score_states_above_100 = sum(
        1
        for row in score_state_rows
        if row["Corrected_prevalence_percent"] > 100
    )

    print("Osmotic group summaries completed.")
    print()
    print(f"MAGs assessed: {len(rows)}")
    print(
        f"Grouping variables assessed: {len(GROUP_COLUMNS)}"
    )
    print(
        f"Unique groups across all variables: {number_of_groups}"
    )
    print(
        f"Groups containing fewer than 5 MAGs: {very_small_groups}"
    )
    print()
    print(
        "Capability prevalence output: "
        f"{CAPABILITY_PREVALENCE_OUTPUT}"
    )
    print(
        "Score-state prevalence output: "
        f"{SCORE_STATE_PREVALENCE_OUTPUT}"
    )
    print(
        "Score summary output: "
        f"{SCORE_SUMMARY_OUTPUT}"
    )
    print(
        "Group audit output: "
        f"{GROUP_AUDIT_OUTPUT}"
    )
    print()
    print(
        "Capability rows with corrected prevalence above 100%: "
        f"{corrected_above_100}"
    )
    print(
        "Exploratory score-state rows with corrected prevalence "
        f"above 100%: {corrected_score_states_above_100}"
    )
    print()
    print(
        "Genome-level score means were not completeness-corrected."
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    rows = read_input()

    groups = build_groups(
        rows
    )

    capability_rows = calculate_capability_prevalence(
        groups
    )

    score_state_rows = calculate_score_state_prevalence(
        groups
    )

    score_summary_rows = calculate_score_summaries(
        groups
    )

    group_audit_rows = calculate_group_audit(
        groups
    )

    write_capability_prevalence(
        capability_rows
    )

    write_score_state_prevalence(
        score_state_rows
    )

    write_score_summaries(
        score_summary_rows
    )

    write_group_audit(
        group_audit_rows
    )

    print_summary(
        rows,
        groups,
        capability_rows,
        score_state_rows,
    )


if __name__ == "__main__":
    main()
