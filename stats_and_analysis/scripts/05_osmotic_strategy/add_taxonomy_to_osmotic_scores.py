#!/usr/bin/env python3

"""
Add GTDB taxonomy to the final osmotic-capability score matrix.

Inputs
------
1. Final osmotic capability matrix:
   osmotic_capability_scores_final.tsv

2. Master gene table:
   final_373_MAGs_gene_master_table.tsv

Outputs
-------
1. New final matrix containing taxonomy:
   osmotic_capability_scores_final_with_taxonomy.tsv

2. Taxonomy merge audit:
   osmotic_taxonomy_merge_audit.tsv

Important
---------
- The existing final score matrix is not overwritten.
- The master gene table contains many rows per MAG, so taxonomy is
  reduced to one consistent MAG-level record before merging.
- Conflicting non-empty taxonomy values for the same MAG are treated
  as an error rather than being silently resolved.
- Only the Python standard library is used.
"""

import csv
import os
import sys
from collections import defaultdict


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = "/path/to/metagenomics_project"

MASTER_GENE_TABLE = os.path.join(
    PROJECT_ROOT,
    "results/functional_annotation/kofamscan/"
    "final_373_MAGs_gene_master_table.tsv",
)

OSMOTIC_SCORE_FILE = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/tables/"
    "osmotic_capability_scores_final.tsv",
)

OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/tables",
)

OUTPUT_MATRIX = os.path.join(
    OUTPUT_DIR,
    "osmotic_capability_scores_final_with_taxonomy.tsv",
)

AUDIT_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_taxonomy_merge_audit.tsv",
)


# ---------------------------------------------------------------------
# Taxonomy columns to recover from the master table
# ---------------------------------------------------------------------

TAXONOMY_COLUMNS = [
    "GTDB_domain",
    "GTDB_classification",
    "classification_method",
    "note",
    "warnings",
    "GTDB_phylum",
    "GTDB_class",
    "GTDB_order",
    "GTDB_family",
    "GTDB_genus",
    "GTDB_species",
]


def clean_value(value):
    """Return a stripped string, converting missing values to empty."""
    if value is None:
        return ""

    return str(value).strip()


def validate_required_columns(fieldnames, required, file_label):
    """Stop if required columns are missing."""
    if fieldnames is None:
        sys.exit(
            f"ERROR: {file_label} does not contain a header."
        )

    missing = sorted(
        set(required).difference(fieldnames)
    )

    if missing:
        sys.exit(
            f"ERROR: {file_label} is missing required columns:\n"
            + "\n".join(missing)
        )


def read_mag_taxonomy():
    """
    Read the gene-level master table and derive one taxonomy record
    per MAG.

    For each MAG and taxonomy field, all non-empty values must agree.
    Conflicting values are reported as errors.
    """
    try:
        handle = open(
            MASTER_GENE_TABLE,
            "r",
            encoding="utf-8-sig",
            newline="",
        )
    except OSError as exc:
        sys.exit(
            "ERROR: Could not open master gene table:\n"
            f"{MASTER_GENE_TABLE}\n"
            f"{exc}"
        )

    # MAG -> field -> set of non-empty values
    observed_values = defaultdict(
        lambda: defaultdict(set)
    )

    gene_row_counts = defaultdict(int)

    with handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        validate_required_columns(
            reader.fieldnames,
            ["Final_MAG"] + TAXONOMY_COLUMNS,
            "Master gene table",
        )

        for row in reader:
            mag = clean_value(
                row.get("Final_MAG")
            )

            if not mag:
                continue

            gene_row_counts[mag] += 1

            for column in TAXONOMY_COLUMNS:
                value = clean_value(
                    row.get(column)
                )

                if value:
                    observed_values[mag][column].add(value)

    if not observed_values:
        sys.exit(
            "ERROR: No MAG taxonomy records were recovered "
            "from the master gene table."
        )

    taxonomy_by_mag = {}
    conflicts = []

    for mag in sorted(observed_values):
        taxonomy_by_mag[mag] = {}

        for column in TAXONOMY_COLUMNS:
            values = observed_values[mag].get(
                column,
                set(),
            )

            if len(values) > 1:
                conflicts.append(
                    {
                        "Final_MAG": mag,
                        "Column": column,
                        "Values": " || ".join(
                            sorted(values)
                        ),
                    }
                )

            if values:
                taxonomy_by_mag[mag][column] = sorted(
                    values
                )[0]
            else:
                taxonomy_by_mag[mag][column] = ""

    if conflicts:
        print(
            "ERROR: Conflicting taxonomy values were found "
            "for one or more MAGs.",
            file=sys.stderr,
        )
        print(
            "The merge was stopped to avoid silently selecting "
            "one conflicting value.",
            file=sys.stderr,
        )
        print(
            file=sys.stderr,
        )

        for conflict in conflicts[:20]:
            print(
                f"{conflict['Final_MAG']}\t"
                f"{conflict['Column']}\t"
                f"{conflict['Values']}",
                file=sys.stderr,
            )

        if len(conflicts) > 20:
            print(
                f"...and {len(conflicts) - 20} additional conflicts.",
                file=sys.stderr,
            )

        sys.exit(1)

    return taxonomy_by_mag, gene_row_counts


def read_osmotic_scores():
    """
    Read the existing final osmotic score matrix.
    """
    try:
        handle = open(
            OSMOTIC_SCORE_FILE,
            "r",
            encoding="utf-8-sig",
            newline="",
        )
    except OSError as exc:
        sys.exit(
            "ERROR: Could not open final osmotic score file:\n"
            f"{OSMOTIC_SCORE_FILE}\n"
            f"{exc}"
        )

    with handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        validate_required_columns(
            reader.fieldnames,
            ["Final_MAG"],
            "Final osmotic score matrix",
        )

        input_columns = list(
            reader.fieldnames
        )

        rows = []
        seen_mags = set()

        for row in reader:
            mag = clean_value(
                row.get("Final_MAG")
            )

            if not mag:
                continue

            if mag in seen_mags:
                sys.exit(
                    "ERROR: Duplicate MAG found in final "
                    "osmotic score matrix:\n"
                    f"{mag}"
                )

            seen_mags.add(mag)
            rows.append(row)

    if not rows:
        sys.exit(
            "ERROR: No MAG rows were found in the final "
            "osmotic score matrix."
        )

    return rows, input_columns


def build_output_columns(input_columns):
    """
    Insert taxonomy after Compartment when possible.

    All original score columns are preserved in their existing order.
    """
    output_columns = []

    inserted_taxonomy = False

    for column in input_columns:
        output_columns.append(column)

        if column == "Compartment":
            for taxonomy_column in TAXONOMY_COLUMNS:
                if taxonomy_column not in output_columns:
                    output_columns.append(
                        taxonomy_column
                    )

            inserted_taxonomy = True

    if not inserted_taxonomy:
        insertion_point = min(
            7,
            len(output_columns),
        )

        output_columns = (
            output_columns[:insertion_point]
            + [
                column
                for column in TAXONOMY_COLUMNS
                if column not in output_columns
            ]
            + output_columns[insertion_point:]
        )

    return output_columns


def merge_taxonomy(score_rows, taxonomy_by_mag, gene_row_counts):
    """
    Add taxonomy to every MAG-level osmotic score row and prepare
    an audit record.
    """
    merged_rows = []
    audit_rows = []

    score_mags = {
        clean_value(row["Final_MAG"])
        for row in score_rows
    }

    taxonomy_mags = set(
        taxonomy_by_mag
    )

    for row in score_rows:
        mag = clean_value(
            row["Final_MAG"]
        )

        taxonomy = taxonomy_by_mag.get(
            mag
        )

        merged_row = dict(row)

        if taxonomy is None:
            for column in TAXONOMY_COLUMNS:
                merged_row[column] = ""

            match_status = "Missing_from_master_gene_table"
            nonempty_taxonomy_fields = 0
            complete_taxonomy = "No"
            gene_rows = 0

        else:
            for column in TAXONOMY_COLUMNS:
                merged_row[column] = taxonomy.get(
                    column,
                    "",
                )

            nonempty_taxonomy_fields = sum(
                1
                for column in TAXONOMY_COLUMNS
                if taxonomy.get(column, "")
            )

            required_rank_columns = [
                "GTDB_domain",
                "GTDB_phylum",
                "GTDB_class",
                "GTDB_order",
                "GTDB_family",
                "GTDB_genus",
                "GTDB_species",
            ]

            complete_taxonomy = (
                "Yes"
                if all(
                    taxonomy.get(column, "")
                    for column in required_rank_columns
                )
                else "No"
            )

            match_status = "Matched"
            gene_rows = gene_row_counts.get(
                mag,
                0,
            )

        merged_rows.append(
            merged_row
        )

        audit_rows.append(
            {
                "Final_MAG": mag,
                "Merge_status": match_status,
                "Master_gene_rows": gene_rows,
                "Nonempty_taxonomy_fields": (
                    nonempty_taxonomy_fields
                ),
                "Complete_domain_to_species": (
                    complete_taxonomy
                ),
                "GTDB_domain": merged_row.get(
                    "GTDB_domain",
                    "",
                ),
                "GTDB_phylum": merged_row.get(
                    "GTDB_phylum",
                    "",
                ),
                "GTDB_class": merged_row.get(
                    "GTDB_class",
                    "",
                ),
                "GTDB_order": merged_row.get(
                    "GTDB_order",
                    "",
                ),
                "GTDB_family": merged_row.get(
                    "GTDB_family",
                    "",
                ),
                "GTDB_genus": merged_row.get(
                    "GTDB_genus",
                    "",
                ),
                "GTDB_species": merged_row.get(
                    "GTDB_species",
                    "",
                ),
            }
        )

    extra_master_mags = sorted(
        taxonomy_mags.difference(score_mags)
    )

    for mag in extra_master_mags:
        taxonomy = taxonomy_by_mag[mag]

        audit_rows.append(
            {
                "Final_MAG": mag,
                "Merge_status": (
                    "Present_in_master_but_missing_from_score_matrix"
                ),
                "Master_gene_rows": gene_row_counts.get(
                    mag,
                    0,
                ),
                "Nonempty_taxonomy_fields": sum(
                    1
                    for column in TAXONOMY_COLUMNS
                    if taxonomy.get(column, "")
                ),
                "Complete_domain_to_species": (
                    "Yes"
                    if all(
                        taxonomy.get(column, "")
                        for column in [
                            "GTDB_domain",
                            "GTDB_phylum",
                            "GTDB_class",
                            "GTDB_order",
                            "GTDB_family",
                            "GTDB_genus",
                            "GTDB_species",
                        ]
                    )
                    else "No"
                ),
                "GTDB_domain": taxonomy.get(
                    "GTDB_domain",
                    "",
                ),
                "GTDB_phylum": taxonomy.get(
                    "GTDB_phylum",
                    "",
                ),
                "GTDB_class": taxonomy.get(
                    "GTDB_class",
                    "",
                ),
                "GTDB_order": taxonomy.get(
                    "GTDB_order",
                    "",
                ),
                "GTDB_family": taxonomy.get(
                    "GTDB_family",
                    "",
                ),
                "GTDB_genus": taxonomy.get(
                    "GTDB_genus",
                    "",
                ),
                "GTDB_species": taxonomy.get(
                    "GTDB_species",
                    "",
                ),
            }
        )

    return merged_rows, audit_rows


def write_merged_matrix(merged_rows, output_columns):
    """
    Write the new taxonomy-enriched osmotic score matrix.
    """
    with open(
        OUTPUT_MATRIX,
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

        for row in merged_rows:
            writer.writerow(row)


def write_audit(audit_rows):
    """
    Write the taxonomy-merge audit.
    """
    audit_columns = [
        "Final_MAG",
        "Merge_status",
        "Master_gene_rows",
        "Nonempty_taxonomy_fields",
        "Complete_domain_to_species",
        "GTDB_domain",
        "GTDB_phylum",
        "GTDB_class",
        "GTDB_order",
        "GTDB_family",
        "GTDB_genus",
        "GTDB_species",
    ]

    with open(
        AUDIT_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=audit_columns,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in audit_rows:
            writer.writerow(row)


def print_summary(merged_rows, audit_rows):
    """
    Print concise merge results.
    """
    matched = sum(
        1
        for row in audit_rows
        if row["Merge_status"] == "Matched"
    )

    missing_from_master = sum(
        1
        for row in audit_rows
        if row["Merge_status"] == (
            "Missing_from_master_gene_table"
        )
    )

    extra_in_master = sum(
        1
        for row in audit_rows
        if row["Merge_status"] == (
            "Present_in_master_but_missing_from_score_matrix"
        )
    )

    complete_taxonomy = sum(
        1
        for row in audit_rows
        if (
            row["Merge_status"] == "Matched"
            and row["Complete_domain_to_species"] == "Yes"
        )
    )

    phylum_present = sum(
        1
        for row in merged_rows
        if clean_value(
            row.get("GTDB_phylum")
        )
    )

    order_present = sum(
        1
        for row in merged_rows
        if clean_value(
            row.get("GTDB_order")
        )
    )

    print("Taxonomy merge completed.")
    print()
    print(
        f"MAGs in osmotic score matrix: {len(merged_rows)}"
    )
    print(
        f"MAGs matched to master taxonomy: {matched}"
    )
    print(
        f"MAGs missing from master taxonomy: {missing_from_master}"
    )
    print(
        f"Master-table MAGs absent from score matrix: {extra_in_master}"
    )
    print(
        "MAGs with complete domain-to-species taxonomy: "
        f"{complete_taxonomy}"
    )
    print(
        f"MAGs with GTDB phylum: {phylum_present}"
    )
    print(
        f"MAGs with GTDB order: {order_present}"
    )
    print()
    print(
        f"Taxonomy-enriched matrix: {OUTPUT_MATRIX}"
    )
    print(
        f"Merge audit: {AUDIT_OUTPUT}"
    )
    print()
    print(
        "The original osmotic_capability_scores_final.tsv "
        "was not modified."
    )


def main():
    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    taxonomy_by_mag, gene_row_counts = (
        read_mag_taxonomy()
    )

    score_rows, input_columns = (
        read_osmotic_scores()
    )

    merged_rows, audit_rows = merge_taxonomy(
        score_rows,
        taxonomy_by_mag,
        gene_row_counts,
    )

    output_columns = build_output_columns(
        input_columns
    )

    write_merged_matrix(
        merged_rows,
        output_columns,
    )

    write_audit(
        audit_rows
    )

    print_summary(
        merged_rows,
        audit_rows,
    )


if __name__ == "__main__":
    main()

