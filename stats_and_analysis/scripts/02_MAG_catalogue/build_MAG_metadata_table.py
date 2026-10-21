#!/usr/bin/env python3

import csv
from pathlib import Path

base = Path("/path/to/metagenomics_project")

input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_metadata.tsv"

metadata_columns = [
    "Final_MAG",
    "Binner",
    "Original_bin",
    "Completeness",
    "Contamination",
    "Quality_class",
    "GTDB_domain",
    "GTDB_classification",
    "classification_method",
    "note",
    "warnings",
    "Site",
    "Treatment",
    "Compartment",
    "GTDB_phylum",
    "GTDB_class",
    "GTDB_order",
    "GTDB_family",
    "GTDB_genus",
    "GTDB_species",
]

mags = {}

with input_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        mag = row["Final_MAG"]
        if mag not in mags:
            mags[mag] = {column: row[column] for column in metadata_columns}

with output_file.open("w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=metadata_columns,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()

    for mag in sorted(mags):
        writer.writerow(mags[mag])

print(f"MAGs: {len(mags)}")
print(f"Output: {output_file}")
