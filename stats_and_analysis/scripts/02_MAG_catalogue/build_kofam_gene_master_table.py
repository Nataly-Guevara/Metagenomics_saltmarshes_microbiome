#!/usr/bin/env python3

import csv
from pathlib import Path

base = Path("/path/to/metagenomics_project")

metadata_file = base / "results/final_summaries/combined_binners/final_373_MAGs_master_table_with_sample_taxonomy.tsv"
kofam_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_kofam_threshold_hits_complete.tsv"
output_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"

# Load MAG metadata
with metadata_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    metadata = {row["Final_MAG"]: row for row in reader}
    metadata_columns = reader.fieldnames

# Longest first prevents accidental partial-prefix matches
mag_ids = sorted(metadata, key=len, reverse=True)

matched = 0
unmatched = 0

with kofam_file.open() as inp, output_file.open("w", newline="") as out:
    reader = csv.DictReader(inp, delimiter="\t")

    fieldnames = (
        ["gene_id", "Final_MAG"]
        + metadata_columns[1:]
        + ["KO", "threshold", "score", "E_value", "KO_definition"]
    )

    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()

    for row in reader:
        gene_id = row["gene_id"]
        mag = next((m for m in mag_ids if gene_id.startswith(m + "__")), None)

        if mag is None:
            unmatched += 1
            continue

        meta = metadata[mag]

        output_row = {
            "gene_id": gene_id,
            "Final_MAG": mag,
            **{col: meta[col] for col in metadata_columns[1:]},
            "KO": row["KO"],
            "threshold": row["threshold"],
            "score": row["score"],
            "E_value": row["E_value"],
            "KO_definition": row["KO_definition"],
        }

        writer.writerow(output_row)
        matched += 1

print(f"Matched KO assignments: {matched}")
print(f"Unmatched KO assignments: {unmatched}")
print(f"Output: {output_file}")
