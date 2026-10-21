#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

INPUT = Path(
    "results/phylogeny_bacteria_352/"
    "tree_taxonomy_summary_352.tsv"
)

OUTPUT = Path(
    "results/phylogeny_bacteria_352/"
    "tree_taxonomy_display_352.tsv"
)

DISPLAY_PHYLA = {
    "Actinomycetota",
    "Pseudomonadota",
    "Bacteroidota",
    "Gemmatimonadota",
    "Chloroflexota",
    "Desulfobacterota",
    "Acidobacteriota",
    "Cyanobacteriota",
    "Myxococcota_A",
    "Planctomycetota",
}

DISPLAY_ORDERS = {
    "Solirubrobacterales",
    "Longimicrobiales",
    "Acidimicrobiales",
    "Nitriliruptorales",
    "Spongiisociales",
    "Promineifilales",
    "Multivoradales",
    "Flavobacteriales",
    "Rhodothermales",
    "Palauibacterales",
    "Cyanobacteriales",
    "Xanthomonadales",
    "Desulfobacterales",
    "Rhizobiales",
    "Rhodobacterales",
    "Pseudomonadales",
}

if not INPUT.exists():
    sys.exit(f"ERROR: input file not found:\n{INPUT}")

with INPUT.open("r", newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    required_columns = {"Final_MAG", "Phylum", "Order"}

    if reader.fieldnames is None:
        sys.exit("ERROR: input table has no header.")

    missing_columns = required_columns - set(reader.fieldnames)

    if missing_columns:
        sys.exit(
            "ERROR: missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    output_rows = []

    for row in reader:
        phylum = row["Phylum"].strip()
        order = row["Order"].strip()

        output_rows.append(
            {
                "Final_MAG": row["Final_MAG"].strip(),
                "Phylum": phylum,
                "Order": order,
                "Phylum_display": (
                    phylum if phylum in DISPLAY_PHYLA else "Other"
                ),
                "Order_display": (
                    order if order in DISPLAY_ORDERS else "Other"
                ),
            }
        )

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

fieldnames = [
    "Final_MAG",
    "Phylum",
    "Order",
    "Phylum_display",
    "Order_display",
]

with OUTPUT.open("w", newline="") as outfile:
    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(output_rows)

print("Finished successfully.")
print(f"MAGs processed: {len(output_rows)}")
print(f"Output: {OUTPUT}")
