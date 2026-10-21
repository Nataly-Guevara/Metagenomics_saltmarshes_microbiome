#!/usr/bin/env python3

import csv
import sys
from pathlib import Path

INPUT = Path(
    "results/phylogeny_bacteria_352/"
    "tree_taxonomy_display_352.tsv"
)

OUTPUT_PHYLA = Path(
    "results/phylogeny_bacteria_352/"
    "itol_phylum_ring_352.txt"
)

OUTPUT_ORDERS = Path(
    "results/phylogeny_bacteria_352/"
    "itol_order_ring_352.txt"
)

PHYLUM_COLORS = {
    "Actinomycetota": "#1f77b4",
    "Pseudomonadota": "#ff7f0e",
    "Bacteroidota": "#2ca02c",
    "Gemmatimonadota": "#d62728",
    "Chloroflexota": "#9467bd",
    "Desulfobacterota": "#8c564b",
    "Acidobacteriota": "#e377c2",
    "Cyanobacteriota": "#17becf",
    "Myxococcota_A": "#bcbd22",
    "Planctomycetota": "#7f7f7f",
    "Other": "#d9d9d9",
}

ORDER_COLORS = {
    "Solirubrobacterales": "#1f77b4",
    "Longimicrobiales": "#ff7f0e",
    "Acidimicrobiales": "#2ca02c",
    "Nitriliruptorales": "#d62728",
    "Spongiisociales": "#9467bd",
    "Promineifilales": "#8c564b",
    "Multivoradales": "#e377c2",
    "Flavobacteriales": "#7f7f7f",
    "Rhodothermales": "#bcbd22",
    "Palauibacterales": "#17becf",
    "Cyanobacteriales": "#393b79",
    "Xanthomonadales": "#637939",
    "Desulfobacterales": "#8c6d31",
    "Rhizobiales": "#843c39",
    "Rhodobacterales": "#7b4173",
    "Pseudomonadales": "#3182bd",
    "Other": "#d9d9d9",
}


def write_colorstrip(
    output_path,
    dataset_label,
    color_map,
    rows,
    display_column,
):
    with output_path.open("w") as outfile:
        outfile.write("DATASET_COLORSTRIP\n")
        outfile.write("SEPARATOR TAB\n")
        outfile.write(f"DATASET_LABEL\t{dataset_label}\n")
        outfile.write("COLOR\t#000000\n")
        outfile.write("STRIP_WIDTH\t25\n")
        outfile.write("MARGIN\t2\n")
        outfile.write("BORDER_WIDTH\t0\n")
        outfile.write("SHOW_INTERNAL\t0\n")
        outfile.write("DATA\n")

        for row in rows:
            mag = row["Final_MAG"].strip()
            category = row[display_column].strip()

            if category not in color_map:
                sys.exit(
                    f"ERROR: no color defined for {category} "
                    f"in {display_column}"
                )

            outfile.write(
                f"{mag}\t{color_map[category]}\t{category}\n"
            )


if not INPUT.exists():
    sys.exit(f"ERROR: input file not found:\n{INPUT}")

with INPUT.open("r", newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    required_columns = {
        "Final_MAG",
        "Phylum_display",
        "Order_display",
    }

    if reader.fieldnames is None:
        sys.exit("ERROR: input table has no header.")

    missing_columns = required_columns - set(reader.fieldnames)

    if missing_columns:
        sys.exit(
            "ERROR: missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    rows = list(reader)

write_colorstrip(
    OUTPUT_PHYLA,
    "Phylum",
    PHYLUM_COLORS,
    rows,
    "Phylum_display",
)

write_colorstrip(
    OUTPUT_ORDERS,
    "Order",
    ORDER_COLORS,
    rows,
    "Order_display",
)

print("Finished successfully.")
print(f"MAGs processed: {len(rows)}")
print(f"Phylum ring: {OUTPUT_PHYLA}")
print(f"Order ring: {OUTPUT_ORDERS}")
