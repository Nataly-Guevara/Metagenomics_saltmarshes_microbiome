#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("/path/to/metagenomics_project")

input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/redox_mechanism_prevalence.tsv"

mechanism_keywords = {
    "Oxidative_stress": [
        "catalase",
        "superoxide dismutase",
        "peroxiredoxin",
        "glutathione peroxidase",
        "thioredoxin",
        "glutaredoxin",
        "methionine sulfoxide reductase",
    ],
    "Aerobic_respiration": [
        "cytochrome c oxidase",
        "cytochrome o ubiquinol oxidase",
        "cytochrome bd ubiquinol oxidase",
    ],
    "Anaerobic_respiration": [
        "nitrate reductase",
        "nitrite reductase",
        "nitric oxide reductase",
        "nitrous oxide reductase",
        "fumarate reductase",
        "dimethyl sulfoxide reductase",
        "trimethylamine-N-oxide reductase",
    ],
    "Electron_transfer": [
        "ferredoxin",
        "flavodoxin",
        "electron-transferring flavoprotein",
        "cytochrome c-type",
        "quinone reductase",
    ],
    "Redox_maintenance": [
        "thioredoxin reductase",
        "glutathione reductase",
        "NADPH:quinone reductase",
        "disulfide reductase",
    ],
}

group_mags = defaultdict(set)
mechanism_mags = defaultdict(lambda: defaultdict(set))
mechanism_kos = defaultdict(set)

with input_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        group = (row["Treatment"], row["Compartment"])
        mag = row["Final_MAG"]
        ko = row["KO"]
        definition = row["KO_definition"].lower()

        group_mags[group].add(mag)

        for mechanism, keywords in mechanism_keywords.items():
            if any(keyword in definition for keyword in keywords):
                mechanism_mags[group][mechanism].add(mag)
                mechanism_kos[mechanism].add(ko)

groups = sorted(group_mags)

with output_file.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ["Treatment", "Compartment", "MAGs"] + list(mechanism_keywords)
    )

    for group in groups:
        n = len(group_mags[group])
        row = [group[0], group[1], n]

        for mechanism in mechanism_keywords:
            prevalence = 100 * len(mechanism_mags[group][mechanism]) / n
            row.append(f"{prevalence:.1f}")

        writer.writerow(row)

print(f"Output: {output_file}")
print()

with output_file.open() as f:
    print(f.read(), end="")

print("\nKOs included per mechanism:")
for mechanism in mechanism_keywords:
    print(f"{mechanism}: {len(mechanism_kos[mechanism])}")
