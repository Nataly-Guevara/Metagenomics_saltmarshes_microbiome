#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("/path/to/metagenomics_project")

input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/specific_redox_pathway_prevalence.tsv"

pathways = {
    "Nitrate_reduction": [
        "nitrate reductase",
    ],
    "Nitrite_reduction": [
        "nitrite reductase",
    ],
    "Nitric_oxide_reduction": [
        "nitric oxide reductase",
    ],
    "Nitrous_oxide_reduction": [
        "nitrous oxide reductase",
    ],
    "Fumarate_respiration": [
        "fumarate reductase",
    ],
    "DMSO_TMAO_respiration": [
        "dimethyl sulfoxide reductase",
        "trimethylamine-N-oxide reductase",
    ],
    "Sulfate_reduction": [
        "adenylylsulfate reductase",
        "dissimilatory sulfite reductase",
    ],
    "Sulfur_sulfide_oxidation": [
        "sulfide:quinone oxidoreductase",
        "sulfide dehydrogenase",
        "thiosulfate dehydrogenase",
        "sulfite dehydrogenase",
        "sulfur oxidation protein",
    ],
    "Low_oxygen_terminal_oxidases": [
        "cytochrome bd ubiquinol oxidase",
        "cytochrome cbb3",
    ],
    "Other_aerobic_terminal_oxidases": [
        "cytochrome o ubiquinol oxidase",
        "cytochrome aa3",
    ],
}

group_mags = defaultdict(set)
pathway_mags = defaultdict(lambda: defaultdict(set))
pathway_kos = defaultdict(dict)

with input_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        group = (row["Treatment"], row["Compartment"])
        mag = row["Final_MAG"]
        ko = row["KO"]
        definition = row["KO_definition"]
        definition_lower = definition.lower()

        group_mags[group].add(mag)

        for pathway, keywords in pathways.items():
            if any(keyword.lower() in definition_lower for keyword in keywords):
                pathway_mags[group][pathway].add(mag)
                pathway_kos[pathway][ko] = definition

with output_file.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ["Treatment", "Compartment", "MAGs"] + list(pathways)
    )

    for group in sorted(group_mags):
        n = len(group_mags[group])
        row = [group[0], group[1], n]

        for pathway in pathways:
            prevalence = 100 * len(pathway_mags[group][pathway]) / n
            row.append(f"{prevalence:.1f}")

        writer.writerow(row)

print(f"Output: {output_file}")
print()

with output_file.open() as f:
    print(f.read(), end="")

print("\nKOs included in each pathway:")
for pathway in pathways:
    print(f"\n### {pathway}")
    for ko, definition in sorted(pathway_kos[pathway].items()):
        print(f"{ko}\t{definition}")
