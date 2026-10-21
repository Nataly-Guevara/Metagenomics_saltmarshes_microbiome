#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("/path/to/metagenomics_project")

input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/curated_carbon_methane_mechanism_prevalence.tsv"

mag_kos = defaultdict(set)
mag_group = {}

with input_file.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        mag = row["Final_MAG"]
        mag_kos[mag].add(row["KO"])
        mag_group[mag] = (row["Treatment"], row["Compartment"])

def detected(kos, mechanism):
    if mechanism == "Calvin_cycle":
        return {"K01601", "K01602"} <= kos

    if mechanism == "Wood_Ljungdahl":
        return {"K00198", "K14138"} <= kos

    if mechanism == "Methanol_oxidation":
        return bool(kos & {"K14028", "K14029", "K23995"})

    if mechanism == "Canonical_methanogenesis":
        return "K00399" in kos

    if mechanism == "Canonical_methane_oxidation":
        return "K10944" in kos

mechanisms = [
    "Calvin_cycle",
    "Wood_Ljungdahl",
    "Methanol_oxidation",
    "Canonical_methanogenesis",
    "Canonical_methane_oxidation",
]

group_mags = defaultdict(set)
mechanism_mags = defaultdict(lambda: defaultdict(set))

for mag, kos in mag_kos.items():
    group = mag_group[mag]
    group_mags[group].add(mag)

    for mechanism in mechanisms:
        if detected(kos, mechanism):
            mechanism_mags[group][mechanism].add(mag)

with output_file.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["Treatment", "Compartment", "MAGs"] + mechanisms)

    for group in sorted(group_mags):
        n = len(group_mags[group])
        row = [group[0], group[1], n]

        for mechanism in mechanisms:
            prevalence = 100 * len(mechanism_mags[group][mechanism]) / n
            row.append(f"{prevalence:.1f}")

        writer.writerow(row)

print(f"Output: {output_file}\n")

with output_file.open() as f:
    print(f.read(), end="")
