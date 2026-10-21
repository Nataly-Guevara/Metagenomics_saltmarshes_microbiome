#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("/path/to/metagenomics_project")
input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/curated_redox_mechanism_prevalence.tsv"

mag_kos = defaultdict(set)
mag_group = {}

with input_file.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        mag = row["Final_MAG"]
        mag_kos[mag].add(row["KO"])
        mag_group[mag] = (row["Treatment"], row["Compartment"])

def detected(kos, mechanism):
    if mechanism == "Nitrate_reduction":
        return bool(kos & {"K00370", "K02567"})
    if mechanism == "Nitrite_to_NO":
        return bool(kos & {"K00368", "K15864"})
    if mechanism == "DNRA":
        return "K03385" in kos
    if mechanism == "NO_to_N2O":
        return "K04561" in kos
    if mechanism == "N2O_to_N2":
        return "K00376" in kos
    if mechanism == "Sulfate_reduction":
        return (
            {"K00394", "K00395"} <= kos
            or {"K11180", "K11181"} <= kos
        )
    if mechanism == "Sulfide_oxidation":
        return bool(kos & {"K17218", "K17229", "K17230"})
    if mechanism == "Sulfite_oxidation":
        return {"K21307", "K21308", "K21309"} <= kos
    if mechanism == "Low_O2_respiration":
        return {"K00425", "K00426"} <= kos

mechanisms = [
    "Nitrate_reduction",
    "Nitrite_to_NO",
    "DNRA",
    "NO_to_N2O",
    "N2O_to_N2",
    "Sulfate_reduction",
    "Sulfide_oxidation",
    "Sulfite_oxidation",
    "Low_O2_respiration",
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
