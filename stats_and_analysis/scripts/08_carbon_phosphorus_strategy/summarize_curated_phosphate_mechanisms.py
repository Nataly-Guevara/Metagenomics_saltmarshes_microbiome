#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("/path/to/metagenomics_project")

input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/curated_phosphate_mechanism_prevalence.tsv"

mag_kos = defaultdict(set)
mag_group = {}

with input_file.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        mag = row["Final_MAG"]
        mag_kos[mag].add(row["KO"])
        mag_group[mag] = (row["Treatment"], row["Compartment"])

def detected(kos, mechanism):
    if mechanism == "High_affinity_P_uptake":
        return {"K02040", "K02037", "K02038", "K02036"} <= kos

    if mechanism == "Alkaline_phosphatase":
        return bool(kos & {"K01077", "K01113"})

    if mechanism == "Polyphosphate_storage":
        return bool(kos & {"K00937", "K22468"})

    if mechanism == "Polyphosphate_degradation":
        return "K01524" in kos

    if mechanism == "Phosphonate_uptake":
        return {"K02041", "K02044", "K02042"} <= kos

mechanisms = [
    "High_affinity_P_uptake",
    "Alkaline_phosphatase",
    "Polyphosphate_storage",
    "Polyphosphate_degradation",
    "Phosphonate_uptake",
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
