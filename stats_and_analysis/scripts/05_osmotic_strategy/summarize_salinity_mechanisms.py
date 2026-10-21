#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("/path/to/metagenomics_project")

input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/salinity_mechanism_prevalence.tsv"

mechanisms = {
    "Trehalose": {
        "K00697", "K01087", "K01236", "K05342",
        "K10236", "K10237", "K10238", "K13057", "K16055"
    },
    "Betaine_proline_transport": {
        "K02000", "K02001", "K02002", "K02168",
        "K03451", "K03762", "K05020",
        "K05845", "K05846", "K05847", "K11928"
    },
    "Ectoine": {
        "K06720"
    },
    "Mechanosensors": {
        "K03282", "K03442", "K05802", "K16053", "K22044"
    },
    "Ion_homeostasis": {
        "K01545", "K01546", "K01547", "K01548",
        "K03498", "K03499", "K03549",
        "K09696", "K09697", "K11105",
        "K15987", "K24160"
    },
}

group_mags = defaultdict(set)
mechanism_mags = defaultdict(lambda: defaultdict(set))

with input_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        group = (row["Treatment"], row["Compartment"])
        mag = row["Final_MAG"]
        ko = row["KO"]

        group_mags[group].add(mag)

        for mechanism, kos in mechanisms.items():
            if ko in kos:
                mechanism_mags[group][mechanism].add(mag)

groups = sorted(group_mags)

with output_file.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ["Treatment", "Compartment", "MAGs"] + list(mechanisms)
    )

    for group in groups:
        n = len(group_mags[group])
        row = [group[0], group[1], n]

        for mechanism in mechanisms:
            prevalence = 100 * len(mechanism_mags[group][mechanism]) / n
            row.append(f"{prevalence:.1f}")

        writer.writerow(row)

print(f"Output: {output_file}")
print()

with output_file.open() as f:
    print(f.read(), end="")
