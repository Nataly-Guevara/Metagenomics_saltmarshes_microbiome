#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

base = Path("/path/to/metagenomics_project")

input_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_gene_master_table.tsv"
output_file = base / "results/functional_annotation/kofamscan/final_373_MAGs_KO_count_matrix.tsv"

counts = defaultdict(lambda: defaultdict(int))
mags = set()
kos = set()

with input_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        mag = row["Final_MAG"]
        ko = row["KO"]

        counts[mag][ko] += 1
        mags.add(mag)
        kos.add(ko)

mags = sorted(mags)
kos = sorted(kos)

with output_file.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["Final_MAG"] + kos)

    for mag in mags:
        writer.writerow([mag] + [counts[mag].get(ko, 0) for ko in kos])

print(f"MAGs: {len(mags)}")
print(f"KOs: {len(kos)}")
print(f"Output: {output_file}")
