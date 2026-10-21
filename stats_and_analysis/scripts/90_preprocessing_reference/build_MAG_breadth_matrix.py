#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

base = Path("mapping/competitive_MAGs")
abundance = base / "MAG_abundance/final_373_MAGs_mean_depth_matrix.tsv"
coverage_dir = base / "contig_depth"

with abundance.open() as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    mags = [row[0] for row in reader]

samples = header[1:]
mag_set = set(mags)

if len(mags) != 373 or len(mag_set) != 373:
    raise ValueError("Expected 373 unique MAGs")

matrix = {mag: {} for mag in mags}

for sample in samples:

    path = coverage_dir / f"{sample}.coverage.tsv"

    total_length = defaultdict(int)
    covered_length = defaultdict(int)
    unmatched = set()

    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            contig = row["#rname"]
            mag = contig.rsplit("__", 1)[0]

            if mag not in mag_set:
                unmatched.add(contig)
                continue

            length = (
                int(row["endpos"])
                - int(row["startpos"]) + 1
            )

            covered = int(row["covbases"])

            if length <= 0 or not 0 <= covered <= length:
                raise ValueError(f"Invalid coverage: {contig}")

            total_length[mag] += length
            covered_length[mag] += covered

    if unmatched:
        raise ValueError(
            f"{sample}: {len(unmatched)} unmatched contigs; "
            f"examples: {list(unmatched)[:3]}"
        )

    values = []

    for mag in mags:

        if total_length[mag] == 0:
            raise ValueError(
                f"Missing reference contigs for {mag} in {sample}"
            )

        breadth = covered_length[mag] / total_length[mag]
        matrix[mag][sample] = breadth
        values.append(breadth)

    print(
        sample.split("_")[1:4],
        ">=10%:", sum(v >= 0.10 for v in values),
        ">=50%:", sum(v >= 0.50 for v in values),
        flush=True
    )

output = base / "MAG_abundance/final_373_MAGs_breadth_matrix.tsv"

with output.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["MAG_ID"] + samples)

    for mag in mags:
        writer.writerow(
            [mag] + [round(matrix[mag][s], 6) for s in samples]
        )

print("\nSUCCESS")
print("MAGs:", len(mags))
print("Samples:", len(samples))
print("Output:", output)
