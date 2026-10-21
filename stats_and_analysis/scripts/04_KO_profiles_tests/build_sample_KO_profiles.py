#!/usr/bin/env python3

import csv
import os
from collections import defaultdict

ABUNDANCE_FILE = (
    "mapping/competitive_MAGs/MAG_abundance/"
    "final_373_MAGs_mean_depth_matrix.tsv"
)

KO_FILE = (
    "results/functional_annotation/kofamscan/"
    "final_373_MAGs_KO_count_matrix.tsv"
)

OUTDIR = "results/functional_annotation/sample_KO_profiles"

RAW_OUT = os.path.join(
    OUTDIR,
    "sample_KO_raw_depth_weighted.tsv"
)

RELATIVE_OUT = os.path.join(
    OUTDIR,
    "sample_KO_relative_MAG_abundance_weighted.tsv"
)

QC_OUT = os.path.join(
    OUTDIR,
    "sample_KO_profiles_QC.txt"
)

os.makedirs(OUTDIR, exist_ok=True)


# ------------------------------------------------------------
# 1. Read MAG abundance matrix
#    Rows: MAGs
#    Columns: samples
# ------------------------------------------------------------

with open(ABUNDANCE_FILE, newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    abundance_header = next(reader)

    if not abundance_header or abundance_header[0] != "MAG_ID":
        raise ValueError(
            "Expected first abundance column to be MAG_ID"
        )

    samples = abundance_header[1:]
    abundance_rows = list(reader)

if len(abundance_rows) != 373:
    raise ValueError(
        f"Expected 373 MAG abundance rows, found {len(abundance_rows)}"
    )

mag_ids = []
abundance_by_mag = {}

for row_number, row in enumerate(abundance_rows, start=2):
    if len(row) != len(abundance_header):
        raise ValueError(
            f"Abundance row {row_number} has {len(row)} columns; "
            f"expected {len(abundance_header)}"
        )

    mag = row[0]

    if mag in abundance_by_mag:
        raise ValueError(f"Duplicate abundance MAG: {mag}")

    values = [float(value) for value in row[1:]]

    if any(value < 0 for value in values):
        raise ValueError(f"Negative abundance detected for MAG: {mag}")

    mag_ids.append(mag)
    abundance_by_mag[mag] = values


# ------------------------------------------------------------
# 2. Read KO copy-number matrix
#    Rows: MAGs
#    Columns: KOs
# ------------------------------------------------------------

with open(KO_FILE, newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    ko_header = next(reader)

    if not ko_header or ko_header[0] != "Final_MAG":
        raise ValueError(
            "Expected first KO matrix column to be Final_MAG"
        )

    ko_ids = ko_header[1:]
    ko_rows = list(reader)

if len(ko_rows) != 373:
    raise ValueError(
        f"Expected 373 KO rows, found {len(ko_rows)}"
    )

if len(set(ko_ids)) != len(ko_ids):
    raise ValueError("Duplicate KO identifiers found in header")

ko_by_mag = {}

for row_number, row in enumerate(ko_rows, start=2):
    if len(row) != len(ko_header):
        raise ValueError(
            f"KO row {row_number} has {len(row)} columns; "
            f"expected {len(ko_header)}"
        )

    mag = row[0]

    if mag in ko_by_mag:
        raise ValueError(f"Duplicate KO MAG: {mag}")

    values = [float(value) for value in row[1:]]

    if any(value < 0 for value in values):
        raise ValueError(f"Negative KO count detected for MAG: {mag}")

    ko_by_mag[mag] = values


# ------------------------------------------------------------
# 3. Validate MAG matching
# ------------------------------------------------------------

abundance_set = set(abundance_by_mag)
ko_set = set(ko_by_mag)

if abundance_set != ko_set:
    only_abundance = sorted(abundance_set - ko_set)
    only_ko = sorted(ko_set - abundance_set)

    raise ValueError(
        "MAG identifiers do not match.\n"
        f"Only in abundance: {len(only_abundance)}\n"
        f"Only in KO matrix: {len(only_ko)}"
    )

number_samples = len(samples)
number_kos = len(ko_ids)

if number_samples != 20:
    raise ValueError(
        f"Expected 20 samples, found {number_samples}"
    )


# ------------------------------------------------------------
# 4. Calculate sample total MAG depth
# ------------------------------------------------------------

sample_total_depth = [0.0] * number_samples

for mag in mag_ids:
    depths = abundance_by_mag[mag]

    for sample_index, depth in enumerate(depths):
        sample_total_depth[sample_index] += depth

zero_total_samples = [
    samples[index]
    for index, total in enumerate(sample_total_depth)
    if total <= 0
]

if zero_total_samples:
    raise ValueError(
        "Samples with zero total MAG depth: "
        + ", ".join(zero_total_samples)
    )


# ------------------------------------------------------------
# 5. Calculate raw and relative KO profiles
#
# Raw:
# sum(MAG mean depth × KO copy number)
#
# Relative:
# sum((MAG depth / total sample MAG depth) × KO copy number)
# ------------------------------------------------------------

raw_profiles = [
    [0.0] * number_kos
    for _ in range(number_samples)
]

relative_profiles = [
    [0.0] * number_kos
    for _ in range(number_samples)
]

for mag_index, mag in enumerate(mag_ids, start=1):
    depths = abundance_by_mag[mag]
    ko_counts = ko_by_mag[mag]

    nonzero_kos = [
        (ko_index, count)
        for ko_index, count in enumerate(ko_counts)
        if count != 0
    ]

    for sample_index, depth in enumerate(depths):
        if depth == 0:
            continue

        relative_depth = depth / sample_total_depth[sample_index]

        for ko_index, count in nonzero_kos:
            raw_profiles[sample_index][ko_index] += depth * count
            relative_profiles[sample_index][ko_index] += (
                relative_depth * count
            )

    if mag_index % 50 == 0 or mag_index == len(mag_ids):
        print(
            f"Processed {mag_index}/{len(mag_ids)} MAGs",
            flush=True
        )


# ------------------------------------------------------------
# 6. Write sample × KO matrices
# ------------------------------------------------------------

def write_profile_matrix(output_file, profiles):
    with open(output_file, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample"] + ko_ids)

        for sample_index, sample in enumerate(samples):
            writer.writerow(
                [sample]
                + [
                    f"{value:.10f}"
                    for value in profiles[sample_index]
                ]
            )


write_profile_matrix(RAW_OUT, raw_profiles)
write_profile_matrix(RELATIVE_OUT, relative_profiles)


# ------------------------------------------------------------
# 7. Write QC report
# ------------------------------------------------------------

raw_nonzero_cells = sum(
    value > 0
    for profile in raw_profiles
    for value in profile
)

relative_nonzero_cells = sum(
    value > 0
    for profile in relative_profiles
    for value in profile
)

kos_present_raw = sum(
    any(raw_profiles[sample_index][ko_index] > 0
        for sample_index in range(number_samples))
    for ko_index in range(number_kos)
)

with open(QC_OUT, "w") as handle:
    handle.write(f"MAGs\t{len(mag_ids)}\n")
    handle.write(f"Samples\t{number_samples}\n")
    handle.write(f"KOs\t{number_kos}\n")
    handle.write(
        f"Output dimensions\t"
        f"{number_samples} samples x {number_kos} KOs\n"
    )
    handle.write(
        f"Raw nonzero cells\t{raw_nonzero_cells}\n"
    )
    handle.write(
        f"Relative nonzero cells\t{relative_nonzero_cells}\n"
    )
    handle.write(
        f"KOs present in at least one sample\t{kos_present_raw}\n"
    )

    handle.write("\nSample\tTotal_MAG_mean_depth\n")

    for sample, total in zip(samples, sample_total_depth):
        handle.write(f"{sample}\t{total:.10f}\n")


print("\nSample-level KO profiles completed successfully.")
print(f"MAGs: {len(mag_ids)}")
print(f"Samples: {number_samples}")
print(f"KOs: {number_kos}")
print(f"Raw profile: {RAW_OUT}")
print(f"Relative profile: {RELATIVE_OUT}")
print(f"QC report: {QC_OUT}")
