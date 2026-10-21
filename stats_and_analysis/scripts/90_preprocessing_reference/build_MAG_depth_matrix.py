#!/usr/bin/env python3

import csv
import glob
import os
import sys
from collections import defaultdict

FAI = "mapping/competitive_MAGs/reference/final_373_MAGs_prefixed.fna.fai"
COVERAGE_GLOB = "mapping/competitive_MAGs/contig_depth/*.coverage.tsv"
OUTDIR = "mapping/competitive_MAGs/MAG_abundance"

WIDE_OUT = os.path.join(OUTDIR, "final_373_MAGs_mean_depth_matrix.tsv")
LONG_OUT = os.path.join(OUTDIR, "final_373_MAGs_mean_depth_long.tsv")
QC_OUT = os.path.join(OUTDIR, "MAG_depth_matrix_QC.txt")


def get_mag_id(contig_id):
    """
    Convert:
    metabat2__sample__bin.15__k127_12345
    into:
    metabat2__sample__bin.15
    """
    if "__k" not in contig_id:
        raise ValueError(
            f"Could not identify MAG prefix in contig: {contig_id}"
        )

    return contig_id.rsplit("__k", 1)[0]


os.makedirs(OUTDIR, exist_ok=True)

# ------------------------------------------------------------
# 1. Read all reference contigs and calculate complete MAG sizes
# ------------------------------------------------------------

contig_to_mag = {}
contig_length = {}
mag_total_length = defaultdict(int)
mag_contig_count = defaultdict(int)

with open(FAI, "r") as handle:
    for line_number, line in enumerate(handle, start=1):
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 2:
            raise ValueError(f"Malformed FAI line {line_number}")

        contig = fields[0]
        length = int(fields[1])
        mag = get_mag_id(contig)

        if contig in contig_to_mag:
            raise ValueError(f"Duplicate contig in FAI: {contig}")

        contig_to_mag[contig] = mag
        contig_length[contig] = length
        mag_total_length[mag] += length
        mag_contig_count[mag] += 1

mags = sorted(mag_total_length)

if len(mags) != 373:
    raise ValueError(
        f"Expected 373 MAGs in reference, but found {len(mags)}"
    )

# ------------------------------------------------------------
# 2. Read each coverage table
# ------------------------------------------------------------

coverage_files = sorted(glob.glob(COVERAGE_GLOB))

if len(coverage_files) != 20:
    raise ValueError(
        f"Expected 20 coverage files, but found {len(coverage_files)}"
    )

# sample -> MAG -> sum(contig length × contig mean depth)
depth_numerator = {}
sample_qc = {}

for coverage_file in coverage_files:

    sample = os.path.basename(coverage_file).replace(".coverage.tsv", "")

    if sample in depth_numerator:
        raise ValueError(f"Duplicate sample name: {sample}")

    numerator = defaultdict(float)
    observed_contigs = set()
    unknown_contigs = []

    with open(coverage_file, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {"#rname", "meandepth"}
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {coverage_file}")

        missing = required_columns.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Missing columns in {coverage_file}: {sorted(missing)}"
            )

        for row in reader:
            contig = row["#rname"]

            if contig not in contig_to_mag:
                unknown_contigs.append(contig)
                continue

            if contig in observed_contigs:
                raise ValueError(
                    f"Duplicate contig row in {coverage_file}: {contig}"
                )

            observed_contigs.add(contig)

            mag = contig_to_mag[contig]
            length = contig_length[contig]
            mean_depth = float(row["meandepth"])

            numerator[mag] += length * mean_depth

    if unknown_contigs:
        raise ValueError(
            f"{len(unknown_contigs)} coverage contigs were absent from "
            f"the reference index in sample {sample}. "
            f"First example: {unknown_contigs[0]}"
        )

    depth_numerator[sample] = numerator
    sample_qc[sample] = {
        "coverage_rows": len(observed_contigs),
        "zero_depth_reference_contigs": (
            len(contig_to_mag) - len(observed_contigs)
        ),
    }

samples = sorted(depth_numerator)

# ------------------------------------------------------------
# 3. Calculate full-length-weighted MAG mean depth
#
# Unmapped contigs have zero depth because the denominator is
# the complete MAG length from the FASTA index.
# ------------------------------------------------------------

mag_depth = {}

for mag in mags:
    mag_depth[mag] = {}

    denominator = mag_total_length[mag]

    if denominator <= 0:
        raise ValueError(f"Invalid total length for MAG: {mag}")

    for sample in samples:
        numerator = depth_numerator[sample].get(mag, 0.0)
        mag_depth[mag][sample] = numerator / denominator

# ------------------------------------------------------------
# 4. Write wide matrix: 373 MAG rows × 20 sample columns
# ------------------------------------------------------------

with open(WIDE_OUT, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["MAG_ID"] + samples)

    for mag in mags:
        writer.writerow(
            [mag] +
            [f"{mag_depth[mag][sample]:.10f}" for sample in samples]
        )

# ------------------------------------------------------------
# 5. Write long-format table
# ------------------------------------------------------------

with open(LONG_OUT, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(
        [
            "MAG_ID",
            "Sample",
            "Mean_depth",
            "MAG_length_bp",
            "MAG_contig_count",
        ]
    )

    for mag in mags:
        for sample in samples:
            writer.writerow(
                [
                    mag,
                    sample,
                    f"{mag_depth[mag][sample]:.10f}",
                    mag_total_length[mag],
                    mag_contig_count[mag],
                ]
            )

# ------------------------------------------------------------
# 6. Write validation summary
# ------------------------------------------------------------

nonzero_cells = sum(
    mag_depth[mag][sample] > 0
    for mag in mags
    for sample in samples
)

total_cells = len(mags) * len(samples)

with open(QC_OUT, "w") as handle:
    handle.write(f"Reference contigs\t{len(contig_to_mag)}\n")
    handle.write(f"MAGs\t{len(mags)}\n")
    handle.write(f"Samples\t{len(samples)}\n")
    handle.write(f"Matrix cells\t{total_cells}\n")
    handle.write(f"Nonzero cells\t{nonzero_cells}\n")
    handle.write(f"Zero cells\t{total_cells - nonzero_cells}\n")
    handle.write("\nSample\tCoverage_rows\tZero_depth_reference_contigs\n")

    for sample in samples:
        handle.write(
            f"{sample}\t"
            f"{sample_qc[sample]['coverage_rows']}\t"
            f"{sample_qc[sample]['zero_depth_reference_contigs']}\n"
        )

print("MAG abundance matrix completed successfully.")
print(f"MAGs: {len(mags)}")
print(f"Samples: {len(samples)}")
print(f"Wide matrix: {WIDE_OUT}")
print(f"Long table: {LONG_OUT}")
print(f"QC report: {QC_OUT}")
