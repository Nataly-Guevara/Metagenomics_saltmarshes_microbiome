#!/usr/bin/env python3

import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("/path/to/metagenomics_project")

INPUT_FILE = (
    BASE
    / "results/functional_annotation/kofamscan/"
      "final_373_MAGs_gene_master_table.tsv"
)

OUTPUT_DIR = (
    BASE
    / "results/ecological_strategies/osmotic_adaptation/tables"
)

MARKER_OUTPUT = OUTPUT_DIR / "osmotic_marker_presence.tsv"
SCORE_OUTPUT = OUTPUT_DIR / "osmotic_strategy_scores.tsv"
SUMMARY_OUTPUT = OUTPUT_DIR / "osmotic_summary.tsv"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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

metadata_columns = [
    "Final_MAG",
    "Completeness",
    "Contamination",
    "Quality_class",
    "Site",
    "Treatment",
    "Compartment",
    "GTDB_domain",
    "GTDB_phylum",
    "GTDB_class",
    "GTDB_order",
    "GTDB_family",
    "GTDB_genus",
    "GTDB_species",
]

all_markers = sorted(
    set().union(*mechanisms.values())
)

mag_metadata = {}
mag_kos = defaultdict(set)

with INPUT_FILE.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    missing_columns = [
        column for column in metadata_columns + ["KO"]
        if column not in reader.fieldnames
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing_columns)
        )

    for row in reader:
        mag = row["Final_MAG"].strip()

        if not mag:
            continue

        if mag not in mag_metadata:
            mag_metadata[mag] = {
                column: row.get(column, "")
                for column in metadata_columns
            }

        ko = row["KO"].strip()

        if ko:
            mag_kos[mag].add(ko)

mags = sorted(mag_metadata)

with MARKER_OUTPUT.open("w", newline="") as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(metadata_columns + all_markers)

    for mag in mags:
        metadata = [
            mag_metadata[mag][column]
            for column in metadata_columns
        ]

        marker_presence = [
            1 if ko in mag_kos[mag] else 0
            for ko in all_markers
        ]

        writer.writerow(metadata + marker_presence)

score_rows = []

with SCORE_OUTPUT.open("w", newline="") as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n"
    )

    mechanism_columns = []

    for mechanism in mechanisms:
        mechanism_columns.extend([
            f"{mechanism}_marker_count",
            f"{mechanism}_presence",
        ])

    writer.writerow(
        metadata_columns
        + mechanism_columns
        + [
            "Mechanisms_detected",
            "Mechanisms_total",
            "Osmotic_strategy_score",
            "Osmotic_strategy_percent",
        ]
    )

    for mag in mags:
        metadata = [
            mag_metadata[mag][column]
            for column in metadata_columns
        ]

        mechanism_values = []
        mechanisms_detected = 0

        for mechanism, markers in mechanisms.items():
            marker_count = len(
                mag_kos[mag].intersection(markers)
            )

            presence = 1 if marker_count > 0 else 0
            mechanisms_detected += presence

            mechanism_values.extend([
                marker_count,
                presence,
            ])

        mechanisms_total = len(mechanisms)
        score = mechanisms_detected / mechanisms_total
        percent = score * 100

        writer.writerow(
            metadata
            + mechanism_values
            + [
                mechanisms_detected,
                mechanisms_total,
                f"{score:.3f}",
                f"{percent:.1f}",
            ]
        )

        score_rows.append({
            "mag": mag,
            "detected": mechanisms_detected,
            "score": score,
        })

score_counts = Counter(
    row["detected"] for row in score_rows
)

scores = [
    row["score"] for row in score_rows
]

with SUMMARY_OUTPUT.open("w", newline="") as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "Summary_type",
        "Category",
        "Value",
    ])

    writer.writerow([
        "Dataset",
        "MAGs",
        len(mags),
    ])

    writer.writerow([
        "Dataset",
        "Osmotic_markers",
        len(all_markers),
    ])

    writer.writerow([
        "Dataset",
        "Osmotic_mechanisms",
        len(mechanisms),
    ])

    writer.writerow([
        "Statistic",
        "Mean_score",
        f"{statistics.mean(scores):.3f}",
    ])

    writer.writerow([
        "Statistic",
        "Median_score",
        f"{statistics.median(scores):.3f}",
    ])

    writer.writerow([
        "Statistic",
        "Minimum_score",
        f"{min(scores):.3f}",
    ])

    writer.writerow([
        "Statistic",
        "Maximum_score",
        f"{max(scores):.3f}",
    ])

    for detected in range(len(mechanisms) + 1):
        writer.writerow([
            "Distribution",
            f"{detected}_of_{len(mechanisms)}_mechanisms",
            score_counts.get(detected, 0),
        ])

print("Osmotic adaptation tables created successfully.")
print()
print(f"MAGs processed: {len(mags)}")
print(f"Curated KO markers: {len(all_markers)}")
print(f"Mechanisms assessed: {len(mechanisms)}")
print()
print(f"Marker table: {MARKER_OUTPUT}")
print(f"Strategy table: {SCORE_OUTPUT}")
print(f"Summary table: {SUMMARY_OUTPUT}")
print()
print("Mechanism-score distribution:")

for detected in range(len(mechanisms) + 1):
    print(
        f"  {detected}/{len(mechanisms)} mechanisms: "
        f"{score_counts.get(detected, 0)} MAGs"
    )
