#!/usr/bin/env python3

from pathlib import Path
import csv

INPUT = Path(
    "results/ecological_strategies/osmotic_adaptation/"
    "tables/osmotic_marker_presence.tsv"
)

OUTPUT = Path(
    "results/ecological_strategies/osmotic_adaptation/"
    "tables/osmotic_marker_prevalence.tsv"
)

with INPUT.open("r", newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")
    rows = list(reader)

if not rows:
    raise ValueError(f"No data rows found in {INPUT}")

ko_columns = [
    column
    for column in reader.fieldnames
    if column.startswith("K") and column[1:].isdigit()
]

summary = []

for ko in ko_columns:
    detected = sum(
        1
        for row in rows
        if float(row.get(ko, 0) or 0) > 0
    )

    total = len(rows)

    summary.append(
        {
            "KO": ko,
            "MAGs_detected": detected,
            "MAGs_total": total,
            "Prevalence_percent": round((detected / total) * 100, 2),
        }
    )

summary.sort(
    key=lambda record: (-record["MAGs_detected"], record["KO"])
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as outfile:
    fieldnames = [
        "KO",
        "MAGs_detected",
        "MAGs_total",
        "Prevalence_percent",
    ]

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(summary)

print("Osmotic marker prevalence table created successfully.")
print()
print(f"MAGs assessed: {len(rows)}")
print(f"KO markers assessed: {len(ko_columns)}")
print(f"Output: {OUTPUT.resolve()}")
print()

header = (
    f"{'KO':<10}"
    f"{'MAGs_detected':>15}"
    f"{'MAGs_total':>12}"
    f"{'Prevalence_%':>15}"
)

print(header)
print("-" * len(header))

for record in summary:
    print(
        f"{record['KO']:<10}"
        f"{record['MAGs_detected']:>15}"
        f"{record['MAGs_total']:>12}"
        f"{record['Prevalence_percent']:>15.2f}"
    )
