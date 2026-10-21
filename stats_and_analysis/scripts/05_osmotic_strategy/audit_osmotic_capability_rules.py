#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import defaultdict

INPUT = Path(
    "results/functional_annotation/kofamscan/"
    "final_373_MAGs_gene_master_table.tsv"
)

OUTPUT_SUMMARY = Path(
    "results/ecological_strategies/osmotic_adaptation/"
    "tables/osmotic_capability_rule_audit.tsv"
)

OUTPUT_MATRIX = Path(
    "results/ecological_strategies/osmotic_adaptation/"
    "tables/osmotic_capability_presence_preliminary.tsv"
)

# Preliminary pathway/system rules.
# "all" means every listed KO must occur in the same MAG.
# "any" means at least one listed KO must occur.
rules = {
    "Trehalose_OtsAB": {
        "mode": "all",
        "kos": ["K00697", "K01087"],
    },
    "Trehalose_K10236_K10237_K10238_system": {
        "mode": "all",
        "kos": ["K10236", "K10237", "K10238"],
    },
    "Trehalose_bifunctional_K16055": {
        "mode": "any",
        "kos": ["K16055"],
    },
    "Osmoprotectant_transporter_K05845_K05846_K05847": {
        "mode": "all",
        "kos": ["K05845", "K05846", "K05847"],
    },
    "ProVWX_transporter_K02000_K02001_K02002": {
        "mode": "all",
        "kos": ["K02000", "K02001", "K02002"],
    },
    "Ectoine_EctABC": {
        "mode": "all",
        "kos": ["K06718", "K00836", "K06720"],
    },
    "Mechanosensitive_channel_any": {
        "mode": "any",
        "kos": ["K03282", "K03442", "K05802", "K16053", "K22044"],
    },
    "Potassium_system_K03498_K03499": {
        "mode": "all",
        "kos": ["K03498", "K03499"],
    },
    "Potassium_system_K09696_K09697": {
        "mode": "all",
        "kos": ["K09696", "K09697"],
    },
    "Ion_system_K01545_K01546_K01547_K01548": {
        "mode": "all",
        "kos": ["K01545", "K01546", "K01547", "K01548"],
    },
    "Sodium_homeostasis_any": {
        "mode": "any",
        "kos": ["K15987", "K24160", "K03549", "K11105"],
    },
}

with INPUT.open("r", newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    if not reader.fieldnames:
        raise ValueError("The input table has no header.")

    mag_column_candidates = [
        "Final_MAG",
        "MAG",
        "Genome",
        "Bin",
    ]

    ko_column_candidates = [
        "KO",
        "KEGG_Orthology",
        "K_number",
    ]

    mag_column = next(
        (column for column in mag_column_candidates if column in reader.fieldnames),
        None,
    )

    ko_column = next(
        (column for column in ko_column_candidates if column in reader.fieldnames),
        None,
    )

    if mag_column is None:
        raise ValueError(
            "Could not identify the MAG column. "
            f"Columns found: {reader.fieldnames}"
        )

    if ko_column is None:
        raise ValueError(
            "Could not identify the KO column. "
            f"Columns found: {reader.fieldnames}"
        )

    mag_kos = defaultdict(set)

    for row in reader:
        mag = row.get(mag_column, "").strip()
        ko_field = row.get(ko_column, "").strip()

        if not mag or not ko_field:
            continue

        # Handles a single KO or multiple KOs separated by common delimiters.
        cleaned = (
            ko_field.replace(",", ";")
            .replace("|", ";")
            .replace(" ", ";")
        )

        for ko in cleaned.split(";"):
            ko = ko.strip()

            if ko.startswith("K") and ko[1:].isdigit():
                mag_kos[mag].add(ko)

mags = sorted(mag_kos)

if not mags:
    raise ValueError(f"No MAG–KO records were found in {INPUT}")

all_rule_kos = sorted(
    {
        ko
        for rule in rules.values()
        for ko in rule["kos"]
    }
)

summary_rows = []
matrix_rows = []

for mag in mags:
    present_kos = mag_kos[mag]
    matrix_row = {"Final_MAG": mag}

    for rule_name, rule in rules.items():
        kos = rule["kos"]
        mode = rule["mode"]

        detected_count = sum(ko in present_kos for ko in kos)

        if mode == "all":
            positive = int(detected_count == len(kos))
        elif mode == "any":
            positive = int(detected_count >= 1)
        else:
            raise ValueError(f"Unsupported rule mode: {mode}")

        matrix_row[rule_name] = positive

    matrix_rows.append(matrix_row)

for rule_name, rule in rules.items():
    positive_mags = sum(
        row[rule_name] == 1
        for row in matrix_rows
    )

    summary_rows.append(
        {
            "Capability_rule": rule_name,
            "Rule_mode": rule["mode"],
            "Required_KOs": ";".join(rule["kos"]),
            "Number_of_KOs": len(rule["kos"]),
            "Positive_MAGs": positive_mags,
            "Total_MAGs": len(mags),
            "Prevalence_percent": round(
                positive_mags / len(mags) * 100,
                2,
            ),
        }
    )

OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT_SUMMARY.open("w", newline="") as outfile:
    fieldnames = [
        "Capability_rule",
        "Rule_mode",
        "Required_KOs",
        "Number_of_KOs",
        "Positive_MAGs",
        "Total_MAGs",
        "Prevalence_percent",
    ]

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(summary_rows)

with OUTPUT_MATRIX.open("w", newline="") as outfile:
    fieldnames = ["Final_MAG"] + list(rules.keys())

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(matrix_rows)

print("Preliminary osmotic capability audit completed.")
print()
print(f"MAGs assessed: {len(mags)}")
print(f"Candidate KOs assessed: {len(all_rule_kos)}")
print(f"Summary: {OUTPUT_SUMMARY.resolve()}")
print(f"Presence matrix: {OUTPUT_MATRIX.resolve()}")
print()

header = (
    f"{'Capability rule':<55}"
    f"{'Positive MAGs':>15}"
    f"{'Prevalence %':>15}"
)

print(header)
print("-" * len(header))

for row in summary_rows:
    print(
        f"{row['Capability_rule']:<55}"
        f"{row['Positive_MAGs']:>15}"
        f"{row['Prevalence_percent']:>15.2f}"
    )
