#!/usr/bin/env python3

import csv
from collections import defaultdict, Counter
from pathlib import Path

BASE = Path("/path/to/metagenomics_project")

INPUT = (
    BASE
    / "results/functional_annotation/kofamscan/"
      "final_373_MAGs_gene_master_table.tsv"
)

OUTPUT_DIR = (
    BASE
    / "results/ecological_strategies/"
      "osmotic_adaptation/tables"
)

MATRIX_OUTPUT = OUTPUT_DIR / "trehalose_route_audit_by_MAG.tsv"
SUMMARY_OUTPUT = OUTPUT_DIR / "trehalose_route_audit_summary.tsv"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------
# Candidate trehalose biosynthesis routes
# ---------------------------------------------------------
#
# OtsAB:
#   K00697 = OtsA
#   K01087 = OtsB
#
# TreYZ:
#   K06044 = TreY
#   K01236 = TreZ
#
# TreS:
#   K13057 = trehalose synthase
#
# Bifunctional route:
#   K16055 = bifunctional trehalose
#            6-phosphate synthase/phosphatase
#
# NOTE:
# K10236/K10237/K10238 are deliberately NOT included here.
# They represent trehalose/maltose transport rather than
# trehalose biosynthesis.

ROUTES = {
    "OtsAB": {
        "mode": "all",
        "kos": {"K00697", "K01087"},
    },

    "TreYZ": {
        "mode": "all",
        "kos": {"K06044", "K01236"},
    },

   "TreS": {
    "mode": "all",
    "kos": {"K05343"},
},

"TreT": {
    "mode": "all",
    "kos": {"K13057"},
},

    "Bifunctional_K16055": {
        "mode": "all",
        "kos": {"K16055"},
    },
}

# ---------------------------------------------------------
# Read gene-level KO annotations
# ---------------------------------------------------------

mag_kos = defaultdict(set)
mag_metadata = {}

with INPUT.open("r", newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    required = {
        "Final_MAG",
        "KO",
        "Completeness",
        "Contamination",
        "Site",
        "Treatment",
        "Compartment",
        "GTDB_phylum",
        "GTDB_order",
    }

    if reader.fieldnames is None:
        raise ValueError("Input file has no header.")

    missing = required - set(reader.fieldnames)

    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(sorted(missing))
        )

    for row in reader:
        mag = row["Final_MAG"].strip()
        ko = row["KO"].strip()

        if not mag:
            continue

        if mag not in mag_metadata:
            mag_metadata[mag] = {
                "Completeness": row["Completeness"],
                "Contamination": row["Contamination"],
                "Site": row["Site"],
                "Treatment": row["Treatment"],
                "Compartment": row["Compartment"],
                "GTDB_phylum": row["GTDB_phylum"],
                "GTDB_order": row["GTDB_order"],
            }

        if ko:
            mag_kos[mag].add(ko)

mags = sorted(mag_metadata)

# ---------------------------------------------------------
# Function to evaluate a route
# ---------------------------------------------------------

def route_present(kos_present, rule):

    required_kos = rule["kos"]
    mode = rule["mode"]

    detected = len(kos_present & required_kos)
    total = len(required_kos)

    # Module evidence score: 0-1
    evidence = detected / total

    if mode == "all":
        present = int(required_kos <= kos_present)

    elif mode == "any":
        present = int(detected >= 1)

    else:
        raise ValueError(
            f"Unknown rule mode: {mode}"
        )

    return detected, total, evidence, present


# ---------------------------------------------------------
# Evaluate every MAG
# ---------------------------------------------------------

results = []

for mag in mags:

    kos_present = mag_kos[mag]

    row = {
        "Final_MAG": mag,
        **mag_metadata[mag],
    }

    positive_routes = []

    for route_name, rule in ROUTES.items():

        detected, total, evidence, present = route_present(
            kos_present,
            rule,
        )

        row[f"{route_name}_markers_detected"] = detected
        row[f"{route_name}_markers_total"] = total
        row[f"{route_name}_evidence"] = f"{evidence:.3f}"
        row[f"{route_name}_present"] = present

        if present:
            positive_routes.append(route_name)

    # -----------------------------------------------------
    # Preliminary trehalose capability
    #
    # Positive if at least one recognised biosynthetic
    # route is complete according to its rule.
    # -----------------------------------------------------

    row["Trehalose_routes_detected"] = len(positive_routes)

    row["Trehalose_route_names"] = (
        ";".join(positive_routes)
        if positive_routes
        else "None"
    )

    row["Trehalose_capability_preliminary"] = int(
        len(positive_routes) >= 1
    )

    results.append(row)


# ---------------------------------------------------------
# Write MAG-level audit table
# ---------------------------------------------------------

fieldnames = [
    "Final_MAG",
    "Completeness",
    "Contamination",
    "Site",
    "Treatment",
    "Compartment",
    "GTDB_phylum",
    "GTDB_order",
]

for route_name in ROUTES:
    fieldnames.extend([
        f"{route_name}_markers_detected",
        f"{route_name}_markers_total",
        f"{route_name}_evidence",
        f"{route_name}_present",
    ])

fieldnames.extend([
    "Trehalose_routes_detected",
    "Trehalose_route_names",
    "Trehalose_capability_preliminary",
])

with MATRIX_OUTPUT.open("w", newline="") as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(results)


# ---------------------------------------------------------
# Summaries
# ---------------------------------------------------------

route_counts = {}

for route_name in ROUTES:

    route_counts[route_name] = sum(
        int(row[f"{route_name}_present"])
        for row in results
    )

trehalose_positive = sum(
    row["Trehalose_capability_preliminary"]
    for row in results
)

route_number_distribution = Counter(
    row["Trehalose_routes_detected"]
    for row in results
)

combination_distribution = Counter(
    row["Trehalose_route_names"]
    for row in results
)


# ---------------------------------------------------------
# Write summary table
# ---------------------------------------------------------

summary_rows = []

for route_name, count in route_counts.items():

    summary_rows.append({
        "Summary_type": "Individual_route",
        "Category": route_name,
        "MAGs": count,
        "Percent": f"{100 * count / len(mags):.2f}",
    })

summary_rows.append({
    "Summary_type": "Combined_capability",
    "Category": "Any_validated_trehalose_route",
    "MAGs": trehalose_positive,
    "Percent": f"{100 * trehalose_positive / len(mags):.2f}",
})

for number in sorted(route_number_distribution):

    count = route_number_distribution[number]

    summary_rows.append({
        "Summary_type": "Number_of_routes",
        "Category": str(number),
        "MAGs": count,
        "Percent": f"{100 * count / len(mags):.2f}",
    })

for combination, count in sorted(
    combination_distribution.items(),
    key=lambda x: (-x[1], x[0]),
):

    summary_rows.append({
        "Summary_type": "Route_combination",
        "Category": combination,
        "MAGs": count,
        "Percent": f"{100 * count / len(mags):.2f}",
    })

with SUMMARY_OUTPUT.open("w", newline="") as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=[
            "Summary_type",
            "Category",
            "MAGs",
            "Percent",
        ],
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(summary_rows)


# ---------------------------------------------------------
# Print useful results
# ---------------------------------------------------------

print("Trehalose route audit completed successfully.")
print()
print(f"MAGs assessed: {len(mags)}")
print()

print("Individual complete routes:")

for route_name, count in route_counts.items():

    percent = 100 * count / len(mags)

    print(
        f"  {route_name:<25} "
        f"{count:>4} MAGs "
        f"({percent:>6.2f}%)"
    )

print()

print(
    "MAGs with >=1 complete trehalose route: "
    f"{trehalose_positive} "
    f"({100 * trehalose_positive / len(mags):.2f}%)"
)

print()
print("Number of complete trehalose routes per MAG:")

for number in sorted(route_number_distribution):

    count = route_number_distribution[number]

    print(
        f"  {number} routes: "
        f"{count} MAGs "
        f"({100 * count / len(mags):.2f}%)"
    )

print()
print("Most common route combinations:")

for combination, count in combination_distribution.most_common(15):

    print(
        f"  {combination:<55} "
        f"{count:>4}"
    )

print()
print(f"MAG-level matrix: {MATRIX_OUTPUT}")
print(f"Summary table:    {SUMMARY_OUTPUT}")
