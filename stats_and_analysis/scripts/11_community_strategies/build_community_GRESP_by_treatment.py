#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict

base = Path("results/ecological_strategies")

files = {
    "Osmotic": base / "osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv",
    "Redox": base / "redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv",
    "Nitrogen": base / "nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv",
    "Carbon_P": base / "carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv",
    "Iron": base / "iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv",
    "Colonisation": base / "colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv"
}

def read_table(path, id_column):
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = reader.fieldnames
        rows = list(reader)

    data = {}
    for row in rows:
        mag = row[id_column]
        if mag in data:
            raise ValueError(f"Duplicate MAG: {mag}")
        data[mag] = row

    return data, columns

# Read all six MAG-level strategies
strategies = {}

for label, path in files.items():
    data, columns = read_table(path, "Final_MAG")

    score_cols = [
        c for c in columns
        if c.lower().endswith("_strategy_score")
    ]

    if not score_cols:
        score_cols = [
            c for c in columns
            if c.lower().endswith("_score")
        ]

    if len(score_cols) != 1:
        raise ValueError(
            f"{label}: cannot identify unique strategy score. "
            f"Candidates: {score_cols}"
        )

    score_col = score_cols[0]
    strategies[label] = {
        mag: float(row[score_col])
        for mag, row in data.items()
    }

    print(label, ":", score_col, len(data), "MAGs")

# Read abundances and breadth
abundance_dir = Path(
    "mapping/competitive_MAGs/MAG_abundance"
)

depth, depth_columns = read_table(
    abundance_dir / "final_373_MAGs_mean_depth_matrix.tsv",
    "MAG_ID"
)

breadth, breadth_columns = read_table(
    abundance_dir / "final_373_MAGs_breadth_matrix.tsv",
    "MAG_ID"
)

samples = depth_columns[1:]

if samples != breadth_columns[1:]:
    raise ValueError("Sample column mismatch")

for label, scores in strategies.items():
    if set(scores) != set(depth):
        raise ValueError(f"MAG mismatch: {label}")

if set(depth) != set(breadth):
    raise ValueError("Breadth MAG mismatch")

# Compute community-weighted strategy scores
results = []

for sample in samples:
    parts = sample.split("_")
    site, treatment, habitat = parts[1:4]

    total_depth = sum(
        float(depth[mag][sample])
        for mag in depth
    )

    for filter_name, threshold in [
        ("Unfiltered", 0),
        ("Breadth10", 0.10),
        ("Breadth50", 0.50)
    ]:
        retained = [
            mag for mag in depth
            if float(breadth[mag][sample]) >= threshold
            and float(depth[mag][sample]) > 0
        ]

        retained_depth = sum(
            float(depth[mag][sample])
            for mag in retained
        )

        if retained_depth <= 0:
            raise ValueError(f"No MAGs retained: {sample}")

        row = {
            "Sample": sample,
            "Site": site,
            "Treatment": treatment,
            "Habitat": habitat,
            "Filter": filter_name,
            "MAGs_retained": len(retained),
            "Retained_depth_pct": 100 * retained_depth / total_depth
        }

        for label, scores in strategies.items():
            weighted_sum = sum(
                float(depth[mag][sample]) * scores[mag]
                for mag in retained
            )

            row[label] = weighted_sum / retained_depth

        results.append(row)

# Save individual sample scores
outdir = base / "community_weighted_GRESP"
outdir.mkdir(parents=True, exist_ok=True)

output = outdir / "community_GRESP_20_samples.tsv"

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

# Summarise by treatment and habitat
groups = defaultdict(list)

for row in results:
    if row["Filter"] == "Breadth10":
        groups[(row["Habitat"], row["Treatment"])].append(row)

summary = outdir / "community_GRESP_treatment_summary.tsv"

with summary.open("w", newline="") as f:
    columns = ["Habitat", "Treatment", "N_sites"] + list(strategies)
    writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
    writer.writeheader()

    for habitat, treatment in sorted(groups):
        group = groups[(habitat, treatment)]
        row = {
            "Habitat": habitat,
            "Treatment": treatment,
            "N_sites": len(group)
        }
        for label in strategies:
            row[label] = sum(r[label] for r in group) / len(group)
        writer.writerow(row)

print("\nSUCCESS")
print("Metagenomic samples:", len(samples))
print("Ecological strategies:", len(strategies))
print("Filter versions:", 3)
print("Sample-level output:", output)
print("Treatment summary:", summary)
