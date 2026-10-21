import csv
from pathlib import Path

BASE = Path("results/ecological_strategies")
AB = Path("mapping/competitive_MAGs/MAG_abundance")

FILES = {
    "Osmotic": "osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv",
    "Redox": "redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv",
    "Nitrogen": "nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv",
    "Carbon_P": "carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv",
    "Iron": "iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv",
    "Colonisation": "colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv"
}

def read_table(path, idcol):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        columns = reader.fieldnames
    data = {r[idcol]: r for r in rows}
    assert len(data) == len(rows), f"Duplicate IDs: {path}"
    return data, columns

depth, cols = read_table(
    AB / "final_373_MAGs_mean_depth_matrix.tsv", "MAG_ID"
)
breadth, bcols = read_table(
    AB / "final_373_MAGs_breadth_matrix.tsv", "MAG_ID"
)

assert cols == bcols
assert set(depth) == set(breadth)
assert len(depth) == 373

samples = cols[1:]
assert len(samples) == 20

# Reuse the exact 37 curated capabilities from Figure 2
source = BASE / "community_weighted_GRESP/community_capabilities_20_samples_long.tsv"

with open(source) as f:
    source_rows = list(csv.DictReader(f, delimiter="\t"))

assert len(source_rows) == 20 * 37 * 3

capabilities = list(dict.fromkeys(
    (r["GRESP"], r["Capability"]) for r in source_rows
))
assert len(capabilities) == 37

reference = {
    (r["Sample"], r["Filter"], r["GRESP"], r["Capability"]):
    float(r["Weighted_prevalence_pct"])
    for r in source_rows
}
assert len(reference) == len(source_rows)

calls = {}

for group, filename in FILES.items():
    data, columns = read_table(BASE / filename, "Final_MAG")
    assert set(data) == set(depth)

    for g, cap in capabilities:
        if g != group:
            continue

        assert cap in columns, (group, cap)

        values = {}
        for mag in depth:
            value = float(data[mag][cap])
            assert value in (0, 1)
            values[mag] = int(value)

        calls[(group, cap)] = values

assert len(calls) == 37

results = []
max_difference = 0

for sample in samples:
    parts = sample.split("_")
    site, treatment, habitat = parts[1:4]

    for label, threshold in [
        ("Breadth10", 0.10),
        ("Breadth50", 0.50)
    ]:
        retained = [
            mag for mag in depth
            if float(depth[mag][sample]) > 0
            and float(breadth[mag][sample]) >= threshold
        ]

        total_depth = sum(
            float(depth[m][sample]) for m in retained
        )
        assert total_depth > 0

        for group, cap in capabilities:
            carriers = [
                mag for mag in retained
                if calls[(group, cap)][mag] == 1
            ]

            weights = [
                float(depth[m][sample]) for m in carriers
            ]

            carrier_depth = sum(weights)

            effective = (
                carrier_depth ** 2 / sum(w*w for w in weights)
                if carrier_depth > 0 else 0.0
            )

            prevalence = 100 * carrier_depth / total_depth

            expected = reference[
                (sample, label, group, cap)
            ]
            difference = abs(prevalence - expected)
            max_difference = max(max_difference, difference)

            if difference > 0.02:
                raise ValueError(
                    f"Prevalence mismatch: {sample}, "
                    f"{label}, {group}, {cap}"
                )

            results.append({
                "Sample": sample,
                "Site": site,
                "Treatment": treatment,
                "Habitat": habitat,
                "Filter": label,
                "GRESP": group,
                "Capability": cap,
                "Carrier_richness": len(carriers),
                "Effective_carriers": round(effective, 6),
                "Weighted_prevalence_pct": round(prevalence, 6)
            })

assert len(results) == 20 * 37 * 2

outdir = BASE / "functional_redundancy"
outdir.mkdir(parents=True, exist_ok=True)
output = outdir / "functional_redundancy_20samples_37capabilities.tsv"

with open(output, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=results[0].keys(), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("SUCCESS")
print("Metagenomes:", len(samples))
print("Capabilities:", len(capabilities))
print("Filters: Breadth10, Breadth50")
print("Results:", len(results))
print("Maximum prevalence QC difference:",
      round(max_difference, 8), "percentage points")
print("Output:", output)
