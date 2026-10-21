#!/usr/bin/env python3

import csv
from pathlib import Path

BASE = Path("results/ecological_strategies")
AB = Path("mapping/competitive_MAGs/MAG_abundance")

GROUPS = {
    "Osmotic": (
        "osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv",
        """Capability_Trehalose
Capability_Compatible_solute_uptake
Capability_Ectoine_biosynthesis
Capability_Mechanosensitive_protection
Capability_Ion_homeostasis""".split()
    ),

    "Redox": (
        "redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv",
        """Capability_Aerobic_respiration
Capability_LowO2_respiration
Capability_Nitrate_respiration
Capability_Nitrite_to_NO
Capability_DNRA
Capability_NO_reduction
Capability_N2O_reduction
Capability_Dissimilatory_sulfate_reduction_core
Capability_Reduced_sulfur_oxidation
Capability_Sulfite_oxidation""".split()
    ),

    "Nitrogen": (
        "nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv",
        """Capability_Nitrogen_fixation
Capability_Ammonium_uptake
Capability_Assimilatory_nitrate_reduction
Capability_Urea_utilisation
Capability_DNRA
Capability_Complete_denitrification
Capability_Nitrification
Capability_Anammox""".split()
    ),

    "Carbon_P": (
        "carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv",
        """Capability_Carbon_fixation
Capability_C1_utilisation
Capability_High_affinity_phosphate_uptake
Capability_Organic_P_mobilisation
Capability_Phosphonate_acquisition
Capability_PolyP_turnover""".split()
    ),

    "Iron": (
        "iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv",
        """Ferrous_iron_uptake
Direct_ferric_iron_uptake
Siderophore_biosynthesis
Iron_complex_uptake""".split()
    ),

    "Colonisation": (
        "colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv",
        """Chemotaxis
Motility
Attachment_biofilm_matrix
Direct_plant_interaction""".split()
    )
}


def read_tsv(path, idcol):
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = reader.fieldnames
        rows = list(reader)

    ids = [r[idcol] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate IDs in {path}")

    return {r[idcol]: r for r in rows}, columns


def binary(value):
    value = value.strip().lower()
    if value in ("1", "1.0", "true"):
        return 1
    if value in ("0", "0.0", "false"):
        return 0
    raise ValueError(f"Unexpected binary capability value: {value}")


depth, depth_cols = read_tsv(
    AB / "final_373_MAGs_mean_depth_matrix.tsv", "MAG_ID"
)

breadth, breadth_cols = read_tsv(
    AB / "final_373_MAGs_breadth_matrix.tsv", "MAG_ID"
)

samples = depth_cols[1:]
mags = list(depth)

assert len(mags) == 373
assert len(samples) == 20
assert depth_cols == breadth_cols
assert set(depth) == set(breadth)

capabilities = []

for group, (filename, fields) in GROUPS.items():

    data, columns = read_tsv(BASE / filename, "Final_MAG")

    if len(data) != 373 or set(data) != set(depth):
        raise ValueError(f"MAG mismatch: {group}")

    for field in fields:
        if field not in columns:
            raise ValueError(f"Missing: {group} / {field}")

        calls = {
            mag: binary(data[mag][field])
            for mag in mags
        }

        capabilities.append((group, field, calls))

assert len(capabilities) == 37

outdir = BASE / "community_weighted_GRESP"
outdir.mkdir(parents=True, exist_ok=True)

output = outdir / "community_capabilities_20_samples_long.tsv"

records = []

for sample in samples:

    parts = sample.split("_")
    site, treatment, habitat = parts[1:4]

    depths = {
        mag: float(depth[mag][sample])
        for mag in mags
    }

    if any(v < 0 for v in depths.values()):
        raise ValueError("Negative MAG depth")

    total_depth = sum(depths.values())
    if total_depth <= 0:
        raise ValueError(f"Zero total depth: {sample}")

    for filter_name, threshold in [
        ("Unfiltered", 0),
        ("Breadth10", 0.10),
        ("Breadth50", 0.50)
    ]:

        retained = [
            mag for mag in mags
            if depths[mag] > 0
            and float(breadth[mag][sample]) >= threshold
        ]

        retained_depth = sum(depths[m] for m in retained)

        if retained_depth <= 0:
            raise ValueError(f"No retained MAGs: {sample}")

        for group, capability, calls in capabilities:

            weighted = sum(
                depths[mag] * calls[mag]
                for mag in retained
            )

            prevalence = 100 * weighted / retained_depth

            if not 0 <= prevalence <= 100.000001:
                raise ValueError("Invalid weighted prevalence")

            records.append({
                "Sample": sample,
                "Site": site,
                "Treatment": treatment,
                "Habitat": habitat,
                "Filter": filter_name,
                "GRESP": group,
                "Capability": capability,
                "Weighted_prevalence_pct": round(prevalence, 6),
                "Retained_depth_pct": round(
                    100 * retained_depth / total_depth, 4
                )
            })

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(records[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(records)

print("SUCCESS")
print("MAGs:", len(mags))
print("Metagenomic samples:", len(samples))
print("GRESP dimensions:", len(GROUPS))
print("Capability entries:", len(capabilities))
print("Filtering approaches:", 3)
print("Output rows:", len(records))
print("Saved:", output)
