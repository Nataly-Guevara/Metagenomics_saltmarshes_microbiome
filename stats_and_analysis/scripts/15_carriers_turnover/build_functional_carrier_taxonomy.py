#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict

BASE = Path("results/ecological_strategies")
AB = Path("mapping/competitive_MAGs/MAG_abundance")
OUT = BASE / "functional_carrier_taxonomy"
OUT.mkdir(parents=True, exist_ok=True)

GROUPS = {
    "Osmotic": (
        "osmotic_adaptation/tables/osmotic_strategy_final_373_MAGs.tsv",
        "Capability_Trehalose Capability_Compatible_solute_uptake Capability_Ectoine_biosynthesis Capability_Mechanosensitive_protection Capability_Ion_homeostasis"
    ),
    "Redox": (
        "redox_flexibility/tables/redox_strategy_final_v2_373_MAGs.tsv",
        "Capability_Aerobic_respiration Capability_LowO2_respiration Capability_Nitrate_respiration Capability_Nitrite_to_NO Capability_DNRA Capability_NO_reduction Capability_N2O_reduction Capability_Dissimilatory_sulfate_reduction_core Capability_Reduced_sulfur_oxidation Capability_Sulfite_oxidation"
    ),
    "Nitrogen": (
        "nitrogen_acquisition_transformation/tables/nitrogen_strategy_final_373_MAGs.tsv",
        "Capability_Nitrogen_fixation Capability_Ammonium_uptake Capability_Assimilatory_nitrate_reduction Capability_Urea_utilisation Capability_DNRA Capability_Complete_denitrification Capability_Nitrification Capability_Anammox"
    ),
    "Carbon_P": (
        "carbon_phosphorus_resource/tables/carbon_phosphorus_strategy_final_373_MAGs.tsv",
        "Capability_Carbon_fixation Capability_C1_utilisation Capability_High_affinity_phosphate_uptake Capability_Organic_P_mobilisation Capability_Phosphonate_acquisition Capability_PolyP_turnover"
    ),
    "Iron": (
        "iron_acquisition/tables/iron_strategy_final_373_MAGs.tsv",
        "Ferrous_iron_uptake Direct_ferric_iron_uptake Siderophore_biosynthesis Iron_complex_uptake"
    ),
    "Colonisation": (
        "colonisation_host_interaction/tables/colonisation_host_strategy_final_373_MAGs.tsv",
        "Chemotaxis Motility Attachment_biofilm_matrix Direct_plant_interaction"
    )
}

def load(path, key):
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == len({r[key] for r in rows})
    return {r[key]: r for r in rows}

depth = load(AB / "final_373_MAGs_mean_depth_matrix.tsv", "MAG_ID")
breadth = load(AB / "final_373_MAGs_breadth_matrix.tsv", "MAG_ID")
taxonomy = load(
    Path("results/final_summaries/combined_binners/final_373_MAGs_GTDB_taxonomy.tsv"),
    "Final_MAG"
)

mags = sorted(depth)
assert len(mags) == 373
assert set(mags) == set(breadth) == set(taxonomy)

samples = [x for x in next(iter(depth.values())) if x != "MAG_ID"]
assert len(samples) == 20

def taxon(mag, level):
    classification = taxonomy[mag]["GTDB_classification"]
    prefix = {"Phylum": "p__", "Order": "o__", "Genus": "g__"}[level]
    for term in classification.split(";"):
        if term.startswith(prefix):
            value = term[len(prefix):]
            return value if value else "Unclassified"
    return "Unclassified"

capabilities = []
for group, (filename, fields) in GROUPS.items():
    table = load(BASE / filename, "Final_MAG")
    assert set(table) == set(mags)
    for field in fields.split():
        assert field in next(iter(table.values()))
        carriers = {
            mag for mag in mags
            if float(table[mag][field]) == 1
        }
        capabilities.append((group, field, carriers))

assert len(capabilities) == 37

carrier_rows = []
lineage_rows = []

for group, capability, carriers in capabilities:
    for mag in sorted(carriers):
        carrier_rows.append({
            "GRESP": group,
            "Capability": capability,
            "MAG_ID": mag,
            "Phylum": taxon(mag, "Phylum"),
            "Order": taxon(mag, "Order"),
            "Genus": taxon(mag, "Genus")
        })

    for sample in samples:
        _, site, treatment, habitat = sample.split("_", 4)[:4]

        detected = [
            mag for mag in mags
            if float(depth[mag][sample]) > 0
            and float(breadth[mag][sample]) >= 0.10
        ]

        total = sum(float(depth[m][sample]) for m in detected)
        assert total > 0

        carriers_here = [m for m in detected if m in carriers]
        carrier_depth = sum(float(depth[m][sample]) for m in carriers_here)

        for level in ("Phylum", "Order"):
            contributions = defaultdict(float)

            for mag in carriers_here:
                contributions[taxon(mag, level)] += float(depth[mag][sample])

            for name, value in contributions.items():
                lineage_rows.append({
                    "Sample": sample,
                    "Site": site,
                    "Treatment": treatment,
                    "Habitat": habitat,
                    "GRESP": group,
                    "Capability": capability,
                    "Taxonomic_level": level,
                    "Taxon": name,
                    "Carrier_MAG_count": sum(
                        taxon(m, level) == name for m in carriers_here
                    ),
                    "Community_abundance_pct": 100 * value / total,
                    "Share_of_capability_carriers_pct":
                        100 * value / carrier_depth if carrier_depth > 0 else 0
                })

def save(rows, filename):
    assert rows
    with (OUT / filename).open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

save(carrier_rows, "MAG_capability_carrier_taxonomy.tsv")
save(lineage_rows, "Sample_capability_lineage_contributions_Breadth10.tsv")

print("SUCCESS: Functional carrier taxonomy analysis")
print("MAGs:", len(mags))
print("Metagenomes:", len(samples))
print("Capabilities:", len(capabilities))
print("Capability-positive MAG records:", len(carrier_rows))
print("Sample-lineage contribution records:", len(lineage_rows))
print("Output directory:", OUT)
