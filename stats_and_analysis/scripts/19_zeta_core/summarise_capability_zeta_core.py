#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict

base = Path("results/ecological_strategies/metagenomic_zeta")
source = base / "Capability_core_occurrence_37_traits.tsv"
output = base / "Capability_core_across_conditions_Breadth10_5pct.tsv"

with source.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

rows = [
    r for r in rows
    if r["Filter"] == "Breadth10"
    and float(r["Threshold_pct"]) == 5.0
    and r["Design"] == "Within_treatment"
]

features = sorted({
    (r["GRESP"], r["Capability"]) for r in rows
})

assert len(features) == 37

conditions = {
    "BS": ["NAT", "BESE", "BARE", "DEGR"],
    "RZ": ["NAT", "BESE", "BARE"]
}

lookup = {
    (r["Habitat"], r["Group"], r["GRESP"], r["Capability"]): r
    for r in rows
}

results = []

for habitat, treatments in conditions.items():
    for group, capability in features:
        record = {
            "Habitat": habitat,
            "GRESP": group,
            "Capability": capability
        }

        core = []

        for treatment in treatments:
            r = lookup[(habitat, treatment, group, capability)]
            present = int(r["Core_all_samples"])
            record[treatment + "_core"] = present
            record[treatment + "_sites_present"] = (
                r["Samples_present"] + "/" + r["N_samples"]
            )
            core.append(present)

        record["Core_across_all_conditions"] = int(all(core))
        record["N_conditions_core"] = sum(core)

        results.append(record)

columns = [
    "Habitat", "GRESP", "Capability",
    "NAT_core", "BESE_core", "BARE_core",
    "DEGR_core", "NAT_sites_present",
    "BESE_sites_present", "BARE_sites_present",
    "DEGR_sites_present",
    "Core_across_all_conditions",
    "N_conditions_core"
]

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=columns, delimiter="\t",
        extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(results)

for habitat, treatments in conditions.items():
    subset = [r for r in results if r["Habitat"] == habitat]
    universal = [
        r for r in subset if r["Core_across_all_conditions"] == 1
    ]
    variable = [
        r for r in subset if 0 < r["N_conditions_core"] < len(treatments)
    ]
    never = [
        r for r in subset if r["N_conditions_core"] == 0
    ]

    print("\n=== {} ===".format(habitat))
    print("Core across all conditions:", len(universal))
    print("Core in some but not all conditions:", len(variable))
    print("Not consistently represented in any condition:", len(never))

    print("\nUNIVERSAL CORE BY ECOLOGICAL STRATEGY")
    for group in ("Osmotic","Redox","Nitrogen",
                  "Carbon_P","Iron","Colonisation"):
        selected = [
            r["Capability"] for r in universal
            if r["GRESP"] == group
        ]
        print("{} ({}): {}".format(
            group, len(selected),
            ", ".join(selected) if selected else "None"
        ))

    print("\nCONDITION-ASSOCIATED CAPABILITIES")
    for r in variable:
        present = [
            t for t in treatments if r[t+"_core"] == 1
        ]
        print("{} / {}: {}".format(
            r["GRESP"],r["Capability"],",".join(present)
        ))

    print("\nNEVER CORE")
    for r in never:
        print("{} / {}".format(r["GRESP"],r["Capability"]))

print("\nSUCCESS")
print("Output:",output)
