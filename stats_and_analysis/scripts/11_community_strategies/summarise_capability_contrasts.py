#!/usr/bin/env python3
import csv
from pathlib import Path

base = Path("results/ecological_strategies/community_weighted_GRESP")
source = base / "community_capabilities_20_samples_long.tsv"

with source.open() as f:
    rows = [
        r for r in csv.DictReader(f, delimiter="\t")
        if r["Filter"] == "Unfiltered"
    ]

data = {}
capabilities = set()

for r in rows:
    key = (r["Habitat"], r["Site"], r["Treatment"],
           r["GRESP"], r["Capability"])

    if key in data:
        raise ValueError(f"Duplicate: {key}")

    data[key] = float(r["Weighted_prevalence_pct"])
    capabilities.add((r["GRESP"], r["Capability"]))

if len(rows) != 740 or len(capabilities) != 37:
    raise ValueError("Unexpected sample or capability count")

contrasts = [
    ("RZ", "BESE", "BARE"),
    ("RZ", "BESE", "NAT"),
    ("RZ", "BARE", "NAT"),
    ("BS", "BESE", "BARE"),
    ("BS", "NAT", "DEGR"),
    ("BS", "BESE", "NAT"),
    ("BS", "BARE", "NAT"),
]

sites = ["ACP", "BASS", "TRUG"]
results = []

for habitat, treatment_a, treatment_b in contrasts:
    for group, capability in sorted(capabilities):
        differences = {}

        for site in sites:
            a = (habitat, site, treatment_a, group, capability)
            b = (habitat, site, treatment_b, group, capability)

            if a in data and b in data:
                differences[site] = data[a] - data[b]

        if not differences:
            continue

        values = list(differences.values())
        mean = sum(values) / len(values)

        consistent = (
            all(v > 0 for v in values)
            or all(v < 0 for v in values)
        )

        results.append({
            "Habitat": habitat,
            "Contrast": f"{treatment_a}_minus_{treatment_b}",
            "GRESP": group,
            "Capability": capability,
            "N_sites": len(values),
            "Mean_difference_pp": round(mean, 4),
            "Consistent_direction": consistent,
            "ACP": differences.get("ACP", ""),
            "BASS": differences.get("BASS", ""),
            "TRUG": differences.get("TRUG", "")
        })

output = base / "capability_site_paired_contrasts.tsv"

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

for habitat, a, b in [
    ("RZ", "BESE", "BARE"),
    ("BS", "BESE", "BARE"),
    ("BS", "NAT", "DEGR")
]:
    subset = [
        r for r in results
        if r["Habitat"] == habitat
        and r["Contrast"] == f"{a}_minus_{b}"
    ]

    subset.sort(
        key=lambda r: abs(r["Mean_difference_pp"]),
        reverse=True
    )

    print(f"\n=== {habitat}: {a} minus {b} ===")
    print("GRESP\tCapability\tN\tDifference_pp\tConsistent")

    for r in subset[:12]:
        print(
            f'{r["GRESP"]}\t{r["Capability"]}\t'
            f'{r["N_sites"]}\t'
            f'{r["Mean_difference_pp"]:+.2f}\t'
            f'{r["Consistent_direction"]}'
        )

print("\nSUCCESS")
print("Capability entries:", len(capabilities))
print("Contrasts saved:", len(results))
print("Output:", output)
