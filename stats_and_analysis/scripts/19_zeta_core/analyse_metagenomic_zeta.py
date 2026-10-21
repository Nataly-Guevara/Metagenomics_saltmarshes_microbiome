#!/usr/bin/env python3
import csv
import itertools
from collections import defaultdict
from pathlib import Path

BASE = Path("results/ecological_strategies")
AB = Path("mapping/competitive_MAGs/MAG_abundance")
KO = Path("results/functional_annotation/sample_KO_profiles")
OUT = BASE / "metagenomic_zeta"
OUT.mkdir(parents=True, exist_ok=True)

def read_tsv(path):
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader), reader.fieldnames

# MAG depth and genome breadth
depth_rows, depth_cols = read_tsv(
    AB / "final_373_MAGs_mean_depth_matrix.tsv"
)
breadth_rows, breadth_cols = read_tsv(
    AB / "final_373_MAGs_breadth_matrix.tsv"
)

assert len(depth_rows) == len(breadth_rows) == 373
assert depth_cols == breadth_cols

depth = {r["MAG_ID"]: r for r in depth_rows}
breadth = {r["MAG_ID"]: r for r in breadth_rows}

assert set(depth) == set(breadth)

samples = depth_cols[1:]
assert len(samples) == 20

# Existing reconstructed KO functional potential
ko_rows, ko_cols = read_tsv(
    KO / "sample_KO_relative_MAG_abundance_weighted.tsv"
)

ko = {r["Sample"]: r for r in ko_rows}
assert len(ko) == 20
assert set(samples) == set(ko)
assert len(ko_cols) - 1 == 6679

metadata = {}
for sample in samples:
    parts = sample.split("_")
    metadata[sample] = {
        "Site": parts[1],
        "Treatment": parts[2],
        "Habitat": parts[3]
    }

# Incidence matrices represented as sets
presence = {
    "MAG": {},
    "KO": {}
}

for sample in samples:
    presence["MAG"][sample] = {
        mag for mag in depth
        if float(depth[mag][sample]) > 0
        and float(breadth[mag][sample]) >= 0.10
    }

    presence["KO"][sample] = {
        feature for feature in ko_cols[1:]
        if float(ko[sample][feature]) > 0
    }

# Define biologically interpretable sample groups.
groups = {}

# Within-site: compare ecological conditions
for site in ("ACP", "BASS", "TRUG"):
    for habitat in ("BS", "RZ"):
        selected = [
            s for s in samples
            if metadata[s]["Site"] == site
            and metadata[s]["Habitat"] == habitat
        ]
        groups[("Within_site", site, habitat)] = selected

# Within-treatment: compare the same condition across sites
for treatment in ("NAT", "BESE", "BARE", "DEGR"):
    for habitat in ("BS", "RZ"):
        selected = [
            s for s in samples
            if metadata[s]["Treatment"] == treatment
            and metadata[s]["Habitat"] == habitat
        ]
        if selected:
            groups[("Within_treatment", treatment, habitat)] = selected

# Habitat-wide curves are descriptive context only:
# They mix sites and restoration conditions.
for habitat in ("BS", "RZ"):
    selected = [
        s for s in samples
        if metadata[s]["Habitat"] == habitat
    ]
    groups[("Habitat_wide", "ALL", habitat)] = selected

results = []

for (design, group, habitat), selected in sorted(groups.items()):
    n = len(selected)

    if n < 1:
        continue

    # All combinations are feasible for this small dataset.
    # Cap habitat-wide curves at order six.
    max_order = min(n, 6) if design == "Habitat_wide" else n

    for feature_type in ("MAG", "KO"):
        sets = presence[feature_type]

        zeta1 = sum(
            len(sets[s]) for s in selected
        ) / float(n)

        for order in range(1, max_order + 1):
            shared_counts = []

            for combination in itertools.combinations(
                selected, order
            ):
                intersection = set(sets[combination[0]])

                for sample in combination[1:]:
                    intersection.intersection_update(
                        sets[sample]
                    )

                shared_counts.append(len(intersection))

            zeta = sum(shared_counts) / float(
                len(shared_counts)
            )

            results.append({
                "Design": design,
                "Group": group,
                "Habitat": habitat,
                "Feature_type": feature_type,
                "N_metagenomes": n,
                "Zeta_order": order,
                "N_combinations": len(shared_counts),
                "Zeta_mean": round(zeta, 6),
                "Zeta_normalised": round(
                    zeta / zeta1 if zeta1 else 0, 6
                )
            })

output = OUT / "Metagenomic_MAG_KO_zeta_diversity.tsv"

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("SUCCESS: Metagenomic zeta diversity")
print("Metagenomes:", len(samples))
print("MAGs:", len(depth))
print("KOs:", len(ko_cols) - 1)
print("Result rows:", len(results))
print("Output:", output)

print("\n=== WITHIN-SITE ZETA ===")
print("Site Habitat Feature n Zeta1 Zeta2 Zeta3 Zeta4")

for site in ("ACP", "BASS", "TRUG"):
    for habitat in ("BS", "RZ"):
        for feature in ("MAG", "KO"):
            subset = {
                r["Zeta_order"]: r
                for r in results
                if r["Design"] == "Within_site"
                and r["Group"] == site
                and r["Habitat"] == habitat
                and r["Feature_type"] == feature
            }

            vals = [
                "{:.1f}".format(subset[k]["Zeta_mean"])
                if k in subset else "-"
                for k in (1, 2, 3, 4)
            ]

            print(
                "{:<5} {:<3} {:<3} {:>2} {}".format(
                    site, habitat, feature,
                    subset[1]["N_metagenomes"],
                    " ".join("{:>8}".format(v) for v in vals)
                )
            )

print("\n=== WITHIN-TREATMENT ZETA ===")
print("Treatment Habitat Feature n Zeta1 Zeta2 Zeta3")

for treatment in ("NAT", "BESE", "BARE", "DEGR"):
    for habitat in ("BS", "RZ"):
        for feature in ("MAG", "KO"):
            subset = {
                r["Zeta_order"]: r
                for r in results
                if r["Design"] == "Within_treatment"
                and r["Group"] == treatment
                and r["Habitat"] == habitat
                and r["Feature_type"] == feature
            }

            if not subset:
                continue

            vals = [
                "{:.1f}".format(subset[k]["Zeta_mean"])
                if k in subset else "-"
                for k in (1, 2, 3)
            ]

            print(
                "{:<9} {:<3} {:<3} {:>2} {}".format(
                    treatment, habitat, feature,
                    subset[1]["N_metagenomes"],
                    " ".join("{:>8}".format(v) for v in vals)
                )
            )
