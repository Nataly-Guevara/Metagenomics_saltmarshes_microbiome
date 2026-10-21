#!/usr/bin/env python3
import csv
import itertools
from pathlib import Path
from collections import defaultdict, Counter

SOURCE = Path(
    "results/ecological_strategies/community_weighted_GRESP/"
    "community_capabilities_20_samples_long.tsv"
)
OUT = Path("results/ecological_strategies/metagenomic_zeta")
OUT.mkdir(parents=True, exist_ok=True)

with SOURCE.open() as f:
    data = list(csv.DictReader(f, delimiter="\t"))

features = sorted(set(
    (r["GRESP"], r["Capability"]) for r in data
))
assert len(features) == 37

samples = {}
profiles = defaultdict(dict)

for r in data:
    sample = r["Sample"]
    filt = r["Filter"]
    feature = (r["GRESP"], r["Capability"])

    samples[sample] = (
        r["Site"], r["Treatment"], r["Habitat"]
    )
    key = (filt, sample)
    assert feature not in profiles[key]
    profiles[key][feature] = float(
        r["Weighted_prevalence_pct"]
    )

assert len(samples) == 20

for filt in ("Breadth10", "Breadth50"):
    for sample in samples:
        assert len(profiles[(filt, sample)]) == 37

groups = {}

# Sharing across ecological conditions within each site
for site in ("ACP", "BASS", "TRUG"):
    for habitat in ("BS", "RZ"):
        selected = [
            s for s, m in samples.items()
            if m[0] == site and m[2] == habitat
        ]
        groups[("Within_site", site, habitat)] = selected

# Sharing across sites within each condition
for treatment in ("NAT", "BESE", "BARE", "DEGR"):
    for habitat in ("BS", "RZ"):
        selected = [
            s for s, m in samples.items()
            if m[1] == treatment and m[2] == habitat
        ]
        if selected:
            groups[("Within_treatment", treatment, habitat)] = selected

thresholds = (1.0, 5.0, 10.0)
zeta_rows = []
feature_rows = []

for filt in ("Breadth10", "Breadth50"):
    for cutoff in thresholds:
        present = {
            sample: {
                feature
                for feature, value in profiles[(filt, sample)].items()
                if value >= cutoff
            }
            for sample in samples
        }

        for (design, group, habitat), selected in sorted(groups.items()):
            n = len(selected)
            if n == 0:
                continue

            z1 = sum(len(present[s]) for s in selected) / float(n)

            for order in range(1, n + 1):
                counts = []
                for combo in itertools.combinations(selected, order):
                    shared = present[combo[0]].copy()
                    for s in combo[1:]:
                        shared.intersection_update(present[s])
                    counts.append(len(shared))

                zeta = sum(counts) / float(len(counts))
                zeta_rows.append({
                    "Filter": filt,
                    "Threshold_pct": cutoff,
                    "Design": design,
                    "Group": group,
                    "Habitat": habitat,
                    "N_samples": n,
                    "Zeta_order": order,
                    "N_combinations": len(counts),
                    "Zeta_mean": round(zeta, 5),
                    "Zeta_normalised": (
                        round(zeta / z1, 5) if z1 else ""
                    )
                })

            # Which capabilities occur in every sample?
            for feature in features:
                n_present = sum(
                    feature in present[s] for s in selected
                )
                feature_rows.append({
                    "Filter": filt,
                    "Threshold_pct": cutoff,
                    "Design": design,
                    "Group": group,
                    "Habitat": habitat,
                    "GRESP": feature[0],
                    "Capability": feature[1],
                    "N_samples": n,
                    "Samples_present": n_present,
                    "Core_all_samples": int(n_present == n)
                })

def save(rows, name):
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

save(zeta_rows, "Capability_zeta_37_traits.tsv")
save(feature_rows, "Capability_core_occurrence_37_traits.tsv")

print("SUCCESS: Capability zeta diversity")
print("Metagenomes:", len(samples))
print("Capabilities:", len(features))
print("Zeta results:", len(zeta_rows))
print("Core capability records:", len(feature_rows))

print("\n=== WITHIN-SITE CAPABILITY ZETA ===")
print("Filter Threshold Site Habitat n Zeta1 Zeta2 Zeta3 Zeta4")

for filt in ("Breadth10", "Breadth50"):
    for cutoff in thresholds:
        for site in ("ACP", "BASS", "TRUG"):
            for habitat in ("BS", "RZ"):
                subset = {
                    r["Zeta_order"]: r
                    for r in zeta_rows
                    if r["Filter"] == filt
                    and r["Threshold_pct"] == cutoff
                    and r["Design"] == "Within_site"
                    and r["Group"] == site
                    and r["Habitat"] == habitat
                }
                if not subset:
                    continue

                values = [
                    "{:.1f}".format(subset[k]["Zeta_mean"])
                    if k in subset else "-"
                    for k in (1, 2, 3, 4)
                ]

                print(
                    "{} {:>4}% {:<5} {:<2} {} {}".format(
                        filt, cutoff, site, habitat,
                        subset[1]["N_samples"],
                        " ".join(values)
                    )
                )

print("\n=== WITHIN-TREATMENT CORE CAPABILITIES ===")
print("Primary filter Breadth10, threshold 5%")
for habitat in ("BS", "RZ"):
    for treatment in ("NAT", "BESE", "BARE", "DEGR"):
        subset = [
            r for r in feature_rows
            if r["Filter"] == "Breadth10"
            and r["Threshold_pct"] == 5.0
            and r["Design"] == "Within_treatment"
            and r["Group"] == treatment
            and r["Habitat"] == habitat
        ]
        if subset:
            core = [r for r in subset if r["Core_all_samples"] == 1]
            print(
                "{} {} | sites={} | shared capabilities={}/37".format(
                    habitat, treatment,
                    subset[0]["N_samples"], len(core)
                )
            )

print("\nOutputs:", OUT)
