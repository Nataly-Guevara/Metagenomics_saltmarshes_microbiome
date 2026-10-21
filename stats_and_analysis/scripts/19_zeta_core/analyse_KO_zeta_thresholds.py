#!/usr/bin/env python3
import csv
import itertools
from pathlib import Path

INPUT = Path(
    "results/functional_annotation/sample_KO_profiles/"
    "sample_KO_relative_MAG_abundance_weighted.tsv"
)

OUT = Path("results/ecological_strategies/metagenomic_zeta")
OUT.mkdir(parents=True, exist_ok=True)

with INPUT.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    features = reader.fieldnames[1:]

assert len(rows) == 20
assert len(features) == 6679

samples = {}

for r in rows:
    sample = r["Sample"]
    parts = sample.split("_")

    values = {
        ko: float(r[ko])
        for ko in features
    }

    total = sum(values.values())
    assert total > 0

    samples[sample] = {
        "Site": parts[1],
        "Treatment": parts[2],
        "Habitat": parts[3],
        "Abundance": {
            ko: 100 * value / total
            for ko, value in values.items()
        }
    }

thresholds = [0.001, 0.01, 0.1]

results = []
detection = []

for threshold in thresholds:

    present = {}

    for sample, info in samples.items():
        present[sample] = {
            ko for ko, value in info["Abundance"].items()
            if value >= threshold
        }

        detection.append({
            "Sample": sample,
            "Site": info["Site"],
            "Treatment": info["Treatment"],
            "Habitat": info["Habitat"],
            "Threshold_pct": threshold,
            "KOs_detected": len(present[sample])
        })

    groups = {}

    for site in ("ACP", "BASS", "TRUG"):
        for habitat in ("BS", "RZ"):
            selected = [
                s for s, info in samples.items()
                if info["Site"] == site
                and info["Habitat"] == habitat
            ]
            groups[("Within_site", site, habitat)] = selected

    for treatment in ("NAT", "BESE", "BARE", "DEGR"):
        for habitat in ("BS", "RZ"):
            selected = [
                s for s, info in samples.items()
                if info["Treatment"] == treatment
                and info["Habitat"] == habitat
            ]
            if selected:
                groups[("Within_treatment", treatment, habitat)] = selected

    for (design, group, habitat), selected in sorted(groups.items()):
        if not selected:
            continue

        zeta1 = sum(len(present[s]) for s in selected) / len(selected)

        for order in range(1, len(selected) + 1):
            counts = []

            for combo in itertools.combinations(selected, order):
                shared = set(present[combo[0]])
                for sample in combo[1:]:
                    shared.intersection_update(present[sample])
                counts.append(len(shared))

            zeta = sum(counts) / len(counts)

            results.append({
                "Threshold_pct": threshold,
                "Design": design,
                "Group": group,
                "Habitat": habitat,
                "N_metagenomes": len(selected),
                "Zeta_order": order,
                "Zeta_mean": zeta,
                "Zeta_normalised": zeta/zeta1 if zeta1 else "",
                "N_combinations": len(counts)
            })

def save(rows, name):
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

save(results, "KO_zeta_abundance_threshold_sensitivity.tsv")
save(detection, "KO_detection_threshold_sensitivity.tsv")

print("SUCCESS: KO zeta sensitivity completed")

print("\n=== KO DETECTION ===")
for threshold in thresholds:
    x = [
        r["KOs_detected"] for r in detection
        if r["Threshold_pct"] == threshold
    ]
    print(
        "Threshold {}% | mean {:.1f} KOs | range {}-{}".format(
            threshold, sum(x)/len(x), min(x), max(x)
        )
    )

print("\n=== WITHIN-SITE KO ZETA ===")
print("Threshold Site Habitat n Zeta1 Zeta2 Zeta3 Zeta4")

for threshold in thresholds:
    for site in ("ACP", "BASS", "TRUG"):
        for habitat in ("BS", "RZ"):
            subset = {
                r["Zeta_order"]: r
                for r in results
                if r["Threshold_pct"] == threshold
                and r["Design"] == "Within_site"
                and r["Group"] == site
                and r["Habitat"] == habitat
            }

            if not subset:
                continue

            vals = [
                "{:.1f}".format(subset[k]["Zeta_mean"])
                if k in subset else "-"
                for k in (1,2,3,4)
            ]

            print(
                "{} {:<5} {:<3} {} {}".format(
                    threshold, site, habitat,
                    subset[1]["N_metagenomes"],
                    " ".join("{:>8}".format(v) for v in vals)
                )
            )

print("\nOutputs:", OUT)
