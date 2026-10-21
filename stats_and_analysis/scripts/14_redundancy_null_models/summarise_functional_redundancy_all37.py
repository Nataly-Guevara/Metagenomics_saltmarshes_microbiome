import csv
from pathlib import Path
from collections import defaultdict
from statistics import mean

BASE = Path(
    "results/ecological_strategies/functional_redundancy"
)

SOURCE = BASE / "functional_redundancy_20samples_37capabilities.tsv"

METRICS = [
    "Carrier_richness",
    "Effective_carriers",
    "Weighted_prevalence_pct"
]

TREATMENTS = {
    "BS": ["NAT", "BESE", "BARE", "DEGR"],
    "RZ": ["NAT", "BESE", "BARE"]
}

CONTRASTS = {
    "BS": [
        ("NAT", "BESE"),
        ("NAT", "BARE"),
        ("NAT", "DEGR"),
        ("BESE", "BARE"),
        ("BESE", "DEGR"),
        ("BARE", "DEGR")
    ],
    "RZ": [
        ("NAT", "BESE"),
        ("NAT", "BARE"),
        ("BESE", "BARE")
    ]
}

# Load all previously calculated redundancy results
with SOURCE.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 1480

data = {}
capabilities = set()
samples = set()

for r in rows:
    site = r["Site"]
    habitat = r["Habitat"]
    treatment = r["Treatment"]
    filter_name = r["Filter"]
    group = r["GRESP"]
    capability = r["Capability"]

    key = (
        site, habitat, treatment,
        filter_name, group, capability
    )

    if key in data:
        raise ValueError(f"Duplicate: {key}")

    data[key] = {
        metric: float(r[metric])
        for metric in METRICS
    }

    capabilities.add((group, capability))
    samples.add(r["Sample"])

assert len(capabilities) == 37
assert len(samples) == 20
assert len(data) == 1480

filters = ["Breadth10", "Breadth50"]
sites = ["ACP", "BASS", "TRUG"]

# -------------------------------------------------
# 1. Treatment-level means across available sites
# -------------------------------------------------

treatment_rows = []

for filt in filters:
    for habitat, treatments in TREATMENTS.items():
        for group, capability in sorted(capabilities):
            for treatment in treatments:

                observations = []

                for site in sites:
                    key = (
                        site, habitat, treatment,
                        filt, group, capability
                    )

                    if key in data:
                        observations.append((site, data[key]))

                if not observations:
                    continue

                record = {
                    "Filter": filt,
                    "Habitat": habitat,
                    "GRESP": group,
                    "Capability": capability,
                    "Treatment": treatment,
                    "N_sites": len(observations),
                    "Sites": ",".join(s for s, _ in observations)
                }

                for metric in METRICS:
                    record["Mean_" + metric] = round(
                        mean(v[metric] for _, v in observations), 6
                    )

                treatment_rows.append(record)

# -------------------------------------------------
# 2. Site-matched treatment contrasts
# -------------------------------------------------

contrast_rows = []

for filt in filters:
    for habitat, contrasts in CONTRASTS.items():
        for group, capability in sorted(capabilities):
            for t1, t2 in contrasts:
                for site in sites:

                    k1 = (
                        site, habitat, t1,
                        filt, group, capability
                    )
                    k2 = (
                        site, habitat, t2,
                        filt, group, capability
                    )

                    # Missing comparisons are not invented
                    if k1 not in data or k2 not in data:
                        continue

                    a = data[k1]
                    b = data[k2]

                    record = {
                        "Filter": filt,
                        "Habitat": habitat,
                        "GRESP": group,
                        "Capability": capability,
                        "Site": site,
                        "Treatment1": t1,
                        "Treatment2": t2,
                        "Contrast": f"{t2}_minus_{t1}"
                    }

                    for metric in METRICS:
                        record["Treatment1_" + metric] = round(
                            a[metric], 6
                        )
                        record["Treatment2_" + metric] = round(
                            b[metric], 6
                        )
                        record["Delta_" + metric] = round(
                            b[metric] - a[metric], 6
                        )

                    contrast_rows.append(record)

# -------------------------------------------------
# 3. Summarise site consistency
# -------------------------------------------------

grouped = defaultdict(list)

for r in contrast_rows:
    key = (
        r["Filter"],
        r["Habitat"],
        r["GRESP"],
        r["Capability"],
        r["Treatment1"],
        r["Treatment2"]
    )
    grouped[key].append(r)

summary_rows = []

for key, observations in sorted(grouped.items()):

    filt, habitat, group, capability, t1, t2 = key

    record = {
        "Filter": filt,
        "Habitat": habitat,
        "GRESP": group,
        "Capability": capability,
        "Contrast": f"{t2}_minus_{t1}",
        "N_sites": len(observations)
    }

    for metric in METRICS:
        field = "Delta_" + metric
        values = [r[field] for r in observations]

        record["Mean_" + field] = round(mean(values), 6)

        positive = sum(v > 0 for v in values)
        negative = sum(v < 0 for v in values)
        zero = sum(v == 0 for v in values)

        record[metric + "_positive_sites"] = positive
        record[metric + "_negative_sites"] = negative
        record[metric + "_zero_sites"] = zero

        for site in sites:
            match = [
                r[field] for r in observations
                if r["Site"] == site
            ]
            record[site + "_" + field] = (
                match[0] if match else ""
            )

    summary_rows.append(record)

# -------------------------------------------------
# 4. Save all results
# -------------------------------------------------

def write_tsv(path, records):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(records[0].keys()),
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(records)

file_means = BASE / "redundancy_all37_treatment_means.tsv"
file_sites = BASE / "redundancy_all37_site_contrasts.tsv"
file_summary = BASE / "redundancy_all37_contrast_summary.tsv"

write_tsv(file_means, treatment_rows)
write_tsv(file_sites, contrast_rows)
write_tsv(file_summary, summary_rows)

assert len(treatment_rows) == 518
assert len(contrast_rows) == 1776
assert len(summary_rows) == 666

print("SUCCESS: Complete functional redundancy analysis")
print("Capabilities:", len(capabilities))
print("Metagenomes:", len(samples))
print("Habitats: BS and RZ")
print("Filters: Breadth10 and Breadth50")
print()
print("Treatment-level rows:", len(treatment_rows))
print("Site-matched contrast rows:", len(contrast_rows))
print("Contrast-summary rows:", len(summary_rows))
print()
print("Outputs:")
print(file_means)
print(file_sites)
print(file_summary)
