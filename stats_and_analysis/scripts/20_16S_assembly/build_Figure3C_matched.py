import csv
import statistics
from pathlib import Path
from collections import defaultdict, Counter

base = Path("betaNTI")
ko_file = Path("results/functional_annotation/figure2_functional_reorganisation/KO_BrayCurtis_distances.tsv")

order = {"NAT": 0, "BESE": 1, "BARE": 2, "DEGR": 3}
aliases = {
    "NAT": "NAT", "NATURAL": "NAT",
    "BESE": "BESE", "BARE": "BARE",
    "DEG": "DEGR", "DEGR": "DEGR", "DEGRADED": "DEGR"
}

def normalise(t):
    return aliases[t.strip().upper()]

def treatment_pair(a, b):
    return tuple(sorted((normalise(a), normalise(b)),
                        key=lambda t: order[t]))

assembly = defaultdict(list)

with (base / "assembly_process_pairwise_with_metadata.csv").open() as f:
    for r in csv.DictReader(f):
        if r["Site1"] != r["Site2"]:
            continue
        if r["Compartment1"] != r["Compartment2"]:
            continue

        pair = treatment_pair(r["Treatment1"], r["Treatment2"])
        if pair[0] == pair[1]:
            continue

        key = (r["Site1"], r["Compartment1"], pair)
        assembly[key].append(r)

ko = {}

with ko_file.open() as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["Site1"] != r["Site2"]:
            continue

        pair = treatment_pair(r["Treatment1"], r["Treatment2"])
        if pair[0] == pair[1]:
            continue

        key = (r["Site1"], r["Habitat"], pair)
        if key in ko:
            raise ValueError(f"Duplicate KO contrast: {key}")

        ko[key] = float(r["BrayCurtis"])

processes = sorted({
    r["Process"]
    for records in assembly.values()
    for r in records
})

rows = []

for key in sorted(ko):
    if key not in assembly:
        raise ValueError(f"Missing assembly contrast: {key}")

    site, habitat, pair = key
    records = assembly[key]
    bnti = [float(r["betaNTI"]) for r in records]
    counts = Counter(r["Process"] for r in records)
    n = len(records)

    row = {
        "Site": site,
        "Habitat": habitat,
        "Contrast": "_vs_".join(pair),
        "N_16S_pairs": n,
        "KO_BrayCurtis": ko[key],
        "Median_betaNTI": statistics.median(bnti),
        "Heterogeneous_selection_pct":
            100 * sum(v > 2 for v in bnti) / n,
        "Homogeneous_selection_pct":
            100 * sum(v < -2 for v in bnti) / n
    }

    for p in processes:
        column = "Process_" + p.replace(" ", "_").replace("-", "_")
        row[column] = 100 * counts[p] / n

    rows.append(row)

assert len(rows) == 24, f"Expected 24 contrasts, found {len(rows)}"

out = base / "Figure3C_matched_assembly_KO_contrasts.tsv"

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=rows[0].keys(), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

print("SUCCESS")
print("Matched contrasts:", len(rows))
print("Bulk soil:", sum(r["Habitat"] == "BS" for r in rows))
print("Rhizosphere:", sum(r["Habitat"] == "RZ" for r in rows))

print("\nSite Habitat Contrast KO_BrayCurtis Median_betaNTI")
for r in rows:
    print(
        r["Site"], r["Habitat"], r["Contrast"],
        f'{r["KO_BrayCurtis"]:.3f}',
        f'{r["Median_betaNTI"]:.2f}'
    )

print("\nOutput:", out)
