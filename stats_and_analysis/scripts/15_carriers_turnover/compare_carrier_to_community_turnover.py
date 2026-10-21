#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict

AB = Path("mapping/competitive_MAGs/MAG_abundance")
BASE = Path(
    "results/ecological_strategies/functional_carrier_taxonomy/carrier_turnover"
)

def load(path):
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {r["MAG_ID"]: r for r in rows}

depth = load(AB / "final_373_MAGs_mean_depth_matrix.tsv")
breadth = load(AB / "final_373_MAGs_breadth_matrix.tsv")

assert len(depth) == len(breadth) == 373
assert set(depth) == set(breadth)

mags = sorted(depth)
samples = [x for x in next(iter(depth.values())) if x != "MAG_ID"]

cells = {}
for sample in samples:
    parts = sample.split("_")
    key = tuple(parts[1:4])
    assert key not in cells
    cells[key] = sample

source = BASE / "Functional_carrier_turnover_site_matched.tsv"

with source.open() as f:
    carrier_rows = list(csv.DictReader(f, delimiter="\t"))

assert len(carrier_rows) == 1776

# Use only unique comparisons rather than repeating
# calculations for each of the 37 capabilities.
comparisons = sorted({
    (
        r["Filter"], r["Site"], r["Habitat"],
        r["Treatment_A"], r["Treatment_B"]
    )
    for r in carrier_rows
})

baseline = {}
records = []

for filt, site, habitat, tr_a, tr_b in comparisons:

    threshold = 0.10 if filt == "Breadth10" else 0.50

    sa = cells[(site, tr_a, habitat)]
    sb = cells[(site, tr_b, habitat)]

    wa = {
        m: float(depth[m][sa])
        for m in mags
        if float(depth[m][sa]) > 0
        and float(breadth[m][sa]) >= threshold
    }

    wb = {
        m: float(depth[m][sb])
        for m in mags
        if float(depth[m][sb]) > 0
        and float(breadth[m][sb]) >= threshold
    }

    total_a = sum(wa.values())
    total_b = sum(wb.values())

    assert total_a > 0 and total_b > 0

    wa = {m: v / total_a for m, v in wa.items()}
    wb = {m: v / total_b for m, v in wb.items()}

    union = set(wa) | set(wb)
    shared = set(wa) & set(wb)

    jaccard = len(shared) / float(len(union))

    bc = 0.5 * sum(
        abs(wa.get(m, 0) - wb.get(m, 0))
        for m in union
    )

    key = (filt, site, habitat, tr_a, tr_b)
    baseline[key] = (jaccard, bc)

    records.append({
        "Filter": filt,
        "Site": site,
        "Habitat": habitat,
        "Contrast": tr_a + "_vs_" + tr_b,
        "MAGs_A": len(wa),
        "MAGs_B": len(wb),
        "Shared_MAGs": len(shared),
        "Community_Jaccard": jaccard,
        "Community_BrayCurtis": bc
    })

# Compare each functional-carrier assemblage against the
# entire recovered MAG assemblage for the same samples.
comparative = []

for r in carrier_rows:
    key = (
        r["Filter"], r["Site"], r["Habitat"],
        r["Treatment_A"], r["Treatment_B"]
    )

    community_j, community_bc = baseline[key]

    carrier_j = (
        float(r["Carrier_Jaccard"])
        if r["Carrier_Jaccard"] else None
    )
    carrier_bc = (
        float(r["Carrier_BrayCurtis"])
        if r["Carrier_BrayCurtis"] else None
    )

    result = dict(r)

    result["Community_Jaccard"] = community_j
    result["Community_BrayCurtis"] = community_bc

    result["Carrier_minus_community_Jaccard"] = (
        carrier_j - community_j
        if carrier_j is not None else ""
    )

    result["Carrier_minus_community_BrayCurtis"] = (
        carrier_bc - community_bc
        if carrier_bc is not None else ""
    )

    comparative.append(result)

def save(rows, filename):
    with (BASE / filename).open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

save(records, "Whole_MAG_community_turnover.tsv")
save(comparative, "Carrier_vs_community_turnover.tsv")

print("SUCCESS: Carrier turnover compared with whole MAG communities")
print("Unique site-matched community comparisons:", len(records))
print("Capability comparisons:", len(comparative))

focal = {
    "Capability_Trehalose",
    "Capability_PolyP_turnover",
    "Iron_complex_uptake",
    "Capability_LowO2_respiration"
}

print("\n=== FOCAL CAPABILITIES: CARRIER VS COMMUNITY ===")

for habitat in ("BS", "RZ"):
    print("\n", habitat)

    for contrast in ("NAT_vs_DEGR", "NAT_vs_BESE",
                     "NAT_vs_BARE", "BESE_vs_BARE"):
        subset = [
            r for r in comparative
            if r["Filter"] == "Breadth10"
            and r["Habitat"] == habitat
            and r["Contrast"] == contrast
            and r["Capability"] in focal
            and r["Carrier_BrayCurtis"]
        ]

        if not subset:
            continue

        community = sum(
            float(r["Community_BrayCurtis"]) for r in subset
        ) / len(subset)

        carrier = sum(
            float(r["Carrier_BrayCurtis"]) for r in subset
        ) / len(subset)

        print(
            "{} | community BC={:.3f} | carrier BC={:.3f} "
            "| carrier-community={:+.3f}".format(
                contrast, community, carrier,
                carrier - community
            )
        )

print("\nOutputs:", BASE)
