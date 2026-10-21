#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict

BASE = Path("results/ecological_strategies/functional_carrier_taxonomy")
AB = Path("mapping/competitive_MAGs/MAG_abundance")
OUT = BASE / "carrier_turnover"
OUT.mkdir(parents=True, exist_ok=True)

def load(path, idcol):
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == len({r[idcol] for r in rows})
    return {r[idcol]: r for r in rows}

depth = load(AB / "final_373_MAGs_mean_depth_matrix.tsv", "MAG_ID")
breadth = load(AB / "final_373_MAGs_breadth_matrix.tsv", "MAG_ID")

with (BASE / "MAG_capability_carrier_taxonomy.tsv").open() as f:
    carrier_rows = list(csv.DictReader(f, delimiter="\t"))

mags = set(depth)
assert len(mags) == 373
assert set(breadth) == mags

capabilities = defaultdict(set)

for r in carrier_rows:
    key = (r["GRESP"], r["Capability"])
    capabilities[key].add(r["MAG_ID"])

# Include capabilities with zero genomic carriers.
# Recover all 37 capability keys from the existing
# functional-redundancy table.
with Path(
    "results/ecological_strategies/functional_redundancy/"
    "functional_redundancy_20samples_37capabilities.tsv"
).open() as f:
    all_capabilities = {
        (r["GRESP"], r["Capability"])
        for r in csv.DictReader(f, delimiter="\t")
    }

assert len(all_capabilities) == 37
assert set(capabilities).issubset(all_capabilities)

samples = [c for c in next(iter(depth.values())) if c != "MAG_ID"]
assert len(samples) == 20

cells = {}
for sample in samples:
    fields = sample.split("_")
    site, treatment, habitat = fields[1:4]
    key = (site, treatment, habitat)
    assert key not in cells
    cells[key] = sample

contrasts = [
    ("NAT", "DEGR"),
    ("NAT", "BESE"),
    ("NAT", "BARE"),
    ("BESE", "BARE"),
    ("BESE", "DEGR"),
    ("BARE", "DEGR")
]

results = []

for threshold_name, threshold in (
    ("Breadth10", 0.10),
    ("Breadth50", 0.50)
):
    # Abundance relative to the complete retained MAG assemblage.
    sample_weights = {}

    for sample in samples:
        retained = {
            mag: float(depth[mag][sample])
            for mag in mags
            if float(depth[mag][sample]) > 0
            and float(breadth[mag][sample]) >= threshold
        }
        total = sum(retained.values())
        assert total > 0
        sample_weights[sample] = {
            mag: value / total
            for mag, value in retained.items()
        }

    for site in ("ACP", "BASS", "TRUG"):
        for habitat in ("BS", "RZ"):
            for treatment_a, treatment_b in contrasts:

                cell_a = (site, treatment_a, habitat)
                cell_b = (site, treatment_b, habitat)

                if cell_a not in cells or cell_b not in cells:
                    continue

                sample_a = cells[cell_a]
                sample_b = cells[cell_b]

                weights_a = sample_weights[sample_a]
                weights_b = sample_weights[sample_b]

                for group, capability in sorted(all_capabilities):
                    carriers = capabilities.get((group, capability), set())

                    a = carriers & set(weights_a)
                    b = carriers & set(weights_b)
                    shared = a & b
                    union = a | b

                    sum_a = sum(weights_a[m] for m in a)
                    sum_b = sum(weights_b[m] for m in b)

                    if union:
                        jaccard = len(shared) / float(len(union))

                        # Conditional carrier Bray-Curtis:
                        # normalise within carriers separately.
                        # Measures distribution, not prevalence.
                        pa = {
                            m: weights_a.get(m, 0) / sum_a
                            if sum_a > 0 else 0
                            for m in union
                        }
                        pb = {
                            m: weights_b.get(m, 0) / sum_b
                            if sum_b > 0 else 0
                            for m in union
                        }

                        if sum_a > 0 and sum_b > 0:
                            bc = sum(
                                abs(pa[m] - pb[m]) for m in union
                            ) / 2
                        else:
                            bc = None
                    else:
                        jaccard = None
                        bc = None

                    shared_a = (
                        100 * sum(weights_a[m] for m in shared) / sum_a
                        if sum_a > 0 else None
                    )
                    shared_b = (
                        100 * sum(weights_b[m] for m in shared) / sum_b
                        if sum_b > 0 else None
                    )

                    results.append({
                        "Filter": threshold_name,
                        "Site": site,
                        "Habitat": habitat,
                        "Treatment_A": treatment_a,
                        "Treatment_B": treatment_b,
                        "Contrast": treatment_a + "_vs_" + treatment_b,
                        "GRESP": group,
                        "Capability": capability,
                        "Carriers_A": len(a),
                        "Carriers_B": len(b),
                        "Shared_carriers": len(shared),
                        "Union_carriers": len(union),
                        "Carrier_Jaccard": jaccard,
                        "Carrier_BrayCurtis": bc,
                        "Prevalence_A_pct": 100 * sum_a,
                        "Prevalence_B_pct": 100 * sum_b,
                        "Shared_abundance_A_pct": shared_a,
                        "Shared_abundance_B_pct": shared_b
                    })

output = OUT / "Functional_carrier_turnover_site_matched.tsv"

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("SUCCESS: Functional-carrier turnover analysis")
print("Capabilities:", len(all_capabilities))
print("Rows:", len(results))
print("Output:", output)

print("\n=== FOCAL RHIZOSPHERE: BESE vs BARE ===")

focal = {
    "Capability_Trehalose",
    "Capability_PolyP_turnover",
    "Iron_complex_uptake",
    "Capability_LowO2_respiration"
}

for capability in sorted(focal):
    subset = [
        r for r in results
        if r["Filter"] == "Breadth10"
        and r["Habitat"] == "RZ"
        and r["Contrast"] == "BESE_vs_BARE"
        and r["Capability"] == capability
    ]

    print("\n", capability)
    for r in subset:
        print(
            "{} | shared={}/{} | Jaccard={} | carrier BC={} "
            "| shared contribution BESE={}%, BARE={}%".format(
                r["Site"],
                r["Shared_carriers"],
                r["Union_carriers"],
                "NA" if r["Carrier_Jaccard"] is None
                else round(r["Carrier_Jaccard"], 3),
                "NA" if r["Carrier_BrayCurtis"] is None
                else round(r["Carrier_BrayCurtis"], 3),
                "NA" if r["Shared_abundance_A_pct"] is None
                else round(r["Shared_abundance_A_pct"], 1),
                "NA" if r["Shared_abundance_B_pct"] is None
                else round(r["Shared_abundance_B_pct"], 1)
            )
        )
