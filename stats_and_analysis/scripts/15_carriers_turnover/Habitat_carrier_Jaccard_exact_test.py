#!/usr/bin/env python3
import csv
import itertools
from pathlib import Path
from collections import defaultdict
from statistics import mean, median

base = Path(
    "results/ecological_strategies/"
    "functional_carrier_taxonomy/carrier_turnover"
)

def load(name):
    with (base / name).open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

carrier = load("Functional_carrier_turnover_site_matched.tsv")
community = load("Whole_MAG_community_turnover.tsv")

C = {
    (r["Filter"], r["Site"], r["Habitat"], r["Contrast"]): r
    for r in community
}

F = defaultdict(dict)

for r in carrier:
    key = (
        r["Filter"], r["Site"], r["Contrast"],
        r["GRESP"], r["Capability"]
    )
    F[key][r["Habitat"]] = r["Carrier_Jaccard"]

contrasts = (
    "NAT_vs_BESE",
    "NAT_vs_BARE",
    "BESE_vs_BARE"
)
sites = ("ACP", "BASS", "TRUG")
filters = ("Breadth10", "Breadth50")

capability_results = []
site_results = []

for filt in filters:
    for contrast in contrasts:
        for site in sites:

            bs_key = (filt, site, "BS", contrast)
            rz_key = (filt, site, "RZ", contrast)

            if bs_key not in C or rz_key not in C:
                continue

            deltas = []

            for key, habitats in F.items():
                if key[:3] != (filt, site, contrast):
                    continue

                bs = habitats.get("BS", "")
                rz = habitats.get("RZ", "")

                if bs == "" or rz == "":
                    continue

                delta = float(rz) - float(bs)
                deltas.append(delta)

                capability_results.append({
                    "Filter": filt,
                    "Site": site,
                    "Contrast": contrast,
                    "GRESP": key[3],
                    "Capability": key[4],
                    "BS_Jaccard": float(bs),
                    "RZ_Jaccard": float(rz),
                    "RZ_minus_BS": delta
                })

            assert deltas, "No evaluable capabilities"

            whole_delta = (
                float(C[rz_key]["Community_Jaccard"])
                - float(C[bs_key]["Community_Jaccard"])
            )

            site_results.append({
                "Filter": filt,
                "Site": site,
                "Contrast": contrast,
                "N_capabilities": len(deltas),
                "Mean_carrier_delta_J": mean(deltas),
                "Median_carrier_delta_J": median(deltas),
                "Positive_capabilities": sum(x > 0 for x in deltas),
                "Negative_capabilities": sum(x < 0 for x in deltas),
                "Whole_MAG_delta_J": whole_delta
            })

def exact_p(values):
    observed = abs(mean(values))
    permutations = [
        abs(mean([x * s for x, s in zip(values, signs)]))
        for signs in itertools.product((-1, 1), repeat=len(values))
    ]
    return (
        sum(x >= observed - 1e-12 for x in permutations)
        / len(permutations)
    )

summary = []

for filt in filters:
    for contrast in contrasts:
        rows = [
            r for r in site_results
            if r["Filter"] == filt and r["Contrast"] == contrast
        ]

        carrier_delta = [r["Mean_carrier_delta_J"] for r in rows]
        whole_delta = [r["Whole_MAG_delta_J"] for r in rows]

        summary.append({
            "Filter": filt,
            "Contrast": contrast,
            "N_sites": len(rows),
            "Sites": ",".join(r["Site"] for r in rows),
            "Mean_carrier_delta_J": mean(carrier_delta),
            "Positive_sites": sum(x > 0 for x in carrier_delta),
            "Exact_two_sided_p": exact_p(carrier_delta),
            "Mean_whole_MAG_delta_J": mean(whole_delta),
            "Whole_MAG_exact_p": exact_p(whole_delta)
        })

def save(rows, filename):
    with (base / filename).open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

save(capability_results, "Habitat_carrier_Jaccard_by_capability.tsv")
save(site_results, "Habitat_carrier_Jaccard_site_effects.tsv")
save(summary, "Habitat_carrier_Jaccard_exact_tests.tsv")

print("\n=== PRIMARY: 10% GENOME BREADTH ===")
print("Positive delta = greater Jaccard sharing in rhizosphere\n")

for r in site_results:
    if r["Filter"] == "Breadth10":
        print(
            "{} {:<13} | Carrier dJ={:+.3f} | MAG dJ={:+.3f} "
            "| capabilities={}".format(
                r["Site"], r["Contrast"],
                r["Mean_carrier_delta_J"],
                r["Whole_MAG_delta_J"],
                r["N_capabilities"]
            )
        )

print("\n=== EXACT SITE-LEVEL TESTS ===")

for r in summary:
    print(
        "{} {:<13} | sites={} | carrier dJ={:+.3f} "
        "| p={:.4f} | MAG dJ={:+.3f}".format(
            r["Filter"], r["Contrast"], r["N_sites"],
            r["Mean_carrier_delta_J"],
            r["Exact_two_sided_p"],
            r["Mean_whole_MAG_delta_J"]
        )
    )

print("\nSUCCESS: Habitat carrier-sharing assessment completed")
print("Outputs:", base)
print("Capability measurements are not independent replicates.")
print("Exact tests swap habitat jointly within each site.")
