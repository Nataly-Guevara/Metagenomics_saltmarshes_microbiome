#!/usr/bin/env python3
import csv
import math
import itertools
import statistics
from collections import defaultdict
from pathlib import Path

base = Path("results/ecological_strategies/functional_redundancy")
source = base / "functional_catalogue_null_999perms_FDR.tsv"
output = base / "Comment06_site_paired_capability_statistics.tsv"

with source.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 666

# Two-sided 95% Student t critical values.
# These intervals are exploratory with only 2-3 sites.
tcrit = {2: 12.706204736, 3: 4.302652730}

results = []

for r in rows:
    vals = []
    for site in ("ACP", "BASS", "TRUG"):
        v = r["{}_delta_Prevalence_pp".format(site)]
        if v not in ("", "NA", "nan"):
            vals.append(float(v))

    n = len(vals)
    assert n == int(r["N_sites"]) and n in (2, 3)

    mean = statistics.mean(vals)
    sd = statistics.stdev(vals)
    se = sd / math.sqrt(n)
    margin = tcrit[n] * se

    # Exact two-sided paired sign-flip test.
    observed = abs(mean)
    possible = [
        abs(statistics.mean(
            [a*b for a,b in zip(vals, signs)]
        ))
        for signs in itertools.product((-1,1), repeat=n)
    ]

    p = sum(x >= observed-1e-10 for x in possible) / len(possible)

    result = dict(r)
    result.update({
        "Paired_mean_difference_pp": mean,
        "Paired_SD_pp": sd,
        "CI95_low_pp": mean-margin,
        "CI95_high_pp": mean+margin,
        "Exact_paired_p": p,
        "Exploratory_status": "Posthoc_capability_screen"
    })
    results.append(result)

# BH correction across 37 capabilities within each
# habitat, contrast and primary/sensitivity threshold.
groups = defaultdict(list)
for r in results:
    groups[(r["Filter"],r["Habitat"],r["Contrast"])].append(r)

for group in groups.values():
    assert len(group) == 37
    ordered = sorted(group,key=lambda r:r["Exact_paired_p"])
    minimum = 1.0
    for rank in range(37,0,-1):
        r = ordered[rank-1]
        minimum = min(minimum,r["Exact_paired_p"]*37/rank)
        r["Paired_FDR_q"] = minimum

with output.open("w",newline="") as f:
    writer = csv.DictWriter(
        f,fieldnames=list(results[0]),delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

focal = {
    "Capability_Trehalose",
    "Capability_PolyP_turnover",
    "Capability_Ammonium_uptake",
    "Iron_complex_uptake",
    "Capability_LowO2_respiration",
    "Capability_Nitrate_respiration"
}

print("\n=== COMMENT 6: FOCAL RHIZOSPHERE CAPABILITIES ===")
print("Differences = BARE minus BESE, percentage points\n")

for threshold in ("Breadth10","Breadth50"):
    print("\n",threshold)
    for r in results:
        if (r["Filter"] == threshold
            and r["Habitat"] == "RZ"
            and r["Contrast"] == "BARE_minus_BESE"
            and r["Capability"] in focal):

            print(
                "{} | mean={:+.2f} | 95% CI=[{:+.2f},{:+.2f}]"
                " | paired p={:.3f} | q={:.3f}".format(
                    r["Capability"],
                    r["Paired_mean_difference_pp"],
                    r["CI95_low_pp"],
                    r["CI95_high_pp"],
                    r["Exact_paired_p"],
                    r["Paired_FDR_q"]
                )
            )

print("\nSUCCESS")
print("Total site-paired contrasts:",len(results))
print("Output:",output)
print("Note: CIs rely on a t approximation with 2-3 sites.")
print("Catalogue-null p-values are separate from paired-test p-values.")
