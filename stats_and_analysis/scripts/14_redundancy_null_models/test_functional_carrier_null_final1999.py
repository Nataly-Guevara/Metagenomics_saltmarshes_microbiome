import csv
import random
import statistics
from pathlib import Path
from collections import Counter

rng = random.Random(20261012)
N_PERM = 1999

BASE = Path("results/ecological_strategies/functional_redundancy")
AB = Path("mapping/competitive_MAGs/MAG_abundance")

def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

depth_rows = read_tsv(
    AB / "final_373_MAGs_mean_depth_matrix.tsv"
)
breadth_rows = read_tsv(
    AB / "final_373_MAGs_breadth_matrix.tsv"
)
redundancy = read_tsv(
    BASE / "functional_redundancy_20samples_37capabilities.tsv"
)

depth = {r["MAG_ID"]: r for r in depth_rows}
breadth = {r["MAG_ID"]: r for r in breadth_rows}

assert len(depth) == len(breadth) == 373
assert set(depth) == set(breadth)
assert len(redundancy) == 1480

# Reuse the exact MAG abundances and detection rules.
communities = {}

for r in redundancy:
    key = (r["Sample"], r["Filter"])

    if key in communities:
        continue

    sample, filt = key
    threshold = {
        "Breadth10": 0.10,
        "Breadth50": 0.50
    }[filt]

    weights = [
        float(depth[m][sample])
        for m in depth
        if float(depth[m][sample]) > 0
        and float(breadth[m][sample]) >= threshold
    ]

    assert len(weights) > 0
    communities[key] = weights

assert len(communities) == 40

def effective_carriers(w):
    total = sum(w)
    return total**2 / sum(x*x for x in w) if total else 0.0

def null_summary(values, observed):
    mu = statistics.mean(values)
    sd = statistics.pstdev(values)

    lower = (1 + sum(v <= observed for v in values)) / (len(values) + 1)
    upper = (1 + sum(v >= observed for v in values)) / (len(values) + 1)

    return {
        "mean": mu,
        "sd": sd,
        "ses": (observed-mu)/sd if sd > 1e-12 else None,
        "p_lower": lower,
        "p_upper": upper,
        "p_two_sided": min(1.0, 2*min(lower, upper))
    }

results = []

for r in redundancy:

    sample = r["Sample"]
    filt = r["Filter"]

    weights = communities[(sample, filt)]
    n = len(weights)
    k = int(r["Carrier_richness"])

    observed_prevalence = float(r["Weighted_prevalence_pct"])
    observed_effective = float(r["Effective_carriers"])

    assert 0 <= k <= n

    record = dict(r)
    record["Detected_MAGs"] = n
    record["Null_randomisations"] = N_PERM

    if k == 0:
        status = "No_carriers"
    elif k == n:
        status = "All_MAGs"
    elif k == 1:
        status = "Single_carrier"
    else:
        status = "Evaluable"

    record["Status"] = status

    # Carrier count is held fixed; each eligible MAG
    # has equal probability of selection.
    prevalence_null = []
    effective_null = []

    if 0 < k < n:

        community_total = sum(weights)

        for _ in range(N_PERM):
            chosen = rng.sample(weights, k)

            prevalence_null.append(
                100 * sum(chosen) / community_total
            )

            if k > 1:
                effective_null.append(
                    effective_carriers(chosen)
                )

    for metric, observed, null_values in [
        ("Prevalence", observed_prevalence, prevalence_null),
        ("Effective_carriers", observed_effective, effective_null)
    ]:

        record["Observed_" + metric] = observed

        if null_values:
            summary = null_summary(null_values, observed)

            record["Null_mean_" + metric] = summary["mean"]
            record["Null_sd_" + metric] = summary["sd"]
            record["Difference_" + metric] = (
                observed - summary["mean"]
            )
            record["SES_" + metric] = summary["ses"]
            record["P_lower_" + metric] = summary["p_lower"]
            record["P_upper_" + metric] = summary["p_upper"]
            record["P_two_sided_" + metric] = summary["p_two_sided"]

        else:
            for prefix in [
                "Null_mean_", "Null_sd_", "Difference_",
                "SES_", "P_lower_", "P_upper_", "P_two_sided_"
            ]:
                record[prefix + metric] = None

    results.append(record)

out = BASE / "functional_carrier_null_v2_final1999.tsv"

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=results[0].keys(), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("SUCCESS: Final functional carrier null model")
print("Metagenomes:", len(communities)//2)
print("Capability observations:", len(results))
print("Randomisations:", N_PERM)

print("\nObservation categories:")
for name, count in sorted(
    Counter(r["Status"] for r in results).items()
):
    print(f"{name}: {count}")

print("\nEvaluable observations:")
print("Prevalence:",
      sum(r["SES_Prevalence"] is not None for r in results))
print("Carrier diversity:",
      sum(r["SES_Effective_carriers"] is not None for r in results))

print("\nOutput:", out)
