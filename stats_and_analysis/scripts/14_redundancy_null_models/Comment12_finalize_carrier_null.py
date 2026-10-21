import csv
import math
from collections import defaultdict
from pathlib import Path

base = Path("results/ecological_strategies/functional_redundancy")
source = base / "functional_carrier_null_v2_final1999.tsv"
output = base / "Comment12_carrier_null_final_FDR.tsv"
summary_out = base / "Comment12_carrier_null_summary.tsv"

with source.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 1480

def number(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None

metrics = ["Prevalence", "Effective_carriers"]

for metric in metrics:
    qcol = "FDR_q_" + metric
    flag = "Null_testable_" + metric

    for r in rows:
        r[qcol] = ""
        r[flag] = "No"

    for filt in ("Breadth10", "Breadth50"):
        valid = []

        for r in rows:
            if r["Filter"] != filt:
                continue

            k = int(r["Carrier_richness"])
            p = number(r["P_two_sided_" + metric])
            sd = number(r["Null_sd_" + metric])

            minimum_k = 1 if metric == "Prevalence" else 2

            if (k >= minimum_k and p is not None
                    and sd is not None and sd > 1e-10):
                valid.append((p, r))
                r[flag] = "Yes"

        valid.sort(key=lambda item: item[0])
        m = len(valid)
        q = 1.0

        for rank in range(m, 0, -1):
            p, r = valid[rank - 1]
            q = min(q, p * m / rank)
            r[qcol] = "{:.6f}".format(q)

        sig = sum(float(r[qcol]) < 0.05 for _, r in valid)

        print("{} {}: {} testable, {} FDR q<0.05".format(
            filt, metric, m, sig
        ))

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

groups = defaultdict(list)
for r in rows:
    groups[(r["Filter"], r["Habitat"],
            r["Treatment"])].append(r)

summary = []
for (filt, habitat, treatment), group in sorted(groups.items()):
    for metric in metrics:
        valid = [
            r for r in group
            if r["Null_testable_" + metric] == "Yes"
        ]
        significant = [
            r for r in valid
            if float(r["FDR_q_" + metric]) < 0.05
        ]

        positive = sum(
            float(r["Difference_" + metric]) > 0
            for r in significant
        )

        summary.append({
            "Filter": filt,
            "Habitat": habitat,
            "Treatment": treatment,
            "Metric": metric,
            "N_testable": len(valid),
            "N_FDR_significant": len(significant),
            "N_above_null": positive,
            "N_below_null": len(significant) - positive
        })

with summary_out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(summary[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(summary)

print("\n=== PRIMARY THRESHOLD: TREATMENT SUMMARY ===")
print("Habitat Treatment Metric                Testable Significant Above Below")

for r in summary:
    if r["Filter"] == "Breadth10":
        print("{:<7} {:<9} {:<21} {:>5} {:>11} {:>5} {:>5}".format(
            r["Habitat"], r["Treatment"], r["Metric"],
            r["N_testable"], r["N_FDR_significant"],
            r["N_above_null"], r["N_below_null"]
        ))

print("\nSUCCESS")
print("Full results:", output)
print("Summary:", summary_out)
print("These are sample-specific null tests, not treatment-effect tests.")
