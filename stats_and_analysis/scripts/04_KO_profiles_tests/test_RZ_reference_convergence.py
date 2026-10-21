import csv
import itertools
import statistics
from pathlib import Path

base = Path("results/functional_annotation/figure2_functional_reorganisation")
source = base / "KO_same_site_distance_to_NAT.tsv"

data = {}

with source.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["Habitat"] != "RZ":
            continue
        key = (row["Site"], row["Treatment"])
        if key in data:
            raise ValueError(f"Duplicate: {key}")
        data[key] = float(row["BrayCurtis_to_NAT"])

differences = []

print("Site\tBESE_to_NAT\tBARE_to_NAT\tBARE_minus_BESE")

for site in ["ACP", "BASS", "TRUG"]:
    bese = data[(site, "BESE")]
    bare = data[(site, "BARE")]
    diff = bare - bese
    differences.append(diff)
    print(f"{site}\t{bese:.4f}\t{bare:.4f}\t{diff:+.4f}")

mean = statistics.mean(differences)
median = statistics.median(differences)

# Exact paired sign-flip permutation test
null = [
    statistics.mean(s * d for s, d in zip(signs, differences))
    for signs in itertools.product([-1, 1], repeat=3)
]

p = sum(abs(x) >= abs(mean) - 1e-12 for x in null) / len(null)

print("\nMean paired difference:", round(mean, 5))
print("Median paired difference:", round(median, 5))
print("Sites where BESE is closer:", sum(d > 0 for d in differences), "/ 3")
print("Exact two-sided permutation p:", p)

output = base / "RZ_BESE_vs_BARE_reference_convergence.tsv"

with output.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["Site", "BARE_minus_BESE"])
    for site, diff in zip(["ACP", "BASS", "TRUG"], differences):
        writer.writerow([site, diff])

print("Saved:", output)
