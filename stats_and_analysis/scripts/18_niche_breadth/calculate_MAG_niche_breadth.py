import csv
from pathlib import Path

base = Path("mapping/competitive_MAGs/MAG_abundance")
depth_file = base / "final_373_MAGs_mean_depth_matrix.tsv"
breadth_file = base / "final_373_MAGs_breadth_matrix.tsv"

def load(path):
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        samples = reader.fieldnames[1:]
        data = {
            r["MAG_ID"]: {s: float(r[s]) for s in samples}
            for r in reader
        }
    return samples, data

samples, depth = load(depth_file)
breadth_samples, breadth = load(breadth_file)

assert samples == breadth_samples
assert len(depth) == len(breadth) == 373
assert set(depth) == set(breadth)

# Relative MAG abundance within each metagenome
totals = {
    s: sum(depth[m][s] for m in depth)
    for s in samples
}
assert all(v > 0 for v in totals.values())

results = []

for mag in sorted(depth):
    row = {"MAG_ID": mag}

    for threshold, suffix in [(0.10, "10"), (0.50, "50")]:
        abundances = {
            s: depth[mag][s] / totals[s]
            if breadth[mag][s] >= threshold else 0.0
            for s in samples
        }

        occupied = [s for s in samples if abundances[s] > 0]
        total = sum(abundances.values())

        if total > 0:
            p = [v / total for v in abundances.values()]
            levins = 1 / sum(x*x for x in p)
        else:
            levins = 0.0

        row[f"Occupancy_{suffix}"] = len(occupied)
        row[f"Levins_{suffix}"] = round(levins, 6)

        # Normalised to the 20 sampled environments
        row[f"Levins_standardised_{suffix}"] = (
            round((levins - 1) / (len(samples) - 1), 6)
            if levins > 0 else 0.0
        )

        row[f"Sites_{suffix}"] = len({
            s.split("_")[1] for s in occupied
        })

        row[f"Treatments_{suffix}"] = len({
            s.split("_")[2] for s in occupied
        })

    results.append(row)

out = Path("results/ecological_strategies/MAG_niche_breadth.tsv")

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=results[0].keys(), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("SUCCESS")
print("MAGs:", len(results))

for suffix in ["10", "50"]:
    values = sorted(r[f"Levins_{suffix}"] for r in results)
    print(f"\nGenome breadth >= {suffix}%")
    print("Levins median:", round(values[len(values)//2], 2))
    print("Levins range:", round(values[0], 2),
          "to", round(values[-1], 2))
    print("MAGs detected across all 3 sites:",
          sum(r[f"Sites_{suffix}"] == 3 for r in results))
    print("MAGs detected across all 4 treatments:",
          sum(r[f"Treatments_{suffix}"] == 4 for r in results))

print("\nOutput:", out)
