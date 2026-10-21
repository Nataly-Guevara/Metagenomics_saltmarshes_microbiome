#!/usr/bin/env python3
import csv
from pathlib import Path
from collections import defaultdict
from statistics import mean

ROOT = Path("results/ecological_strategies/functional_carrier_taxonomy/carrier_turnover")
AB = Path("mapping/competitive_MAGs/MAG_abundance")

def read_rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

depth_rows = read_rows(AB / "final_373_MAGs_mean_depth_matrix.tsv")
breadth_rows = read_rows(AB / "final_373_MAGs_breadth_matrix.tsv")
decomp = read_rows(ROOT / "Carrier_prevalence_decomposition_all.tsv")

depth = {r["MAG_ID"]: r for r in depth_rows}
breadth = {r["MAG_ID"]: r for r in breadth_rows}

assert len(depth) == 373
assert set(depth) == set(breadth)

samples = [c for c in depth_rows[0] if c != "MAG_ID"]
assert len(samples) == 20

cells = {}
for sample in samples:
    parts = sample.split("_")
    key = tuple(parts[1:4])
    assert len(key) == 3 and key not in cells
    cells[key] = sample

def save(name, rows):
    with (ROOT / name).open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

grouped = defaultdict(list)
for r in decomp:
    key = (
        r["Filter"], r["Site"],
        r["Habitat"], r["Contrast"]
    )
    grouped[key].append(r)

contrasts = (
    ("NAT", "BESE"),
    ("NAT", "BARE"),
    ("BESE", "BARE")
)

results = []

for filter_name, threshold in (
    ("Breadth10", 0.10),
    ("Breadth50", 0.50)
):
    weights = {}

    for sample in samples:
        retained = {
            mag: float(depth[mag][sample])
            for mag in depth
            if float(depth[mag][sample]) > 0
            and float(breadth[mag][sample]) >= threshold
        }

        total = sum(retained.values())
        assert total > 0

        weights[sample] = {
            mag: value / total
            for mag, value in retained.items()
        }

    for site in ("ACP", "BASS", "TRUG"):
        for habitat in ("BS", "RZ"):
            for a, b in contrasts:

                ka = (site, a, habitat)
                kb = (site, b, habitat)

                if ka not in cells or kb not in cells:
                    continue

                wa = weights[cells[ka]]
                wb = weights[cells[kb]]

                shared = set(wa) & set(wb)

                shared_a = 100 * sum(wa[m] for m in shared)
                shared_b = 100 * sum(wb[m] for m in shared)

                # Both communities contribute equally.
                community_pct = (shared_a + shared_b) / 2

                contrast = a + "_vs_" + b
                key = (filter_name, site, habitat, contrast)

                rows = grouped[key]
                assert len(rows) == 37

                numerator = sum(
                    float(r["Shared_community_A_pp"])
                    + float(r["Shared_community_B_pp"])
                    for r in rows
                )
                denominator = sum(
                    float(r["Prevalence_A_pct"])
                    + float(r["Prevalence_B_pct"])
                    for r in rows
                )

                capability_pct = (
                    100 * numerator / denominator
                    if denominator > 0 else float("nan")
                )

                results.append({
                    "Filter": filter_name,
                    "Site": site,
                    "Habitat": habitat,
                    "Contrast": contrast,
                    "MAGs_A": len(wa),
                    "MAGs_B": len(wb),
                    "Shared_MAGs": len(shared),
                    "Whole_MAG_shared_A_pct": shared_a,
                    "Whole_MAG_shared_B_pct": shared_b,
                    "Whole_MAG_shared_mean_pct": community_pct,
                    "Capability_weighted_shared_pct": capability_pct,
                    "Capability_minus_whole_MAG_pp":
                        capability_pct - community_pct
                })

save("Whole_MAG_vs_capability_shared_abundance.tsv", results)

paired = []

for filter_name in ("Breadth10", "Breadth50"):
    for a, b in contrasts:
        contrast = a + "_vs_" + b

        for site in ("ACP", "BASS", "TRUG"):
            found = {
                r["Habitat"]: r
                for r in results
                if r["Filter"] == filter_name
                and r["Contrast"] == contrast
                and r["Site"] == site
            }

            if not {"BS", "RZ"}.issubset(found):
                continue

            bs, rz = found["BS"], found["RZ"]

            paired.append({
                "Filter": filter_name,
                "Site": site,
                "Contrast": contrast,
                "BS_whole_MAG_pct":
                    bs["Whole_MAG_shared_mean_pct"],
                "RZ_whole_MAG_pct":
                    rz["Whole_MAG_shared_mean_pct"],
                "RZ_minus_BS_whole_MAG_pp":
                    rz["Whole_MAG_shared_mean_pct"]
                    - bs["Whole_MAG_shared_mean_pct"],
                "BS_capability_pct":
                    bs["Capability_weighted_shared_pct"],
                "RZ_capability_pct":
                    rz["Capability_weighted_shared_pct"],
                "RZ_minus_BS_capability_pp":
                    rz["Capability_weighted_shared_pct"]
                    - bs["Capability_weighted_shared_pct"],
                "BS_capability_minus_whole_pp":
                    bs["Capability_minus_whole_MAG_pp"],
                "RZ_capability_minus_whole_pp":
                    rz["Capability_minus_whole_MAG_pp"]
            })

save("Whole_MAG_vs_capability_paired_habitats.tsv", paired)

print("\n=== SITE-MATCHED COMPARISON: BREADTH10 ===")
print("Percent abundance from MAGs detected in both communities\n")

for r in paired:
    if r["Filter"] != "Breadth10":
        continue

    print(
        f'{r["Site"]} {r["Contrast"]}: '
        f'whole MAG BS={r["BS_whole_MAG_pct"]:.1f}%, '
        f'RZ={r["RZ_whole_MAG_pct"]:.1f}% | '
        f'capability BS={r["BS_capability_pct"]:.1f}%, '
        f'RZ={r["RZ_capability_pct"]:.1f}%'
    )

print("\n=== MEAN PAIRED HABITAT DIFFERENCES ===")

for filter_name in ("Breadth10", "Breadth50"):
    for a, b in contrasts:
        contrast = a + "_vs_" + b

        rows = [
            r for r in paired
            if r["Filter"] == filter_name
            and r["Contrast"] == contrast
        ]

        if not rows:
            continue

        whole = mean(
            r["RZ_minus_BS_whole_MAG_pp"] for r in rows
        )
        capability = mean(
            r["RZ_minus_BS_capability_pp"] for r in rows
        )

        print(
            f"{filter_name} {contrast} | sites={len(rows)} | "
            f"whole MAG delta={whole:+.1f} pp | "
            f"capability delta={capability:+.1f} pp | "
            f"difference={capability-whole:+.1f} pp"
        )

print("\nSUCCESS")
print("Outputs:", ROOT)
print("Shared abundance uses retained MAG depth.")
print("Capability-weighted abundance counts MAGs across traits.")
print("Descriptive comparisons, not independent capability tests.")
