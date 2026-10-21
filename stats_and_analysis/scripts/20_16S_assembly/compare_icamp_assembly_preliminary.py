import csv
from collections import defaultdict
from pathlib import Path

base = Path("betaNTI")

with (base / "metadata_16S.csv").open(
    encoding="utf-8-sig"
) as f:
    metadata = {
        r["SampleID"]: r
        for r in csv.DictReader(f)
    }

with (base / "icamp_scale_103_rand19" /
      "icamp_scale_103_rand19.process.CbMPDiCBraya.csv").open() as f:
    icamp = list(csv.DictReader(f))

with (base / "assembly_process_pairwise_with_metadata.csv").open() as f:
    whole = {
        tuple(sorted((r["Sample1"], r["Sample2"]))): r["Process"]
        for r in csv.DictReader(f)
    }

processes = [
    "Heterogeneous.Selection",
    "Homogeneous.Selection",
    "Dispersal.Limitation",
    "Homogenizing.Dispersal",
    "Drift.and.Others"
]

groups = defaultdict(list)

for r in icamp:
    a, b = r["sample1"], r["sample2"]
    ma, mb = metadata[a], metadata[b]

    if ma["Site"] != mb["Site"]:
        continue
    if ma["SoilPortion"] != mb["SoilPortion"]:
        continue

    t1, t2 = ma["Treatment"], mb["Treatment"]

    if t1 == t2:
        if t1 not in ("BESE", "BARE"):
            continue
        contrast = t1 + "_within"
    else:
        ts = set((t1, t2))
        if ts == {"NATURAL", "BESE"}:
            contrast = "NAT_BESE"
        elif ts == {"NATURAL", "BARE"}:
            contrast = "NAT_BARE"
        elif ts == {"BESE", "BARE"}:
            contrast = "BESE_BARE"
        else:
            continue

    values = [float(r[p]) for p in processes]
    assert abs(sum(values) - 1) < 0.001

    key = tuple(sorted((a, b)))
    assert key in whole

    groups[
        (ma["Site"], ma["SoilPortion"], contrast)
    ].append((values, whole[key]))

results = []

for (site, habitat, contrast), records in sorted(groups.items()):
    n = len(records)

    means = [
        sum(v[0][j] for v in records) / n * 100
        for j in range(5)
    ]

    whole_hom = sum(
        label == "Homogeneous selection"
        for _, label in records
    ) / n * 100

    results.append({
        "Site": site,
        "Habitat": habitat,
        "Contrast": contrast,
        "N_pairs": n,
        "Whole_homogeneous_pct": round(whole_hom, 2),
        "iCAMP_heterogeneous_pct": round(means[0], 2),
        "iCAMP_homogeneous_pct": round(means[1], 2),
        "iCAMP_dispersal_limitation_pct": round(means[2], 2),
        "iCAMP_homogenizing_dispersal_pct": round(means[3], 2),
        "iCAMP_drift_other_pct": round(means[4], 2)
    })

out = base / "iCAMP_19rand_preliminary_treatment_comparison.tsv"

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("=== RHIZOSPHERE: PRELIMINARY iCAMP ===")
print("Site  Contrast       N  Whole-Hom%  iCAMP-Hom%  iCAMP-DispLim%")

for r in results:
    if r["Habitat"] == "RZ":
        print("{:<5} {:<13} {:>2} {:>10.1f} {:>11.1f} {:>15.1f}".format(
            r["Site"], r["Contrast"], r["N_pairs"],
            r["Whole_homogeneous_pct"],
            r["iCAMP_homogeneous_pct"],
            r["iCAMP_dispersal_limitation_pct"]
        ))

print("\nSUCCESS")
print("Site-level summary rows:", len(results))
print("Output:", out)
print("Exploratory only: iCAMP used 19 randomisations.")
