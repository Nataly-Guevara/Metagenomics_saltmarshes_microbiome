import csv
from pathlib import Path
from collections import Counter

base = Path("results/functional_annotation/figure2_functional_reorganisation")
source = base / "KO_BrayCurtis_distances.tsv"
output = base / "Comment01_matched_habitat_contrasts.tsv"

metadata = {}
distances = {}

with source.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        habitat = r["Habitat"]

        for side in (1, 2):
            sample = r["Sample{}".format(side)]
            info = (
                r["Site{}".format(side)],
                r["Treatment{}".format(side)],
                habitat
            )
            if sample in metadata:
                assert metadata[sample] == info
            metadata[sample] = info

        key = tuple(sorted((r["Sample1"], r["Sample2"])))
        distances[key] = float(r["BrayCurtis"])

assert len(metadata) == 20

cells = {}
for sample, info in metadata.items():
    assert info not in cells
    cells[info] = sample

assert Counter(x[2] for x in cells) == {"BS": 11, "RZ": 9}

sites = ["ACP", "BASS", "TRUG"]
contrasts = [
    ("NAT", "BESE"),
    ("NAT", "BARE"),
    ("BESE", "BARE")
]

print("\n=== METAGENOME DESIGN ===")
for site in sites:
    for habitat in ["BS", "RZ"]:
        present = [
            t for t in ["NAT", "BESE", "BARE", "DEGR"]
            if (site, t, habitat) in cells
        ]
        print(site, habitat, ":", ", ".join(present))

print("\n=== MATCHED HABITAT CONTRASTS ===")
print("Site  Contrast       BS       RZ    RZ-minus-BS")

results = []

for a, b in contrasts:
    for site in sites:
        keys = [
            (site, a, "BS"), (site, b, "BS"),
            (site, a, "RZ"), (site, b, "RZ")
        ]
        if not all(k in cells for k in keys):
            continue

        s1, s2, s3, s4 = [cells[k] for k in keys]

        bs = distances[tuple(sorted((s1, s2)))]
        rz = distances[tuple(sorted((s3, s4)))]

        row = {
            "Site": site,
            "Contrast": "{}_vs_{}".format(a, b),
            "BS_BrayCurtis": bs,
            "RZ_BrayCurtis": rz,
            "RZ_minus_BS": rz - bs
        }
        results.append(row)

        print("{:<5} {:<13} {:8.5f} {:8.5f} {:+11.5f}".format(
            site, row["Contrast"], bs, rz, rz-bs
        ))

assert len(results) == 7

with output.open("w") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

print("\nSUCCESS")
print("Matched BS-RZ comparisons:", len(results))
print("Output:", output)
