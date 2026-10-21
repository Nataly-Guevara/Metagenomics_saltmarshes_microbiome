#!/usr/bin/env python3
import csv
import itertools
import math
from collections import Counter, defaultdict
from pathlib import Path

base = Path("betaNTI")
out = base / "Comment03_assembly_block_permutation.tsv"
siteout = base / "Comment03_assembly_site_effects.tsv"

with (base / "metadata_16S.csv").open(encoding="utf-8-sig") as f:
    metadata = list(csv.DictReader(f))

with (base / "assembly_process_pairwise_with_metadata.csv").open() as f:
    pairrows = list(csv.DictReader(f))

assert len(metadata) == 103
assert len(pairrows) == 5253

meta = {r["SampleID"]: r for r in metadata}
assert len(meta) == 103

process = {}
for r in pairrows:
    a, b = r["Sample1"], r["Sample2"]
    assert a in meta and b in meta and a != b
    key = tuple(sorted((a, b)))
    assert key not in process
    process[key] = r["Process"]

categories = sorted(set(process.values()))
assert len(categories) == 5

def get_process(a, b):
    return process[tuple(sorted((a, b)))]

def proportions(group_a, group_b=None):
    if group_b is None:
        pairs = itertools.combinations(group_a, 2)
    else:
        pairs = itertools.product(group_a, group_b)

    counts = Counter(get_process(a, b) for a, b in pairs)
    n = sum(counts.values())
    assert n > 0
    return [counts[c] / float(n) for c in categories]

def difference(a, b):
    return tuple(x - y for x, y in zip(a, b))

def analyse_site(site, habitat):
    records = [
        r for r in metadata
        if r["Site"] == site and r["SoilPortion"] == habitat
    ]

    natural = [
        r["SampleID"] for r in records
        if r["Treatment"] == "NATURAL"
    ]

    by_treatment = {}
    for tr in ("BESE", "BARE"):
        groups = defaultdict(list)
        for r in records:
            if r["Treatment"] == tr:
                groups[r["Block"]].append(r["SampleID"])
        by_treatment[tr] = groups

    blocks = sorted(
        set(by_treatment["BESE"]) &
        set(by_treatment["BARE"])
    )

    assert natural, "No Natural samples: {} {}".format(site, habitat)
    assert len(blocks) >= 2, "Insufficient matched blocks"

    matched = []
    for block in blocks:
        a = by_treatment["BESE"][block]
        b = by_treatment["BARE"][block]
        assert len(a) == len(b) == 1
        matched.append((a[0], b[0]))

    n_bits = len(matched)
    outcomes = {"reference": [], "within": []}

    for flags in itertools.product((0, 1), repeat=n_bits):
        bese, bare = [], []

        for (a, b), flip in zip(matched, flags):
            bese.append(b if flip else a)
            bare.append(a if flip else b)

        reference = difference(
            proportions(natural, bese),
            proportions(natural, bare)
        )

        within = difference(
            proportions(bese),
            proportions(bare)
        )

        outcomes["reference"].append(reference)
        outcomes["within"].append(within)

    return {
        "site": site,
        "habitat": habitat,
        "blocks": len(blocks),
        "natural": len(natural),
        "outcomes": outcomes
    }

results = []
site_results = []

hom = categories.index("Homogeneous selection")

for habitat in ("RZ", "BS"):
    sites = [
        analyse_site(site, habitat)
        for site in ("ACP", "BASS", "TRUG")
    ]

    print("\n=== {} ===".format(habitat))
    for s in sites:
        print(
            "{}: {} matched BESE-BARE blocks; {} Natural samples".format(
                s["site"], s["blocks"], s["natural"]
            )
        )

    for comparison in ("reference", "within"):
        arrays = [
            s["outcomes"][comparison]
            for s in sites
        ]

        observed_vectors = [a[0] for a in arrays]

        def average_vectors(vectors):
            return tuple(
                sum(v[j] for v in vectors) / len(vectors)
                for j in range(len(categories))
            )

        observed = average_vectors(observed_vectors)
        observed_omnibus = sum(x*x for x in observed)
        observed_hom = observed[hom]

        ge_omnibus = 0
        ge_hom = 0
        nperm = 0

        for vectors in itertools.product(*arrays):
            mean = average_vectors(vectors)
            statistic = sum(x*x for x in mean)

            if statistic >= observed_omnibus - 1e-12:
                ge_omnibus += 1

            if abs(mean[hom]) >= abs(observed_hom) - 1e-12:
                ge_hom += 1

            nperm += 1

        result = {
            "Habitat": habitat,
            "Comparison": comparison,
            "N_sites": len(sites),
            "Matched_blocks": sum(s["blocks"] for s in sites),
            "N_permutations": nperm,
            "Omnibus_statistic": observed_omnibus,
            "Omnibus_exact_p": ge_omnibus / float(nperm),
            "Homogeneous_selection_BESE_minus_BARE_pp":
                100 * observed_hom,
            "Homogeneous_selection_exact_p":
                ge_hom / float(nperm)
        }
        results.append(result)

        for s in sites:
            vec = s["outcomes"][comparison][0]
            for j, category in enumerate(categories):
                site_results.append({
                    "Site": s["site"],
                    "Habitat": habitat,
                    "Comparison": comparison,
                    "Process": category,
                    "BESE_minus_BARE_pp": 100 * vec[j],
                    "Matched_blocks": s["blocks"]
                })

        print(
            "{} | omnibus p={:.6f} | homogeneous-selection "
            "difference={:+.2f} pp | p={:.6f} | permutations={}".format(
                comparison,
                result["Omnibus_exact_p"],
                result["Homogeneous_selection_BESE_minus_BARE_pp"],
                result["Homogeneous_selection_exact_p"],
                nperm
            )
        )

with out.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(results)

with siteout.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(site_results[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(site_results)

print("\nSUCCESS: Comment 3 treatment comparisons completed")
print("Results:", out)
print("Site effects:", siteout)
