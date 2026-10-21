import csv
import itertools
from pathlib import Path
from collections import defaultdict, Counter

base = Path("betaNTI")

processes = [
    "Homogeneous selection",
    "Heterogeneous selection",
    "Dispersal limitation",
    "Homogenizing dispersal",
    "Undominated"
]

sites = ["ACP", "BASS", "TRUG"]

def read_csv(path):
    with path.open(encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

metadata = read_csv(base / "metadata_16S.csv")
pairrows = read_csv(
    base / "assembly_process_pairwise_with_metadata.csv"
)

meta = {r["SampleID"]: r for r in metadata}

assert len(meta) == 103
assert len(pairrows) == 5253

# Check corrected BASS assignments.
assert meta["9_run2_BASS_BARE_BS"]["SoilPortion"] == "BS"
assert meta["69_run2_BASS_BESE_RZ"]["Treatment"] == "BESE"

# Index the previously calculated assembly classifications.
pair_process = {}

for r in pairrows:
    key = tuple(sorted((r["Sample1"], r["Sample2"])))
    assert key not in pair_process
    assert r["Process"] in processes
    pair_process[key] = r["Process"]

def process_between(sample_a, sample_b):
    key = tuple(sorted((sample_a, sample_b)))
    return pair_process[key]

# Each vector represents the difference in process
# percentages between one BESE and one BARE
# sample when each is compared with Degraded samples.
block_vectors = []
site_results = []

for site in sites:

    def samples(treatment):
        return [
            r for r in metadata
            if r["Site"] == site
            and r["SoilPortion"] == "BS"
            and r["Treatment"] == treatment
        ]

    bese = samples("BESE")
    bare = samples("BARE")
    degraded = samples("DEGRADED")

    assert len(bese) == 5
    assert len(bare) == 5
    assert len(degraded) >= 4

    # BESE and BARE must have matching field blocks.
    bese_by_block = {r["Block"]: r["SampleID"] for r in bese}
    bare_by_block = {r["Block"]: r["SampleID"] for r in bare}

    assert len(bese_by_block) == 5
    assert len(bare_by_block) == 5
    assert set(bese_by_block) == set(bare_by_block)

    deg_ids = [r["SampleID"] for r in degraded]

    counts_bese = Counter()
    counts_bare = Counter()

    site_vectors = []

    for block in sorted(bese_by_block):

        b = bese_by_block[block]
        a = bare_by_block[block]

        cb = Counter(
            process_between(b, deg)
            for deg in deg_ids
        )

        ca = Counter(
            process_between(a, deg)
            for deg in deg_ids
        )

        counts_bese.update(cb)
        counts_bare.update(ca)

        vector = [
            100.0 * (cb[p] - ca[p]) / len(deg_ids)
            for p in processes
        ]

        assert abs(sum(vector)) < 1e-8

        site_vectors.append(vector)
        block_vectors.append(vector)

    expected = 5 * len(deg_ids)

    assert sum(counts_bese.values()) == expected
    assert sum(counts_bare.values()) == expected

    result = {
        "Site": site,
        "Habitat": "BS",
        "N_BESE": len(bese),
        "N_BARE": len(bare),
        "N_Degraded": len(degraded),
        "Matched_blocks": 5,
        "N_pairs_per_comparison": expected
    }

    for p in processes:
        short = p.replace(" ", "_")

        pb = 100.0 * counts_bese[p] / expected
        pa = 100.0 * counts_bare[p] / expected

        result[short + "_BESE_DEGR_pct"] = pb
        result[short + "_BARE_DEGR_pct"] = pa
        result[short + "_difference_pp"] = pb - pa

    site_results.append(result)

assert len(block_vectors) == 15

def mean_vector(vectors, signs):
    n = float(len(vectors))
    return [
        sum(
            signs[i] * vectors[i][k]
            for i in range(len(vectors))
        ) / n
        for k in range(len(processes))
    ]

def omnibus_stat(v):
    return sum(x * x for x in v)

observed = mean_vector(
    block_vectors,
    [1] * len(block_vectors)
)

observed_omnibus = omnibus_stat(observed)
observed_hom = observed[0]

ge_omnibus = 0
ge_hom = 0
nperm = 0

# Exact matched-block label permutations.
for signs in itertools.product(
    (1, -1), repeat=len(block_vectors)
):

    permuted = mean_vector(block_vectors, signs)

    if (
        omnibus_stat(permuted)
        >= observed_omnibus - 1e-10
    ):
        ge_omnibus += 1

    if abs(permuted[0]) >= abs(observed_hom) - 1e-10:
        ge_hom += 1

    nperm += 1

assert nperm == 32768

main_result = {
    "Habitat": "BS",
    "Comparison": "BESE_DEGR_vs_BARE_DEGR",
    "N_sites": 3,
    "Matched_blocks": 15,
    "N_permutations": nperm,
    "Omnibus_statistic": observed_omnibus,
    "Omnibus_exact_p": ge_omnibus / float(nperm),
    "Homogeneous_selection_BESE_minus_BARE_pp": observed_hom,
    "Homogeneous_selection_exact_p": ge_hom / float(nperm)
}

for k, p in enumerate(processes):
    main_result[
        p.replace(" ", "_") + "_difference_pp"
    ] = observed[k]

out1 = base / (
    "Comment03B_Degraded_BS_permutation_20261019.tsv"
)
out2 = base / (
    "Comment03B_Degraded_BS_site_effects_20261019.tsv"
)

def save(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0]),
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

save(out1, [main_result])
save(out2, site_results)

print("\n=== BESE-DEGRADED VS BARE-DEGRADED ===")
print("Habitat: Bulk soil")
print("Sites:", len(sites))
print("Matched field blocks:", len(block_vectors))
print("Exact permutations:", nperm)

print("\n=== SITE-LEVEL HOMOGENEOUS SELECTION ===")
for r in site_results:
    print(
        r["Site"],
        "BESE-DEGR:",
        round(r["Homogeneous_selection_BESE_DEGR_pct"], 2),
        "BARE-DEGR:",
        round(r["Homogeneous_selection_BARE_DEGR_pct"], 2),
        "Difference:",
        round(r["Homogeneous_selection_difference_pp"], 2)
    )

print("\n=== FIVE-PROCESS DIFFERENCES ===")
for p, value in zip(processes, observed):
    print(p, ":", round(value, 3), "pp")

print("\n=== EXACT PERMUTATION TESTS ===")
print(
    "Omnibus p =",
    round(main_result["Omnibus_exact_p"], 6)
)
print(
    "Homogeneous selection difference =",
    round(observed_hom, 3), "pp"
)
print(
    "Homogeneous-selection p =",
    round(main_result["Homogeneous_selection_exact_p"], 6)
)

print("\n=== OUTPUT FILES ===")
print(out1)
print(out2)
print("\nSUCCESS: Degraded-reference test completed")
