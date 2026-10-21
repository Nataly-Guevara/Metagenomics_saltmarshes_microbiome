#!/usr/bin/env python3
import csv
import itertools
from collections import defaultdict
from pathlib import Path
from statistics import mean

BASE = Path(
    "results/ecological_strategies/"
    "functional_carrier_taxonomy/carrier_turnover"
)
INPUT = BASE / "Functional_carrier_turnover_site_matched.tsv"

def write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def shared_contribution(prevalence, conditional_share):
    # Original Shared_abundance_*_pct uses carrier abundance
    # as its denominator, NOT the whole MAG community.
    if conditional_share in ("", "NA", "None"):
        if abs(prevalence) > 1e-8:
            raise ValueError("Missing shared fraction with positive prevalence")
        return 0.0

    fraction = float(conditional_share)
    if not -1e-8 <= fraction <= 100 + 1e-8:
        raise ValueError("Invalid conditional shared percentage")

    return prevalence * fraction / 100.0

with INPUT.open() as f:
    original = list(csv.DictReader(f, delimiter="\t"))

out = []
maximum_error = 0.0
empty_unions = defaultdict(int)
seen = set()

for r in original:
    key = (
        r["Filter"], r["Site"], r["Habitat"],
        r["Contrast"], r["GRESP"], r["Capability"]
    )
    if key in seen:
        raise ValueError(f"Duplicated comparison: {key}")
    seen.add(key)

    pa = float(r["Prevalence_A_pct"])
    pb = float(r["Prevalence_B_pct"])

    sa = shared_contribution(
        pa, r["Shared_abundance_A_pct"]
    )
    sb = shared_contribution(
        pb, r["Shared_abundance_B_pct"]
    )

    ua = pa - sa
    ub = pb - sb

    # B minus A:
    shared_change = sb - sa
    unique_gain = ub
    unique_loss = -ua
    net_change = pb - pa

    reconstructed = (
        shared_change + unique_gain + unique_loss
    )
    error = abs(reconstructed - net_change)
    maximum_error = max(maximum_error, error)

    if error > 1e-7 or min(sa, sb, ua, ub) < -1e-7:
        raise ValueError(f"Decomposition failed: {key}")

    union = int(r["Union_carriers"])
    j = r["Carrier_Jaccard"]

    if union == 0:
        if j not in ("", "NA", "None"):
            raise ValueError("Empty union has defined Jaccard")
        if max(abs(pa), abs(pb)) > 1e-7:
            raise ValueError("Empty union with positive prevalence")
        empty_unions[(r["Filter"], r["Contrast"])] += 1
    elif j in ("", "NA", "None"):
        raise ValueError("Nonempty union has undefined Jaccard")

    row = dict(r)
    row.update({
        "Shared_community_A_pp": sa,
        "Shared_community_B_pp": sb,
        "A_only_community_pp": ua,
        "B_only_community_pp": ub,
        "Shared_carrier_change_pp": shared_change,
        "B_only_gain_pp": unique_gain,
        "A_only_loss_pp": unique_loss,
        "Net_prevalence_change_pp": net_change,
        "Reconstruction_error": error
    })
    out.append(row)

write_tsv(
    BASE / "Carrier_prevalence_decomposition_all.tsv",
    out
)

# Summary within each site first.
by_site = defaultdict(list)
for r in out:
    key = (
        r["Filter"], r["Site"],
        r["Habitat"], r["Contrast"]
    )
    by_site[key].append(r)

site_summary = []
for (filt, site, habitat, contrast), rows in sorted(by_site.items()):
    if len(rows) != 37:
        raise ValueError(
            f"Expected 37 capabilities, got {len(rows)}: "
            f"{(filt, site, habitat, contrast)}"
        )
    site_summary.append({
        "Filter": filt,
        "Site": site,
        "Habitat": habitat,
        "Contrast": contrast,
        "N_capabilities": len(rows),
        "Mean_abs_net_change_pp": mean(
            abs(r["Net_prevalence_change_pp"]) for r in rows
        ),
        "Mean_abs_shared_change_pp": mean(
            abs(r["Shared_carrier_change_pp"]) for r in rows
        ),
        "Mean_B_only_gain_pp": mean(
            r["B_only_gain_pp"] for r in rows
        ),
        "Mean_A_only_loss_magnitude_pp": mean(
            -r["A_only_loss_pp"] for r in rows
        )
    })

write_tsv(
    BASE / "Carrier_prevalence_decomposition_site_summary.tsv",
    site_summary
)

# Focal traits, retaining individual sites (no pseudoreplication).
focal = {
    "Capability_Trehalose",
    "Capability_PolyP_turnover",
    "Iron_complex_uptake",
    "Capability_N2O_reduction"
}
focal_rows = [r for r in out if r["Capability"] in focal]
write_tsv(
    BASE / "Carrier_prevalence_decomposition_focal.tsv",
    focal_rows
)

# Sensitivity: a common set of capabilities must have defined
# carrier Jaccard in both BS and RZ at ALL available matched sites.
lookup = {
    (r["Filter"], r["Contrast"], r["Site"],
     r["Habitat"], r["GRESP"], r["Capability"]): r
    for r in out
}

contrasts = ("NAT_vs_BESE", "NAT_vs_BARE", "BESE_vs_BARE")
sites = ("ACP", "BASS", "TRUG")
filters = ("Breadth10", "Breadth50")

common_results = []

for filt in filters:
    for contrast in contrasts:
        matched_sites = []

        for site in sites:
            has_bs = (filt, site, "BS", contrast) in by_site
            has_rz = (filt, site, "RZ", contrast) in by_site
            if has_bs and has_rz:
                matched_sites.append(site)

        eligible_sets = []

        for site in matched_sites:
            for habitat in ("BS", "RZ"):
                eligible_sets.append({
                    (r["GRESP"], r["Capability"])
                    for r in by_site[(filt, site, habitat, contrast)]
                    if r["Carrier_Jaccard"] not in ("", "NA", "None")
                })

        common = set.intersection(*eligible_sets)
        if not common:
            raise ValueError(f"No common eligible capabilities: {contrast}")

        deltas = []
        for site in matched_sites:
            values = []
            for group, capability in sorted(common):
                bs = lookup[(
                    filt, contrast, site, "BS", group, capability
                )]
                rz = lookup[(
                    filt, contrast, site, "RZ", group, capability
                )]
                values.append(
                    float(rz["Carrier_Jaccard"])
                    - float(bs["Carrier_Jaccard"])
                )
            deltas.append(mean(values))

        # Joint habitat swap within each site.
        observed = abs(mean(deltas))
        permuted = [
            abs(mean([a*s for a, s in zip(deltas, signs)]))
            for signs in itertools.product(
                (-1, 1), repeat=len(deltas)
            )
        ]
        p = sum(
            x >= observed - 1e-12 for x in permuted
        ) / len(permuted)

        common_results.append({
            "Filter": filt,
            "Contrast": contrast,
            "Matched_sites": ",".join(matched_sites),
            "N_sites": len(matched_sites),
            "N_common_eligible_capabilities": len(common),
            "Site_deltas": ",".join(
                f"{site}:{delta:+.4f}"
                for site, delta in zip(matched_sites, deltas)
            ),
            "Mean_RZ_minus_BS_Jaccard": mean(deltas),
            "Exact_site_p": p
        })

write_tsv(
    BASE / "Habitat_carrier_Jaccard_common_eligibility.tsv",
    common_results
)

print("\n=== DECOMPOSITION VALIDATION ===")
print("Input comparisons:", len(original))
print("Maximum reconstruction error:", f"{maximum_error:.3g}")
print("Empty carrier unions are retained as zero prevalence")
print("but excluded from Jaccard means.")

print("\n=== COMMON-ELIGIBILITY SENSITIVITY ===")
for r in common_results:
    print(
        f'{r["Filter"]} {r["Contrast"]} | '
        f'sites={r["N_sites"]} | '
        f'common capabilities={r["N_common_eligible_capabilities"]} | '
        f'delta J={r["Mean_RZ_minus_BS_Jaccard"]:+.3f} | '
        f'p={r["Exact_site_p"]:.3f}'
    )
    print("  ", r["Site_deltas"])

print("\n=== FOCAL BESE MINUS BARE, RHIZOSPHERE, BREADTH10 ===")
for r in focal_rows:
    if (
        r["Filter"] == "Breadth10"
        and r["Habitat"] == "RZ"
        and r["Contrast"] == "BESE_vs_BARE"
    ):
        print(
            f'{r["Site"]} {r["Capability"]} | '
            f'net={r["Net_prevalence_change_pp"]:+.2f} | '
            f'shared={r["Shared_carrier_change_pp"]:+.2f} | '
            f'B-only={r["B_only_gain_pp"]:+.2f} | '
            f'A-only={r["A_only_loss_pp"]:+.2f}'
        )

print("\nSUCCESS")
print("Outputs:", BASE)
print("All contributions are percentage points of reconstructed")
print("MAG-community prevalence; values are genomic potential.")
print("Capabilities are not independent biological replicates.")
