#!/usr/bin/env python3

"""
Build final REDOX FLEXIBILITY ecological-strategy tables.

Framework:
    KO annotations
        ->
    (1) module/pathway evidence
        ->
    (2) biological capability calls (0/1)
        ->
    (3) standardized redox-flexibility breadth (0-1)

Nine equally weighted redox capabilities:

    1. Nitrate respiration
    2. Nitrite -> NO reduction
    3. DNRA / nitrite -> ammonium
    4. Nitric oxide reduction
    5. Nitrous oxide reduction
    6. Dissimilatory sulfate-reduction core
    7. Reduced-sulfur oxidation
    8. Sulfite oxidation
    9. Low-oxygen respiration

Important:
- Individual KOs are NOT counted as independent ecological traits.
- Alternative biochemical routes collapse into one capability.
- Module evidence is retained for audit/QC.
- The strategy score measures breadth of encoded redox capabilities,
  not measured physiological activity.
- Individual MAG scores are NOT corrected for genome completeness.
"""

import csv
import re
from collections import defaultdict, Counter
from pathlib import Path


# ============================================================
# PATHS
# ============================================================

BASE = Path(
    "/path/to/metagenomics_project"
)

INPUT = (
    BASE
    / "results/functional_annotation/kofamscan/"
      "final_373_MAGs_gene_master_table.tsv"
)

TREE = (
    BASE
    / "results/phylogeny_bacteria_352/"
      "bacteria_352_bac120.tree"
)

OUTPUT_DIR = (
    BASE
    / "results/ecological_strategies/"
      "redox_flexibility"
)

TABLE_DIR = OUTPUT_DIR / "tables"
ITOL_DIR = OUTPUT_DIR / "itol"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
ITOL_DIR.mkdir(parents=True, exist_ok=True)

FINAL_TABLE = (
    TABLE_DIR
    / "redox_strategy_final_373_MAGs.tsv"
)

MODULE_TABLE = (
    TABLE_DIR
    / "redox_module_evidence_373_MAGs.tsv"
)

SUMMARY_TABLE = (
    TABLE_DIR
    / "redox_strategy_final_summary.tsv"
)

ITOL_FILE = (
    ITOL_DIR
    / "itol_redox_strategy_score_352.txt"
)


# ============================================================
# HELPERS
# ============================================================

def parse_ko_field(value):
    """
    Extract KEGG K numbers from a KO field.
    """

    if not value:
        return set()

    return set(
        re.findall(r"K\d{5}", value)
    )


def evaluate_all(kos_present, required):
    """
    All markers are required for the module.
    """

    required = set(required)

    detected = len(
        kos_present & required
    )

    total = len(required)

    evidence = (
        detected / total
        if total
        else 0
    )

    complete = int(
        required <= kos_present
    )

    return (
        detected,
        total,
        evidence,
        complete,
    )


def evaluate_single(kos_present, ko):
    """
    Single diagnostic marker.
    """

    present = int(
        ko in kos_present
    )

    return (
        present,
        1,
        float(present),
        present,
    )


def add_module(
    row,
    kos_present,
    prefix,
    required,
    mode="all",
):
    """
    Evaluate one module and add its evidence columns.
    """

    if mode == "all":

        detected, total, evidence, complete = (
            evaluate_all(
                kos_present,
                required,
            )
        )

    elif mode == "single":

        detected, total, evidence, complete = (
            evaluate_single(
                kos_present,
                required[0],
            )
        )

    else:

        raise ValueError(
            f"Unknown mode: {mode}"
        )

    row[
        f"{prefix}_markers_detected"
    ] = detected

    row[
        f"{prefix}_markers_total"
    ] = total

    row[
        f"{prefix}_module_evidence"
    ] = round(
        evidence,
        3,
    )

    row[
        f"{prefix}_complete"
    ] = complete

    return complete


# ============================================================
# READ MASTER TABLE
# ============================================================

mag_kos = defaultdict(set)
mag_metadata = {}

with INPUT.open(
    "r",
    newline=""
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t",
    )

    if reader.fieldnames is None:
        raise ValueError(
            "Input table has no header."
        )

    required_columns = {
        "Final_MAG",
        "KO",
        "Completeness",
        "Contamination",
        "Site",
        "Treatment",
        "Compartment",
        "GTDB_phylum",
        "GTDB_order",
    }

    missing = (
        required_columns
        - set(reader.fieldnames)
    )

    if missing:

        raise ValueError(
            "Missing required columns: "
            + ", ".join(
                sorted(missing)
            )
        )

    for row in reader:

        mag = row[
            "Final_MAG"
        ].strip()

        if not mag:
            continue

        if mag not in mag_metadata:

            mag_metadata[mag] = {
                "Completeness":
                    row["Completeness"],
                "Contamination":
                    row["Contamination"],
                "Site":
                    row["Site"],
                "Treatment":
                    row["Treatment"],
                "Compartment":
                    row["Compartment"],
                "GTDB_phylum":
                    row["GTDB_phylum"],
                "GTDB_order":
                    row["GTDB_order"],
            }

        for ko in parse_ko_field(
            row["KO"]
        ):
            mag_kos[mag].add(
                ko
            )


mags = sorted(
    mag_metadata
)

if len(mags) != 373:

    print(
        "WARNING: expected 373 MAGs, "
        f"but found {len(mags)}."
    )


# ============================================================
# MODULE DEFINITIONS
# ============================================================

#
# 1. NITRATE RESPIRATION
#
# Alternative routes:
#
# NarGHI:
# K00370 NarG
# K00371 NarH
# K00374 NarI
#
# NapAB:
# K02567 NapA
# K02568 NapB
#

NAR_GHI = [
    "K00370",
    "K00371",
    "K00374",
]

NAP_AB = [
    "K02567",
    "K02568",
]


#
# 2. NITRITE -> NO
#
# Alternative catalytic enzymes:
#
# K00368 NirK
# K15864 NirS
#

NIRK = "K00368"
NIRS = "K15864"


#
# 3. DNRA / NITRITE -> AMMONIUM
#
# Alternative systems:
#
# NirBD:
# K00362 NirB
# K00363 NirD
#
# NrfAH:
# K03385 NrfA
# K15876 NrfH
#

NIR_BD = [
    "K00362",
    "K00363",
]

NRF_AH = [
    "K03385",
    "K15876",
]


#
# 4. NITRIC OXIDE REDUCTION
#
# K04561 NorB
#
# NorB is the catalytic membrane component and includes
# qNor-type nitric oxide reductases in KEGG annotations.
#
# K02305 NorC is retained as supporting evidence but is
# not required because qNor systems can function without NorC.
#

NORB = "K04561"
NORC = "K02305"


#
# 5. NITROUS OXIDE REDUCTION
#
# K00376 NosZ
#

NOSZ = "K00376"


#
# 6. DISSIMILATORY SULFATE-REDUCTION CORE
#
# Conservative core:
#
# K00958 Sat
# K00394 AprA
# K00395 AprB
# K11180 DsrA
# K11181 DsrB
#
# This captures the principal Sat-AprAB-DsrAB core.
#
# Directionality can require additional phylogenetic/contextual
# interpretation because Dsr-related systems can also function
# in sulfur oxidation. Therefore all module evidence is retained.
#

SULFATE_REDUCTION_CORE = [
    "K00958",
    "K00394",
    "K00395",
    "K11180",
    "K11181",
]


#
# 7. REDUCED-SULFUR OXIDATION
#
# Alternative routes:
#
# Sox core:
# K17222 SoxA
# K17223 SoxX
# K17224 SoxB
# K17226 SoxY
# K17227 SoxZ
#
# SoxC/D are not required for the core call because KEGG M00595
# treats SoxCD as optional components in the module definition.
#
# Alternative:
# K17218 Sqr, sulfide:quinone oxidoreductase
#

SOX_CORE = [
    "K17222",
    "K17223",
    "K17224",
    "K17226",
    "K17227",
]

SQR = "K17218"


#
# 8. SULFITE OXIDATION
#
# Alternative systems:
#
# SoeABC:
# K21307 SoeA
# K21308 SoeB
# K21309 SoeC
#
# SorAB:
# K05301 SorA
# K00386 SorB
#

SOE_ABC = [
    "K21307",
    "K21308",
    "K21309",
]

SOR_AB = [
    "K05301",
    "K00386",
]


#
# 9. LOW-OXYGEN RESPIRATION
#
# Alternative terminal oxidase systems.
#
# cbb3-type route A:
# K00404 + K00405 + K00406 + K00407
#
# cbb3-type fused route:
# K15862 + K00406 + K00407
#
# Cytochrome bd core:
# K00425 + K00426
#
# Small cytochrome-bd subunits K00424/K22501 are retained
# as supporting evidence but not required in the primary call,
# because small ORFs can be missed in MAG annotations.
#

CBB3_ROUTE_A = [
    "K00404",
    "K00405",
    "K00406",
    "K00407",
]

CBB3_ROUTE_FUSED = [
    "K15862",
    "K00406",
    "K00407",
]

CYD_AB = [
    "K00425",
    "K00426",
]

CYD_SMALL_SUBUNITS = [
    "K00424",
    "K22501",
]


# ============================================================
# SCORE EACH MAG
# ============================================================

results = []

for mag in mags:

    kos_present = mag_kos[
        mag
    ]

    row = {
        "Final_MAG": mag,
        **mag_metadata[mag],
    }


    # ========================================================
    # 1. NITRATE RESPIRATION
    # ========================================================

    nar = add_module(
        row,
        kos_present,
        "Nitrate_NarGHI",
        NAR_GHI,
    )

    nap = add_module(
        row,
        kos_present,
        "Nitrate_NapAB",
        NAP_AB,
    )

    nitrate_routes = []

    if nar:
        nitrate_routes.append(
            "NarGHI"
        )

    if nap:
        nitrate_routes.append(
            "NapAB"
        )

    row[
        "Nitrate_respiration_routes"
    ] = (
        ";".join(nitrate_routes)
        if nitrate_routes
        else "None"
    )

    row[
        "Capability_Nitrate_respiration"
    ] = int(
        nar or nap
    )


    # ========================================================
    # 2. NITRITE -> NO REDUCTION
    # ========================================================

    nirk = add_module(
        row,
        kos_present,
        "Nitrite_NO_NirK",
        [NIRK],
        mode="single",
    )

    nirs = add_module(
        row,
        kos_present,
        "Nitrite_NO_NirS",
        [NIRS],
        mode="single",
    )

    nitrite_no_routes = []

    if nirk:
        nitrite_no_routes.append(
            "NirK"
        )

    if nirs:
        nitrite_no_routes.append(
            "NirS"
        )

    row[
        "Nitrite_NO_routes"
    ] = (
        ";".join(
            nitrite_no_routes
        )
        if nitrite_no_routes
        else "None"
    )

    row[
        "Capability_Nitrite_to_NO"
    ] = int(
        nirk or nirs
    )


    # ========================================================
    # 3. DNRA
    # ========================================================

    nirbd = add_module(
        row,
        kos_present,
        "DNRA_NirBD",
        NIR_BD,
    )

    nrfah = add_module(
        row,
        kos_present,
        "DNRA_NrfAH",
        NRF_AH,
    )

    dnra_routes = []

    if nirbd:
        dnra_routes.append(
            "NirBD"
        )

    if nrfah:
        dnra_routes.append(
            "NrfAH"
        )

    row[
        "DNRA_routes"
    ] = (
        ";".join(
            dnra_routes
        )
        if dnra_routes
        else "None"
    )

    row[
        "Capability_DNRA"
    ] = int(
        nirbd or nrfah
    )


    # ========================================================
    # 4. NITRIC OXIDE REDUCTION
    # ========================================================

    norb = add_module(
        row,
        kos_present,
        "NO_reduction_NorB",
        [NORB],
        mode="single",
    )

    norc = add_module(
        row,
        kos_present,
        "NO_reduction_NorC_support",
        [NORC],
        mode="single",
    )

    row[
        "NO_reduction_NorBC_both"
    ] = int(
        norb and norc
    )

    row[
        "Capability_NO_reduction"
    ] = norb


    # ========================================================
    # 5. NITROUS OXIDE REDUCTION
    # ========================================================

    nosz = add_module(
        row,
        kos_present,
        "N2O_reduction_NosZ",
        [NOSZ],
        mode="single",
    )

    row[
        "Capability_N2O_reduction"
    ] = nosz


    # ========================================================
    # 6. DISSIMILATORY SULFATE-REDUCTION CORE
    # ========================================================

    sulfate_core = add_module(
        row,
        kos_present,
        "Sulfate_reduction_Sat_AprAB_DsrAB",
        SULFATE_REDUCTION_CORE,
    )

    row[
        "Capability_Dissimilatory_sulfate_reduction_core"
    ] = sulfate_core


    # ========================================================
    # 7. REDUCED-SULFUR OXIDATION
    # ========================================================

    sox = add_module(
        row,
        kos_present,
        "Sulfur_oxidation_Sox_core",
        SOX_CORE,
    )

    sqr = add_module(
        row,
        kos_present,
        "Sulfur_oxidation_SQR",
        [SQR],
        mode="single",
    )

    sulfur_oxidation_routes = []

    if sox:
        sulfur_oxidation_routes.append(
            "Sox"
        )

    if sqr:
        sulfur_oxidation_routes.append(
            "SQR"
        )

    row[
        "Reduced_sulfur_oxidation_routes"
    ] = (
        ";".join(
            sulfur_oxidation_routes
        )
        if sulfur_oxidation_routes
        else "None"
    )

    row[
        "Capability_Reduced_sulfur_oxidation"
    ] = int(
        sox or sqr
    )


    # ========================================================
    # 8. SULFITE OXIDATION
    # ========================================================

    soe = add_module(
        row,
        kos_present,
        "Sulfite_oxidation_SoeABC",
        SOE_ABC,
    )

    sor = add_module(
        row,
        kos_present,
        "Sulfite_oxidation_SorAB",
        SOR_AB,
    )

    sulfite_routes = []

    if soe:
        sulfite_routes.append(
            "SoeABC"
        )

    if sor:
        sulfite_routes.append(
            "SorAB"
        )

    row[
        "Sulfite_oxidation_routes"
    ] = (
        ";".join(
            sulfite_routes
        )
        if sulfite_routes
        else "None"
    )

    row[
        "Capability_Sulfite_oxidation"
    ] = int(
        soe or sor
    )


    # ========================================================
    # 9. LOW-OXYGEN RESPIRATION
    # ========================================================

    cbb3_a = add_module(
        row,
        kos_present,
        "LowO2_cbb3_standard",
        CBB3_ROUTE_A,
    )

    cbb3_fused = add_module(
        row,
        kos_present,
        "LowO2_cbb3_fused",
        CBB3_ROUTE_FUSED,
    )

    cydab = add_module(
        row,
        kos_present,
        "LowO2_CydAB",
        CYD_AB,
    )

    #
    # Supporting small cytochrome-bd subunit evidence.
    #

    cyd_small_present = [
        ko
        for ko in CYD_SMALL_SUBUNITS
        if ko in kos_present
    ]

    row[
        "LowO2_Cyd_small_subunits_detected"
    ] = len(
        cyd_small_present
    )

    row[
        "LowO2_Cyd_small_subunit_names"
    ] = (
        ";".join(
            cyd_small_present
        )
        if cyd_small_present
        else "None"
    )

    low_o2_routes = []

    if cbb3_a:
        low_o2_routes.append(
            "cbb3_standard"
        )

    if cbb3_fused:
        low_o2_routes.append(
            "cbb3_fused"
        )

    if cydab:
        low_o2_routes.append(
            "cytochrome_bd"
        )

    row[
        "LowO2_respiration_routes"
    ] = (
        ";".join(
            low_o2_routes
        )
        if low_o2_routes
        else "None"
    )

    row[
        "Capability_LowO2_respiration"
    ] = int(
        cbb3_a
        or cbb3_fused
        or cydab
    )


    # ========================================================
    # FINAL REDOX STRATEGY
    # ========================================================

    capability_columns = [
        "Capability_Nitrate_respiration",
        "Capability_Nitrite_to_NO",
        "Capability_DNRA",
        "Capability_NO_reduction",
        "Capability_N2O_reduction",
        "Capability_Dissimilatory_sulfate_reduction_core",
        "Capability_Reduced_sulfur_oxidation",
        "Capability_Sulfite_oxidation",
        "Capability_LowO2_respiration",
    ]

    capability_count = sum(
        int(row[col])
        for col in capability_columns
    )

    row[
        "Redox_capabilities_detected"
    ] = capability_count

    row[
        "Redox_capabilities_total"
    ] = len(
        capability_columns
    )

    row[
        "Redox_strategy_score"
    ] = round(
        capability_count
        / len(capability_columns),
        3,
    )

    results.append(row)


# ============================================================
# WRITE FINAL MAG TABLE
# ============================================================

if not results:

    raise ValueError(
        "No MAG results generated."
    )

fieldnames = list(
    results[0].keys()
)

with FINAL_TABLE.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(
        results
    )


# ============================================================
# WRITE MODULE-EVIDENCE TABLE
# ============================================================

module_columns = [
    col
    for col in fieldnames
    if (
        col == "Final_MAG"
        or col in {
            "Completeness",
            "Contamination",
            "Site",
            "Treatment",
            "Compartment",
            "GTDB_phylum",
            "GTDB_order",
        }
        or "markers_detected" in col
        or "markers_total" in col
        or "module_evidence" in col
        or col.endswith("_complete")
        or col.endswith("_routes")
        or col.endswith("_both")
        or "small_subunit" in col
    )
]

with MODULE_TABLE.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=module_columns,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()

    for row in results:
        writer.writerow(row)


# ============================================================
# SUMMARY
# ============================================================

capability_summary = [
    (
        "Nitrate_respiration",
        "Capability_Nitrate_respiration",
    ),
    (
        "Nitrite_to_NO",
        "Capability_Nitrite_to_NO",
    ),
    (
        "DNRA",
        "Capability_DNRA",
    ),
    (
        "NO_reduction",
        "Capability_NO_reduction",
    ),
    (
        "N2O_reduction",
        "Capability_N2O_reduction",
    ),
    (
        "Dissimilatory_sulfate_reduction_core",
        "Capability_Dissimilatory_sulfate_reduction_core",
    ),
    (
        "Reduced_sulfur_oxidation",
        "Capability_Reduced_sulfur_oxidation",
    ),
    (
        "Sulfite_oxidation",
        "Capability_Sulfite_oxidation",
    ),
    (
        "LowO2_respiration",
        "Capability_LowO2_respiration",
    ),
]

summary_rows = []

for label, column in capability_summary:

    count = sum(
        int(row[column])
        for row in results
    )

    summary_rows.append({
        "Summary_type":
            "Capability_prevalence",
        "Category":
            label,
        "MAGs":
            count,
        "Percent":
            f"{100 * count / len(results):.2f}",
    })


# ------------------------------------------------------------
# Individual module prevalence
# ------------------------------------------------------------

complete_columns = [
    col
    for col in fieldnames
    if col.endswith(
        "_complete"
    )
]

for column in complete_columns:

    count = sum(
        int(row[column])
        for row in results
    )

    summary_rows.append({
        "Summary_type":
            "Module_prevalence",
        "Category":
            column.replace(
                "_complete",
                "",
            ),
        "MAGs":
            count,
        "Percent":
            f"{100 * count / len(results):.2f}",
    })


# ------------------------------------------------------------
# Strategy breadth distribution
# ------------------------------------------------------------

score_distribution = Counter(
    row[
        "Redox_capabilities_detected"
    ]
    for row in results
)

for count_value in sorted(
    score_distribution
):

    n = score_distribution[
        count_value
    ]

    summary_rows.append({
        "Summary_type":
            "Capability_count_distribution",
        "Category":
            f"{count_value}_of_9",
        "MAGs":
            n,
        "Percent":
            f"{100 * n / len(results):.2f}",
    })


with SUMMARY_TABLE.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=[
            "Summary_type",
            "Category",
            "MAGs",
            "Percent",
        ],
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(
        summary_rows
    )


# ============================================================
# TREE TIP MATCHING
# ============================================================

tree_tips = set()

if TREE.exists():

    newick = (
        TREE
        .read_text()
        .strip()
    )

    matches = re.findall(
        r"(?<=[(,])([^():,]+)(?=:)",
        newick
    )

    for item in matches:

        item = (
            item.strip()
            .strip("'")
            .strip('"')
        )

        if item:
            tree_tips.add(
                item
            )

else:

    print(
        f"WARNING: tree not found: "
        f"{TREE}"
    )


tree_results = [
    row
    for row in results
    if row[
        "Final_MAG"
    ] in tree_tips
]


# ============================================================
# iTOL DATASET
# ============================================================

with ITOL_FILE.open(
    "w"
) as outfile:

    outfile.write(
        "DATASET_GRADIENT\n"
    )

    outfile.write(
        "SEPARATOR TAB\n"
    )

    outfile.write(
        "DATASET_LABEL\t"
        "Redox flexibility\n"
    )

    outfile.write(
        "COLOR\t#000000\n"
    )

    outfile.write(
        "COLOR_MIN\t#FFFFFF\n"
    )

    outfile.write(
        "COLOR_MAX\t#000000\n"
    )

    outfile.write(
        "LEGEND_TITLE\t"
        "Redox capability breadth\n"
    )

    outfile.write(
        "DATA\n"
    )

    for row in tree_results:

        outfile.write(
            f"{row['Final_MAG']}\t"
            f"{row['Redox_strategy_score']}\n"
        )


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print(
    "=" * 72
)

print(
    "FINAL REDOX-FLEXIBILITY STRATEGY ANALYSIS"
)

print(
    "=" * 72
)

print()
print(
    f"MAGs assessed: "
    f"{len(results)}"
)

print(
    f"Bacterial tree MAGs matched: "
    f"{len(tree_results)}"
)

if tree_tips:

    print(
        f"Tree tips detected: "
        f"{len(tree_tips)}"
    )


print()
print(
    "Capability prevalence:"
)

print()

for label, column in capability_summary:

    count = sum(
        int(row[column])
        for row in results
    )

    percent = (
        100
        * count
        / len(results)
    )

    print(
        f"  {label:<42}"
        f"{count:>4} MAGs "
        f"({percent:6.2f}%)"
    )


print()
print(
    "Redox strategy distribution:"
)

for capability_count in sorted(
    score_distribution
):

    n = score_distribution[
        capability_count
    ]

    print(
        f"  {capability_count}/9 : "
        f"{n:>4} MAGs "
        f"({100*n/len(results):6.2f}%)"
    )


print()
print(
    "Selected pathway/module prevalence:"
)

for column in complete_columns:

    count = sum(
        int(row[column])
        for row in results
    )

    if count > 0:

        label = (
            column
            .replace(
                "_complete",
                "",
            )
        )

        print(
            f"  {label:<42}"
            f"{count:>4}"
        )


print()
print(
    "Outputs:"
)

print(
    f"  Final MAG table: "
    f"{FINAL_TABLE}"
)

print(
    f"  Module evidence: "
    f"{MODULE_TABLE}"
)

print(
    f"  Summary:         "
    f"{SUMMARY_TABLE}"
)

print(
    f"  iTOL:            "
    f"{ITOL_FILE}"
)

print()
print(
    "Interpretation:"
)

print(
    "  Redox_strategy_score = "
    "number of detected redox capabilities / 9"
)

print(
    "  Alternative molecular routes contribute "
    "only once to the corresponding ecological capability."
)

print(
    "  The score represents genome-encoded redox breadth, "
    "not measured physiological flexibility."
)

print(
    "  Sulfate-reduction core calls retain pathway evidence "
    "because Dsr-related pathway direction can require "
    "additional biological interpretation."
)

print(
    "=" * 72
)
