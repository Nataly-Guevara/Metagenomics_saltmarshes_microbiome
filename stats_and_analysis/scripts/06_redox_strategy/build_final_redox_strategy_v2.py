#!/usr/bin/env python3

"""
FINAL REDOX FLEXIBILITY STRATEGY - VERSION 2

Framework:
    KOs
      ->
    (1) module/pathway evidence
      ->
    (2) binary ecological capabilities
      ->
    (3) standardized redox-strategy breadth

Ten equally weighted capabilities:

 1. Standard aerobic respiration
 2. Low-oxygen respiration
 3. Nitrate respiration
 4. Nitrite -> NO reduction
 5. DNRA / nitrite -> ammonium
 6. Nitric oxide reduction
 7. Nitrous oxide reduction
 8. Dissimilatory sulfate-reduction core
 9. Reduced-sulfur oxidation
10. Sulfite oxidation

Redox_strategy_score =
    number of detected capabilities / 10

Important:
- Individual KOs are NOT counted as independent strategies.
- Alternative molecular routes collapse into one capability.
- Fractional module evidence is retained for QC.
- The main ecological score uses binary capabilities.
- Genome completeness is retained but does not inflate
  individual MAG scores.
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

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

ITOL_DIR.mkdir(
    parents=True,
    exist_ok=True
)

FINAL_TABLE = (
    TABLE_DIR
    / "redox_strategy_final_v2_373_MAGs.tsv"
)

MODULE_TABLE = (
    TABLE_DIR
    / "redox_module_evidence_v2_373_MAGs.tsv"
)

SUMMARY_TABLE = (
    TABLE_DIR
    / "redox_strategy_final_v2_summary.tsv"
)

ITOL_FILE = (
    ITOL_DIR
    / "itol_redox_strategy_score_v2_352.txt"
)


# ============================================================
# BASIC FUNCTIONS
# ============================================================

def parse_ko_field(value):

    if not value:
        return set()

    return set(
        re.findall(
            r"K\d{5}",
            value
        )
    )


def fraction_present(
    kos_present,
    required
):

    required = set(required)

    if not required:
        return 0.0

    return (
        len(
            kos_present
            & required
        )
        / len(required)
    )


def all_present(
    kos_present,
    required
):

    return int(
        set(required)
        <= kos_present
    )


def any_present(
    kos_present,
    alternatives
):

    return int(
        bool(
            kos_present
            & set(alternatives)
        )
    )


def record_simple_module(
    row,
    kos_present,
    name,
    required
):

    detected = len(
        kos_present
        & set(required)
    )

    total = len(
        required
    )

    evidence = (
        detected / total
        if total
        else 0
    )

    complete = int(
        set(required)
        <= kos_present
    )

    row[
        f"{name}_markers_detected"
    ] = detected

    row[
        f"{name}_markers_total"
    ] = total

    row[
        f"{name}_module_evidence"
    ] = round(
        evidence,
        3
    )

    row[
        f"{name}_complete"
    ] = complete

    return complete


def record_single(
    row,
    kos_present,
    name,
    ko
):

    present = int(
        ko in kos_present
    )

    row[
        f"{name}_markers_detected"
    ] = present

    row[
        f"{name}_markers_total"
    ] = 1

    row[
        f"{name}_module_evidence"
    ] = float(present)

    row[
        f"{name}_complete"
    ] = present

    return present


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
        delimiter="\t"
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


# ------------------------------------------------------------
# 1. STANDARD AEROBIC RESPIRATION
# ------------------------------------------------------------
#
# KEGG M00155
#
# K02275 +
#     (
#       K02274 + K02276
#       OR
#       K15408
#     )
#
# K02277 / K27109 are optional subunit-IV alternatives
# and are retained as supporting evidence only.
#
# Alternative aerobic terminal oxidase:
#
# KEGG M00417
# K02297 + K02298 + K02299 + K02300
#

COX_ROUTE_STANDARD = [
    "K02275",
    "K02274",
    "K02276",
]

COX_ROUTE_FUSED = [
    "K02275",
    "K15408",
]

COX_OPTIONAL = [
    "K02277",
    "K27109",
]

CYTOCHROME_O = [
    "K02297",
    "K02298",
    "K02299",
    "K02300",
]


# ------------------------------------------------------------
# 2. LOW-OXYGEN RESPIRATION
# ------------------------------------------------------------
#
# KEGG M00156 cbb3-type cytochrome c oxidase
#
# standard:
# K00404 + K00405 + K00407 + K00406
#
# fused:
# K15862 + K00407 + K00406
#
# KEGG M00153 cytochrome bd:
#
# K00425 + K00426 +
#   (K00424 OR K22501)
#

CBB3_STANDARD = [
    "K00404",
    "K00405",
    "K00407",
    "K00406",
]

CBB3_FUSED = [
    "K15862",
    "K00407",
    "K00406",
]

CYD_ROUTE_K00424 = [
    "K00425",
    "K00426",
    "K00424",
]

CYD_ROUTE_K22501 = [
    "K00425",
    "K00426",
    "K22501",
]

CYD_CORE = [
    "K00425",
    "K00426",
]


# ------------------------------------------------------------
# 3. NITRATE RESPIRATION
# ------------------------------------------------------------

NAR_GHI = [
    "K00370",
    "K00371",
    "K00374",
]

NAP_AB = [
    "K02567",
    "K02568",
]


# ------------------------------------------------------------
# 4. NITRITE -> NO
# ------------------------------------------------------------

NIRK = "K00368"
NIRS = "K15864"


# ------------------------------------------------------------
# 5. DNRA
# ------------------------------------------------------------

NIR_BD = [
    "K00362",
    "K00363",
]

NRF_AH = [
    "K03385",
    "K15876",
]


# ------------------------------------------------------------
# 6. NITRIC OXIDE REDUCTION
# ------------------------------------------------------------
#
# K04561 = NorB
#
# NorB is retained as the catalytic marker because K04561
# can include qNor systems that do not require NorC.
#
# K02305 NorC is retained as supporting evidence.
#

NORB = "K04561"
NORC = "K02305"


# ------------------------------------------------------------
# 7. NITROUS OXIDE REDUCTION
# ------------------------------------------------------------

NOSZ = "K00376"


# ------------------------------------------------------------
# 8. DISSIMILATORY SULFATE-REDUCTION CORE
# ------------------------------------------------------------
#
# Conservative Sat-AprAB-DsrAB core.
#
# This is deliberately called a "core" rather than claiming
# the complete KEGG M00596 module.
#

SULFATE_REDUCTION_CORE = [
    "K00958",
    "K00394",
    "K00395",
    "K11180",
    "K11181",
]


# ------------------------------------------------------------
# 9. REDUCED-SULFUR OXIDATION
# ------------------------------------------------------------
#
# KEGG Sox core:
#
# K17222 SoxA
# K17223 SoxX
# K17224 SoxB
# K17226 SoxY
# K17227 SoxZ
#
# SoxC/SoxD are optional in KEGG M00595.
#
# Alternative:
# K17218 SQR
#

SOX_CORE = [
    "K17222",
    "K17223",
    "K17224",
    "K17226",
    "K17227",
]

SQR = "K17218"


# ------------------------------------------------------------
# 10. SULFITE OXIDATION
# ------------------------------------------------------------

SOE_ABC = [
    "K21307",
    "K21308",
    "K21309",
]

SOR_AB = [
    "K05301",
    "K00386",
]


# ============================================================
# SCORE MAGs
# ============================================================

results = []

for mag in mags:

    kos = mag_kos[
        mag
    ]

    row = {
        "Final_MAG":
            mag,
        **mag_metadata[mag],
    }


    # ========================================================
    # 1. STANDARD AEROBIC RESPIRATION
    # ========================================================

    cox_standard_evidence = (
        fraction_present(
            kos,
            COX_ROUTE_STANDARD
        )
    )

    cox_fused_evidence = (
        fraction_present(
            kos,
            COX_ROUTE_FUSED
        )
    )

    cox_standard = all_present(
        kos,
        COX_ROUTE_STANDARD
    )

    cox_fused = all_present(
        kos,
        COX_ROUTE_FUSED
    )

    cox_complete = int(
        cox_standard
        or cox_fused
    )

    cox_evidence = max(
        cox_standard_evidence,
        cox_fused_evidence
    )

    row[
        "Aerobic_M00155_module_evidence"
    ] = round(
        cox_evidence,
        3
    )

    row[
        "Aerobic_M00155_complete"
    ] = cox_complete

    row[
        "Aerobic_M00155_route"
    ] = (
        "standard"
        if cox_standard
        else (
            "fused"
            if cox_fused
            else "None"
        )
    )

    optional_cox = [
        ko
        for ko in COX_OPTIONAL
        if ko in kos
    ]

    row[
        "Aerobic_M00155_optional_subunits"
    ] = (
        ";".join(
            optional_cox
        )
        if optional_cox
        else "None"
    )


    cyt_o = record_simple_module(
        row,
        kos,
        "Aerobic_M00417_cytochrome_o",
        CYTOCHROME_O,
    )


    aerobic_routes = []

    if cox_complete:
        aerobic_routes.append(
            "cytochrome_c_oxidase"
        )

    if cyt_o:
        aerobic_routes.append(
            "cytochrome_o"
        )

    row[
        "Aerobic_respiration_routes"
    ] = (
        ";".join(
            aerobic_routes
        )
        if aerobic_routes
        else "None"
    )

    row[
        "Capability_Aerobic_respiration"
    ] = int(
        cox_complete
        or cyt_o
    )


    # ========================================================
    # 2. LOW-OXYGEN RESPIRATION
    # ========================================================

    cbb3_standard_evidence = (
        fraction_present(
            kos,
            CBB3_STANDARD
        )
    )

    cbb3_fused_evidence = (
        fraction_present(
            kos,
            CBB3_FUSED
        )
    )

    cbb3_standard = all_present(
        kos,
        CBB3_STANDARD
    )

    cbb3_fused = all_present(
        kos,
        CBB3_FUSED
    )

    cbb3_complete = int(
        cbb3_standard
        or cbb3_fused
    )

    row[
        "LowO2_M00156_cbb3_module_evidence"
    ] = round(
        max(
            cbb3_standard_evidence,
            cbb3_fused_evidence
        ),
        3
    )

    row[
        "LowO2_M00156_cbb3_complete"
    ] = cbb3_complete

    row[
        "LowO2_M00156_cbb3_route"
    ] = (
        "standard"
        if cbb3_standard
        else (
            "fused"
            if cbb3_fused
            else "None"
        )
    )


    cyd_route1_evidence = (
        fraction_present(
            kos,
            CYD_ROUTE_K00424
        )
    )

    cyd_route2_evidence = (
        fraction_present(
            kos,
            CYD_ROUTE_K22501
        )
    )

    cyd_route1 = all_present(
        kos,
        CYD_ROUTE_K00424
    )

    cyd_route2 = all_present(
        kos,
        CYD_ROUTE_K22501
    )

    cyd_exact = int(
        cyd_route1
        or cyd_route2
    )

    cyd_core = all_present(
        kos,
        CYD_CORE
    )

    row[
        "LowO2_M00153_bd_module_evidence"
    ] = round(
        max(
            cyd_route1_evidence,
            cyd_route2_evidence
        ),
        3
    )

    row[
        "LowO2_M00153_bd_complete"
    ] = cyd_exact

    row[
        "LowO2_CydAB_core_present"
    ] = cyd_core

    row[
        "LowO2_M00153_bd_route"
    ] = (
        "K00424"
        if cyd_route1
        else (
            "K22501"
            if cyd_route2
            else "None"
        )
    )


    low_o2_routes = []

    if cbb3_complete:
        low_o2_routes.append(
            "cbb3"
        )

    if cyd_exact:
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
        cbb3_complete
        or cyd_exact
    )


    # Sensitivity evidence:
    # CydAB core without required small subunit.

    row[
        "Capability_LowO2_respiration_CydAB_core_sensitivity"
    ] = int(
        cbb3_complete
        or cyd_core
    )


    # ========================================================
    # 3. NITRATE RESPIRATION
    # ========================================================

    nar = record_simple_module(
        row,
        kos,
        "Nitrate_NarGHI",
        NAR_GHI,
    )

    nap = record_simple_module(
        row,
        kos,
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
        ";".join(
            nitrate_routes
        )
        if nitrate_routes
        else "None"
    )

    row[
        "Capability_Nitrate_respiration"
    ] = int(
        nar
        or nap
    )


    # ========================================================
    # 4. NITRITE -> NO
    # ========================================================

    nirk = record_single(
        row,
        kos,
        "Nitrite_NO_NirK",
        NIRK,
    )

    nirs = record_single(
        row,
        kos,
        "Nitrite_NO_NirS",
        NIRS,
    )

    nitrite_routes = []

    if nirk:
        nitrite_routes.append(
            "NirK"
        )

    if nirs:
        nitrite_routes.append(
            "NirS"
        )

    row[
        "Nitrite_NO_routes"
    ] = (
        ";".join(
            nitrite_routes
        )
        if nitrite_routes
        else "None"
    )

    row[
        "Capability_Nitrite_to_NO"
    ] = int(
        nirk
        or nirs
    )


    # ========================================================
    # 5. DNRA
    # ========================================================

    nirbd = record_simple_module(
        row,
        kos,
        "DNRA_NirBD",
        NIR_BD,
    )

    nrfah = record_simple_module(
        row,
        kos,
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
        nirbd
        or nrfah
    )


    # ========================================================
    # 6. NITRIC OXIDE REDUCTION
    # ========================================================

    norb = record_single(
        row,
        kos,
        "NO_reduction_NorB",
        NORB,
    )

    norc = record_single(
        row,
        kos,
        "NO_reduction_NorC_support",
        NORC,
    )

    row[
        "NO_reduction_NorBC_both"
    ] = int(
        norb
        and norc
    )

    row[
        "Capability_NO_reduction"
    ] = norb


    # ========================================================
    # 7. NITROUS OXIDE REDUCTION
    # ========================================================

    nosz = record_single(
        row,
        kos,
        "N2O_reduction_NosZ",
        NOSZ,
    )

    row[
        "Capability_N2O_reduction"
    ] = nosz


    # ========================================================
    # 8. DISSIMILATORY SULFATE REDUCTION CORE
    # ========================================================

    sulfate_core = (
        record_simple_module(
            row,
            kos,
            "Sulfate_reduction_Sat_AprAB_DsrAB",
            SULFATE_REDUCTION_CORE,
        )
    )

    row[
        "Capability_Dissimilatory_sulfate_reduction_core"
    ] = sulfate_core


    # ========================================================
    # 9. REDUCED-SULFUR OXIDATION
    # ========================================================

    sox = record_simple_module(
        row,
        kos,
        "Sulfur_oxidation_Sox_core",
        SOX_CORE,
    )

    sqr = record_single(
        row,
        kos,
        "Sulfur_oxidation_SQR",
        SQR,
    )

    sulfur_routes = []

    if sox:
        sulfur_routes.append(
            "Sox"
        )

    if sqr:
        sulfur_routes.append(
            "SQR"
        )

    row[
        "Reduced_sulfur_oxidation_routes"
    ] = (
        ";".join(
            sulfur_routes
        )
        if sulfur_routes
        else "None"
    )

    row[
        "Capability_Reduced_sulfur_oxidation"
    ] = int(
        sox
        or sqr
    )


    # ========================================================
    # 10. SULFITE OXIDATION
    # ========================================================

    soe = record_simple_module(
        row,
        kos,
        "Sulfite_oxidation_SoeABC",
        SOE_ABC,
    )

    sor = record_simple_module(
        row,
        kos,
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
        soe
        or sor
    )


    # ========================================================
    # FINAL REDOX STRATEGY SCORE
    # ========================================================

    capability_columns = [
        "Capability_Aerobic_respiration",
        "Capability_LowO2_respiration",
        "Capability_Nitrate_respiration",
        "Capability_Nitrite_to_NO",
        "Capability_DNRA",
        "Capability_NO_reduction",
        "Capability_N2O_reduction",
        "Capability_Dissimilatory_sulfate_reduction_core",
        "Capability_Reduced_sulfur_oxidation",
        "Capability_Sulfite_oxidation",
    ]

    capability_count = sum(
        int(
            row[column]
        )
        for column
        in capability_columns
    )

    row[
        "Redox_capabilities_detected"
    ] = capability_count

    row[
        "Redox_capabilities_total"
    ] = 10

    row[
        "Redox_strategy_score"
    ] = round(
        capability_count
        / 10,
        3
    )


    # ========================================================
    # LOW-O2 SENSITIVITY SCORE
    # ========================================================

    sensitivity_count = (
        capability_count
        - row[
            "Capability_LowO2_respiration"
        ]
        + row[
            "Capability_LowO2_respiration_CydAB_core_sensitivity"
        ]
    )

    row[
        "Redox_capabilities_detected_CydAB_sensitivity"
    ] = sensitivity_count

    row[
        "Redox_strategy_score_CydAB_sensitivity"
    ] = round(
        sensitivity_count
        / 10,
        3
    )

    results.append(
        row
    )


# ============================================================
# WRITE FINAL TABLE
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
# MODULE-EVIDENCE TABLE
# ============================================================

module_columns = [
    column
    for column in fieldnames
    if (
        column
        == "Final_MAG"
        or column
        in {
            "Completeness",
            "Contamination",
            "Site",
            "Treatment",
            "Compartment",
            "GTDB_phylum",
            "GTDB_order",
        }
        or "module_evidence"
        in column
        or "markers_detected"
        in column
        or "markers_total"
        in column
        or column.endswith(
            "_complete"
        )
        or column.endswith(
            "_routes"
        )
        or column.endswith(
            "_route"
        )
        or "support"
        in column
        or "optional"
        in column
        or "core_present"
        in column
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
        writer.writerow(
            row
        )


# ============================================================
# SUMMARY
# ============================================================

capability_summary = [
    (
        "Aerobic_respiration",
        "Capability_Aerobic_respiration",
    ),
    (
        "LowO2_respiration",
        "Capability_LowO2_respiration",
    ),
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
]

summary_rows = []

for label, column in capability_summary:

    count = sum(
        int(
            row[column]
        )
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


# Low-O2 sensitivity prevalence

low_o2_sensitivity = sum(
    row[
        "Capability_LowO2_respiration_CydAB_core_sensitivity"
    ]
    for row in results
)

summary_rows.append({
    "Summary_type":
        "Sensitivity",
    "Category":
        "LowO2_with_CydAB_core",
    "MAGs":
        low_o2_sensitivity,
    "Percent":
        f"{100 * low_o2_sensitivity / len(results):.2f}",
})


# Score distribution

score_distribution = Counter(
    row[
        "Redox_capabilities_detected"
    ]
    for row in results
)

for number in sorted(
    score_distribution
):

    count = score_distribution[
        number
    ]

    summary_rows.append({
        "Summary_type":
            "Capability_count_distribution",
        "Category":
            f"{number}_of_10",
        "MAGs":
            count,
        "Percent":
            f"{100 * count / len(results):.2f}",
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
    "FINAL REDOX-FLEXIBILITY STRATEGY ANALYSIS - V2"
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
        int(
            row[column]
        )
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
    "Low-O2 sensitivity:"
)

print(
    "  Primary KEGG-defined low-O2 capability: "
    f"{sum(row['Capability_LowO2_respiration'] for row in results)}"
)

print(
    "  Including CydAB core sensitivity:       "
    f"{low_o2_sensitivity}"
)


print()
print(
    "Redox strategy distribution:"
)

for number in sorted(
    score_distribution
):

    count = score_distribution[
        number
    ]

    print(
        f"  {number}/10 : "
        f"{count:>4} MAGs "
        f"({100 * count / len(results):6.2f}%)"
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
    "number of detected biological capabilities / 10"
)

print(
    "  Aerobic and low-O2 respiration are treated "
    "as distinct ecological capabilities."
)

print(
    "  Alternative molecular routes contribute only "
    "once to each capability."
)

print(
    "  The score measures genome-encoded redox breadth, "
    "not experimentally measured physiological flexibility."
)

print(
    "  Genome completeness is retained but individual "
    "MAG scores are not completeness-corrected."
)

print(
    "=" * 72
)
