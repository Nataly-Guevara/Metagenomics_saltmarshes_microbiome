#!/usr/bin/env python3

"""
FINAL NITROGEN ACQUISITION AND TRANSFORMATION STRATEGY

Framework:
    KO annotations
        ->
    (1) module/pathway evidence
        ->
    (2) binary biological capability calls
        ->
    (3) standardized nitrogen-strategy breadth (0-1)

Eight equally weighted nitrogen capabilities:

    1. Nitrogen fixation
    2. Ammonium uptake
    3. Assimilatory nitrate reduction
    4. Urea utilisation
    5. DNRA
    6. Complete denitrification
    7. Nitrification
    8. Anammox

Nitrogen_strategy_score =
    number of capabilities detected / 8

Important:
- Individual KOs are NOT counted as independent ecological traits.
- Alternative molecular routes collapse into one capability.
- Complete multigene processes require their defined components.
- Module/pathway evidence is retained for QC.
- Scores represent genome-encoded functional breadth,
  not measured pathway activity.
- Genome completeness is retained but individual MAG scores
  are not completeness-corrected.
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
      "nitrogen_acquisition_transformation"
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
    / "nitrogen_strategy_final_373_MAGs.tsv"
)

MODULE_TABLE = (
    TABLE_DIR
    / "nitrogen_module_evidence_373_MAGs.tsv"
)

SUMMARY_TABLE = (
    TABLE_DIR
    / "nitrogen_strategy_final_summary.tsv"
)

ITOL_FILE = (
    ITOL_DIR
    / "itol_nitrogen_strategy_score_352.txt"
)


# ============================================================
# HELPERS
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

    required = set(
        required
    )

    if not required:
        return 0.0

    return (
        len(
            kos_present
            & required
        )
        / len(
            required
        )
    )


def all_present(
    kos_present,
    required
):

    return int(
        set(
            required
        )
        <= kos_present
    )


def any_present(
    kos_present,
    alternatives
):

    return int(
        bool(
            kos_present
            & set(
                alternatives
            )
        )
    )


def record_module(
    row,
    kos_present,
    name,
    required
):

    detected = len(
        kos_present
        & set(
            required
        )
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
        set(
            required
        )
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
    ] = float(
        present
    )

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
        - set(
            reader.fieldnames
        )
    )

    if missing:

        raise ValueError(
            "Missing required columns: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

    for row in reader:

        mag = row[
            "Final_MAG"
        ].strip()

        if not mag:
            continue

        if mag not in mag_metadata:

            mag_metadata[
                mag
            ] = {
                "Completeness":
                    row[
                        "Completeness"
                    ],
                "Contamination":
                    row[
                        "Contamination"
                    ],
                "Site":
                    row[
                        "Site"
                    ],
                "Treatment":
                    row[
                        "Treatment"
                    ],
                "Compartment":
                    row[
                        "Compartment"
                    ],
                "GTDB_phylum":
                    row[
                        "GTDB_phylum"
                    ],
                "GTDB_order":
                    row[
                        "GTDB_order"
                    ],
            }

        for ko in parse_ko_field(
            row[
                "KO"
            ]
        ):

            mag_kos[
                mag
            ].add(
                ko
            )


mags = sorted(
    mag_metadata
)

if len(
    mags
) != 373:

    print(
        "WARNING: expected 373 MAGs, "
        f"but found {len(mags)}."
    )


# ============================================================
# MODULE DEFINITIONS
# ============================================================


# ------------------------------------------------------------
# 1. NITROGEN FIXATION
# KEGG M00175
# ------------------------------------------------------------
#
# Canonical nitrogenase:
# K02588 NifH
# K02586 NifD
# K02591 NifK
#
# K00531 is optional in KEGG M00175 and is therefore not
# required here.
#
# Alternative nitrogenase route in KEGG:
# K22896 + K22897 + K22898 + K22899
#

NIF_CANONICAL = [
    "K02588",
    "K02586",
    "K02591",
]

NIF_ALTERNATIVE = [
    "K22896",
    "K22897",
    "K22898",
    "K22899",
]


# ------------------------------------------------------------
# 2. AMMONIUM UPTAKE
# ------------------------------------------------------------

AMT = "K03320"


# ------------------------------------------------------------
# 3. ASSIMILATORY NITRATE REDUCTION
# KEGG M00531
# ------------------------------------------------------------
#
# Step 1: nitrate -> nitrite
#
# Alternative routes:
#   K00367
#   K10534
#   K00372
#
# K00360 is optional with K00372 in the KEGG module.
#
# Step 2: nitrite -> ammonium
#
# Alternative routes:
#   K00366
#   K17877
#   K26139 + K26138
#   K00361
#

ASSIM_NITRATE_SINGLE_ROUTES = [
    "K00367",
    "K10534",
    "K00372",
]

ASSIM_NITRITE_SINGLE_ROUTES = [
    "K00366",
    "K17877",
    "K00361",
]

ASSIM_NITRITE_TWO_COMPONENT = [
    "K26139",
    "K26138",
]


# ------------------------------------------------------------
# 4. UREA UTILISATION
# ------------------------------------------------------------
#
# Recognised urease routes:
#
# Single-chain urease:
#   K01427
#
# Three-subunit urease:
#   K01428 + K01429 + K01430
#
# Fused gamma/beta form:
#   K01428 + K14048
#

UREASE_SINGLE = "K01427"

UREASE_ABC = [
    "K01428",
    "K01429",
    "K01430",
]

UREASE_FUSED = [
    "K01428",
    "K14048",
]


# ------------------------------------------------------------
# SHARED RESPIRATORY NITROGEN MODULES
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

NIR_BD = [
    "K00362",
    "K00363",
]

NRF_AH = [
    "K03385",
    "K15876",
]

NIRK = "K00368"
NIRS = "K15864"

NORB = "K04561"
NORC = "K02305"

NOSZ = "K00376"


# ------------------------------------------------------------
# 5. DNRA
# KEGG M00530
# ------------------------------------------------------------
#
# Complete DNRA requires:
#
# nitrate -> nitrite:
#   NarGHI OR NapAB
#
# AND
#
# nitrite -> ammonium:
#   NirBD OR NrfAH
#


# ------------------------------------------------------------
# 6. COMPLETE DENITRIFICATION
# KEGG M00529
# ------------------------------------------------------------
#
# nitrate -> nitrite:
#   NarGHI OR NapAB
#
# nitrite -> NO:
#   NirK OR NirS
#
# NO -> N2O:
#   NorB + NorC
#
# N2O -> N2:
#   NosZ
#
# This deliberately represents COMPLETE denitrification.
# Partial denitrification steps are retained in the evidence table
# but do not establish the full capability.
#


# ------------------------------------------------------------
# 7. NITRIFICATION
# KEGG M00528
# ------------------------------------------------------------
#
# ammonia -> hydroxylamine:
#   K28504 + K10945 + K10946
#
# hydroxylamine -> nitrite:
#   K10535
#
# Primary rule requires all four.
#
# An ammonia-monooxygenase-core sensitivity call is also retained.
#

AMO_CORE = [
    "K28504",
    "K10945",
    "K10946",
]

HAO = "K10535"


# ------------------------------------------------------------
# 8. ANAMMOX
# KEGG M00973
# ------------------------------------------------------------
#
# nitrite -> NO:
#   NirK OR NirS
#
# hydrazine synthesis:
#   K20932 + K20933 + K20934
#
# hydrazine -> N2:
#   K20935
#

HZS_ABC = [
    "K20932",
    "K20933",
    "K20934",
]

HDH = "K20935"


# ============================================================
# SCORE EACH MAG
# ============================================================

results = []

for mag in mags:

    kos = mag_kos[
        mag
    ]

    row = {
        "Final_MAG":
            mag,
        **mag_metadata[
            mag
        ],
    }


    # ========================================================
    # 1. NITROGEN FIXATION
    # ========================================================

    nif_canonical = record_module(
        row,
        kos,
        "Nfix_NifHDK",
        NIF_CANONICAL,
    )

    nif_alternative = record_module(
        row,
        kos,
        "Nfix_alternative_nitrogenase",
        NIF_ALTERNATIVE,
    )

    nfix_routes = []

    if nif_canonical:

        nfix_routes.append(
            "NifHDK"
        )

    if nif_alternative:

        nfix_routes.append(
            "alternative_nitrogenase"
        )

    row[
        "Nitrogen_fixation_routes"
    ] = (
        ";".join(
            nfix_routes
        )
        if nfix_routes
        else "None"
    )

    row[
        "Capability_Nitrogen_fixation"
    ] = int(
        nif_canonical
        or nif_alternative
    )


    # ========================================================
    # 2. AMMONIUM UPTAKE
    # ========================================================

    amt = record_single(
        row,
        kos,
        "Ammonium_uptake_Amt",
        AMT,
    )

    row[
        "Capability_Ammonium_uptake"
    ] = amt


    # ========================================================
    # 3. ASSIMILATORY NITRATE REDUCTION
    # ========================================================

    # nitrate -> nitrite
    nitrate_assim_single_hits = [
        ko
        for ko in ASSIM_NITRATE_SINGLE_ROUTES
        if ko in kos
    ]

    nitrate_assim = int(
        len(
            nitrate_assim_single_hits
        ) >= 1
    )

    row[
        "Assim_nitrate_reduction_marker_names"
    ] = (
        ";".join(
            nitrate_assim_single_hits
        )
        if nitrate_assim_single_hits
        else "None"
    )

    row[
        "Assim_nitrate_to_nitrite_complete"
    ] = nitrate_assim


    # nitrite -> ammonium

    nitrite_assim_single_hits = [
        ko
        for ko in ASSIM_NITRITE_SINGLE_ROUTES
        if ko in kos
    ]

    nitrite_assim_two_component = (
        all_present(
            kos,
            ASSIM_NITRITE_TWO_COMPONENT
        )
    )

    nitrite_assim = int(
        bool(
            nitrite_assim_single_hits
        )
        or
        nitrite_assim_two_component
    )

    row[
        "Assim_nitrite_single_marker_names"
    ] = (
        ";".join(
            nitrite_assim_single_hits
        )
        if nitrite_assim_single_hits
        else "None"
    )

    row[
        "Assim_nitrite_K26139_K26138_complete"
    ] = nitrite_assim_two_component

    row[
        "Assim_nitrite_to_ammonium_complete"
    ] = nitrite_assim


    row[
        "Capability_Assimilatory_nitrate_reduction"
    ] = int(
        nitrate_assim
        and
        nitrite_assim
    )


    # ========================================================
    # 4. UREA UTILISATION
    # ========================================================

    urease_single = record_single(
        row,
        kos,
        "Urea_single_chain_urease",
        UREASE_SINGLE,
    )

    urease_abc = record_module(
        row,
        kos,
        "Urea_UreABC",
        UREASE_ABC,
    )

    urease_fused = record_module(
        row,
        kos,
        "Urea_fused_gamma_beta",
        UREASE_FUSED,
    )

    urea_routes = []

    if urease_single:

        urea_routes.append(
            "single_chain"
        )

    if urease_abc:

        urea_routes.append(
            "UreABC"
        )

    if urease_fused:

        urea_routes.append(
            "fused_gamma_beta"
        )

    row[
        "Urea_utilisation_routes"
    ] = (
        ";".join(
            urea_routes
        )
        if urea_routes
        else "None"
    )

    row[
        "Capability_Urea_utilisation"
    ] = int(
        urease_single
        or urease_abc
        or urease_fused
    )


    # ========================================================
    # SHARED RESPIRATORY MODULE EVIDENCE
    # ========================================================

    nar = record_module(
        row,
        kos,
        "Resp_NarGHI",
        NAR_GHI,
    )

    nap = record_module(
        row,
        kos,
        "Resp_NapAB",
        NAP_AB,
    )

    respiratory_nitrate = int(
        nar
        or nap
    )


    nirbd = record_module(
        row,
        kos,
        "DNRA_NirBD",
        NIR_BD,
    )

    nrfah = record_module(
        row,
        kos,
        "DNRA_NrfAH",
        NRF_AH,
    )

    dnra_nitrite_to_ammonium = int(
        nirbd
        or nrfah
    )


    nirk = record_single(
        row,
        kos,
        "Denitrification_NirK",
        NIRK,
    )

    nirs = record_single(
        row,
        kos,
        "Denitrification_NirS",
        NIRS,
    )

    denit_nitrite_to_no = int(
        nirk
        or nirs
    )


    norb = record_single(
        row,
        kos,
        "Denitrification_NorB",
        NORB,
    )

    norc = record_single(
        row,
        kos,
        "Denitrification_NorC",
        NORC,
    )

    denit_no_to_n2o = int(
        norb
        and norc
    )


    nosz = record_single(
        row,
        kos,
        "Denitrification_NosZ",
        NOSZ,
    )


    # ========================================================
    # 5. DNRA
    # ========================================================

    row[
        "DNRA_nitrate_to_nitrite"
    ] = respiratory_nitrate

    row[
        "DNRA_nitrite_to_ammonium"
    ] = dnra_nitrite_to_ammonium

    row[
        "Capability_DNRA"
    ] = int(
        respiratory_nitrate
        and
        dnra_nitrite_to_ammonium
    )


    # ========================================================
    # 6. COMPLETE DENITRIFICATION
    # ========================================================

    denitrification_steps = (
        respiratory_nitrate
        +
        denit_nitrite_to_no
        +
        denit_no_to_n2o
        +
        nosz
    )

    row[
        "Denitrification_steps_detected"
    ] = denitrification_steps

    row[
        "Denitrification_steps_total"
    ] = 4

    row[
        "Denitrification_module_evidence"
    ] = round(
        denitrification_steps
        / 4,
        3
    )

    row[
        "Capability_Complete_denitrification"
    ] = int(
        denitrification_steps
        == 4
    )


    # ========================================================
    # 7. NITRIFICATION
    # ========================================================

    amo_core = record_module(
        row,
        kos,
        "Nitrification_AmoABC",
        AMO_CORE,
    )

    hao = record_single(
        row,
        kos,
        "Nitrification_HAO",
        HAO,
    )

    row[
        "Capability_Nitrification"
    ] = int(
        amo_core
        and
        hao
    )

    # Sensitivity / supporting evidence:
    row[
        "Nitrification_AmoABC_core_only"
    ] = amo_core


    # ========================================================
    # 8. ANAMMOX
    # ========================================================

    anammox_nitrite_to_no = int(
        nirk
        or nirs
    )

    hzs = record_module(
        row,
        kos,
        "Anammox_HzsABC",
        HZS_ABC,
    )

    hdh = record_single(
        row,
        kos,
        "Anammox_Hdh",
        HDH,
    )

    anammox_steps = (
        anammox_nitrite_to_no
        +
        hzs
        +
        hdh
    )

    row[
        "Anammox_steps_detected"
    ] = anammox_steps

    row[
        "Anammox_steps_total"
    ] = 3

    row[
        "Anammox_module_evidence"
    ] = round(
        anammox_steps
        / 3,
        3
    )

    row[
        "Capability_Anammox"
    ] = int(
        anammox_steps
        == 3
    )


    # ========================================================
    # FINAL NITROGEN STRATEGY
    # ========================================================

    capability_columns = [
        "Capability_Nitrogen_fixation",
        "Capability_Ammonium_uptake",
        "Capability_Assimilatory_nitrate_reduction",
        "Capability_Urea_utilisation",
        "Capability_DNRA",
        "Capability_Complete_denitrification",
        "Capability_Nitrification",
        "Capability_Anammox",
    ]

    capability_count = sum(
        int(
            row[
                column
            ]
        )
        for column
        in capability_columns
    )

    row[
        "Nitrogen_capabilities_detected"
    ] = capability_count

    row[
        "Nitrogen_capabilities_total"
    ] = len(
        capability_columns
    )

    row[
        "Nitrogen_strategy_score"
    ] = round(
        capability_count
        / len(
            capability_columns
        ),
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
            "_route_names"
        )
        or "steps_detected"
        in column
        or "steps_total"
        in column
        or "core_only"
        in column
        or "marker_names"
        in column
        or "nitrate_to_nitrite"
        in column
        or "nitrite_to_ammonium"
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
        "Nitrogen_fixation",
        "Capability_Nitrogen_fixation",
    ),
    (
        "Ammonium_uptake",
        "Capability_Ammonium_uptake",
    ),
    (
        "Assimilatory_nitrate_reduction",
        "Capability_Assimilatory_nitrate_reduction",
    ),
    (
        "Urea_utilisation",
        "Capability_Urea_utilisation",
    ),
    (
        "DNRA",
        "Capability_DNRA",
    ),
    (
        "Complete_denitrification",
        "Capability_Complete_denitrification",
    ),
    (
        "Nitrification",
        "Capability_Nitrification",
    ),
    (
        "Anammox",
        "Capability_Anammox",
    ),
]

summary_rows = []

for label, column in capability_summary:

    count = sum(
        int(
            row[
                column
            ]
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


# ------------------------------------------------------------
# Supporting / sensitivity summaries
# ------------------------------------------------------------

amo_core_count = sum(
    int(
        row[
            "Nitrification_AmoABC_core_only"
        ]
    )
    for row in results
)

summary_rows.append({
    "Summary_type":
        "Supporting_evidence",

    "Category":
        "AmoABC_core_only",

    "MAGs":
        amo_core_count,

    "Percent":
        f"{100 * amo_core_count / len(results):.2f}",
})


# ------------------------------------------------------------
# Strategy breadth distribution
# ------------------------------------------------------------

score_distribution = Counter(
    row[
        "Nitrogen_capabilities_detected"
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
            f"{number}_of_8",

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
            item
            .strip()
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
        "Nitrogen acquisition and transformation\n"
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
        "Nitrogen capability breadth\n"
    )

    outfile.write(
        "DATA\n"
    )

    for row in tree_results:

        outfile.write(
            f"{row['Final_MAG']}\t"
            f"{row['Nitrogen_strategy_score']}\n"
        )


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print(
    "=" * 72
)

print(
    "FINAL NITROGEN ACQUISITION AND TRANSFORMATION ANALYSIS"
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
            row[
                column
            ]
        )
        for row in results
    )

    percent = (
        100
        * count
        / len(
            results
        )
    )

    print(
        f"  {label:<40}"
        f"{count:>4} MAGs "
        f"({percent:6.2f}%)"
    )


print()
print(
    "Supporting evidence:"
)

print(
    f"  AmoABC core detected: "
    f"{amo_core_count} MAGs "
    f"({100 * amo_core_count / len(results):.2f}%)"
)


print()
print(
    "Nitrogen strategy distribution:"
)

for number in sorted(
    score_distribution
):

    count = score_distribution[
        number
    ]

    print(
        f"  {number}/8 : "
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
    "  Nitrogen_strategy_score = "
    "number of detected biological capabilities / 8"
)

print(
    "  The score represents genome-encoded breadth of "
    "nitrogen acquisition and transformation potential."
)

print(
    "  Complete DNRA, denitrification, nitrification and "
    "anammox require their predefined multistep pathway rules."
)

print(
    "  Partial pathway evidence is retained in the module "
    "table but does not automatically establish a capability."
)

print(
    "  The nitrogen and redox axes may share individual "
    "pathways because they answer different ecological questions."
)

print(
    "  The six ecological-strategy axes should therefore not "
    "be summed into one global ecological score."
)

print(
    "=" * 72
)
