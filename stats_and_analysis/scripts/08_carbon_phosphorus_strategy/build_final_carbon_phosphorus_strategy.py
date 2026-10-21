#!/usr/bin/env python3

"""
FINAL CARBON + PHOSPHORUS RESOURCE STRATEGY

Framework
---------
KO annotations
    ->
(1) module/pathway evidence
    ->
(2) binary biological capabilities
    ->
(3) carbon and phosphorus subdomain scores
    ->
(4) combined C+P resource strategy score

CARBON SUBDOMAIN
----------------
1. Autotrophic carbon fixation
   Positive when at least one recognised carbon-fixation
   route is supported:

      - Calvin cycle diagnostic core
      - Wood-Ljungdahl diagnostic core
      - reductive TCA diagnostic core
      - 3-hydroxypropionate bicycle diagnostic core

2. C1 / methylotrophic carbon utilisation
   Requires:
      methanol oxidation
      AND
      a recognised formaldehyde-assimilation pathway
      (RuMP OR serine pathway)

Carbon score:
    S_C = (carbon fixation + C1 utilisation) / 2


PHOSPHORUS SUBDOMAIN
--------------------
1. High-affinity phosphate uptake
   PstSCAB

2. Organic-P mobilisation
   alkaline phosphatase

3. Phosphonate acquisition
   PhnCDE

4. Polyphosphate turnover
   synthesis/storage AND mobilisation

Phosphorus score:
    S_P = sum(P capabilities) / 4


COMBINED SCORE
--------------
    S_C+P = (S_C + S_P) / 2


Important
---------
- Individual KOs are NOT counted as ecological traits.
- Alternative pathways collapse into one capability.
- Carbon and phosphorus receive equal final weight.
- Module evidence is retained for QC.
- CAZymes / Salicornia-derived polysaccharide utilisation
  are NOT included here. They will be analysed separately.
- Scores represent genome-encoded potential, not measured activity.
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
      "carbon_phosphorus_resource"
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
    / "carbon_phosphorus_strategy_final_373_MAGs.tsv"
)

MODULE_TABLE = (
    TABLE_DIR
    / "carbon_phosphorus_module_evidence_373_MAGs.tsv"
)

SUMMARY_TABLE = (
    TABLE_DIR
    / "carbon_phosphorus_strategy_final_summary.tsv"
)

ITOL_FILE = (
    ITOL_DIR
    / "itol_carbon_phosphorus_strategy_score_352.txt"
)


# ============================================================
# BASIC HELPERS
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
        / len(required)
    )


def record_module(
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
        else 0.0
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

            mag_metadata[
                mag
            ] = {
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

            mag_kos[
                mag
            ].add(
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
# CARBON MODULE DEFINITIONS
# ============================================================


# ------------------------------------------------------------
# 1A. CALVIN CYCLE DIAGNOSTIC CORE
# ------------------------------------------------------------
#
# Require:
#
# phosphoribulokinase:
#   K00855 OR K28462
#
# AND
#
# RuBisCO large subunit:
#   K01601
#
# K01602 (small subunit) is retained as supporting evidence
# but is not required because not all bacterial RuBisCO
# forms use the canonical small subunit architecture.
#

CALVIN_PRK = [
    "K00855",
    "K28462",
]

RUBISCO_LARGE = "K01601"
RUBISCO_SMALL = "K01602"


# ------------------------------------------------------------
# 1B. WOOD-LJUNGDAHL DIAGNOSTIC CORE
# ------------------------------------------------------------
#
# Conservative CODH/ACS core.
#
# K00198
# K14138
# K00197
# K00194
#
# This represents strong genomic evidence for the central
# carbonyl branch / acetyl-CoA synthase complex.
#
# Because this machinery can participate in reversible
# acetyl-CoA metabolism, we retain the exact pathway evidence
# and describe this as Wood-Ljungdahl core potential.
#

WOOD_LJUNGDAHL_CORE = [
    "K00198",
    "K14138",
    "K00197",
    "K00194",
]


# ------------------------------------------------------------
# 1C. REDUCTIVE TCA DIAGNOSTIC CORE
# ------------------------------------------------------------
#
# ATP citrate lyase:
#
# canonical:
#   K15230 + K15231
#
# alternative citryl-CoA route:
#   K15232 + K15233 + K15234
#

RTCA_ACL = [
    "K15230",
    "K15231",
]

RTCA_CITRYL = [
    "K15232",
    "K15233",
    "K15234",
]


# ------------------------------------------------------------
# 1D. 3-HYDROXYPROPIONATE BICYCLE
# ------------------------------------------------------------
#
# Use the highly diagnostic core enzymes:
#
# K14468 = malonyl-CoA reductase
# K14469 = multifunctional 3HP enzyme
#
# The full KEGG pathway contains many additional enzymes
# shared with central metabolism.
#
# We therefore retain a strict diagnostic-core call and
# also store marker evidence.
#

HP_BICYCLE_CORE = [
    "K14468",
    "K14469",
]


# ------------------------------------------------------------
# 2. C1 / METHYLOTROPHIC UTILISATION
# ------------------------------------------------------------

#
# Methanol oxidation
#
# MxaFI:
#   K14028 + K14029
#
# OR
#
# XoxF:
#   K23995
#

METHANOL_MXAFI = [
    "K14028",
    "K14029",
]

METHANOL_XOXF = "K23995"


# ------------------------------------------------------------
# RuMP formaldehyde-assimilation pathway
#
# KEGG M00345 logic simplified into its functional steps:
#
# Step 1:
# HPS/PHI route:
#
#   (K08093 OR K13812) + K08094
#   OR
#   K13831
#
# Step 2:
#   K00850 OR K16370 OR K21071
#
# Step 3:
#   K01624
#

RUMP_HPS_ALTERNATIVES = [
    "K08093",
    "K13812",
]

RUMP_PHI = "K08094"

RUMP_FUSED = "K13831"

RUMP_STEP2 = [
    "K00850",
    "K16370",
    "K21071",
]

RUMP_STEP3 = "K01624"


# ------------------------------------------------------------
# Serine formaldehyde-assimilation pathway
# KEGG M00346
# ------------------------------------------------------------

SERINE_PATHWAY = [
    "K00600",
    "K00830",
    "K00018",
    "K11529",
    "K01689",
    "K01595",
    "K00024",
    "K08691",
]

# Malate-CoA ligase is represented by:
# K08692 + K14067
SERINE_MALATE_COA_LIGASE = [
    "K08692",
    "K14067",
]


# ============================================================
# PHOSPHORUS MODULE DEFINITIONS
# ============================================================


# ------------------------------------------------------------
# 1. HIGH-AFFINITY PHOSPHATE UPTAKE
# PstSCAB
# ------------------------------------------------------------

PST_SCAB = [
    "K02040",
    "K02037",
    "K02038",
    "K02036",
]


# ------------------------------------------------------------
# 2. ORGANIC-P MOBILISATION
# alkaline phosphatases
# ------------------------------------------------------------

ALKALINE_PHOSPHATASES = [
    "K01077",
    "K01113",
]


# ------------------------------------------------------------
# 3. PHOSPHONATE ACQUISITION
# PhnCDE
# ------------------------------------------------------------

PHN_CDE = [
    "K02041",
    "K02042",
    "K02044",
]


# ------------------------------------------------------------
# 4. POLYPHOSPHATE TURNOVER
# ------------------------------------------------------------
#
# Storage/synthesis:
#
# K00937 = PPK1
# K22468 = PPK2
#
# Mobilisation:
#
# K01524 = exopolyphosphatase / GppA
#
# Primary capability requires evidence for BOTH
# storage/synthesis and mobilisation.
#
# Individual components remain separately available.
#

POLYP_SYNTHESIS = [
    "K00937",
    "K22468",
]

POLYP_MOBILISATION = "K01524"


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
    # CARBON 1:
    # AUTOTROPHIC CARBON FIXATION
    # ========================================================


    # --------------------------------------------------------
    # Calvin diagnostic core
    # --------------------------------------------------------

    prk_hits = [
        ko
        for ko in CALVIN_PRK
        if ko in kos
    ]

    prk = int(
        bool(
            prk_hits
        )
    )

    rubisco_large = int(
        RUBISCO_LARGE
        in kos
    )

    rubisco_small = int(
        RUBISCO_SMALL
        in kos
    )

    calvin = int(
        prk
        and
        rubisco_large
    )

    row[
        "Carbon_Calvin_PRK_marker_names"
    ] = (
        ";".join(
            prk_hits
        )
        if prk_hits
        else "None"
    )

    row[
        "Carbon_Calvin_PRK_present"
    ] = prk

    row[
        "Carbon_Calvin_RuBisCO_large_present"
    ] = rubisco_large

    row[
        "Carbon_Calvin_RuBisCO_small_support"
    ] = rubisco_small

    row[
        "Carbon_Calvin_diagnostic_core_complete"
    ] = calvin


    # --------------------------------------------------------
    # Wood-Ljungdahl core
    # --------------------------------------------------------

    wlp = record_module(
        row,
        kos,
        "Carbon_Wood_Ljungdahl_core",
        WOOD_LJUNGDAHL_CORE,
    )


    # --------------------------------------------------------
    # reductive TCA
    # --------------------------------------------------------

    rtca_acl = record_module(
        row,
        kos,
        "Carbon_rTCA_ATP_citrate_lyase",
        RTCA_ACL,
    )

    rtca_citryl = record_module(
        row,
        kos,
        "Carbon_rTCA_citryl_CoA_route",
        RTCA_CITRYL,
    )

    rtca = int(
        rtca_acl
        or
        rtca_citryl
    )

    row[
        "Carbon_rTCA_diagnostic_core_complete"
    ] = rtca


    # --------------------------------------------------------
    # 3HP bicycle
    # --------------------------------------------------------

    hp_bicycle = record_module(
        row,
        kos,
        "Carbon_3HP_bicycle_core",
        HP_BICYCLE_CORE,
    )


    # --------------------------------------------------------
    # Alternative carbon-fixation routes collapse
    # into ONE carbon-fixation capability.
    # --------------------------------------------------------

    fixation_routes = []

    if calvin:
        fixation_routes.append(
            "Calvin"
        )

    if wlp:
        fixation_routes.append(
            "Wood_Ljungdahl_core"
        )

    if rtca:
        fixation_routes.append(
            "reductive_TCA"
        )

    if hp_bicycle:
        fixation_routes.append(
            "3HP_bicycle_core"
        )

    row[
        "Carbon_fixation_routes_detected"
    ] = len(
        fixation_routes
    )

    row[
        "Carbon_fixation_route_names"
    ] = (
        ";".join(
            fixation_routes
        )
        if fixation_routes
        else "None"
    )

    row[
        "Capability_Carbon_fixation"
    ] = int(
        bool(
            fixation_routes
        )
    )


    # ========================================================
    # CARBON 2:
    # C1 / METHYLOTROPHIC UTILISATION
    # ========================================================


    # --------------------------------------------------------
    # Methanol oxidation
    # --------------------------------------------------------

    mxafi = record_module(
        row,
        kos,
        "C1_MxaFI",
        METHANOL_MXAFI,
    )

    xoxf = record_single(
        row,
        kos,
        "C1_XoxF",
        METHANOL_XOXF,
    )

    methanol_oxidation = int(
        mxafi
        or
        xoxf
    )

    methanol_routes = []

    if mxafi:
        methanol_routes.append(
            "MxaFI"
        )

    if xoxf:
        methanol_routes.append(
            "XoxF"
        )

    row[
        "C1_methanol_oxidation_routes"
    ] = (
        ";".join(
            methanol_routes
        )
        if methanol_routes
        else "None"
    )

    row[
        "C1_methanol_oxidation_complete"
    ] = methanol_oxidation


    # --------------------------------------------------------
    # RuMP assimilation
    # --------------------------------------------------------

    rump_hps = any_present(
        kos,
        RUMP_HPS_ALTERNATIVES
    )

    rump_phi = int(
        RUMP_PHI
        in kos
    )

    rump_fused = int(
        RUMP_FUSED
        in kos
    )

    rump_first_step = int(
        (
            rump_hps
            and
            rump_phi
        )
        or
        rump_fused
    )

    rump_step2 = any_present(
        kos,
        RUMP_STEP2
    )

    rump_step3 = int(
        RUMP_STEP3
        in kos
    )

    rump = int(
        rump_first_step
        and
        rump_step2
        and
        rump_step3
    )

    row[
        "C1_RuMP_initial_step_complete"
    ] = rump_first_step

    row[
        "C1_RuMP_step2_complete"
    ] = rump_step2

    row[
        "C1_RuMP_step3_complete"
    ] = rump_step3

    row[
        "C1_RuMP_pathway_complete"
    ] = rump


    # --------------------------------------------------------
    # Serine pathway
    # --------------------------------------------------------

    serine_core_detected = len(
        kos
        & set(
            SERINE_PATHWAY
        )
    )

    serine_core_fraction = (
        serine_core_detected
        / len(
            SERINE_PATHWAY
        )
    )

    serine_core = all_present(
        kos,
        SERINE_PATHWAY
    )

    serine_malate_coa = all_present(
        kos,
        SERINE_MALATE_COA_LIGASE
    )

    serine_pathway = int(
        serine_core
        and
        serine_malate_coa
    )

    row[
        "C1_Serine_core_markers_detected"
    ] = serine_core_detected

    row[
        "C1_Serine_core_markers_total"
    ] = len(
        SERINE_PATHWAY
    )

    row[
        "C1_Serine_core_module_evidence"
    ] = round(
        serine_core_fraction,
        3
    )

    row[
        "C1_Serine_malate_CoA_ligase_complete"
    ] = serine_malate_coa

    row[
        "C1_Serine_pathway_complete"
    ] = serine_pathway


    # --------------------------------------------------------
    # C1 utilisation capability
    #
    # Requires:
    # methanol oxidation
    # AND
    # RuMP or serine assimilation
    # --------------------------------------------------------

    c1_assimilation = int(
        rump
        or
        serine_pathway
    )

    c1_utilisation = int(
        methanol_oxidation
        and
        c1_assimilation
    )

    row[
        "C1_assimilation_route"
    ] = (
        "RuMP;Serine"
        if (
            rump
            and
            serine_pathway
        )
        else (
            "RuMP"
            if rump
            else (
                "Serine"
                if serine_pathway
                else "None"
            )
        )
    )

    row[
        "Capability_C1_utilisation"
    ] = c1_utilisation


    # ========================================================
    # CARBON SUBDOMAIN SCORE
    # ========================================================

    carbon_capability_count = (
        row[
            "Capability_Carbon_fixation"
        ]
        +
        row[
            "Capability_C1_utilisation"
        ]
    )

    row[
        "Carbon_capabilities_detected"
    ] = carbon_capability_count

    row[
        "Carbon_capabilities_total"
    ] = 2

    row[
        "Carbon_subdomain_score"
    ] = round(
        carbon_capability_count
        / 2,
        3
    )


    # ========================================================
    # PHOSPHORUS 1:
    # HIGH-AFFINITY PHOSPHATE UPTAKE
    # ========================================================

    pst = record_module(
        row,
        kos,
        "P_PstSCAB",
        PST_SCAB,
    )

    row[
        "Capability_High_affinity_phosphate_uptake"
    ] = pst


    # ========================================================
    # PHOSPHORUS 2:
    # ORGANIC-P MOBILISATION
    # ========================================================

    alkaline_hits = [
        ko
        for ko in ALKALINE_PHOSPHATASES
        if ko in kos
    ]

    alkaline = int(
        bool(
            alkaline_hits
        )
    )

    row[
        "P_Alkaline_phosphatase_marker_names"
    ] = (
        ";".join(
            alkaline_hits
        )
        if alkaline_hits
        else "None"
    )

    row[
        "Capability_Organic_P_mobilisation"
    ] = alkaline


    # ========================================================
    # PHOSPHORUS 3:
    # PHOSPHONATE ACQUISITION
    # ========================================================

    phn = record_module(
        row,
        kos,
        "P_PhnCDE",
        PHN_CDE,
    )

    row[
        "Capability_Phosphonate_acquisition"
    ] = phn


    # ========================================================
    # PHOSPHORUS 4:
    # POLYPHOSPHATE TURNOVER
    # ========================================================

    polyP_synthesis_hits = [
        ko
        for ko in POLYP_SYNTHESIS
        if ko in kos
    ]

    polyP_synthesis = int(
        bool(
            polyP_synthesis_hits
        )
    )

    polyP_mobilisation = int(
        POLYP_MOBILISATION
        in kos
    )

    polyP_turnover = int(
        polyP_synthesis
        and
        polyP_mobilisation
    )

    row[
        "P_PolyP_synthesis_marker_names"
    ] = (
        ";".join(
            polyP_synthesis_hits
        )
        if polyP_synthesis_hits
        else "None"
    )

    row[
        "P_PolyP_synthesis_present"
    ] = polyP_synthesis

    row[
        "P_PolyP_mobilisation_K01524_present"
    ] = polyP_mobilisation

    row[
        "Capability_PolyP_turnover"
    ] = polyP_turnover


    # ========================================================
    # PHOSPHORUS SUBDOMAIN SCORE
    # ========================================================

    phosphorus_capability_columns = [
        "Capability_High_affinity_phosphate_uptake",
        "Capability_Organic_P_mobilisation",
        "Capability_Phosphonate_acquisition",
        "Capability_PolyP_turnover",
    ]

    phosphorus_capability_count = sum(
        int(
            row[
                column
            ]
        )
        for column
        in phosphorus_capability_columns
    )

    row[
        "Phosphorus_capabilities_detected"
    ] = phosphorus_capability_count

    row[
        "Phosphorus_capabilities_total"
    ] = 4

    row[
        "Phosphorus_subdomain_score"
    ] = round(
        phosphorus_capability_count
        / 4,
        3
    )


    # ========================================================
    # FINAL C + P RESOURCE STRATEGY
    # ========================================================
    #
    # Carbon and phosphorus receive equal weight.
    #
    # We do NOT simply average all 6 capabilities because
    # that would give phosphorus twice the influence of carbon.
    # ========================================================

    combined_score = (
        row[
            "Carbon_subdomain_score"
        ]
        +
        row[
            "Phosphorus_subdomain_score"
        ]
    ) / 2

    row[
        "Carbon_Phosphorus_strategy_score"
    ] = round(
        combined_score,
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
# WRITE MODULE-EVIDENCE TABLE
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

        or "marker_names"
        in column

        or "route_names"
        in column

        or "routes_detected"
        in column

        or "_present"
        in column

        or "assimilation_route"
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
        "Carbon_fixation",
        "Capability_Carbon_fixation",
    ),

    (
        "C1_utilisation",
        "Capability_C1_utilisation",
    ),

    (
        "High_affinity_phosphate_uptake",
        "Capability_High_affinity_phosphate_uptake",
    ),

    (
        "Organic_P_mobilisation",
        "Capability_Organic_P_mobilisation",
    ),

    (
        "Phosphonate_acquisition",
        "Capability_Phosphonate_acquisition",
    ),

    (
        "PolyP_turnover",
        "Capability_PolyP_turnover",
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
# Carbon route prevalence
# ------------------------------------------------------------

carbon_route_columns = [
    (
        "Calvin_core",
        "Carbon_Calvin_diagnostic_core_complete",
    ),

    (
        "Wood_Ljungdahl_core",
        "Carbon_Wood_Ljungdahl_core_complete",
    ),

    (
        "reductive_TCA_core",
        "Carbon_rTCA_diagnostic_core_complete",
    ),

    (
        "3HP_bicycle_core",
        "Carbon_3HP_bicycle_core_complete",
    ),

    (
        "RuMP_assimilation",
        "C1_RuMP_pathway_complete",
    ),

    (
        "Serine_assimilation",
        "C1_Serine_pathway_complete",
    ),

    (
        "Methanol_oxidation",
        "C1_methanol_oxidation_complete",
    ),
]

for label, column in carbon_route_columns:

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
            "Carbon_route_prevalence",

        "Category":
            label,

        "MAGs":
            count,

        "Percent":
            f"{100 * count / len(results):.2f}",
    })


# ------------------------------------------------------------
# Score distributions
# ------------------------------------------------------------

carbon_distribution = Counter(
    row[
        "Carbon_capabilities_detected"
    ]
    for row in results
)

for number in sorted(
    carbon_distribution
):

    count = carbon_distribution[
        number
    ]

    summary_rows.append({
        "Summary_type":
            "Carbon_capability_distribution",

        "Category":
            f"{number}_of_2",

        "MAGs":
            count,

        "Percent":
            f"{100 * count / len(results):.2f}",
    })


phosphorus_distribution = Counter(
    row[
        "Phosphorus_capabilities_detected"
    ]
    for row in results
)

for number in sorted(
    phosphorus_distribution
):

    count = phosphorus_distribution[
        number
    ]

    summary_rows.append({
        "Summary_type":
            "Phosphorus_capability_distribution",

        "Category":
            f"{number}_of_4",

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
        "Carbon + phosphorus resource potential\n"
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
        "C+P resource strategy breadth\n"
    )

    outfile.write(
        "DATA\n"
    )

    for row in tree_results:

        outfile.write(
            f"{row['Final_MAG']}\t"
            f"{row['Carbon_Phosphorus_strategy_score']}\n"
        )


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print(
    "=" * 72
)

print(
    "FINAL CARBON + PHOSPHORUS RESOURCE STRATEGY ANALYSIS"
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
    "Carbon pathway evidence:"
)

for label, column in carbon_route_columns:

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
    "Carbon capability distribution:"
)

for number in sorted(
    carbon_distribution
):

    count = carbon_distribution[
        number
    ]

    print(
        f"  {number}/2 : "
        f"{count:>4} MAGs "
        f"({100 * count / len(results):6.2f}%)"
    )


print()
print(
    "Phosphorus capability distribution:"
)

for number in sorted(
    phosphorus_distribution
):

    count = phosphorus_distribution[
        number
    ]

    print(
        f"  {number}/4 : "
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
    "  Carbon and phosphorus are scored as separate "
    "resource subdomains and then given equal weight."
)

print(
    "  Carbon fixation represents presence of at least one "
    "recognised diagnostic carbon-fixation route."
)

print(
    "  C1 utilisation requires methanol oxidation plus "
    "a recognised formaldehyde-assimilation pathway."
)

print(
    "  Phosphorus capabilities represent high-affinity "
    "inorganic-P uptake, organic-P mobilisation, "
    "phosphonate acquisition and polyP turnover."
)

print(
    "  CAZyme-based Salicornia polysaccharide utilisation "
    "is intentionally excluded from this primary score "
    "and will be analysed separately."
)

print(
    "=" * 72
)
