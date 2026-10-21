#!/usr/bin/env python3

"""
Build final Colonisation and Plant-Interaction ecological strategy.

Framework
---------
genes -> defined functional systems -> binary capabilities
      -> standardized ecological capability breadth

FINAL CAPABILITIES
==================

1. Chemotaxis
   Core markers:
       K03406  MCP
       K03407  CheA
       K03408  CheW
       K03413  CheY

   Final rule:
       >= 3 of 4 core markers detected

   CheR/CheB are retained as supporting evidence only.


2. Motility
   Alternative routes collapse into ONE capability.

   A. Flagellar motility

      Four submodules:

      Motor:
          K02556 + K02557

      Switch:
          K02409 + K02410 + K02416 + K02417

      Export:
          K02400 + K02401 + K02419 + K02420 + K02421

      Hook/filament:
          K02390 + K02406 + K02407

      Final flagellar rule:
          >= 3 of the 4 submodules complete

   OR

   B. Type-IV-pilus / twitching motility

      Assembly:
          K02662 K02663 K02664 K02665 K02666

      Retraction:
          PilT K02669

      Final twitching rule:
          >= 4 of 5 assembly markers
          AND PilT present

      PilU and twitching signalling genes are retained as
      supporting evidence only.


3. Attachment / biofilm-matrix potential

   Any ONE recognized attachment/matrix route qualifies:

   A. PGA matrix
      >= 3 of 4:
          K11935 K11931 K11936 K11937

   B. Curli
      complete:
          K04334 + K04335

   C. LapA adhesion system
      LapA K12549
      AND complete transport:
          K12541 + K12542 + K12543

   D. Tad/Flp attachment
      >= 6 of 7:
          K02651 K02278 K02279 K02280
          K02283 K12510 K12511

   E. Succinoglycan/EPS
      complete:
          K16566 + K16567 + K16568

   Alternative routes collapse into ONE capability.


4. Direct plant-interaction potential

   Final rule:
       ACC deaminase K01505

   IAA IAM evidence is retained only as supplementary evidence:
       K00466 + K21801

   Broad amidase K01426 is also retained as support only and
   NEVER used by itself to infer IAA production.


FINAL SCORE
===========

Colonisation_host_strategy_score =
    (
        Chemotaxis
        + Motility
        + Attachment_biofilm_matrix
        + Direct_plant_interaction
    ) / 4

Possible scores:
    0.00
    0.25
    0.50
    0.75
    1.00

The score represents standardized ecological capability breadth.
It does NOT represent physiological performance or percentage
plant-growth promotion.


OUTPUTS
=======

results/ecological_strategies/colonisation_host_interaction/tables/

    colonisation_host_strategy_final_373_MAGs.tsv
    colonisation_host_module_evidence_373_MAGs.tsv
    colonisation_host_strategy_final_summary.tsv

results/ecological_strategies/colonisation_host_interaction/itol/

    itol_colonisation_host_strategy_score_352.txt
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

TABLE_DIR = (
    BASE
    / "results/ecological_strategies/"
      "colonisation_host_interaction/tables"
)

ITOL_DIR = (
    BASE
    / "results/ecological_strategies/"
      "colonisation_host_interaction/itol"
)

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
    / "colonisation_host_strategy_final_373_MAGs.tsv"
)

EVIDENCE_TABLE = (
    TABLE_DIR
    / "colonisation_host_module_evidence_373_MAGs.tsv"
)

SUMMARY_TABLE = (
    TABLE_DIR
    / "colonisation_host_strategy_final_summary.tsv"
)

ITOL_FILE = (
    ITOL_DIR
    / "itol_colonisation_host_strategy_score_352.txt"
)


# ============================================================
# CAPABILITY 1: CHEMOTAXIS
# ============================================================

CHEMOTAXIS_CORE = {
    "K03406",  # MCP
    "K03407",  # CheA
    "K03408",  # CheW
    "K03413",  # CheY
}

CHEMOTAXIS_SUPPORT = {
    "K00575",  # CheR
    "K03412",  # CheB
}


# ============================================================
# CAPABILITY 2A: FLAGELLAR MOTILITY
# ============================================================

FLAGELLAR_MODULES = {

    "Motor": {
        "K02556",
        "K02557",
    },

    "Switch": {
        "K02409",
        "K02410",
        "K02416",
        "K02417",
    },

    "Export": {
        "K02400",
        "K02401",
        "K02419",
        "K02420",
        "K02421",
    },

    "Hook_filament": {
        "K02390",
        "K02406",
        "K02407",
    },
}


# ============================================================
# CAPABILITY 2B: TYPE-IV-PILUS / TWITCHING
# ============================================================

TYPEIV_ASSEMBLY = {
    "K02662",
    "K02663",
    "K02664",
    "K02665",
    "K02666",
}

PILT = "K02669"
PILU = "K02670"

TWITCHING_SIGNAL = {
    "K02657",
    "K02658",
    "K02659",
    "K02660",
}


# ============================================================
# CAPABILITY 3: ATTACHMENT / BIOFILM-MATRIX
# ============================================================

PGA_MATRIX = {
    "K11935",
    "K11931",
    "K11936",
    "K11937",
}

CURLI = {
    "K04334",
    "K04335",
}

LAPA = "K12549"

LAPA_TRANSPORT = {
    "K12541",
    "K12542",
    "K12543",
}

TAD_FLP = {
    "K02651",
    "K02278",
    "K02279",
    "K02280",
    "K02283",
    "K12510",
    "K12511",
}

SUCCINOGLYCAN_EPS = {
    "K16566",
    "K16567",
    "K16568",
}


# ============================================================
# CAPABILITY 4: DIRECT PLANT INTERACTION
# ============================================================

ACC_DEAMINASE = "K01505"


# ============================================================
# SUPPLEMENTARY IAA EVIDENCE
# ============================================================

IAA_M = "K00466"
IAA_H_SPECIFIC = "K21801"
IAA_BROAD_AMIDASE = "K01426"


# ============================================================
# METADATA
# ============================================================

METADATA_COLUMNS = [

    "Completeness",
    "Contamination",
    "Quality_class",

    "GTDB_domain",
    "GTDB_phylum",
    "GTDB_class",
    "GTDB_order",
    "GTDB_family",
    "GTDB_genus",
    "GTDB_species",

    "Site",
    "Treatment",
    "Compartment",
]


# ============================================================
# FUNCTIONS
# ============================================================

def extract_kos(value):

    if not value:
        return []

    return re.findall(
        r"K\d{5}",
        value
    )


def detect_count(
    mag_kos,
    expected
):

    return len(
        set(expected)
        & mag_kos
    )


def fraction_detected(
    mag_kos,
    expected
):

    expected = set(
        expected
    )

    if not expected:
        return 0.0

    return (
        len(
            expected
            & mag_kos
        )
        / len(
            expected
        )
    )


def complete(
    mag_kos,
    expected
):

    return int(
        set(
            expected
        ).issubset(
            mag_kos
        )
    )


def breadth_class(score):

    if score == 0:
        return "Absent"

    elif score <= 0.25:
        return "Low"

    elif score <= 0.50:
        return "Intermediate"

    else:
        return "High"


# ============================================================
# V4 iTOL DISPLAY
# ============================================================

HOST_COLOUR = "#3A923A"

TRIANGLE_SYMBOL = 4

LOW_SIZE = 0.8
INTERMEDIATE_SIZE = 2.2

HIGH_MIN_SIZE = 3.5
HIGH_MAX_SIZE = 5.0

ITOL_MAXIMUM_SIZE = 16
SYMBOL_SPACING = 5


def score_to_itol_size(score):

    if score <= 0:

        return None


    if score <= 0.25:

        return LOW_SIZE


    if score <= 0.50:

        return INTERMEDIATE_SIZE


    high_fraction = (
        score
        - 0.50
    ) / 0.50


    return (
        HIGH_MIN_SIZE
        +
        high_fraction
        * (
            HIGH_MAX_SIZE
            - HIGH_MIN_SIZE
        )
    )


# ============================================================
# CHECK INPUTS
# ============================================================

if not INPUT.exists():

    raise FileNotFoundError(
        f"Input not found:\n{INPUT}"
    )


if not TREE.exists():

    raise FileNotFoundError(
        f"Tree not found:\n{TREE}"
    )


# ============================================================
# READ GENE MASTER TABLE
# ============================================================

mag_kos = defaultdict(
    set
)

metadata = {}


print()

print(
    "=" * 72
)

print(
    "FINAL COLONISATION AND PLANT-INTERACTION STRATEGY"
)

print(
    "=" * 72
)

print()

print(
    f"Reading:\n  {INPUT}"
)

print()


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


    required = {
        "Final_MAG",
        "KO",
    }


    missing = (
        required
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


    rows_read = 0


    for row in reader:

        rows_read += 1


        mag = row[
            "Final_MAG"
        ].strip()


        if not mag:

            continue


        if mag not in metadata:

            metadata[
                mag
            ] = {

                column:
                    row.get(
                        column,
                        ""
                    ).strip()

                for column
                in METADATA_COLUMNS
            }


        for ko in extract_kos(
            row.get(
                "KO",
                ""
            )
        ):

            mag_kos[
                mag
            ].add(
                ko
            )


# Keep MAGs without KO calls.

for mag in metadata:

    if mag not in mag_kos:

        mag_kos[
            mag
        ] = set()


mags = sorted(
    metadata.keys()
)


N = len(
    mags
)


print(
    f"Gene rows read: "
    f"{rows_read:,}"
)

print(
    f"Unique MAGs: "
    f"{N}"
)

print()


if N != 373:

    print(
        "WARNING: expected 373 MAGs, "
        f"found {N}."
    )

    print()


# ============================================================
# BUILD FINAL CAPABILITIES
# ============================================================

final_rows = []

evidence_rows = []


for mag in mags:

    kos = mag_kos[
        mag
    ]


    # ========================================================
    # 1. CHEMOTAXIS
    # ========================================================

    chem_detected = sorted(
        CHEMOTAXIS_CORE
        & kos
    )

    chem_n = len(
        chem_detected
    )

    chem_fraction = (
        chem_n
        / len(
            CHEMOTAXIS_CORE
        )
    )


    # >=3 of 4

    chemotaxis = int(
        chem_n >= 3
    )


    chem_support_detected = sorted(
        CHEMOTAXIS_SUPPORT
        & kos
    )


    # ========================================================
    # 2A. FLAGELLAR MOTILITY
    # ========================================================

    flagellar_module_calls = {}

    flagellar_module_fractions = {}

    flagellar_module_detected = {}


    for (
        module_name,
        module_kos
    ) in FLAGELLAR_MODULES.items():

        detected = sorted(
            module_kos
            & kos
        )

        fraction = (
            len(
                detected
            )
            / len(
                module_kos
            )
        )


        module_complete = int(
            len(
                detected
            )
            == len(
                module_kos
            )
        )


        flagellar_module_calls[
            module_name
        ] = module_complete


        flagellar_module_fractions[
            module_name
        ] = fraction


        flagellar_module_detected[
            module_name
        ] = detected


    n_flagellar_modules_complete = sum(
        flagellar_module_calls.values()
    )


    flagellar_motility = int(
        n_flagellar_modules_complete
        >= 3
    )


    # ========================================================
    # 2B. TYPE-IV-PILUS / TWITCHING MOTILITY
    # ========================================================

    typeiv_detected = sorted(
        TYPEIV_ASSEMBLY
        & kos
    )

    typeiv_n = len(
        typeiv_detected
    )

    typeiv_fraction = (
        typeiv_n
        / len(
            TYPEIV_ASSEMBLY
        )
    )


    pilT_present = int(
        PILT in kos
    )

    pilU_present = int(
        PILU in kos
    )


    twitch_signal_detected = sorted(
        TWITCHING_SIGNAL
        & kos
    )


    # >=4/5 assembly markers + PilT

    twitching_motility = int(

        typeiv_n >= 4

        and

        pilT_present == 1
    )


    # ========================================================
    # FINAL MOTILITY CAPABILITY
    # ========================================================

    motility = int(

        flagellar_motility

        or

        twitching_motility
    )


    # ========================================================
    # 3. ATTACHMENT / BIOFILM-MATRIX
    # ========================================================

    # --------------------------------------------------------
    # PGA
    # >=3/4
    # --------------------------------------------------------

    pga_detected = sorted(
        PGA_MATRIX
        & kos
    )

    pga_n = len(
        pga_detected
    )

    pga_route = int(
        pga_n >= 3
    )


    # --------------------------------------------------------
    # Curli
    # complete 2/2
    # --------------------------------------------------------

    curli_detected = sorted(
        CURLI
        & kos
    )

    curli_route = int(
        len(
            curli_detected
        )
        == len(
            CURLI
        )
    )


    # --------------------------------------------------------
    # LapA
    # adhesin + complete transport
    # --------------------------------------------------------

    lapA_present = int(
        LAPA in kos
    )


    lap_transport_detected = sorted(
        LAPA_TRANSPORT
        & kos
    )


    lap_transport_complete = int(

        len(
            lap_transport_detected
        )

        ==

        len(
            LAPA_TRANSPORT
        )
    )


    lap_route = int(

        lapA_present

        and

        lap_transport_complete
    )


    # --------------------------------------------------------
    # Tad / Flp
    # >=6/7
    # --------------------------------------------------------

    tad_detected = sorted(
        TAD_FLP
        & kos
    )

    tad_n = len(
        tad_detected
    )


    tad_route = int(
        tad_n >= 6
    )


    # --------------------------------------------------------
    # Succinoglycan / EPS
    # complete 3/3
    # --------------------------------------------------------

    eps_detected = sorted(
        SUCCINOGLYCAN_EPS
        & kos
    )


    eps_route = int(

        len(
            eps_detected
        )

        ==

        len(
            SUCCINOGLYCAN_EPS
        )
    )


    # --------------------------------------------------------
    # Collapse alternative attachment routes
    # --------------------------------------------------------

    positive_attachment_routes = []


    if pga_route:

        positive_attachment_routes.append(
            "PGA"
        )


    if curli_route:

        positive_attachment_routes.append(
            "Curli"
        )


    if lap_route:

        positive_attachment_routes.append(
            "LapA"
        )


    if tad_route:

        positive_attachment_routes.append(
            "Tad_Flp"
        )


    if eps_route:

        positive_attachment_routes.append(
            "Succinoglycan_EPS"
        )


    attachment_matrix = int(
        len(
            positive_attachment_routes
        )
        > 0
    )


    # ========================================================
    # 4. DIRECT PLANT-INTERACTION POTENTIAL
    # ========================================================

    acc_deaminase = int(
        ACC_DEAMINASE
        in kos
    )


    direct_plant_interaction = (
        acc_deaminase
    )


    # ========================================================
    # SUPPLEMENTARY IAA EVIDENCE
    # ========================================================

    iaaM_present = int(
        IAA_M
        in kos
    )


    iaaH_specific_present = int(
        IAA_H_SPECIFIC
        in kos
    )


    iaa_IAM_complete = int(

        iaaM_present

        and

        iaaH_specific_present
    )


    broad_amidase_present = int(
        IAA_BROAD_AMIDASE
        in kos
    )


    # ========================================================
    # FINAL STRATEGY SCORE
    # ========================================================

    capability_count = (

        chemotaxis

        + motility

        + attachment_matrix

        + direct_plant_interaction
    )


    strategy_score = (
        capability_count
        / 4
    )


    strategy_class = breadth_class(
        strategy_score
    )


    # ========================================================
    # FINAL TABLE
    # ========================================================

    final_row = {

        "Final_MAG":
            mag,

        "Chemotaxis":
            chemotaxis,

        "Motility":
            motility,

        "Attachment_biofilm_matrix":
            attachment_matrix,

        "Direct_plant_interaction":
            direct_plant_interaction,

        "Colonisation_host_capability_count":
            capability_count,

        "Colonisation_host_strategy_score":
            f"{strategy_score:.4f}",

        "Colonisation_host_strategy_class":
            strategy_class,

        "Supplementary_IAM_IAA_complete":
            iaa_IAM_complete,
    }


    for column in METADATA_COLUMNS:

        final_row[
            column
        ] = metadata[
            mag
        ].get(
            column,
            ""
        )


    final_rows.append(
        final_row
    )


    # ========================================================
    # DETAILED EVIDENCE TABLE
    # ========================================================

    evidence_row = {

        "Final_MAG":
            mag,


        # ----------------------------------------------------
        # Chemotaxis
        # ----------------------------------------------------

        "Chemotaxis_core_detected_n":
            chem_n,

        "Chemotaxis_core_expected_n":
            len(
                CHEMOTAXIS_CORE
            ),

        "Chemotaxis_core_fraction":
            f"{chem_fraction:.4f}",

        "Chemotaxis_core_detected_KOs":
            ";".join(
                chem_detected
            ),

        "Chemotaxis_support_detected_KOs":
            ";".join(
                chem_support_detected
            ),

        "Chemotaxis_final":
            chemotaxis,


        # ----------------------------------------------------
        # Flagellar motility
        # ----------------------------------------------------

        "Flagellar_modules_complete_n":
            n_flagellar_modules_complete,

        "Flagellar_motility_final_route":
            flagellar_motility,


        # ----------------------------------------------------
        # Twitching motility
        # ----------------------------------------------------

        "TypeIV_assembly_detected_n":
            typeiv_n,

        "TypeIV_assembly_fraction":
            f"{typeiv_fraction:.4f}",

        "TypeIV_assembly_detected_KOs":
            ";".join(
                typeiv_detected
            ),

        "PilT_K02669":
            pilT_present,

        "PilU_K02670_support":
            pilU_present,

        "Twitching_signalling_detected_n":
            len(
                twitch_signal_detected
            ),

        "Twitching_signalling_detected_KOs":
            ";".join(
                twitch_signal_detected
            ),

        "Twitching_motility_final_route":
            twitching_motility,

        "Motility_final":
            motility,


        # ----------------------------------------------------
        # Attachment
        # ----------------------------------------------------

        "PGA_detected_n":
            pga_n,

        "PGA_detected_KOs":
            ";".join(
                pga_detected
            ),

        "PGA_route":
            pga_route,


        "Curli_detected_n":
            len(
                curli_detected
            ),

        "Curli_detected_KOs":
            ";".join(
                curli_detected
            ),

        "Curli_route":
            curli_route,


        "LapA_K12549":
            lapA_present,

        "Lap_transport_detected_n":
            len(
                lap_transport_detected
            ),

        "Lap_transport_detected_KOs":
            ";".join(
                lap_transport_detected
            ),

        "Lap_transport_complete":
            lap_transport_complete,

        "LapA_route":
            lap_route,


        "Tad_Flp_detected_n":
            tad_n,

        "Tad_Flp_detected_KOs":
            ";".join(
                tad_detected
            ),

        "Tad_Flp_route":
            tad_route,


        "Succinoglycan_EPS_detected_n":
            len(
                eps_detected
            ),

        "Succinoglycan_EPS_detected_KOs":
            ";".join(
                eps_detected
            ),

        "Succinoglycan_EPS_route":
            eps_route,


        "N_positive_attachment_routes":
            len(
                positive_attachment_routes
            ),

        "Positive_attachment_routes":
            ";".join(
                positive_attachment_routes
            ),

        "Attachment_biofilm_matrix_final":
            attachment_matrix,


        # ----------------------------------------------------
        # Plant interaction
        # ----------------------------------------------------

        "ACC_deaminase_K01505":
            acc_deaminase,

        "Direct_plant_interaction_final":
            direct_plant_interaction,


        # ----------------------------------------------------
        # Supplementary IAA evidence
        # ----------------------------------------------------

        "IaaM_K00466":
            iaaM_present,

        "Specific_IaaH_K21801":
            iaaH_specific_present,

        "IAM_IAA_complete":
            iaa_IAM_complete,

        "Broad_amidase_K01426_support":
            broad_amidase_present,


        # ----------------------------------------------------
        # Final strategy
        # ----------------------------------------------------

        "Colonisation_host_capability_count":
            capability_count,

        "Colonisation_host_strategy_score":
            f"{strategy_score:.4f}",
    }


    # Add detailed flagellar submodule evidence

    for module_name in FLAGELLAR_MODULES:

        evidence_row[
            f"Flagellar_{module_name}_complete"
        ] = flagellar_module_calls[
            module_name
        ]

        evidence_row[
            f"Flagellar_{module_name}_fraction"
        ] = (
            f"{flagellar_module_fractions[module_name]:.4f}"
        )

        evidence_row[
            f"Flagellar_{module_name}_detected_KOs"
        ] = ";".join(
            flagellar_module_detected[
                module_name
            ]
        )


    evidence_rows.append(
        evidence_row
    )


# ============================================================
# WRITE FINAL TABLE
# ============================================================

final_fields = list(
    final_rows[
        0
    ].keys()
)


with FINAL_TABLE.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=final_fields,
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        final_rows
    )


# ============================================================
# WRITE EVIDENCE TABLE
# ============================================================

evidence_fields = list(
    evidence_rows[
        0
    ].keys()
)


with EVIDENCE_TABLE.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=evidence_fields,
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        evidence_rows
    )


# ============================================================
# SUMMARY COUNTS
# ============================================================

def count_final(
    column
):

    return sum(

        int(
            row[
                column
            ]
        )

        for row
        in final_rows
    )


chemotaxis_n = count_final(
    "Chemotaxis"
)

motility_n = count_final(
    "Motility"
)

attachment_n = count_final(
    "Attachment_biofilm_matrix"
)

plant_n = count_final(
    "Direct_plant_interaction"
)


flagellar_route_n = sum(

    int(
        row[
            "Flagellar_motility_final_route"
        ]
    )

    for row
    in evidence_rows
)


twitching_route_n = sum(

    int(
        row[
            "Twitching_motility_final_route"
        ]
    )

    for row
    in evidence_rows
)


pga_n = sum(
    int(
        row[
            "PGA_route"
        ]
    )
    for row
    in evidence_rows
)

curli_n = sum(
    int(
        row[
            "Curli_route"
        ]
    )
    for row
    in evidence_rows
)

lap_n = sum(
    int(
        row[
            "LapA_route"
        ]
    )
    for row
    in evidence_rows
)

tad_n = sum(
    int(
        row[
            "Tad_Flp_route"
        ]
    )
    for row
    in evidence_rows
)

eps_n = sum(
    int(
        row[
            "Succinoglycan_EPS_route"
        ]
    )
    for row
    in evidence_rows
)


iaa_n = sum(
    int(
        row[
            "IAM_IAA_complete"
        ]
    )
    for row
    in evidence_rows
)


score_counter = Counter(

    row[
        "Colonisation_host_strategy_score"
    ]

    for row
    in final_rows
)


class_counter = Counter(

    row[
        "Colonisation_host_strategy_class"
    ]

    for row
    in final_rows
)


# ============================================================
# WRITE SUMMARY TABLE
# ============================================================

summary_rows = []


def add_summary(
    section,
    item,
    count,
    denominator=N
):

    percentage = (
        100
        * count
        / denominator
        if denominator
        else 0
    )


    summary_rows.append({

        "Section":
            section,

        "Item":
            item,

        "Count":
            count,

        "Denominator":
            denominator,

        "Percent":
            f"{percentage:.2f}",
    })


# Main capabilities

add_summary(
    "Final capability",
    "Chemotaxis",
    chemotaxis_n
)

add_summary(
    "Final capability",
    "Motility",
    motility_n
)

add_summary(
    "Final capability",
    "Attachment / biofilm-matrix",
    attachment_n
)

add_summary(
    "Final capability",
    "Direct plant interaction",
    plant_n
)


# Motility routes

add_summary(
    "Motility route",
    "Flagellar",
    flagellar_route_n
)

add_summary(
    "Motility route",
    "Type IV pilus / twitching",
    twitching_route_n
)


# Attachment routes

add_summary(
    "Attachment route",
    "PGA",
    pga_n
)

add_summary(
    "Attachment route",
    "Curli",
    curli_n
)

add_summary(
    "Attachment route",
    "LapA",
    lap_n
)

add_summary(
    "Attachment route",
    "Tad / Flp",
    tad_n
)

add_summary(
    "Attachment route",
    "Succinoglycan / EPS",
    eps_n
)


# Supplementary

add_summary(
    "Supplementary",
    "Complete IAM IAA route",
    iaa_n
)


# Score distribution

for score in [

    "0.0000",
    "0.2500",
    "0.5000",
    "0.7500",
    "1.0000",

]:

    add_summary(
        "Strategy score",
        score,
        score_counter[
            score
        ]
    )


# Visual classes

for category in [

    "Absent",
    "Low",
    "Intermediate",
    "High",

]:

    add_summary(
        "Visual breadth class",
        category,
        class_counter[
            category
        ]
    )


with SUMMARY_TABLE.open(
    "w",
    newline=""
) as outfile:

    fields = [

        "Section",
        "Item",
        "Count",
        "Denominator",
        "Percent",
    ]


    writer = csv.DictWriter(
        outfile,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        summary_rows
    )


# ============================================================
# READ 352-TIP BACTERIAL TREE
# ============================================================

newick = (
    TREE
    .read_text()
    .strip()
)


matches = re.findall(
    r"(?<=[(,])([^():,]+)(?=:)",
    newick
)


tree_tips = []

seen = set()


for item in matches:

    tip = (
        item
        .strip()
        .strip("'")
        .strip('"')
    )


    if (
        tip
        and tip not in seen
    ):

        tree_tips.append(
            tip
        )

        seen.add(
            tip
        )


score_by_mag = {

    row[
        "Final_MAG"
    ]:
        float(
            row[
                "Colonisation_host_strategy_score"
            ]
        )

    for row
    in final_rows
}


missing_tree_scores = [

    tip
    for tip in tree_tips
    if tip not in score_by_mag
]


if missing_tree_scores:

    raise ValueError(

        f"{len(missing_tree_scores)} "
        "tree MAGs missing Colonisation/Host scores."
    )


# ============================================================
# WRITE iTOL GREEN TRIANGLE RING
# ============================================================

with ITOL_FILE.open(
    "w"
) as outfile:

    outfile.write(
        "DATASET_SYMBOL\n"
    )

    outfile.write(
        "SEPARATOR TAB\n"
    )

    outfile.write(
        "DATASET_LABEL\t"
        "Colonisation + plant interaction\n"
    )

    outfile.write(
        f"COLOR\t"
        f"{HOST_COLOUR}\n"
    )

    outfile.write(
        f"MAXIMUM_SIZE\t"
        f"{ITOL_MAXIMUM_SIZE}\n"
    )

    outfile.write(
        f"SYMBOL_SPACING\t"
        f"{SYMBOL_SPACING}\n"
    )

    outfile.write(
        "SHOW_LABELS\t0\n"
    )

    outfile.write(
        "\nDATA\n"
    )


    for mag in tree_tips:

        score = score_by_mag[
            mag
        ]


        size = score_to_itol_size(
            score
        )


        # True absence -> no symbol.

        if size is None:

            continue


        category = breadth_class(
            score
        )


        tooltip = (

            "Colonisation + plant interaction: "
            f"{score:.3f} "
            f"({category} breadth)"
        )


        outfile.write(

            f"{mag}\t"

            f"{TRIANGLE_SYMBOL}\t"

            f"{size:.4f}\t"

            f"{HOST_COLOUR}\t"

            f"1\t"

            f"-1\t"

            f"{tooltip}\n"
        )


# ============================================================
# TERMINAL SUMMARY
# ============================================================

print(
    "=" * 72
)

print(
    "FINAL COLONISATION / HOST STRATEGY CREATED SUCCESSFULLY"
)

print(
    "=" * 72
)

print()


print(
    f"MAGs scored: {N}"
)

print(
    f"Bacterial tree tips: "
    f"{len(tree_tips)}"
)

print()


print(
    "FINAL CAPABILITIES"
)

print(
    "-" * 72
)


for label, count in [

    (
        "Chemotaxis",
        chemotaxis_n
    ),

    (
        "Motility",
        motility_n
    ),

    (
        "Attachment / biofilm-matrix",
        attachment_n
    ),

    (
        "Direct plant interaction",
        plant_n
    ),

]:

    print(

        f"{label}: "
        f"{count}/{N} "
        f"({100 * count / N:.2f}%)"
    )


print()

print(
    "MOTILITY ROUTES"
)

print(
    "-" * 72
)

print(
    f"Flagellar route: "
    f"{flagellar_route_n}/{N} "
    f"({100 * flagellar_route_n / N:.2f}%)"
)

print(
    f"Type-IV-pilus/twitching route: "
    f"{twitching_route_n}/{N} "
    f"({100 * twitching_route_n / N:.2f}%)"
)


print()

print(
    "ATTACHMENT ROUTES"
)

print(
    "-" * 72
)


for label, count in [

    (
        "PGA",
        pga_n
    ),

    (
        "Curli",
        curli_n
    ),

    (
        "LapA",
        lap_n
    ),

    (
        "Tad / Flp",
        tad_n
    ),

    (
        "Succinoglycan / EPS",
        eps_n
    ),

]:

    print(
        f"{label}: "
        f"{count}/{N} "
        f"({100 * count / N:.2f}%)"
    )


print()

print(
    "SUPPLEMENTARY"
)

print(
    "-" * 72
)

print(
    f"Complete IAM IAA route: "
    f"{iaa_n}/{N} "
    f"({100 * iaa_n / N:.2f}%)"
)


print()

print(
    "STRATEGY SCORE DISTRIBUTION"
)

print(
    "-" * 72
)


for score in [

    "0.0000",
    "0.2500",
    "0.5000",
    "0.7500",
    "1.0000",

]:

    count = score_counter[
        score
    ]

    print(

        f"{score}: "
        f"{count} "
        f"({100 * count / N:.2f}%)"
    )


print()

print(
    "VISUAL BREADTH CLASSES"
)

print(
    "-" * 72
)


for category in [

    "Absent",
    "Low",
    "Intermediate",
    "High",

]:

    count = class_counter[
        category
    ]

    print(

        f"{category}: "
        f"{count} "
        f"({100 * count / N:.2f}%)"
    )


print()

print(
    "FINAL RULE"
)

print(
    "-" * 72
)

print(
    "Colonisation/Host score = "
    "(chemotaxis + motility + "
    "attachment/biofilm-matrix + "
    "direct plant interaction) / 4"
)

print()

print(
    "Chemotaxis:"
)

print(
    "  >=3/4 MCP-CheA-CheW-CheY core markers"
)

print()

print(
    "Motility:"
)

print(
    "  >=3/4 complete flagellar submodules"
)

print(
    "  OR"
)

print(
    "  >=4/5 Type-IV-pilus assembly markers + PilT"
)

print()

print(
    "Attachment / biofilm-matrix:"
)

print(
    "  PGA >=3/4"
)

print(
    "  OR Curli complete"
)

print(
    "  OR LapA + complete transport"
)

print(
    "  OR Tad/Flp >=6/7"
)

print(
    "  OR succinoglycan/EPS complete"
)

print()

print(
    "Direct plant interaction:"
)

print(
    "  ACC deaminase K01505"
)

print()

print(
    "OUTPUTS"
)

print(
    "-" * 72
)

print(
    f"Final strategy:\n"
    f"  {FINAL_TABLE}"
)

print()

print(
    f"Module evidence:\n"
    f"  {EVIDENCE_TABLE}"
)

print()

print(
    f"Summary:\n"
    f"  {SUMMARY_TABLE}"
)

print()

print(
    f"iTOL ring:\n"
    f"  {ITOL_FILE}"
)

print()

print(
    "iTOL:"
)

print(
    "  triangle = Host interaction"
)

print(
    f"  colour = {HOST_COLOUR}"
)

print(
    "  symbol size = standardized capability breadth"
)

print(
    "=" * 72
)
