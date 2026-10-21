#!/usr/bin/env python3

"""
Build final Iron Acquisition ecological strategy.

Framework
---------
genes -> functional evidence -> binary capabilities
      -> standardized ecological capability breadth

FINAL IRON CAPABILITIES
=======================

1. Ferrous iron uptake
   Main diagnostic:
       FeoB K04759

   FeoA K04758 is retained as supporting evidence but is not
   required for the main binary capability call.

2. Direct ferric iron uptake
   Complete Afu/Fbp Fe(III) ABC transporter:
       K02010 + K02011 + K02012

3. Siderophore biosynthesis
   Presence of a sufficiently complete recognized siderophore
   biosynthetic route.

   Rule:
       - routes containing <=3 markers require ALL markers
       - routes containing >=4 markers require >=75% of markers
         AND at least 3 markers

   Alternative routes collapse into ONE capability.

4. Iron-complex uptake
   Presence of either:

       A. complete iron-complex ABC transporter:
          K02013 + K02015 + K02016

       OR

       B. TonB-dependent iron-complex uptake:
          K02014 + K03832

   These alternatives collapse into ONE capability.

SUPPLEMENTARY
=============

K16087 haem-associated iron receptor is retained separately and
does NOT contribute to the four-capability Iron strategy score.

FINAL SCORE
===========

Iron_strategy_score =
    (
        Ferrous_iron_uptake
        + Direct_ferric_iron_uptake
        + Siderophore_biosynthesis
        + Iron_complex_uptake
    ) / 4

Therefore possible scores are:

0.00
0.25
0.50
0.75
1.00

The score represents standardized ecological capability breadth.
It does NOT represent percentage physiological capacity.

Outputs
-------
results/ecological_strategies/iron_acquisition/tables/

    iron_strategy_final_373_MAGs.tsv
    iron_module_evidence_373_MAGs.tsv
    iron_strategy_final_summary.tsv

results/ecological_strategies/iron_acquisition/itol/

    itol_iron_strategy_score_352.txt
"""

import csv
import math
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
      "iron_acquisition/tables"
)

ITOL_DIR = (
    BASE
    / "results/ecological_strategies/"
      "iron_acquisition/itol"
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
    / "iron_strategy_final_373_MAGs.tsv"
)

EVIDENCE_TABLE = (
    TABLE_DIR
    / "iron_module_evidence_373_MAGs.tsv"
)

SUMMARY_TABLE = (
    TABLE_DIR
    / "iron_strategy_final_summary.tsv"
)

ITOL_FILE = (
    ITOL_DIR
    / "itol_iron_strategy_score_352.txt"
)


# ============================================================
# FINAL CAPABILITY RULES
# ============================================================

# ------------------------------------------------------------
# Capability 1: Ferrous iron uptake
# ------------------------------------------------------------

FEOA = "K04758"
FEOB = "K04759"


# ------------------------------------------------------------
# Capability 2: Direct ferric iron uptake
# ------------------------------------------------------------

AFU_FBP = {
    "K02010",
    "K02011",
    "K02012",
}


# ------------------------------------------------------------
# Capability 4: Iron-complex uptake
# ------------------------------------------------------------

IRON_COMPLEX_ABC = {
    "K02013",
    "K02015",
    "K02016",
}

TONB_IRON_COMPLEX = {
    "K02014",
    "K03832",
}


# ------------------------------------------------------------
# Supplementary haem-associated evidence
# ------------------------------------------------------------

HEME_RECEPTOR = "K16087"


# ============================================================
# SIDEROPHORE BIOSYNTHESIS ROUTES
# ============================================================
#
# These are the same candidate routes audited previously.
# Alternative pathways are collapsed into one final
# siderophore-biosynthesis capability.
#
# ============================================================

SIDEROPHORE_ROUTES = {

    "Enterobactin": [
        "K02361",
        "K01252",
        "K00216",
        "K02363",
        "K02364",
    ],

    "Aerobactin_M00918": [
        "K03897",
        "K03896",
        "K03894",
        "K03895",
    ],

    "Pyochelin": [
        "K04781",
        "K04782",
        "K12238",
        "K12239",
        "K12240",
        "K12241",
    ],

    "Yersiniabactin": [
        "K04781",
        "K04783",
        "K04784",
        "K04785",
        "K04786",
    ],

    "Mycobactin": [
        "K04781",
        "K04787",
        "K04788",
        "K04789",
        "K04790",
        "K04791",
        "K04792",
        "K04793",
    ],

    "Staphyloferrin_B_M00875": [
        "K23371",
        "K21949",
        "K21721",
        "K23372",
        "K23373",
        "K23374",
        "K23375",
    ],

    "Staphyloferrin_A_M00876": [
        "K21898",
        "K23446",
        "K23447",
    ],

    "Deferoxamine_M01053": [
        "K28682",
        "K27501",
        "K28683",
        "K28684",
    ],

    "Schizokinen_M01054": [
        "K00836",
        "K13745",
        "K27501",
        "K28694",
        "K28695",
        "K28696",
    ],

    "Putrebactin_M01055": [
        "K27501",
        "K28697",
        "K28698",
    ],

    "Alcaligin_M01056": [
        "K27501",
        "K28708",
        "K28699",
        "K28700",
    ],

    "Rhizoferrin_M01057": [
        "K28702",
        "K28703",
    ],

    "Achromobactin_M01058": [
        "K28704",
        "K28705",
        "K23374",
        "K23375",
    ],

    "Vibrioferrin_M01059": [
        "K28704",
        "K28705",
        "K28706",
        "K28707",
    ],

    "Petrobactin_M01060": [
        "K15652",
        "K24110",
        "K24111",
        "K24108",
        "K24109",
        "K24112",
    ],
}


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


def all_present(
    mag_kos,
    required
):

    return set(
        required
    ).issubset(
        mag_kos
    )


def siderophore_route_evidence(
    mag_kos,
    route_kos
):

    """
    Evaluate one siderophore route.

    For routes with <=3 markers:
        require all markers.

    For routes with >=4 markers:
        require >=75% and at least 3 markers.
    """

    route_set = set(
        route_kos
    )

    detected = (
        route_set
        & mag_kos
    )

    n_expected = len(
        route_set
    )

    n_detected = len(
        detected
    )

    fraction = (
        n_detected
        / n_expected
        if n_expected
        else 0
    )


    if n_expected <= 3:

        required_n = n_expected

    else:

        required_n = max(
            3,
            math.ceil(
                0.75
                * n_expected
            )
        )


    call = int(
        n_detected
        >= required_n
    )


    return {

        "call":
            call,

        "n_detected":
            n_detected,

        "n_expected":
            n_expected,

        "required_n":
            required_n,

        "fraction":
            fraction,

        "detected":
            sorted(
                detected
            ),

        "missing":
            sorted(
                route_set
                - detected
            ),
    }


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

IRON_COLOUR = "#8E5EA2"

SQUARE_SYMBOL = 1

LOW_SIZE = 0.8

INTERMEDIATE_SIZE = 2.2

HIGH_MIN_SIZE = 3.5

HIGH_MAX_SIZE = 5.0

ITOL_MAXIMUM_SIZE = 16

SYMBOL_SPACING = 5


def score_to_itol_size(
    score
):

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
# CHECK INPUT
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
# READ GENE-LEVEL MASTER TABLE
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
    "FINAL IRON ACQUISITION STRATEGY"
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


    required_columns = {
        "Final_MAG",
        "KO",
    }


    missing_columns = (
        required_columns
        - set(
            reader.fieldnames
        )
    )


    if missing_columns:

        raise ValueError(
            "Missing required columns: "
            + ", ".join(
                sorted(
                    missing_columns
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


# Ensure MAGs without detected KOs remain represented.

for mag in metadata:

    if mag not in mag_kos:

        mag_kos[
            mag
        ] = set()


mags = sorted(
    metadata.keys()
)


print(
    f"Gene rows read: "
    f"{rows_read:,}"
)

print(
    f"Unique MAGs: "
    f"{len(mags)}"
)

print()


if len(mags) != 373:

    print(
        "WARNING: expected 373 MAGs, "
        f"found {len(mags)}."
    )

    print()


# ============================================================
# BUILD FINAL CAPABILITY CALLS
# ============================================================

final_rows = []

evidence_rows = []


for mag in mags:

    kos = mag_kos[
        mag
    ]


    # --------------------------------------------------------
    # 1. FERROUS IRON UPTAKE
    # --------------------------------------------------------

    feoA = int(
        FEOA in kos
    )

    feoB = int(
        FEOB in kos
    )


    ferrous_iron_uptake = int(
        feoB == 1
    )


    feoAB_strict = int(
        feoA == 1
        and feoB == 1
    )


    # --------------------------------------------------------
    # 2. DIRECT FERRIC IRON UPTAKE
    # --------------------------------------------------------

    afu_detected = sorted(
        AFU_FBP
        & kos
    )


    afu_fraction = (
        len(
            afu_detected
        )
        / len(
            AFU_FBP
        )
    )


    direct_ferric_iron_uptake = int(
        all_present(
            kos,
            AFU_FBP
        )
    )


    # --------------------------------------------------------
    # 3. SIDEROPHORE BIOSYNTHESIS
    # --------------------------------------------------------

    siderophore_positive_routes = []

    best_route = ""

    best_route_fraction = 0.0

    best_route_n = 0

    route_evidence = {}


    for (
        route_name,
        route_kos
    ) in SIDEROPHORE_ROUTES.items():

        evidence = (
            siderophore_route_evidence(
                kos,
                route_kos
            )
        )


        route_evidence[
            route_name
        ] = evidence


        if evidence[
            "call"
        ]:

            siderophore_positive_routes.append(
                route_name
            )


        if (
            evidence[
                "fraction"
            ]
            > best_route_fraction
        ):

            best_route = (
                route_name
            )

            best_route_fraction = (
                evidence[
                    "fraction"
                ]
            )

            best_route_n = (
                evidence[
                    "n_detected"
                ]
            )


    siderophore_biosynthesis = int(
        len(
            siderophore_positive_routes
        )
        > 0
    )


    # --------------------------------------------------------
    # 4. IRON-COMPLEX UPTAKE
    # --------------------------------------------------------

    complex_abc_detected = sorted(
        IRON_COMPLEX_ABC
        & kos
    )


    complex_abc_fraction = (
        len(
            complex_abc_detected
        )
        / len(
            IRON_COMPLEX_ABC
        )
    )


    iron_complex_abc = int(
        all_present(
            kos,
            IRON_COMPLEX_ABC
        )
    )


    tonb_detected = sorted(
        TONB_IRON_COMPLEX
        & kos
    )


    tonb_fraction = (
        len(
            tonb_detected
        )
        / len(
            TONB_IRON_COMPLEX
        )
    )


    tonb_iron_complex = int(
        all_present(
            kos,
            TONB_IRON_COMPLEX
        )
    )


    iron_complex_uptake = int(
        iron_complex_abc
        or tonb_iron_complex
    )


    # --------------------------------------------------------
    # Supplementary haem-associated acquisition
    # --------------------------------------------------------

    haem_receptor = int(
        HEME_RECEPTOR
        in kos
    )


    # --------------------------------------------------------
    # FINAL IRON SCORE
    # --------------------------------------------------------

    capability_count = (

        ferrous_iron_uptake

        + direct_ferric_iron_uptake

        + siderophore_biosynthesis

        + iron_complex_uptake
    )


    iron_strategy_score = (
        capability_count
        / 4
    )


    strategy_class = breadth_class(
        iron_strategy_score
    )


    # --------------------------------------------------------
    # FINAL STRATEGY TABLE
    # --------------------------------------------------------

    final_row = {

        "Final_MAG":
            mag,

        "Ferrous_iron_uptake":
            ferrous_iron_uptake,

        "Direct_ferric_iron_uptake":
            direct_ferric_iron_uptake,

        "Siderophore_biosynthesis":
            siderophore_biosynthesis,

        "Iron_complex_uptake":
            iron_complex_uptake,

        "Iron_capability_count":
            capability_count,

        "Iron_strategy_score":
            f"{iron_strategy_score:.4f}",

        "Iron_strategy_class":
            strategy_class,

        "Supplementary_haem_receptor":
            haem_receptor,
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


    # --------------------------------------------------------
    # DETAILED EVIDENCE TABLE
    # --------------------------------------------------------

    evidence_row = {

        "Final_MAG":
            mag,

        # Ferrous iron
        "FeoA_K04758":
            feoA,

        "FeoB_K04759":
            feoB,

        "FeoAB_strict":
            feoAB_strict,

        "Ferrous_iron_uptake_final":
            ferrous_iron_uptake,

        # Direct ferric iron
        "Afu_Fbp_detected_n":
            len(
                afu_detected
            ),

        "Afu_Fbp_expected_n":
            len(
                AFU_FBP
            ),

        "Afu_Fbp_fraction":
            f"{afu_fraction:.4f}",

        "Afu_Fbp_detected_KOs":
            ";".join(
                afu_detected
            ),

        "Direct_ferric_iron_uptake_final":
            direct_ferric_iron_uptake,

        # Siderophore
        "N_positive_siderophore_routes":
            len(
                siderophore_positive_routes
            ),

        "Positive_siderophore_routes":
            ";".join(
                siderophore_positive_routes
            ),

        "Best_siderophore_route":
            best_route,

        "Best_siderophore_detected_n":
            best_route_n,

        "Best_siderophore_fraction":
            f"{best_route_fraction:.4f}",

        "Siderophore_biosynthesis_final":
            siderophore_biosynthesis,

        # Iron-complex ABC
        "Iron_complex_ABC_detected_n":
            len(
                complex_abc_detected
            ),

        "Iron_complex_ABC_fraction":
            f"{complex_abc_fraction:.4f}",

        "Iron_complex_ABC_detected_KOs":
            ";".join(
                complex_abc_detected
            ),

        "Iron_complex_ABC_complete":
            iron_complex_abc,

        # TonB-dependent
        "TonB_iron_complex_detected_n":
            len(
                tonb_detected
            ),

        "TonB_iron_complex_fraction":
            f"{tonb_fraction:.4f}",

        "TonB_iron_complex_detected_KOs":
            ";".join(
                tonb_detected
            ),

        "TonB_iron_complex_complete":
            tonb_iron_complex,

        "Iron_complex_uptake_final":
            iron_complex_uptake,

        # Supplementary
        "Haem_receptor_K16087":
            haem_receptor,

        # Final score
        "Iron_capability_count":
            capability_count,

        "Iron_strategy_score":
            f"{iron_strategy_score:.4f}",
    }


    # --------------------------------------------------------
    # Add route-specific siderophore evidence
    # --------------------------------------------------------

    for (
        route_name,
        evidence
    ) in route_evidence.items():

        evidence_row[
            f"{route_name}_detected_n"
        ] = evidence[
            "n_detected"
        ]

        evidence_row[
            f"{route_name}_expected_n"
        ] = evidence[
            "n_expected"
        ]

        evidence_row[
            f"{route_name}_required_n"
        ] = evidence[
            "required_n"
        ]

        evidence_row[
            f"{route_name}_fraction"
        ] = (
            f"{evidence['fraction']:.4f}"
        )

        evidence_row[
            f"{route_name}_call"
        ] = evidence[
            "call"
        ]

        evidence_row[
            f"{route_name}_detected_KOs"
        ] = ";".join(
            evidence[
                "detected"
            ]
        )


    evidence_rows.append(
        evidence_row
    )


# ============================================================
# WRITE FINAL STRATEGY TABLE
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
# FINAL SUMMARY COUNTS
# ============================================================

N = len(
    final_rows
)


def capability_count(
    column
):

    return sum(
        int(
            row[
                column
            ]
        )
        for row in final_rows
    )


ferrous_n = capability_count(
    "Ferrous_iron_uptake"
)

ferric_n = capability_count(
    "Direct_ferric_iron_uptake"
)

siderophore_n = capability_count(
    "Siderophore_biosynthesis"
)

complex_n = capability_count(
    "Iron_complex_uptake"
)

heme_n = capability_count(
    "Supplementary_haem_receptor"
)


score_counter = Counter(

    row[
        "Iron_strategy_score"
    ]

    for row in final_rows
)


class_counter = Counter(

    row[
        "Iron_strategy_class"
    ]

    for row in final_rows
)


# ============================================================
# SIDEROPHORE ROUTE POSITIVE COUNTS
# ============================================================

route_positive_counts = Counter()


for row in evidence_rows:

    for route_name in SIDEROPHORE_ROUTES:

        if int(
            row[
                f"{route_name}_call"
            ]
        ) == 1:

            route_positive_counts[
                route_name
            ] += 1


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


add_summary(
    "Final capability",
    "Ferrous iron uptake",
    ferrous_n
)

add_summary(
    "Final capability",
    "Direct ferric iron uptake",
    ferric_n
)

add_summary(
    "Final capability",
    "Siderophore biosynthesis",
    siderophore_n
)

add_summary(
    "Final capability",
    "Iron-complex uptake",
    complex_n
)

add_summary(
    "Supplementary",
    "Haem-associated receptor K16087",
    heme_n
)


for route_name in SIDEROPHORE_ROUTES:

    add_summary(
        "Siderophore route",
        route_name,
        route_positive_counts[
            route_name
        ]
    )


for score in [
    "0.0000",
    "0.2500",
    "0.5000",
    "0.7500",
    "1.0000",
]:

    add_summary(
        "Iron strategy score",
        score,
        score_counter[
            score
        ]
    )


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
# READ 352 BACTERIAL TREE TIPS
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
                "Iron_strategy_score"
            ]
        )

    for row in final_rows
}


missing_tree_scores = [

    tip
    for tip in tree_tips
    if tip not in score_by_mag
]


if missing_tree_scores:

    raise ValueError(
        f"{len(missing_tree_scores)} "
        "tree MAGs missing Iron scores."
    )


# ============================================================
# WRITE iTOL IRON RING
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
        "Iron acquisition\n"
    )

    outfile.write(
        f"COLOR\t"
        f"{IRON_COLOUR}\n"
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


        # score = 0 -> true absence -> no symbol

        if size is None:
            continue


        category = breadth_class(
            score
        )


        tooltip = (
            f"Iron acquisition: "
            f"{score:.3f} "
            f"({category} breadth)"
        )


        outfile.write(
            f"{mag}\t"
            f"{SQUARE_SYMBOL}\t"
            f"{size:.4f}\t"
            f"{IRON_COLOUR}\t"
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
    "FINAL IRON STRATEGY CREATED SUCCESSFULLY"
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


for (
    label,
    count
) in [

    (
        "Ferrous iron uptake",
        ferrous_n
    ),

    (
        "Direct ferric iron uptake",
        ferric_n
    ),

    (
        "Siderophore biosynthesis",
        siderophore_n
    ),

    (
        "Iron-complex uptake",
        complex_n
    ),

]:

    print(
        f"{label}: "
        f"{count}/{N} "
        f"({100 * count / N:.2f}%)"
    )


print()

print(
    "SIDEROPHORE ROUTES PASSING FINAL RULE"
)

print(
    "-" * 72
)


for route_name in SIDEROPHORE_ROUTES:

    print(
        f"{route_name}: "
        f"{route_positive_counts[route_name]}"
    )


print()

print(
    "SUPPLEMENTARY"
)

print(
    "-" * 72
)

print(
    f"Haem-associated receptor K16087: "
    f"{heme_n}/{N} "
    f"({100 * heme_n / N:.2f}%)"
)

print()


print(
    "IRON STRATEGY SCORE DISTRIBUTION"
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
    "Iron score = "
    "(Fe2+ uptake + direct Fe3+ uptake + "
    "siderophore biosynthesis + "
    "iron-complex uptake) / 4"
)

print()

print(
    "Fe2+ uptake:"
)

print(
    "  FeoB K04759"
)

print()

print(
    "Direct Fe3+ uptake:"
)

print(
    "  K02010 + K02011 + K02012"
)

print()

print(
    "Siderophore biosynthesis:"
)

print(
    "  >=75% of a recognized route for routes "
    "with >=4 markers;"
)

print(
    "  complete route required for routes "
    "with <=3 markers."
)

print()

print(
    "Iron-complex uptake:"
)

print(
    "  K02013 + K02015 + K02016"
)

print(
    "  OR"
)

print(
    "  K02014 + K03832"
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
    f"iTOL Iron ring:\n"
    f"  {ITOL_FILE}"
)

print()

print(
    "iTOL:"
)

print(
    f"  square = Resource acquisition"
)

print(
    f"  colour = {IRON_COLOUR}"
)

print(
    "  symbol size = standardized Iron capability breadth"
)

print(
    "=" * 72
)
