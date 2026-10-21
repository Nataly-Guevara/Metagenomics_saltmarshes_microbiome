#!/usr/bin/env python3

"""
Audit candidate genomic markers for the Iron Acquisition ecological strategy.

IMPORTANT
---------
This script DOES NOT create the final iron strategy score.

Its purpose is to examine, across all 373 MAGs:

1. Ferrous iron uptake evidence
2. Direct ferric iron uptake evidence
3. Iron-complex / siderophore uptake evidence
4. Siderophore biosynthesis evidence across multiple recognised routes
5. Supplementary haem-iron evidence

Outputs are intended to help define biologically defensible binary
capability rules before building the final Iron strategy.

Input
-----
Gene-level KOfam master table:
results/functional_annotation/kofamscan/
final_373_MAGs_gene_master_table.tsv

Outputs
-------
results/ecological_strategies/iron_acquisition/audit/

    iron_candidate_KO_prevalence.tsv
    iron_candidate_route_summary.tsv
    iron_strategy_audit_373_MAGs.tsv
    iron_strategy_audit_summary.txt
"""

import csv
import re
from collections import defaultdict
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

OUTPUT_DIR = (
    BASE
    / "results/ecological_strategies/"
      "iron_acquisition/audit"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

KO_PREVALENCE_OUT = (
    OUTPUT_DIR
    / "iron_candidate_KO_prevalence.tsv"
)

ROUTE_SUMMARY_OUT = (
    OUTPUT_DIR
    / "iron_candidate_route_summary.tsv"
)

MAG_AUDIT_OUT = (
    OUTPUT_DIR
    / "iron_strategy_audit_373_MAGs.tsv"
)

SUMMARY_OUT = (
    OUTPUT_DIR
    / "iron_strategy_audit_summary.txt"
)


# ============================================================
# CORE IRON ACQUISITION SYSTEMS TO AUDIT
# ============================================================

CORE_SYSTEMS = {

    # --------------------------------------------------------
    # Ferrous iron uptake
    # --------------------------------------------------------

    "FeoAB_ferrous_uptake": {
        "category": "Ferrous iron uptake",
        "kos": [
            "K04758",  # FeoA
            "K04759",  # FeoB
        ],
    },


    # --------------------------------------------------------
    # Direct ferric iron ABC uptake
    # afuC/fbpC + afuB/fbpB + afuA/fbpA
    # --------------------------------------------------------

    "Afu_Fbp_ferric_ABC": {
        "category": "Direct ferric iron uptake",
        "kos": [
            "K02010",
            "K02011",
            "K02012",
        ],
    },


    # --------------------------------------------------------
    # Iron-complex ABC transporter
    # Not assumed here to be exclusively siderophore-specific.
    # --------------------------------------------------------

    "Iron_complex_ABC": {
        "category": "Iron-complex uptake",
        "kos": [
            "K02013",  # ATP-binding
            "K02015",  # permease
            "K02016",  # substrate-binding
        ],
    },


    # --------------------------------------------------------
    # TonB-dependent outer-membrane iron-complex uptake
    #
    # K02014 = iron-complex outer membrane receptor
    # K03832 = TonB
    #
    # We audit these together but do NOT yet make a final
    # capability call.
    # --------------------------------------------------------

    "TonB_iron_complex_receptor": {
        "category": "Iron-complex uptake",
        "kos": [
            "K02014",
            "K03832",
        ],
    },


    # --------------------------------------------------------
    # Expanded TonB energy-transduction evidence
    #
    # ExbB + TonB.
    #
    # This is supportive evidence only because TonB systems
    # can transport substrates other than siderophores.
    # --------------------------------------------------------

    "TonB_ExbB_support": {
        "category": "Iron-complex uptake supporting evidence",
        "kos": [
            "K03561",  # ExbB
            "K03832",  # TonB
        ],
    },
}


# ============================================================
# SIDEROPHORE BIOSYNTHESIS ROUTES
# ============================================================
#
# These are AUDIT marker sets.
#
# For every route we report:
#
#   number of expected KOs detected
#   fraction of marker set detected
#   whether ALL listed KOs occur in the same MAG
#
# We do NOT yet assume that "all markers present" is the final
# criterion for every biological route.
#
# ============================================================

SIDEROPHORE_ROUTES = {

    # --------------------------------------------------------
    # Enterobactin-like NRPS siderophore pathway
    # --------------------------------------------------------

    "Enterobactin": [
        "K02361",  # EntC
        "K01252",  # EntB
        "K00216",  # EntA
        "K02363",  # EntE
        "K02364",  # EntF
    ],


    # --------------------------------------------------------
    # Aerobactin
    # KEGG M00918
    # --------------------------------------------------------

    "Aerobactin_M00918": [
        "K03897",
        "K03896",
        "K03894",
        "K03895",
    ],


    # --------------------------------------------------------
    # Pyochelin
    # --------------------------------------------------------

    "Pyochelin": [
        "K04781",
        "K04782",
        "K12238",
        "K12239",
        "K12240",
        "K12241",
    ],


    # --------------------------------------------------------
    # Yersiniabactin
    # --------------------------------------------------------

    "Yersiniabactin": [
        "K04781",
        "K04783",
        "K04784",
        "K04785",
        "K04786",
    ],


    # --------------------------------------------------------
    # Mycobactin
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Staphyloferrin B
    # KEGG M00875
    # --------------------------------------------------------

    "Staphyloferrin_B_M00875": [
        "K23371",
        "K21949",
        "K21721",
        "K23372",
        "K23373",
        "K23374",
        "K23375",
    ],


    # --------------------------------------------------------
    # Staphyloferrin A
    # KEGG M00876
    # --------------------------------------------------------

    "Staphyloferrin_A_M00876": [
        "K21898",
        "K23446",
        "K23447",
    ],


    # --------------------------------------------------------
    # Deferoxamine
    # KEGG M01053
    # --------------------------------------------------------

    "Deferoxamine_M01053": [
        "K28682",
        "K27501",
        "K28683",
        "K28684",
    ],


    # --------------------------------------------------------
    # Schizokinen
    # KEGG M01054
    # --------------------------------------------------------

    "Schizokinen_M01054": [
        "K00836",
        "K13745",
        "K27501",
        "K28694",
        "K28695",
        "K28696",
    ],


    # --------------------------------------------------------
    # Putrebactin
    # KEGG M01055
    # --------------------------------------------------------

    "Putrebactin_M01055": [
        "K27501",
        "K28697",
        "K28698",
    ],


    # --------------------------------------------------------
    # Alcaligin
    # KEGG M01056
    # --------------------------------------------------------

    "Alcaligin_M01056": [
        "K27501",
        "K28708",
        "K28699",
        "K28700",
    ],


    # --------------------------------------------------------
    # Rhizoferrin
    # KEGG M01057
    # --------------------------------------------------------

    "Rhizoferrin_M01057": [
        "K28702",
        "K28703",
    ],


    # --------------------------------------------------------
    # Achromobactin
    # KEGG M01058
    # --------------------------------------------------------

    "Achromobactin_M01058": [
        "K28704",
        "K28705",
        "K23374",
        "K23375",
    ],


    # --------------------------------------------------------
    # Vibrioferrin
    # KEGG M01059
    # --------------------------------------------------------

    "Vibrioferrin_M01059": [
        "K28704",
        "K28705",
        "K28706",
        "K28707",
    ],


    # --------------------------------------------------------
    # Petrobactin
    # KEGG M01060
    # --------------------------------------------------------

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
# SUPPLEMENTARY IRON MARKERS
# ============================================================

SUPPLEMENTARY_MARKERS = {

    # Haem receptor retained separately from the main strategy
    # until its ecological interpretation is decided.

    "K16087": "Haem/haemoglobin-associated iron receptor",

}


# ============================================================
# BUILD COMPLETE CANDIDATE KO SET
# ============================================================

candidate_kos = set()


for system in CORE_SYSTEMS.values():

    candidate_kos.update(
        system["kos"]
    )


for kos in SIDEROPHORE_ROUTES.values():

    candidate_kos.update(
        kos
    )


candidate_kos.update(
    SUPPLEMENTARY_MARKERS.keys()
)


# ============================================================
# HUMAN-READABLE KO LABELS
# ============================================================

KO_LABELS = {

    "K04758": "FeoA ferrous iron transport protein A",
    "K04759": "FeoB ferrous iron transport protein B",

    "K02010": "AfuC/FbpC Fe(III) transporter ATP-binding protein",
    "K02011": "AfuB/FbpB Fe(III) transporter permease",
    "K02012": "AfuA/FbpA Fe(III) transporter substrate-binding protein",

    "K02013": "Iron-complex transporter ATP-binding protein",
    "K02015": "Iron-complex transporter permease",
    "K02016": "Iron-complex transporter substrate-binding protein",
    "K02014": "Iron-complex outer-membrane receptor",

    "K03832": "TonB",
    "K03561": "ExbB",

    "K16087": "Haem/haemoglobin-associated iron receptor",

    "K02361": "EntC",
    "K01252": "EntB",
    "K00216": "EntA",
    "K02363": "EntE",
    "K02364": "EntF",

    "K03897": "Aerobactin biosynthesis marker",
    "K03896": "Aerobactin biosynthesis marker",
    "K03894": "Aerobactin biosynthesis marker",
    "K03895": "Aerobactin biosynthesis marker",

    "K04781": "Salicylate synthetase family siderophore marker",
    "K04782": "PchB",
    "K12238": "PchD",
    "K12239": "PchE",
    "K12240": "PchF",
    "K12241": "PchG",

    "K04783": "Yersiniabactin YbtE/Irp5",
    "K04784": "Yersiniabactin Irp2/HMWP2",
    "K04785": "Yersiniabactin YbtU/Irp3",
    "K04786": "Yersiniabactin Irp1/HMWP1",

    "K04787": "Mycobactin MbtA",
    "K04788": "Mycobactin MbtB",
    "K04789": "Mycobactin MbtE",
    "K04790": "Mycobactin MbtC",
    "K04791": "Mycobactin MbtD",
    "K04792": "Mycobactin MbtF",
    "K04793": "Mycobactin MbtG",
}


# ============================================================
# FUNCTIONS
# ============================================================

def extract_kos(value):

    """
    Extract KEGG K numbers robustly from a KO field.
    """

    if not value:

        return []


    return re.findall(
        r"K\d{5}",
        value
    )


def route_stats(
    mag_kos,
    expected_kos
):

    """
    Return evidence statistics for one marker set in one MAG.
    """

    expected = set(
        expected_kos
    )

    detected = (
        expected
        & mag_kos
    )

    n_expected = len(
        expected
    )

    n_detected = len(
        detected
    )

    fraction = (
        n_detected
        / n_expected
        if n_expected > 0
        else 0
    )

    complete = int(
        n_detected
        == n_expected
        and n_expected > 0
    )

    return {
        "n_detected": n_detected,
        "n_expected": n_expected,
        "fraction": fraction,
        "complete": complete,
        "detected_kos": sorted(
            detected
        ),
        "missing_kos": sorted(
            expected
            - detected
        ),
    }


# ============================================================
# CHECK INPUT
# ============================================================

if not INPUT.exists():

    raise FileNotFoundError(
        f"Input file not found:\n{INPUT}"
    )


# ============================================================
# READ MASTER TABLE
# ============================================================

mag_kos = defaultdict(
    set
)

metadata = {}

ko_gene_counts = defaultdict(
    lambda: defaultdict(int)
)


metadata_columns = [

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


print()
print("=" * 72)
print("IRON ACQUISITION STRATEGY — CANDIDATE MARKER AUDIT")
print("=" * 72)
print()
print(f"Reading:")
print(f"  {INPUT}")
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
            "Input file has no header."
        )


    required = {
        "Final_MAG",
        "KO",
    }


    missing_required = (
        required
        - set(reader.fieldnames)
    )


    if missing_required:

        raise ValueError(
            "Missing required columns: "
            + ", ".join(
                sorted(
                    missing_required
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


        # ----------------------------------------------------
        # Retain MAG metadata from first observed gene row
        # ----------------------------------------------------

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
                in metadata_columns
            }


        # ----------------------------------------------------
        # Extract KOs
        # ----------------------------------------------------

        kos = extract_kos(
            row.get(
                "KO",
                ""
            )
        )


        for ko in kos:

            mag_kos[
                mag
            ].add(
                ko
            )


            if ko in candidate_kos:

                ko_gene_counts[
                    ko
                ][
                    mag
                ] += 1


# Ensure MAGs that had no KO remain represented.

for mag in metadata:

    if mag not in mag_kos:

        mag_kos[
            mag
        ] = set()


mags = sorted(
    metadata.keys()
)


N_MAGS = len(
    mags
)


print(
    f"Gene rows read: "
    f"{rows_read:,}"
)

print(
    f"Unique MAGs: "
    f"{N_MAGS}"
)

print(
    f"Candidate iron KOs audited: "
    f"{len(candidate_kos)}"
)

print()


if N_MAGS != 373:

    print(
        "WARNING: expected 373 MAGs "
        f"but found {N_MAGS}."
    )

    print()


# ============================================================
# INDIVIDUAL KO PREVALENCE
# ============================================================

ko_prevalence_rows = []


for ko in sorted(
    candidate_kos
):

    mags_with = [

        mag
        for mag in mags

        if ko in mag_kos[
            mag
        ]
    ]


    n_mags = len(
        mags_with
    )


    prevalence = (
        100
        * n_mags
        / N_MAGS
        if N_MAGS
        else 0
    )


    total_gene_hits = sum(

        ko_gene_counts[
            ko
        ].values()
    )


    ko_prevalence_rows.append({

        "KO":
            ko,

        "Description":
            KO_LABELS.get(
                ko,
                ""
            ),

        "MAGs_with_KO":
            n_mags,

        "Prevalence_percent":
            f"{prevalence:.2f}",

        "Total_gene_hits":
            total_gene_hits,
    })


with KO_PREVALENCE_OUT.open(
    "w",
    newline=""
) as outfile:

    fields = [

        "KO",
        "Description",
        "MAGs_with_KO",
        "Prevalence_percent",
        "Total_gene_hits",
    ]


    writer = csv.DictWriter(
        outfile,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        ko_prevalence_rows
    )


# ============================================================
# ROUTE / SYSTEM SUMMARY
# ============================================================

route_summary_rows = []


all_routes = {}


for name, info in CORE_SYSTEMS.items():

    all_routes[
        name
    ] = {

        "type":
            "Core iron acquisition system",

        "category":
            info[
                "category"
            ],

        "kos":
            info[
                "kos"
            ],
    }


for name, kos in SIDEROPHORE_ROUTES.items():

    all_routes[
        name
    ] = {

        "type":
            "Siderophore biosynthesis route",

        "category":
            "Siderophore biosynthesis",

        "kos":
            kos,
    }


for route_name, info in all_routes.items():

    expected = info[
        "kos"
    ]


    complete_count = 0

    at_least_75_count = 0

    at_least_50_count = 0

    any_marker_count = 0


    fractions = []


    for mag in mags:

        stats = route_stats(
            mag_kos[
                mag
            ],
            expected
        )


        fraction = stats[
            "fraction"
        ]


        fractions.append(
            fraction
        )


        if stats[
            "n_detected"
        ] > 0:

            any_marker_count += 1


        if fraction >= 0.50:

            at_least_50_count += 1


        if fraction >= 0.75:

            at_least_75_count += 1


        if stats[
            "complete"
        ]:

            complete_count += 1


    mean_fraction = (

        sum(
            fractions
        )
        / N_MAGS

        if N_MAGS
        else 0
    )


    route_summary_rows.append({

        "Route":
            route_name,

        "Type":
            info[
                "type"
            ],

        "Category":
            info[
                "category"
            ],

        "Expected_KOs":
            ";".join(
                expected
            ),

        "N_expected_KOs":
            len(
                expected
            ),

        "MAGs_any_marker":
            any_marker_count,

        "MAGs_at_least_50pct":
            at_least_50_count,

        "MAGs_at_least_75pct":
            at_least_75_count,

        "MAGs_all_markers":
            complete_count,

        "All_markers_percent":
            (
                f"{100 * complete_count / N_MAGS:.2f}"
                if N_MAGS
                else "0.00"
            ),

        "Mean_marker_fraction":
            f"{mean_fraction:.4f}",
    })


with ROUTE_SUMMARY_OUT.open(
    "w",
    newline=""
) as outfile:

    fields = [

        "Route",
        "Type",
        "Category",
        "Expected_KOs",
        "N_expected_KOs",

        "MAGs_any_marker",
        "MAGs_at_least_50pct",
        "MAGs_at_least_75pct",
        "MAGs_all_markers",

        "All_markers_percent",
        "Mean_marker_fraction",
    ]


    writer = csv.DictWriter(
        outfile,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        route_summary_rows
    )


# ============================================================
# PER-MAG AUDIT TABLE
# ============================================================

mag_rows = []


for mag in mags:

    row = {

        "Final_MAG":
            mag,
    }


    # --------------------------------------------------------
    # Metadata
    # --------------------------------------------------------

    for column in metadata_columns:

        row[
            column
        ] = metadata[
            mag
        ].get(
            column,
            ""
        )


    # --------------------------------------------------------
    # Core system evidence
    # --------------------------------------------------------

    for system_name, info in CORE_SYSTEMS.items():

        stats = route_stats(
            mag_kos[
                mag
            ],
            info[
                "kos"
            ]
        )


        row[
            f"{system_name}_n"
        ] = stats[
            "n_detected"
        ]


        row[
            f"{system_name}_fraction"
        ] = (
            f"{stats['fraction']:.4f}"
        )


        row[
            f"{system_name}_all"
        ] = stats[
            "complete"
        ]


        row[
            f"{system_name}_detected_KOs"
        ] = ";".join(
            stats[
                "detected_kos"
            ]
        )


    # --------------------------------------------------------
    # Individual Feo evidence
    # --------------------------------------------------------

    row[
        "FeoA_K04758"
    ] = int(
        "K04758"
        in mag_kos[
            mag
        ]
    )


    row[
        "FeoB_K04759"
    ] = int(
        "K04759"
        in mag_kos[
            mag
        ]
    )


    # --------------------------------------------------------
    # Individual uptake evidence
    # --------------------------------------------------------

    for ko in [

        "K02010",
        "K02011",
        "K02012",

        "K02013",
        "K02014",
        "K02015",
        "K02016",

        "K03561",
        "K03832",

    ]:

        row[
            ko
        ] = int(
            ko
            in mag_kos[
                mag
            ]
        )


    # --------------------------------------------------------
    # Siderophore biosynthesis routes
    # --------------------------------------------------------

    complete_siderophore_routes = []


    best_siderophore_fraction = 0.0


    best_siderophore_route = ""


    for route_name, kos in SIDEROPHORE_ROUTES.items():

        stats = route_stats(
            mag_kos[
                mag
            ],
            kos
        )


        row[
            f"{route_name}_n"
        ] = stats[
            "n_detected"
        ]


        row[
            f"{route_name}_fraction"
        ] = (
            f"{stats['fraction']:.4f}"
        )


        row[
            f"{route_name}_all"
        ] = stats[
            "complete"
        ]


        row[
            f"{route_name}_detected_KOs"
        ] = ";".join(
            stats[
                "detected_kos"
            ]
        )


        if stats[
            "complete"
        ]:

            complete_siderophore_routes.append(
                route_name
            )


        if stats[
            "fraction"
        ] > best_siderophore_fraction:

            best_siderophore_fraction = stats[
                "fraction"
            ]

            best_siderophore_route = route_name


    row[
        "N_complete_siderophore_routes"
    ] = len(
        complete_siderophore_routes
    )


    row[
        "Complete_siderophore_routes"
    ] = ";".join(
        complete_siderophore_routes
    )


    row[
        "Best_siderophore_route"
    ] = best_siderophore_route


    row[
        "Best_siderophore_marker_fraction"
    ] = (
        f"{best_siderophore_fraction:.4f}"
    )


    # --------------------------------------------------------
    # Supplementary haem evidence
    # --------------------------------------------------------

    for ko, description in SUPPLEMENTARY_MARKERS.items():

        row[
            f"Supplementary_{ko}"
        ] = int(
            ko
            in mag_kos[
                mag
            ]
        )


    mag_rows.append(
        row
    )


# Build field order dynamically from first row.

if not mag_rows:

    raise ValueError(
        "No MAG rows were generated."
    )


fieldnames = list(
    mag_rows[
        0
    ].keys()
)


with MAG_AUDIT_OUT.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        mag_rows
    )


# ============================================================
# CONVENIENCE COUNTS
# ============================================================

def count_mags_with_ko(
    ko
):

    return sum(

        1

        for mag in mags

        if ko in mag_kos[
            mag
        ]
    )


def count_complete(
    kos
):

    return sum(

        1

        for mag in mags

        if set(
            kos
        ).issubset(
            mag_kos[
                mag
            ]
        )
    )


# Core system counts

feoA_n = count_mags_with_ko(
    "K04758"
)

feoB_n = count_mags_with_ko(
    "K04759"
)

feoAB_n = count_complete(
    CORE_SYSTEMS[
        "FeoAB_ferrous_uptake"
    ][
        "kos"
    ]
)


afu_n = count_complete(
    CORE_SYSTEMS[
        "Afu_Fbp_ferric_ABC"
    ][
        "kos"
    ]
)


iron_complex_abc_n = count_complete(
    CORE_SYSTEMS[
        "Iron_complex_ABC"
    ][
        "kos"
    ]
)


tonB_receptor_n = count_complete(
    CORE_SYSTEMS[
        "TonB_iron_complex_receptor"
    ][
        "kos"
    ]
)


exbB_tonB_n = count_complete(
    CORE_SYSTEMS[
        "TonB_ExbB_support"
    ][
        "kos"
    ]
)


heme_n = count_mags_with_ko(
    "K16087"
)


# Any complete siderophore route

any_complete_siderophore = 0


for mag in mags:

    has_complete_route = False


    for kos in SIDEROPHORE_ROUTES.values():

        if set(
            kos
        ).issubset(
            mag_kos[
                mag
            ]
        ):

            has_complete_route = True

            break


    if has_complete_route:

        any_complete_siderophore += 1


# ============================================================
# WRITE TEXT SUMMARY
# ============================================================

summary_lines = []


summary_lines.append(
    "=" * 72
)

summary_lines.append(
    "IRON ACQUISITION STRATEGY — AUDIT SUMMARY"
)

summary_lines.append(
    "=" * 72
)

summary_lines.append(
    ""
)

summary_lines.append(
    f"MAGs audited: {N_MAGS}"
)

summary_lines.append(
    f"Candidate KOs audited: {len(candidate_kos)}"
)

summary_lines.append(
    ""
)

summary_lines.append(
    "CORE IRON ACQUISITION EVIDENCE"
)

summary_lines.append(
    "-" * 72
)

summary_lines.append(
    f"FeoA K04758: {feoA_n} MAGs"
)

summary_lines.append(
    f"FeoB K04759: {feoB_n} MAGs"
)

summary_lines.append(
    f"FeoA + FeoB: {feoAB_n} MAGs"
)

summary_lines.append(
    ""
)

summary_lines.append(
    "Direct Fe(III) Afu/Fbp ABC "
    f"(K02010+K02011+K02012): {afu_n} MAGs"
)

summary_lines.append(
    "Iron-complex ABC "
    f"(K02013+K02015+K02016): "
    f"{iron_complex_abc_n} MAGs"
)

summary_lines.append(
    "Iron-complex receptor + TonB "
    f"(K02014+K03832): "
    f"{tonB_receptor_n} MAGs"
)

summary_lines.append(
    "ExbB + TonB support "
    f"(K03561+K03832): "
    f"{exbB_tonB_n} MAGs"
)

summary_lines.append(
    ""
)

summary_lines.append(
    "SIDEROPHORE BIOSYNTHESIS ROUTES"
)

summary_lines.append(
    "-" * 72
)


for route_name, kos in SIDEROPHORE_ROUTES.items():

    complete_n = count_complete(
        kos
    )


    summary_lines.append(

        f"{route_name}: "
        f"{complete_n}/{N_MAGS} "
        f"with all audited markers"
    )


summary_lines.append(
    ""
)

summary_lines.append(
    "MAGs with >=1 complete audited "
    f"siderophore route: "
    f"{any_complete_siderophore}"
)

summary_lines.append(
    ""
)

summary_lines.append(
    "SUPPLEMENTARY"
)

summary_lines.append(
    "-" * 72
)

summary_lines.append(
    f"K16087 haem-associated receptor: "
    f"{heme_n} MAGs"
)

summary_lines.append(
    ""
)

summary_lines.append(
    "IMPORTANT:"
)

summary_lines.append(
    "These are audit results only."
)

summary_lines.append(
    "No final Iron capability calls have yet been assigned."
)

summary_lines.append(
    "Marker fractions are evidence summaries, "
    "not physiological percentages."
)

summary_lines.append(
    "=" * 72
)


with SUMMARY_OUT.open(
    "w"
) as outfile:

    outfile.write(
        "\n".join(
            summary_lines
        )
        + "\n"
    )


# ============================================================
# PRINT SUMMARY TO TERMINAL
# ============================================================

print(
    "\n".join(
        summary_lines
    )
)

print()

print(
    "OUTPUT FILES"
)

print(
    "-" * 72
)

print(
    f"KO prevalence:\n"
    f"  {KO_PREVALENCE_OUT}"
)

print()

print(
    f"Route summary:\n"
    f"  {ROUTE_SUMMARY_OUT}"
)

print()

print(
    f"Per-MAG audit:\n"
    f"  {MAG_AUDIT_OUT}"
)

print()

print(
    f"Text summary:\n"
    f"  {SUMMARY_OUT}"
)

print(
    "=" * 72
)

