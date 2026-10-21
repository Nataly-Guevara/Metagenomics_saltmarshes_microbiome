#!/usr/bin/env python3

"""
Audit candidate markers for the Colonisation and Plant-Interaction
ecological strategy across the 373 saltmarsh MAGs.

PURPOSE
=======
This script DOES NOT create the final Strategy 6 score.

It audits four biological blocks:

1. Chemotaxis
   - environmental sensing and directional response

2. Motility
   - flagellar motility
   - type IV pilus / twitching motility

3. Attachment / biofilm-matrix potential
   - PGA biofilm polysaccharide
   - curli
   - LapA adhesion system
   - Tad/Flp attachment pili
   - succinoglycan/EPS markers

4. Direct plant-interaction potential
   - ACC deaminase
   - indole-3-acetamide (IAM) IAA biosynthesis route

The final binary capability rules will be decided AFTER examining
the prevalence and co-occurrence patterns produced by this audit.

INPUT
=====
results/functional_annotation/kofamscan/
final_373_MAGs_gene_master_table.tsv

OUTPUTS
=======
results/ecological_strategies/colonisation_host_interaction/audit/

    colonisation_candidate_KO_prevalence.tsv
    colonisation_candidate_system_summary.tsv
    colonisation_host_audit_373_MAGs.tsv
    colonisation_host_audit_summary.txt
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
      "colonisation_host_interaction/audit"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

KO_PREVALENCE_OUT = (
    OUTPUT_DIR
    / "colonisation_candidate_KO_prevalence.tsv"
)

SYSTEM_SUMMARY_OUT = (
    OUTPUT_DIR
    / "colonisation_candidate_system_summary.tsv"
)

MAG_AUDIT_OUT = (
    OUTPUT_DIR
    / "colonisation_host_audit_373_MAGs.tsv"
)

SUMMARY_OUT = (
    OUTPUT_DIR
    / "colonisation_host_audit_summary.txt"
)


# ============================================================
# CANDIDATE SYSTEMS
# ============================================================

SYSTEMS = {

    # ========================================================
    # 1. CHEMOTAXIS
    # ========================================================

    "Chemotaxis_core": {

        "category":
            "Chemotaxis",

        "interpretation":
            "Core sensing and directional-response machinery",

        "kos": [
            "K03406",  # MCP
            "K03407",  # CheA
            "K03408",  # CheW
            "K03413",  # CheY
        ],
    },


    "Chemotaxis_adaptation": {

        "category":
            "Chemotaxis",

        "interpretation":
            "Chemotaxis receptor adaptation machinery",

        "kos": [
            "K00575",  # CheR
            "K03412",  # CheB
        ],
    },


    # ========================================================
    # 2A. FLAGELLAR MOTILITY
    # ========================================================

    "Flagellar_motor": {

        "category":
            "Motility_flagellar",

        "interpretation":
            "Flagellar stator/motor machinery",

        "kos": [
            "K02556",  # MotA/PomA
            "K02557",  # MotB/PomB
        ],
    },


    "Flagellar_switch": {

        "category":
            "Motility_flagellar",

        "interpretation":
            "Flagellar basal-body/switch machinery",

        "kos": [
            "K02409",  # FliF
            "K02410",  # FliG
            "K02416",  # FliM
            "K02417",  # FliN
        ],
    },


    "Flagellar_export": {

        "category":
            "Motility_flagellar",

        "interpretation":
            "Flagellar protein-export machinery",

        "kos": [
            "K02400",  # FlhA
            "K02401",  # FlhB
            "K02419",  # FliP
            "K02420",  # FliQ
            "K02421",  # FliR
        ],
    },


    "Flagellar_hook_filament": {

        "category":
            "Motility_flagellar",

        "interpretation":
            "Flagellar hook/filament evidence",

        "kos": [
            "K02390",  # FlgE hook
            "K02406",  # FliC flagellin
            "K02407",  # FliD filament cap
        ],
    },


    # ========================================================
    # 2B. TYPE IV PILUS / TWITCHING MOTILITY
    # ========================================================

    "Twitching_signalling": {

        "category":
            "Motility_twitching",

        "interpretation":
            "Type-IV-pilus twitching chemosensory machinery",

        "kos": [
            "K02657",  # PilG
            "K02658",  # PilH
            "K02659",  # PilI
            "K02660",  # PilJ
        ],
    },


    "Twitching_retraction": {

        "category":
            "Motility_twitching",

        "interpretation":
            "Type-IV-pilus twitching retraction ATPases",

        "kos": [
            "K02669",  # PilT
            "K02670",  # PilU
        ],
    },


    "TypeIV_pilus_assembly": {

        "category":
            "Motility_twitching",

        "interpretation":
            "Type IV pilus assembly machinery",

        "kos": [
            "K02662",  # PilM
            "K02663",  # PilN
            "K02664",  # PilO
            "K02665",  # PilP
            "K02666",  # PilQ
        ],
    },


    # ========================================================
    # 3. ATTACHMENT / BIOFILM-MATRIX
    # ========================================================

    "PGA_biofilm_matrix": {

        "category":
            "Attachment_biofilm",

        "interpretation":
            "Poly-beta-1,6-N-acetyl-D-glucosamine biofilm system",

        "kos": [
            "K11935",  # PgaA
            "K11931",  # PgaB
            "K11936",  # PgaC
            "K11937",  # PgaD
        ],
    },


    "Curli_attachment": {

        "category":
            "Attachment_biofilm",

        "interpretation":
            "Curli extracellular attachment fibres",

        "kos": [
            "K04334",  # CsgA
            "K04335",  # CsgB
        ],
    },


    "LapA_adhesin": {

        "category":
            "Attachment_biofilm",

        "interpretation":
            "Large surface adhesin LapA",

        "kos": [
            "K12549",  # LapA
        ],
    },


    "LapA_transport": {

        "category":
            "Attachment_biofilm",

        "interpretation":
            "Type I adhesin export machinery",

        "kos": [
            "K12541",  # LapB / ABC component
            "K12542",  # membrane fusion protein
            "K12543",  # outer membrane component LapE
        ],
    },


    "Tad_Flp_attachment": {

        "category":
            "Attachment_biofilm",

        "interpretation":
            "Tad/Flp pili associated with stable surface attachment",

        "kos": [
            "K02651",  # Flp pilin
            "K02278",  # CpaA/TadV
            "K02279",  # CpaB
            "K02280",  # CpaC
            "K02283",  # CpaF
            "K12510",  # TadB
            "K12511",  # TadC
        ],
    },


    "Succinoglycan_EPS_markers": {

        "category":
            "Attachment_biofilm",

        "interpretation":
            "Rhizobial-type EPS/succinoglycan-associated markers",

        "kos": [
            "K16566",  # ExoY
            "K16567",  # ExoQ
            "K16568",  # ExoZ
        ],
    },


    # ========================================================
    # 4. DIRECT PLANT-INTERACTION POTENTIAL
    # ========================================================

    "ACC_deaminase": {

        "category":
            "Plant_interaction",

        "interpretation":
            "1-aminocyclopropane-1-carboxylate deaminase",

        "kos": [
            "K01505",
        ],
    },


    "IAA_IAM_route": {

        "category":
            "Plant_interaction",

        "interpretation":
            "Indole-3-acetamide pathway to IAA",

        "kos": [
            "K00466",  # IaaM / tryptophan 2-monooxygenase
            "K21801",  # IaaH / indoleacetamide hydrolase
        ],
    },


    # --------------------------------------------------------
    # SUPPORTING ONLY:
    #
    # K01426 is a broad amidase. It can catalyse IAM
    # hydrolysis but is not sufficiently specific to score
    # IAA production alone.
    # --------------------------------------------------------

    "IAA_amidase_support": {

        "category":
            "Plant_interaction_support",

        "interpretation":
            "Broad amidase capable of IAM hydrolysis; support only",

        "kos": [
            "K01426",
        ],
    },
}


# ============================================================
# METADATA COLUMNS
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
# BUILD CANDIDATE KO SET
# ============================================================

candidate_kos = set()


for info in SYSTEMS.values():

    candidate_kos.update(
        info["kos"]
    )


# ============================================================
# FUNCTIONS
# ============================================================

def extract_kos(value):

    """
    Robustly extract K numbers from KO field.
    """

    if not value:

        return []

    return re.findall(
        r"K\d{5}",
        value
    )


def system_stats(
    mag_kos,
    expected_kos
):

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
        if n_expected
        else 0
    )

    complete = int(
        n_detected
        == n_expected
        and n_expected > 0
    )

    return {

        "n_detected":
            n_detected,

        "n_expected":
            n_expected,

        "fraction":
            fraction,

        "complete":
            complete,

        "detected_kos":
            sorted(
                detected
            ),

        "missing_kos":
            sorted(
                expected - detected
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

ko_definitions = {}


print()

print(
    "=" * 72
)

print(
    "COLONISATION / HOST-INTERACTION MARKER AUDIT"
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
            "Input file has no header."
        )


    required = {
        "Final_MAG",
        "KO",
    }


    missing_required = (
        required
        - set(
            reader.fieldnames
        )
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


                definition = (
                    row.get(
                        "KO_definition",
                        ""
                    )
                    .strip()
                )


                if (
                    definition
                    and ko not in ko_definitions
                ):

                    ko_definitions[
                        ko
                    ] = definition


# Ensure MAGs without KOs remain represented.

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
    f"Candidate KOs audited: "
    f"{len(candidate_kos)}"
)

print(
    f"Candidate systems audited: "
    f"{len(SYSTEMS)}"
)

print()


if N_MAGS != 373:

    print(
        "WARNING: expected 373 MAGs, "
        f"found {N_MAGS}."
    )

    print()


# ============================================================
# MAP KO -> SYSTEMS
# ============================================================

ko_to_systems = defaultdict(
    list
)


for system_name, info in SYSTEMS.items():

    for ko in info[
        "kos"
    ]:

        ko_to_systems[
            ko
        ].append(
            system_name
        )


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

        "KO_definition":
            ko_definitions.get(
                ko,
                ""
            ),

        "Candidate_systems":
            ";".join(
                ko_to_systems[
                    ko
                ]
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
        "KO_definition",
        "Candidate_systems",
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
# SYSTEM SUMMARY
# ============================================================

system_summary_rows = []


for system_name, info in SYSTEMS.items():

    expected = info[
        "kos"
    ]


    any_count = 0

    half_count = 0

    three_quarter_count = 0

    complete_count = 0

    fractions = []


    for mag in mags:

        stats = system_stats(
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

            any_count += 1


        if fraction >= 0.50:

            half_count += 1


        if fraction >= 0.75:

            three_quarter_count += 1


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


    system_summary_rows.append({

        "System":
            system_name,

        "Category":
            info[
                "category"
            ],

        "Interpretation":
            info[
                "interpretation"
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
            any_count,

        "MAGs_at_least_50pct":
            half_count,

        "MAGs_at_least_75pct":
            three_quarter_count,

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


with SYSTEM_SUMMARY_OUT.open(
    "w",
    newline=""
) as outfile:

    fields = [

        "System",
        "Category",
        "Interpretation",
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
        system_summary_rows
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


    for column in METADATA_COLUMNS:

        row[
            column
        ] = metadata[
            mag
        ].get(
            column,
            ""
        )


    # --------------------------------------------------------
    # Store evidence for every candidate system
    # --------------------------------------------------------

    system_calls = {}


    for system_name, info in SYSTEMS.items():

        stats = system_stats(
            mag_kos[
                mag
            ],
            info[
                "kos"
            ]
        )


        system_calls[
            system_name
        ] = stats


        row[
            f"{system_name}_n"
        ] = stats[
            "n_detected"
        ]


        row[
            f"{system_name}_expected"
        ] = stats[
            "n_expected"
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


    # ========================================================
    # EXPLORATORY SUMMARIES
    #
    # These are NOT final capability calls.
    # They are included to help us decide final rules.
    # ========================================================


    # --------------------------------------------------------
    # Chemotaxis
    # --------------------------------------------------------

    row[
        "EXPLORE_Chemotaxis_core_complete"
    ] = system_calls[
        "Chemotaxis_core"
    ][
        "complete"
    ]


    row[
        "EXPLORE_Chemotaxis_core_75pct"
    ] = int(
        system_calls[
            "Chemotaxis_core"
        ][
            "fraction"
        ]
        >= 0.75
    )


    # --------------------------------------------------------
    # Flagellar motility:
    # count how many of four flagellar submodules are complete.
    # --------------------------------------------------------

    flagellar_modules = [

        "Flagellar_motor",

        "Flagellar_switch",

        "Flagellar_export",

        "Flagellar_hook_filament",
    ]


    n_flagellar_modules_complete = sum(

        system_calls[
            module
        ][
            "complete"
        ]

        for module
        in flagellar_modules
    )


    row[
        "EXPLORE_N_flagellar_modules_complete"
    ] = n_flagellar_modules_complete


    row[
        "EXPLORE_Flagellar_3of4_modules"
    ] = int(
        n_flagellar_modules_complete
        >= 3
    )


    # --------------------------------------------------------
    # Twitching/type-IV-pilus evidence
    # --------------------------------------------------------

    pilus_assembly_complete = (
        system_calls[
            "TypeIV_pilus_assembly"
        ][
            "complete"
        ]
    )


    pilT_present = int(
        "K02669"
        in mag_kos[
            mag
        ]
    )


    pilU_present = int(
        "K02670"
        in mag_kos[
            mag
        ]
    )


    row[
        "EXPLORE_PilT_present"
    ] = pilT_present


    row[
        "EXPLORE_PilU_present"
    ] = pilU_present


    row[
        "EXPLORE_TypeIV_assembly_complete"
    ] = pilus_assembly_complete


    row[
        "EXPLORE_TypeIV_assembly_plus_PilT"
    ] = int(
        pilus_assembly_complete
        and pilT_present
    )


    # --------------------------------------------------------
    # Attachment systems
    # --------------------------------------------------------

    lapA_present = (
        system_calls[
            "LapA_adhesin"
        ][
            "complete"
        ]
    )


    lap_transport_complete = (
        system_calls[
            "LapA_transport"
        ][
            "complete"
        ]
    )


    row[
        "EXPLORE_LapA_plus_transport"
    ] = int(
        lapA_present
        and lap_transport_complete
    )


    attachment_complete_routes = []


    for route in [

        "PGA_biofilm_matrix",

        "Curli_attachment",

        "Tad_Flp_attachment",

        "Succinoglycan_EPS_markers",
    ]:

        if system_calls[
            route
        ][
            "complete"
        ]:

            attachment_complete_routes.append(
                route
            )


    if (
        lapA_present
        and lap_transport_complete
    ):

        attachment_complete_routes.append(
            "LapA_complete_system"
        )


    row[
        "EXPLORE_N_complete_attachment_routes"
    ] = len(
        attachment_complete_routes
    )


    row[
        "EXPLORE_complete_attachment_routes"
    ] = ";".join(
        attachment_complete_routes
    )


    # --------------------------------------------------------
    # Direct plant-interaction evidence
    # --------------------------------------------------------

    acc_present = (
        system_calls[
            "ACC_deaminase"
        ][
            "complete"
        ]
    )


    iaa_iam_complete = (
        system_calls[
            "IAA_IAM_route"
        ][
            "complete"
        ]
    )


    iaaM_present = int(
        "K00466"
        in mag_kos[
            mag
        ]
    )


    iaaH_specific_present = int(
        "K21801"
        in mag_kos[
            mag
        ]
    )


    broad_amidase_present = int(
        "K01426"
        in mag_kos[
            mag
        ]
    )


    row[
        "EXPLORE_ACC_deaminase"
    ] = acc_present


    row[
        "EXPLORE_IaaM_K00466"
    ] = iaaM_present


    row[
        "EXPLORE_specific_IaaH_K21801"
    ] = iaaH_specific_present


    row[
        "EXPLORE_broad_amidase_K01426"
    ] = broad_amidase_present


    row[
        "EXPLORE_complete_IAM_IAA_route"
    ] = iaa_iam_complete


    row[
        "EXPLORE_any_direct_plant_marker"
    ] = int(
        acc_present
        or iaa_iam_complete
    )


    mag_rows.append(
        row
    )


if not mag_rows:

    raise ValueError(
        "No MAG rows generated."
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

def count_mag_rows(
    column
):

    return sum(

        int(
            row[
                column
            ]
        )

        for row
        in mag_rows
    )


chemotaxis_strict_n = count_mag_rows(
    "EXPLORE_Chemotaxis_core_complete"
)

chemotaxis_75_n = count_mag_rows(
    "EXPLORE_Chemotaxis_core_75pct"
)

flagellar_3of4_n = count_mag_rows(
    "EXPLORE_Flagellar_3of4_modules"
)

pilus_plus_pilT_n = count_mag_rows(
    "EXPLORE_TypeIV_assembly_plus_PilT"
)

lap_complete_n = count_mag_rows(
    "EXPLORE_LapA_plus_transport"
)

acc_n = count_mag_rows(
    "EXPLORE_ACC_deaminase"
)

iaa_complete_n = count_mag_rows(
    "EXPLORE_complete_IAM_IAA_route"
)

iaaM_n = count_mag_rows(
    "EXPLORE_IaaM_K00466"
)

iaaH_n = count_mag_rows(
    "EXPLORE_specific_IaaH_K21801"
)

broad_amidase_n = count_mag_rows(
    "EXPLORE_broad_amidase_K01426"
)

any_direct_plant_n = count_mag_rows(
    "EXPLORE_any_direct_plant_marker"
)


# ============================================================
# SYSTEM COMPLETE COUNTS
# ============================================================

complete_counts = {}


for system_name, info in SYSTEMS.items():

    complete_counts[
        system_name
    ] = sum(

        1

        for mag in mags

        if system_stats(
            mag_kos[
                mag
            ],
            info[
                "kos"
            ]
        )[
            "complete"
        ]
    )


# ============================================================
# WRITE TEXT SUMMARY
# ============================================================

summary_lines = []


summary_lines.append(
    "=" * 72
)

summary_lines.append(
    "COLONISATION / HOST-INTERACTION AUDIT SUMMARY"
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
    f"Candidate systems audited: {len(SYSTEMS)}"
)

summary_lines.append(
    ""
)


# ------------------------------------------------------------
# Chemotaxis
# ------------------------------------------------------------

summary_lines.append(
    "1. CHEMOTAXIS"
)

summary_lines.append(
    "-" * 72
)

summary_lines.append(
    f"Core MCP+CheA+CheW+CheY complete: "
    f"{chemotaxis_strict_n}/{N_MAGS}"
)

summary_lines.append(
    f"Core chemotaxis >=75%: "
    f"{chemotaxis_75_n}/{N_MAGS}"
)

summary_lines.append(
    f"CheR+CheB adaptation complete: "
    f"{complete_counts['Chemotaxis_adaptation']}/{N_MAGS}"
)

summary_lines.append(
    ""
)


# ------------------------------------------------------------
# Motility
# ------------------------------------------------------------

summary_lines.append(
    "2. MOTILITY"
)

summary_lines.append(
    "-" * 72
)

summary_lines.append(
    f"Flagellar motor complete: "
    f"{complete_counts['Flagellar_motor']}/{N_MAGS}"
)

summary_lines.append(
    f"Flagellar switch complete: "
    f"{complete_counts['Flagellar_switch']}/{N_MAGS}"
)

summary_lines.append(
    f"Flagellar export complete: "
    f"{complete_counts['Flagellar_export']}/{N_MAGS}"
)

summary_lines.append(
    f"Flagellar hook/filament complete: "
    f"{complete_counts['Flagellar_hook_filament']}/{N_MAGS}"
)

summary_lines.append(
    f"At least 3/4 flagellar submodules complete: "
    f"{flagellar_3of4_n}/{N_MAGS}"
)

summary_lines.append(
    ""
)

summary_lines.append(
    f"Twitching signalling complete: "
    f"{complete_counts['Twitching_signalling']}/{N_MAGS}"
)

summary_lines.append(
    f"Twitching PilT+PilU complete: "
    f"{complete_counts['Twitching_retraction']}/{N_MAGS}"
)

summary_lines.append(
    f"Type-IV-pilus assembly complete: "
    f"{complete_counts['TypeIV_pilus_assembly']}/{N_MAGS}"
)

summary_lines.append(
    f"Type-IV-pilus assembly + PilT: "
    f"{pilus_plus_pilT_n}/{N_MAGS}"
)

summary_lines.append(
    ""
)


# ------------------------------------------------------------
# Attachment
# ------------------------------------------------------------

summary_lines.append(
    "3. ATTACHMENT / BIOFILM-MATRIX"
)

summary_lines.append(
    "-" * 72
)

summary_lines.append(
    f"PGA system complete: "
    f"{complete_counts['PGA_biofilm_matrix']}/{N_MAGS}"
)

summary_lines.append(
    f"Curli CsgA+CsgB complete: "
    f"{complete_counts['Curli_attachment']}/{N_MAGS}"
)

summary_lines.append(
    f"LapA detected: "
    f"{complete_counts['LapA_adhesin']}/{N_MAGS}"
)

summary_lines.append(
    f"Lap transport complete: "
    f"{complete_counts['LapA_transport']}/{N_MAGS}"
)

summary_lines.append(
    f"LapA + complete transport: "
    f"{lap_complete_n}/{N_MAGS}"
)

summary_lines.append(
    f"Tad/Flp marker set complete: "
    f"{complete_counts['Tad_Flp_attachment']}/{N_MAGS}"
)

summary_lines.append(
    f"Succinoglycan/EPS marker set complete: "
    f"{complete_counts['Succinoglycan_EPS_markers']}/{N_MAGS}"
)

summary_lines.append(
    ""
)


# ------------------------------------------------------------
# Plant interaction
# ------------------------------------------------------------

summary_lines.append(
    "4. DIRECT PLANT-INTERACTION POTENTIAL"
)

summary_lines.append(
    "-" * 72
)

summary_lines.append(
    f"ACC deaminase K01505: "
    f"{acc_n}/{N_MAGS}"
)

summary_lines.append(
    f"IaaM K00466: "
    f"{iaaM_n}/{N_MAGS}"
)

summary_lines.append(
    f"Specific IaaH K21801: "
    f"{iaaH_n}/{N_MAGS}"
)

summary_lines.append(
    f"Complete IaaM + specific IaaH IAM route: "
    f"{iaa_complete_n}/{N_MAGS}"
)

summary_lines.append(
    f"Broad amidase K01426 "
    f"(support only): "
    f"{broad_amidase_n}/{N_MAGS}"
)

summary_lines.append(
    f"ACC deaminase OR complete IAM route: "
    f"{any_direct_plant_n}/{N_MAGS}"
)

summary_lines.append(
    ""
)


# ------------------------------------------------------------
# Important caveat
# ------------------------------------------------------------

summary_lines.append(
    "IMPORTANT"
)

summary_lines.append(
    "-" * 72
)

summary_lines.append(
    "These are audit results only."
)

summary_lines.append(
    "No final Colonisation/Host capability calls "
    "have yet been assigned."
)

summary_lines.append(
    "Absence of a complete marker set does not by itself "
    "prove absence of the phenotype."
)

summary_lines.append(
    "Marker fractions are genomic evidence summaries, "
    "not physiological percentages."
)

summary_lines.append(
    "Generic amidase K01426 is supplementary evidence only "
    "and must not be used alone to infer IAA production."
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
# PRINT SUMMARY
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
    f"System summary:\n"
    f"  {SYSTEM_SUMMARY_OUT}"
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
