#!/usr/bin/env python3

"""
Audit genome-resolved nitrogen acquisition and transformation rules.

The intended Ricci-style nitrogen strategy contains nine candidate
binary capabilities:

Nitrogen acquisition
--------------------
1. Nitrogen fixation
2. Ammonium acquisition
3. Nitrate/nitrite acquisition
4. Urea utilisation

Assimilation
------------
5. Assimilatory nitrate reduction

Nitrogen transformation
-----------------------
6. Nitrification
7. DNRA
8. Denitrification
9. Anammox

This script does not define the final score. It audits strict KEGG
modules, MAG-tolerant alternatives and individual pathway components.

Important methodological point
------------------------------
DNRA and denitrification markers also occur in the redox-flexibility
axis. This is intentional:

- Redox flexibility asks whether a MAG can function across alternative
  electron-acceptor and redox conditions.
- Nitrogen acquisition and transformation asks how a MAG acquires,
  retains, oxidises or removes nitrogen.

Shared genes may therefore contribute to both axes, but are counted
only once within each individual strategy axis.

Outputs
-------
1. nitrogen_capability_rule_audit_preliminary.tsv
2. nitrogen_capability_presence_preliminary.tsv
3. nitrogen_candidate_KO_prevalence.tsv
4. nitrogen_multigene_system_diagnostics.tsv
5. nitrogen_rule_definitions_preliminary.tsv
6. nitrogen_candidate_KO_definition_audit.tsv

No existing osmotic or redox files are modified.
Only the Python standard library is used.
"""

import csv
import itertools
import os
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = "/path/to/metagenomics_project"

INPUT_FILE = os.path.join(
    PROJECT_ROOT,
    "results/functional_annotation/kofamscan/"
    "final_373_MAGs_gene_master_table.tsv",
)

OUTPUT_ROOT = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/"
    "nitrogen_acquisition_transformation",
)

TABLE_OUTPUT_DIR = os.path.join(
    OUTPUT_ROOT,
    "tables",
)

RULE_AUDIT_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "nitrogen_capability_rule_audit_preliminary.tsv",
)

PRESENCE_MATRIX_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "nitrogen_capability_presence_preliminary.tsv",
)

KO_PREVALENCE_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "nitrogen_candidate_KO_prevalence.tsv",
)

SYSTEM_DIAGNOSTICS_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "nitrogen_multigene_system_diagnostics.tsv",
)

RULE_DEFINITIONS_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "nitrogen_rule_definitions_preliminary.tsv",
)

KO_DEFINITION_AUDIT_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "nitrogen_candidate_KO_definition_audit.tsv",
)


# ---------------------------------------------------------------------
# Metadata preserved per MAG
# ---------------------------------------------------------------------

METADATA_COLUMNS = [
    "Completeness",
    "Contamination",
    "Quality_class",
    "Site",
    "Treatment",
    "Compartment",
    "GTDB_domain",
    "GTDB_classification",
    "GTDB_phylum",
    "GTDB_class",
    "GTDB_order",
    "GTDB_family",
    "GTDB_genus",
    "GTDB_species",
]


KO_PATTERN = re.compile(r"K\d{5}")


# ---------------------------------------------------------------------
# Candidate KO descriptions
#
# Actual KO_definition values from the master table are exported
# separately and can be compared with these expected functions.
# ---------------------------------------------------------------------

CANDIDATE_KO_EXPECTATIONS = {
    # Nitrogen fixation
    "K02588": "nifH; nitrogenase iron protein",
    "K02586": "nifD; molybdenum nitrogenase alpha subunit",
    "K02591": "nifK; molybdenum nitrogenase beta subunit",
    "K00531": "nitrogenase delta/accessory subunit",
    "K22896": "vnfD; vanadium nitrogenase alpha subunit",
    "K22897": "vnfK; vanadium nitrogenase beta subunit",
    "K22898": "vnfG; vanadium nitrogenase delta subunit",
    "K22899": "vnfH; vanadium nitrogenase iron protein",

    # Ammonium acquisition
    "K03320": "amt; Amt-family ammonium transporter",

    # Nitrate and nitrite acquisition
    "K02575": "narK/nrtP/nasA; nitrate/nitrite transporter",
    "K15576": "nrtA; nitrate/nitrite transporter binding protein",
    "K15577": "nrtB; nitrate/nitrite transporter permease",
    "K15578": "nrtC; nitrate/nitrite transporter ATPase",
    "K15579": "nitrate/nitrite transporter accessory component",

    # Urea transport and hydrolysis
    "K11959": "urtA; urea transporter binding protein",
    "K11960": "urtB; urea transporter permease",
    "K11961": "urtC; urea transporter permease",
    "K11962": "urtD; urea transporter ATPase",
    "K11963": "urtE; urea transporter ATPase",
    "K01427": "single-chain urease",
    "K01428": "ureC; urease alpha subunit",
    "K01429": "ureB; urease beta subunit",
    "K01430": "ureA; urease gamma subunit",
    "K14048": "ureAB; fused urease gamma/beta subunit",

    # Assimilatory nitrate reduction
    "K00367": "nasA; assimilatory nitrate reductase",
    "K10534": "assimilatory nitrate reductase alternative",
    "K00372": "assimilatory nitrate reductase catalytic component",
    "K00360": "assimilatory nitrate reductase accessory component",
    "K00366": "nirA; ferredoxin-nitrite reductase",
    "K17877": "assimilatory nitrite reductase alternative",
    "K26139": "assimilatory nitrite reductase subunit",
    "K26138": "assimilatory nitrite reductase subunit",
    "K00361": "NAD(P)H-dependent assimilatory nitrite reductase",

    # Nitrification
    "K28504": "amoA; ammonia monooxygenase subunit A",
    "K10944": "amoA/pmoA; ammonia/methane monooxygenase subunit A",
    "K10945": "amoB/pmoB; ammonia/methane monooxygenase subunit B",
    "K10946": "amoC/pmoC; ammonia/methane monooxygenase subunit C",
    "K10535": "hao; hydroxylamine oxidoreductase",

    # Respiratory nitrate reduction shared by DNRA/denitrification
    "K00370": "narG; respiratory nitrate reductase alpha subunit",
    "K00371": "narH; respiratory nitrate reductase beta subunit",
    "K00374": "narI; respiratory nitrate reductase membrane subunit",
    "K02567": "napA; periplasmic nitrate reductase catalytic subunit",
    "K02568": "napB; periplasmic nitrate reductase cytochrome subunit",

    # DNRA
    "K00362": "nirB; nitrite reductase large subunit",
    "K00363": "nirD; nitrite reductase small subunit",
    "K03385": "nrfA; cytochrome-c nitrite reductase",
    "K15876": "nrfH; NrfA electron-transfer subunit",

    # Denitrification
    "K00368": "nirK; copper-containing nitrite reductase",
    "K15864": "nirS; cytochrome-cd1 nitrite reductase",
    "K04561": "norB; nitric-oxide reductase subunit B",
    "K02305": "norC; nitric-oxide reductase subunit C",
    "K00376": "nosZ; nitrous-oxide reductase",

    # Anammox
    "K20932": "hzsA; hydrazine synthase subunit A",
    "K20933": "hzsB; hydrazine synthase subunit B",
    "K20934": "hzsC; hydrazine synthase subunit C",
    "K20935": "hdh; hydrazine dehydrogenase",
}


# ---------------------------------------------------------------------
# Rule construction helpers
# ---------------------------------------------------------------------

def route(*required_kos):
    """Create one route requiring every supplied KO."""
    return tuple(required_kos)


def combine_stages(*stage_route_lists):
    """
    Combine alternative stages into complete pathway routes.

    Example:
        stage 1: NarGHI OR NapAB
        stage 2: NirBD OR NrfAH

    produces:
        NarGHI + NirBD
        NarGHI + NrfAH
        NapAB + NirBD
        NapAB + NrfAH
    """
    combined_routes = []

    for combination in itertools.product(*stage_route_lists):
        merged = []

        for stage_route in combination:
            for ko in stage_route:
                if ko not in merged:
                    merged.append(ko)

        combined_routes.append(
            tuple(merged)
        )

    return combined_routes


def rule_is_positive(ko_set, routes):
    """A rule is positive when at least one complete route is present."""
    return any(
        all(ko in ko_set for ko in pathway_route)
        for pathway_route in routes
    )


def format_rule_logic(routes):
    """Create a readable AND/OR rule description."""
    formatted_routes = []

    for pathway_route in routes:
        if len(pathway_route) == 1:
            formatted_routes.append(
                pathway_route[0]
            )
        else:
            formatted_routes.append(
                "(" + " + ".join(pathway_route) + ")"
            )

    return " OR ".join(
        formatted_routes
    )


# ---------------------------------------------------------------------
# Reusable pathway stages
# ---------------------------------------------------------------------

NITROGEN_FIXATION_ROUTES = [
    route("K02588", "K02586", "K02591"),
    route("K22896", "K22897", "K22898", "K22899"),
]

NITRATE_RESPIRATION_ROUTES = [
    route("K00370", "K00371", "K00374"),
    route("K02567", "K02568"),
]

DNRA_NITRITE_TO_AMMONIUM_ROUTES = [
    route("K00362", "K00363"),
    route("K03385", "K15876"),
]

DENITRIFICATION_NITRITE_ROUTES = [
    route("K00368"),
    route("K15864"),
]

DENITRIFICATION_NO_ROUTES = [
    route("K04561", "K02305"),
]

DENITRIFICATION_N2O_ROUTES = [
    route("K00376"),
]

ASSIMILATORY_NITRATE_ROUTES = [
    route("K00367"),
    route("K10534"),
    route("K00372"),
]

ASSIMILATORY_NITRITE_ROUTES = [
    route("K00366"),
    route("K17877"),
    route("K26139", "K26138"),
    route("K00361"),
]

UREASE_ROUTES = [
    route("K01427"),
    route("K01428", "K01429", "K01430"),
    route("K01428", "K14048"),
]

NITRATE_ACQUISITION_ROUTES = [
    route("K02575"),
    route("K15576", "K15577", "K15578"),
]

NITRIFICATION_CURRENT_ROUTES = [
    route("K28504", "K10945", "K10946", "K10535"),
]

NITRIFICATION_LEGACY_ROUTES = [
    route("K10944", "K10945", "K10946", "K10535"),
]

ANAMMOX_ROUTES = combine_stages(
    DENITRIFICATION_NITRITE_ROUTES,
    [
        route("K20932", "K20933", "K20934"),
    ],
    [
        route("K20935"),
    ],
)

DNRA_COMPLETE_ROUTES = combine_stages(
    NITRATE_RESPIRATION_ROUTES,
    DNRA_NITRITE_TO_AMMONIUM_ROUTES,
)

DENITRIFICATION_COMPLETE_ROUTES = combine_stages(
    NITRATE_RESPIRATION_ROUTES,
    DENITRIFICATION_NITRITE_ROUTES,
    DENITRIFICATION_NO_ROUTES,
    DENITRIFICATION_N2O_ROUTES,
)

ASSIMILATORY_NITRATE_REDUCTION_ROUTES = combine_stages(
    ASSIMILATORY_NITRATE_ROUTES,
    ASSIMILATORY_NITRITE_ROUTES,
)


# ---------------------------------------------------------------------
# Candidate rules
# ---------------------------------------------------------------------

RULES = [
    # ================================================================
    # 1. Nitrogen fixation
    # ================================================================
    {
        "capability": "Nitrogen fixation",
        "rule_name": "Nitrogen_fixation_NifHDK",
        "routes": [
            route("K02588", "K02586", "K02591"),
        ],
        "status": "Retain_candidate_complete_route",
        "evidence": "KEGG M00175",
        "rationale": (
            "Requires the conventional molybdenum nitrogenase "
            "NifH, NifD and NifK machinery."
        ),
    },
    {
        "capability": "Nitrogen fixation",
        "rule_name": "Nitrogen_fixation_VnfDGKH",
        "routes": [
            route("K22896", "K22897", "K22898", "K22899"),
        ],
        "status": "Retain_candidate_alternative_route",
        "evidence": "KEGG M00175",
        "rationale": (
            "Requires the complete vanadium-dependent nitrogenase "
            "alternative defined by KEGG."
        ),
    },
    {
        "capability": "Nitrogen fixation",
        "rule_name": "Nitrogen_fixation_preliminary",
        "routes": NITROGEN_FIXATION_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG M00175",
        "rationale": (
            "Positive when either conventional NifHDK or complete "
            "vanadium-dependent nitrogenase is detected."
        ),
    },

    # ================================================================
    # 2. Ammonium acquisition
    # ================================================================
    {
        "capability": "Ammonium acquisition",
        "rule_name": "Ammonium_acquisition_Amt_K03320",
        "routes": [
            route("K03320"),
        ],
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG K03320",
        "rationale": (
            "Amt-family ammonium transporters independently establish "
            "a recognised ammonium-acquisition mechanism."
        ),
    },

    # ================================================================
    # 3. Nitrate/nitrite acquisition
    # ================================================================
    {
        "capability": "Nitrate/nitrite acquisition",
        "rule_name": "Nitrate_nitrite_acquisition_NarK",
        "routes": [
            route("K02575"),
        ],
        "status": "Retain_diagnostic_single_transporter",
        "evidence": "KEGG K02575 and M00615",
        "rationale": (
            "A NarK/NrtP/NasA-family nitrate/nitrite transporter may "
            "independently establish uptake capability."
        ),
    },
    {
        "capability": "Nitrate/nitrite acquisition",
        "rule_name": "Nitrate_nitrite_acquisition_NrtABC",
        "routes": [
            route("K15576", "K15577", "K15578"),
        ],
        "status": "Retain_candidate_complete_transporter",
        "evidence": "KEGG M00615",
        "rationale": (
            "Requires substrate-binding, permease and ATP-binding "
            "components of the NrtABC uptake system."
        ),
    },
    {
        "capability": "Nitrate/nitrite acquisition",
        "rule_name": "Nitrate_nitrite_acquisition_preliminary",
        "routes": NITRATE_ACQUISITION_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG M00615",
        "rationale": (
            "Positive when NarK-family transport or complete NrtABC "
            "is detected."
        ),
    },

    # ================================================================
    # 4. Urea utilisation
    # ================================================================
    {
        "capability": "Urea utilisation",
        "rule_name": "Urea_hydrolysis_single_chain_K01427",
        "routes": [
            route("K01427"),
        ],
        "status": "Retain_candidate_single_enzyme",
        "evidence": "KEGG urease reaction R00131",
        "rationale": (
            "Single-chain urease can independently catalyse urea "
            "hydrolysis."
        ),
    },
    {
        "capability": "Urea utilisation",
        "rule_name": "Urea_hydrolysis_UreABC",
        "routes": [
            route("K01428", "K01429", "K01430"),
        ],
        "status": "Retain_candidate_complete_enzyme",
        "evidence": "KEGG urease reaction R00131",
        "rationale": (
            "Requires the alpha, beta and gamma urease subunits."
        ),
    },
    {
        "capability": "Urea utilisation",
        "rule_name": "Urea_hydrolysis_UreC_UreAB_fusion",
        "routes": [
            route("K01428", "K14048"),
        ],
        "status": "Retain_candidate_complete_enzyme",
        "evidence": "KEGG urease reaction R00131",
        "rationale": (
            "Alternative complete urease architecture containing "
            "UreC and the fused UreAB gamma/beta subunit."
        ),
    },
    {
        "capability": "Urea utilisation",
        "rule_name": "Urea_hydrolysis_preliminary",
        "routes": UREASE_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG urease reaction R00131",
        "rationale": (
            "Positive when any recognised complete urease architecture "
            "is detected."
        ),
    },
    {
        "capability": "Urea utilisation",
        "rule_name": "Urea_transport_UrtABCDE",
        "routes": [
            route(
                "K11959",
                "K11960",
                "K11961",
                "K11962",
                "K11963",
            ),
        ],
        "status": "Secondary_acquisition_trait",
        "evidence": "KEGG transporter hierarchy",
        "rationale": (
            "Complete five-component urea ABC transporter. It is "
            "audited separately because urea hydrolysis can occur "
            "without recovery of the complete transporter."
        ),
    },

    # ================================================================
    # 5. Assimilatory nitrate reduction
    # ================================================================
    {
        "capability": "Assimilatory nitrate reduction",
        "rule_name": "Assimilatory_nitrate_reductase_stage",
        "routes": ASSIMILATORY_NITRATE_ROUTES,
        "status": "Stage_diagnostic",
        "evidence": "KEGG M00531",
        "rationale": (
            "Presence of one accepted assimilatory nitrate-to-nitrite "
            "reduction route."
        ),
    },
    {
        "capability": "Assimilatory nitrate reduction",
        "rule_name": "Assimilatory_nitrite_reductase_stage",
        "routes": ASSIMILATORY_NITRITE_ROUTES,
        "status": "Stage_diagnostic",
        "evidence": "KEGG M00531",
        "rationale": (
            "Presence of one accepted assimilatory nitrite-to-ammonium "
            "reduction route."
        ),
    },
    {
        "capability": "Assimilatory nitrate reduction",
        "rule_name": "Assimilatory_nitrate_reduction_preliminary",
        "routes": ASSIMILATORY_NITRATE_REDUCTION_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG M00531",
        "rationale": (
            "Requires both nitrate-to-nitrite and nitrite-to-ammonium "
            "assimilation stages."
        ),
    },

    # ================================================================
    # 6. Nitrification
    # ================================================================
    {
        "capability": "Nitrification",
        "rule_name": "Nitrification_current_KEGG_M00528",
        "routes": NITRIFICATION_CURRENT_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "Current KEGG M00528",
        "rationale": (
            "Requires current ammonia monooxygenase subunits and "
            "hydroxylamine oxidoreductase."
        ),
    },
    {
        "capability": "Nitrification",
        "rule_name": "Nitrification_legacy_AmoABC_Hao",
        "routes": NITRIFICATION_LEGACY_ROUTES,
        "status": "Sensitivity_rule",
        "evidence": "Legacy amo/pmo KO architecture",
        "rationale": (
            "Uses the older shared AmoA/PmoA marker. This may not "
            "distinguish ammonia oxidation from methane oxidation and "
            "must therefore remain a sensitivity rule."
        ),
    },

    # ================================================================
    # 7. DNRA
    # ================================================================
    {
        "capability": "DNRA",
        "rule_name": "DNRA_nitrite_to_ammonium_core",
        "routes": DNRA_NITRITE_TO_AMMONIUM_ROUTES,
        "status": "MAG_tolerant_stage_diagnostic",
        "evidence": "KEGG M00530 second stage",
        "rationale": (
            "Detects NirBD or NrfAH but does not by itself establish "
            "complete nitrate-to-ammonium DNRA."
        ),
    },
    {
        "capability": "DNRA",
        "rule_name": "DNRA_complete_preliminary",
        "routes": DNRA_COMPLETE_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG M00530",
        "rationale": (
            "Requires respiratory nitrate reduction plus a complete "
            "nitrite-to-ammonium route."
        ),
    },

    # ================================================================
    # 8. Denitrification
    # ================================================================
    {
        "capability": "Denitrification",
        "rule_name": "Denitrification_nitrate_reduction_stage",
        "routes": NITRATE_RESPIRATION_ROUTES,
        "status": "Stage_diagnostic",
        "evidence": "KEGG M00529 stage 1",
        "rationale": "Respiratory nitrate reduction to nitrite.",
    },
    {
        "capability": "Denitrification",
        "rule_name": "Denitrification_nitrite_reduction_stage",
        "routes": DENITRIFICATION_NITRITE_ROUTES,
        "status": "Stage_diagnostic",
        "evidence": "KEGG M00529 stage 2",
        "rationale": "Nitrite reduction to nitric oxide by NirK or NirS.",
    },
    {
        "capability": "Denitrification",
        "rule_name": "Denitrification_NO_reduction_stage",
        "routes": DENITRIFICATION_NO_ROUTES,
        "status": "Stage_diagnostic",
        "evidence": "KEGG M00529 stage 3",
        "rationale": "Nitric-oxide reduction by complete NorBC.",
    },
    {
        "capability": "Denitrification",
        "rule_name": "Denitrification_N2O_reduction_stage",
        "routes": DENITRIFICATION_N2O_ROUTES,
        "status": "Stage_diagnostic",
        "evidence": "KEGG M00529 stage 4",
        "rationale": "Nitrous-oxide reduction by NosZ.",
    },
    {
        "capability": "Denitrification",
        "rule_name": "Denitrification_complete_preliminary",
        "routes": DENITRIFICATION_COMPLETE_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG M00529",
        "rationale": (
            "Requires all four denitrification stages from nitrate "
            "reduction through nitrous-oxide reduction."
        ),
    },

    # ================================================================
    # 9. Anammox
    # ================================================================
    {
        "capability": "Anammox",
        "rule_name": "Anammox_complete_preliminary",
        "routes": ANAMMOX_ROUTES,
        "status": "Preliminary_capability_rule",
        "evidence": "KEGG M00973",
        "rationale": (
            "Requires NirK or NirS, complete hydrazine synthase "
            "HzsABC and hydrazine dehydrogenase."
        ),
    },
]


# ---------------------------------------------------------------------
# Nine candidate capability rules used only for preliminary scoring
# ---------------------------------------------------------------------

PRELIMINARY_CAPABILITY_RULES = [
    "Nitrogen_fixation_preliminary",
    "Ammonium_acquisition_Amt_K03320",
    "Nitrate_nitrite_acquisition_preliminary",
    "Urea_hydrolysis_preliminary",
    "Assimilatory_nitrate_reduction_preliminary",
    "Nitrification_current_KEGG_M00528",
    "DNRA_complete_preliminary",
    "Denitrification_complete_preliminary",
    "Anammox_complete_preliminary",
]


DENITRIFICATION_STAGE_RULES = [
    "Denitrification_nitrate_reduction_stage",
    "Denitrification_nitrite_reduction_stage",
    "Denitrification_NO_reduction_stage",
    "Denitrification_N2O_reduction_stage",
]


# ---------------------------------------------------------------------
# Multigene component diagnostics
# ---------------------------------------------------------------------

MULTIGENE_SYSTEMS = [
    {
        "system": "Nitrogenase_NifHDK",
        "components": ("K02588", "K02586", "K02591"),
    },
    {
        "system": "Nitrogenase_VnfDGKH",
        "components": ("K22896", "K22897", "K22898", "K22899"),
    },
    {
        "system": "Nitrate_transporter_NrtABC",
        "components": ("K15576", "K15577", "K15578"),
    },
    {
        "system": "Urea_transporter_UrtABCDE",
        "components": (
            "K11959",
            "K11960",
            "K11961",
            "K11962",
            "K11963",
        ),
    },
    {
        "system": "Urease_UreABC",
        "components": ("K01428", "K01429", "K01430"),
    },
    {
        "system": "Current_ammonia_monooxygenase",
        "components": ("K28504", "K10945", "K10946"),
    },
    {
        "system": "Legacy_AmoABC_PmoABC",
        "components": ("K10944", "K10945", "K10946"),
    },
    {
        "system": "Respiratory_nitrate_reductase_NarGHI",
        "components": ("K00370", "K00371", "K00374"),
    },
    {
        "system": "Periplasmic_nitrate_reductase_NapAB",
        "components": ("K02567", "K02568"),
    },
    {
        "system": "DNRA_NirBD",
        "components": ("K00362", "K00363"),
    },
    {
        "system": "DNRA_NrfAH",
        "components": ("K03385", "K15876"),
    },
    {
        "system": "Nitric_oxide_reductase_NorBC",
        "components": ("K04561", "K02305"),
    },
    {
        "system": "Anammox_HzsABC",
        "components": ("K20932", "K20933", "K20934"),
    },
]


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def clean_text(value):
    if value is None:
        return ""

    return str(value).strip()


def safe_percent(numerator, denominator):
    if denominator == 0:
        return 0.0

    return numerator / denominator * 100.0


def validate_input_columns(fieldnames):
    if fieldnames is None:
        sys.exit(
            "ERROR: Input table does not contain a header."
        )

    required = {
        "Final_MAG",
        "KO",
        "KO_definition",
    }

    missing = sorted(
        required.difference(fieldnames)
    )

    if missing:
        sys.exit(
            "ERROR: Input table is missing required columns:\n"
            + "\n".join(missing)
        )


# ---------------------------------------------------------------------
# Read master gene table
# ---------------------------------------------------------------------

def read_master_table():
    try:
        handle = open(
            INPUT_FILE,
            "r",
            encoding="utf-8-sig",
            newline="",
        )
    except OSError as exc:
        sys.exit(
            "ERROR: Could not open master gene table:\n"
            f"{INPUT_FILE}\n"
            f"{exc}"
        )

    mag_kos = defaultdict(set)
    mag_metadata = {}
    gene_rows_per_mag = defaultdict(int)
    observed_ko_definitions = defaultdict(set)

    with handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        validate_input_columns(
            reader.fieldnames
        )

        metadata_columns = [
            column
            for column in METADATA_COLUMNS
            if column in reader.fieldnames
        ]

        for row in reader:
            mag = clean_text(
                row.get("Final_MAG")
            )

            if not mag:
                continue

            # Ensure MAG exists even if a row has no KO.
            mag_kos[mag]

            gene_rows_per_mag[mag] += 1

            if mag not in mag_metadata:
                mag_metadata[mag] = {
                    column: clean_text(
                        row.get(column)
                    )
                    for column in metadata_columns
                }

            ko_text = clean_text(
                row.get("KO")
            )

            ko_definition = clean_text(
                row.get("KO_definition")
            )

            for ko in KO_PATTERN.findall(ko_text):
                mag_kos[mag].add(
                    ko
                )

                if (
                    ko in CANDIDATE_KO_EXPECTATIONS
                    and ko_definition
                ):
                    observed_ko_definitions[ko].add(
                        ko_definition
                    )

    if not mag_kos:
        sys.exit(
            "ERROR: No MAGs were recovered from the master table."
        )

    return (
        mag_kos,
        mag_metadata,
        metadata_columns,
        gene_rows_per_mag,
        observed_ko_definitions,
    )


# ---------------------------------------------------------------------
# Calculate rule calls
# ---------------------------------------------------------------------

def calculate_rule_calls(mag_kos):
    calls = {}

    for mag, ko_set in mag_kos.items():
        calls[mag] = {}

        for rule_definition in RULES:
            rule_name = rule_definition["rule_name"]

            calls[mag][rule_name] = int(
                rule_is_positive(
                    ko_set,
                    rule_definition["routes"],
                )
            )

        denitrification_step_count = sum(
            calls[mag][stage_rule]
            for stage_rule in DENITRIFICATION_STAGE_RULES
        )

        calls[mag][
            "Denitrification_step_count_0_4"
        ] = denitrification_step_count

        calls[mag][
            "Denitrification_at_least_2_steps"
        ] = int(
            denitrification_step_count >= 2
        )

        preliminary_count = sum(
            calls[mag][rule_name]
            for rule_name in PRELIMINARY_CAPABILITY_RULES
        )

        calls[mag][
            "Nitrogen_preliminary_capability_count_0_9"
        ] = preliminary_count

        calls[mag][
            "Nitrogen_preliminary_strategy_score_0_1"
        ] = preliminary_count / 9.0

    return calls


# ---------------------------------------------------------------------
# Write rule definitions
# ---------------------------------------------------------------------

def write_rule_definitions():
    columns = [
        "Capability",
        "Rule_name",
        "Rule_logic",
        "Decision_status",
        "Evidence_basis",
        "Biological_rationale",
        "Included_in_preliminary_9_capability_score",
        "Overlap_with_redox_axis",
    ]

    with open(
        RULE_DEFINITIONS_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )

        writer.writeheader()

        for rule_definition in RULES:
            rule_name = rule_definition["rule_name"]

            overlap = (
                "Yes"
                if rule_definition["capability"]
                in {"DNRA", "Denitrification"}
                else "No"
            )

            writer.writerow(
                {
                    "Capability": rule_definition["capability"],
                    "Rule_name": rule_name,
                    "Rule_logic": format_rule_logic(
                        rule_definition["routes"]
                    ),
                    "Decision_status": rule_definition["status"],
                    "Evidence_basis": rule_definition["evidence"],
                    "Biological_rationale": (
                        rule_definition["rationale"]
                    ),
                    "Included_in_preliminary_9_capability_score": (
                        "Yes"
                        if rule_name in PRELIMINARY_CAPABILITY_RULES
                        else "No"
                    ),
                    "Overlap_with_redox_axis": overlap,
                }
            )


# ---------------------------------------------------------------------
# Write rule audit
# ---------------------------------------------------------------------

def write_rule_audit(mag_kos, calls):
    total_mags = len(mag_kos)

    columns = [
        "Capability",
        "Rule_name",
        "Rule_logic",
        "Positive_MAGs",
        "Negative_MAGs",
        "Total_MAGs",
        "Prevalence_percent",
        "Decision_status",
        "Evidence_basis",
        "Biological_rationale",
        "Included_in_preliminary_9_capability_score",
        "Overlap_with_redox_axis",
    ]

    with open(
        RULE_AUDIT_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )

        writer.writeheader()

        for rule_definition in RULES:
            rule_name = rule_definition["rule_name"]

            positive = sum(
                calls[mag][rule_name]
                for mag in mag_kos
            )

            overlap = (
                "Yes"
                if rule_definition["capability"]
                in {"DNRA", "Denitrification"}
                else "No"
            )

            writer.writerow(
                {
                    "Capability": rule_definition["capability"],
                    "Rule_name": rule_name,
                    "Rule_logic": format_rule_logic(
                        rule_definition["routes"]
                    ),
                    "Positive_MAGs": positive,
                    "Negative_MAGs": total_mags - positive,
                    "Total_MAGs": total_mags,
                    "Prevalence_percent": (
                        f"{safe_percent(positive, total_mags):.2f}"
                    ),
                    "Decision_status": rule_definition["status"],
                    "Evidence_basis": rule_definition["evidence"],
                    "Biological_rationale": (
                        rule_definition["rationale"]
                    ),
                    "Included_in_preliminary_9_capability_score": (
                        "Yes"
                        if rule_name in PRELIMINARY_CAPABILITY_RULES
                        else "No"
                    ),
                    "Overlap_with_redox_axis": overlap,
                }
            )


# ---------------------------------------------------------------------
# Write MAG-level presence matrix
# ---------------------------------------------------------------------

def write_presence_matrix(
    mag_kos,
    mag_metadata,
    metadata_columns,
    gene_rows_per_mag,
    calls,
):
    rule_columns = [
        rule_definition["rule_name"]
        for rule_definition in RULES
    ]

    diagnostic_columns = [
        "Denitrification_step_count_0_4",
        "Denitrification_at_least_2_steps",
        "Nitrogen_preliminary_capability_count_0_9",
        "Nitrogen_preliminary_strategy_score_0_1",
        "Genome_level_completeness_correction",
    ]

    output_columns = (
        ["Final_MAG"]
        + metadata_columns
        + ["Master_gene_rows"]
        + rule_columns
        + diagnostic_columns
    )

    with open(
        PRESENCE_MATRIX_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_columns,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for mag in sorted(mag_kos):
            output_row = {
                "Final_MAG": mag,
                "Master_gene_rows": gene_rows_per_mag.get(
                    mag,
                    0,
                ),
                "Genome_level_completeness_correction": (
                    "Not_applied"
                ),
            }

            output_row.update(
                mag_metadata.get(
                    mag,
                    {},
                )
            )

            output_row.update(
                calls[mag]
            )

            output_row[
                "Nitrogen_preliminary_strategy_score_0_1"
            ] = (
                f"{calls[mag]['Nitrogen_preliminary_strategy_score_0_1']:.4f}"
            )

            writer.writerow(
                output_row
            )


# ---------------------------------------------------------------------
# Write KO prevalence
# ---------------------------------------------------------------------

def write_ko_prevalence(mag_kos):
    total_mags = len(mag_kos)

    columns = [
        "KO",
        "Expected_annotation",
        "MAGs_detected",
        "MAGs_total",
        "Prevalence_percent",
    ]

    output_rows = []

    for ko in sorted(CANDIDATE_KO_EXPECTATIONS):
        positive = sum(
            1
            for ko_set in mag_kos.values()
            if ko in ko_set
        )

        output_rows.append(
            {
                "KO": ko,
                "Expected_annotation": (
                    CANDIDATE_KO_EXPECTATIONS[ko]
                ),
                "MAGs_detected": positive,
                "MAGs_total": total_mags,
                "Prevalence_percent": safe_percent(
                    positive,
                    total_mags,
                ),
            }
        )

    output_rows.sort(
        key=lambda row: (
            -row["MAGs_detected"],
            row["KO"],
        )
    )

    with open(
        KO_PREVALENCE_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )

        writer.writeheader()

        for row in output_rows:
            output_row = dict(row)

            output_row[
                "Prevalence_percent"
            ] = (
                f"{row['Prevalence_percent']:.2f}"
            )

            writer.writerow(
                output_row
            )


# ---------------------------------------------------------------------
# Write KO-definition audit
# ---------------------------------------------------------------------

def write_ko_definition_audit(
    mag_kos,
    observed_ko_definitions,
):
    total_mags = len(mag_kos)

    columns = [
        "KO",
        "Expected_annotation",
        "Observed_master_table_definitions",
        "Number_of_distinct_observed_definitions",
        "MAGs_detected",
        "MAGs_total",
        "Definition_present_in_master_table",
        "Manual_review_required",
    ]

    with open(
        KO_DEFINITION_AUDIT_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )

        writer.writeheader()

        for ko in sorted(CANDIDATE_KO_EXPECTATIONS):
            definitions = sorted(
                observed_ko_definitions.get(
                    ko,
                    set(),
                )
            )

            positive = sum(
                1
                for ko_set in mag_kos.values()
                if ko in ko_set
            )

            writer.writerow(
                {
                    "KO": ko,
                    "Expected_annotation": (
                        CANDIDATE_KO_EXPECTATIONS[ko]
                    ),
                    "Observed_master_table_definitions": (
                        " || ".join(definitions)
                    ),
                    "Number_of_distinct_observed_definitions": (
                        len(definitions)
                    ),
                    "MAGs_detected": positive,
                    "MAGs_total": total_mags,
                    "Definition_present_in_master_table": (
                        "Yes" if definitions else "No"
                    ),
                    "Manual_review_required": (
                        "Yes"
                        if positive > 0 and not definitions
                        else "No"
                    ),
                }
            )


# ---------------------------------------------------------------------
# Write component-recovery diagnostics
# ---------------------------------------------------------------------

def write_system_diagnostics(mag_kos):
    total_mags = len(mag_kos)

    columns = [
        "System",
        "Required_components",
        "Number_required",
        "Components_detected",
        "MAGs",
        "MAGs_total",
        "Prevalence_percent",
        "Recovery_class",
    ]

    with open(
        SYSTEM_DIAGNOSTICS_OUTPUT,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )

        writer.writeheader()

        for system_definition in MULTIGENE_SYSTEMS:
            system_name = system_definition["system"]
            components = system_definition["components"]
            number_required = len(components)

            distribution = defaultdict(int)

            for ko_set in mag_kos.values():
                number_detected = sum(
                    1
                    for ko in components
                    if ko in ko_set
                )

                distribution[number_detected] += 1

            for number_detected in range(
                number_required + 1
            ):
                mag_count = distribution[number_detected]

                if number_detected == 0:
                    recovery_class = "No_components"
                elif number_detected == number_required:
                    recovery_class = "Complete"
                else:
                    recovery_class = "Partial"

                writer.writerow(
                    {
                        "System": system_name,
                        "Required_components": (
                            " + ".join(components)
                        ),
                        "Number_required": number_required,
                        "Components_detected": number_detected,
                        "MAGs": mag_count,
                        "MAGs_total": total_mags,
                        "Prevalence_percent": (
                            f"{safe_percent(mag_count, total_mags):.2f}"
                        ),
                        "Recovery_class": recovery_class,
                    }
                )


# ---------------------------------------------------------------------
# Terminal summary
# ---------------------------------------------------------------------

def print_summary(mag_kos, calls):
    total_mags = len(mag_kos)

    print(
        "Preliminary nitrogen acquisition and transformation "
        "audit completed."
    )
    print()
    print(f"MAGs assessed: {total_mags}")
    print(
        "Candidate KOs assessed: "
        f"{len(CANDIDATE_KO_EXPECTATIONS)}"
    )
    print(
        f"Candidate rules assessed: {len(RULES)}"
    )
    print()
    print(f"Rule audit: {RULE_AUDIT_OUTPUT}")
    print(f"Presence matrix: {PRESENCE_MATRIX_OUTPUT}")
    print(f"KO prevalence: {KO_PREVALENCE_OUTPUT}")
    print(
        f"System diagnostics: {SYSTEM_DIAGNOSTICS_OUTPUT}"
    )
    print(
        f"Rule definitions: {RULE_DEFINITIONS_OUTPUT}"
    )
    print(
        f"KO-definition audit: {KO_DEFINITION_AUDIT_OUTPUT}"
    )
    print()

    print("Nine preliminary nitrogen capabilities")
    print("---------------------------------------")

    for rule_name in PRELIMINARY_CAPABILITY_RULES:
        positive = sum(
            calls[mag][rule_name]
            for mag in mag_kos
        )

        prevalence = safe_percent(
            positive,
            total_mags,
        )

        print(
            f"{rule_name:<58}"
            f"{positive:>5} MAGs "
            f"({prevalence:>6.2f}%)"
        )

    score_distribution = defaultdict(int)

    for mag in mag_kos:
        capability_count = calls[mag][
            "Nitrogen_preliminary_capability_count_0_9"
        ]

        score_distribution[
            capability_count
        ] += 1

    print()
    print("Preliminary nitrogen capability-count distribution")
    print("--------------------------------------------------")

    for capability_count in range(0, 10):
        mag_count = score_distribution[
            capability_count
        ]

        prevalence = safe_percent(
            mag_count,
            total_mags,
        )

        print(
            f"{capability_count}/9: "
            f"{mag_count} MAGs "
            f"({prevalence:.2f}%)"
        )

    denitrification_distribution = defaultdict(int)

    for mag in mag_kos:
        step_count = calls[mag][
            "Denitrification_step_count_0_4"
        ]

        denitrification_distribution[
            step_count
        ] += 1

    print()
    print("Denitrification stage-count distribution")
    print("----------------------------------------")

    for step_count in range(0, 5):
        mag_count = denitrification_distribution[
            step_count
        ]

        prevalence = safe_percent(
            mag_count,
            total_mags,
        )

        print(
            f"{step_count}/4 stages: "
            f"{mag_count} MAGs "
            f"({prevalence:.2f}%)"
        )

    print()
    print(
        "The 0-9 score is preliminary and must not yet be used "
        "as the final nitrogen strategy."
    )
    print(
        "No individual-MAG completeness correction was applied."
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    os.makedirs(
        TABLE_OUTPUT_DIR,
        exist_ok=True,
    )

    (
        mag_kos,
        mag_metadata,
        metadata_columns,
        gene_rows_per_mag,
        observed_ko_definitions,
    ) = read_master_table()

    calls = calculate_rule_calls(
        mag_kos
    )

    write_rule_definitions()

    write_rule_audit(
        mag_kos,
        calls,
    )

    write_presence_matrix(
        mag_kos,
        mag_metadata,
        metadata_columns,
        gene_rows_per_mag,
        calls,
    )

    write_ko_prevalence(
        mag_kos
    )

    write_ko_definition_audit(
        mag_kos,
        observed_ko_definitions,
    )

    write_system_diagnostics(
        mag_kos
    )

    print_summary(
        mag_kos,
        calls,
    )


if __name__ == "__main__":
    main()
