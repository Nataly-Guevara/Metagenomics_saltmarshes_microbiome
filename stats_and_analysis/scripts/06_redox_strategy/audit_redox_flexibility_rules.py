#!/usr/bin/env python3

"""
Audit candidate rules for genome-resolved redox flexibility.

The intended Ricci-style redox-flexibility axis contains eight
binary biological capabilities:

1. Nitrate respiration
2. Nitrite and/or nitric-oxide reduction
3. Dissimilatory nitrate reduction to ammonium (DNRA)
4. Nitrous-oxide reduction
5. Dissimilatory sulfate reduction
6. Sulfide oxidation
7. Sulfite oxidation
8. Low-oxygen respiration

This script does NOT calculate the final redox-flexibility score.

Instead, it:

1. Reads the final 373-MAG gene master table.
2. Collects all detected KOs within each MAG.
3. Tests strict pathway rules and alternative diagnostic rules.
4. Calculates exact within-MAG component co-occurrence.
5. Preserves MAG metadata and GTDB taxonomy.
6. Produces KO-definition and rule-definition audits.
7. Creates a preliminary eight-capability presence matrix.

All outputs receive new filenames. Existing osmotic and functional
analysis files are not modified.

Only the Python standard library is used.
"""

import csv
import os
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------

PROJECT_ROOT = "/path/to/metagenomics_project"

INPUT_FILE = os.path.join(
    PROJECT_ROOT,
    "results/functional_annotation/kofamscan/"
    "final_373_MAGs_gene_master_table.tsv",
)

OUTPUT_ROOT = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/redox_flexibility",
)

TABLE_OUTPUT_DIR = os.path.join(
    OUTPUT_ROOT,
    "tables",
)

RULE_AUDIT_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "redox_capability_rule_audit_preliminary.tsv",
)

PRESENCE_MATRIX_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "redox_capability_presence_preliminary.tsv",
)

KO_PREVALENCE_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "redox_candidate_KO_prevalence.tsv",
)

SYSTEM_DIAGNOSTICS_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "redox_multigene_system_diagnostics.tsv",
)

RULE_DEFINITIONS_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "redox_rule_definitions_preliminary.tsv",
)

KO_DEFINITION_AUDIT_OUTPUT = os.path.join(
    TABLE_OUTPUT_DIR,
    "redox_candidate_KO_definition_audit.tsv",
)


# ---------------------------------------------------------------------
# Metadata preserved for every MAG
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


# ---------------------------------------------------------------------
# Candidate KO annotations
#
# These expected descriptions are for auditing only. The script also
# extracts the actual KO_definition values present in the master table.
# ---------------------------------------------------------------------

CANDIDATE_KO_EXPECTATIONS = {
    # Nitrate respiration
    "K00370": "narG/narZ/nxrA; nitrate reductase alpha subunit",
    "K00371": "narH/narY/nxrB; nitrate reductase beta subunit",
    "K00374": "narI/narV; nitrate reductase gamma subunit",
    "K02567": "napA; periplasmic nitrate reductase catalytic subunit",
    "K02568": "napB; periplasmic nitrate reductase cytochrome subunit",

    # Nitrite and nitric-oxide reduction
    "K00368": "nirK; copper-containing nitrite reductase",
    "K15864": "nirS; cytochrome cd1 nitrite reductase",
    "K04561": "norB; nitric-oxide reductase subunit B",
    "K02305": "norC; nitric-oxide reductase subunit C",

    # DNRA
    "K00362": "nirB; NADH-dependent nitrite reductase large subunit",
    "K00363": "nirD; NADH-dependent nitrite reductase small subunit",
    "K03385": "nrfA; cytochrome-c nitrite reductase",
    "K15876": "nrfH; cytochrome-c nitrite reductase electron-transfer subunit",

    # Nitrous-oxide reduction
    "K00376": "nosZ; nitrous-oxide reductase",

    # Dissimilatory sulfate reduction
    "K00958": "sat; sulfate adenylyltransferase",
    "K00394": "aprA; adenylylsulfate reductase alpha subunit",
    "K00395": "aprB; adenylylsulfate reductase beta subunit",
    "K11180": "dsrA; dissimilatory sulfite reductase alpha subunit",
    "K11181": "dsrB; dissimilatory sulfite reductase beta subunit",
    "K27187": "DsrMKJOP complex subunit M",
    "K27188": "DsrMKJOP complex subunit K",
    "K27189": "DsrMKJOP complex subunit J",
    "K27190": "DsrMKJOP complex subunit O",
    "K27191": "DsrMKJOP complex subunit P",
    "K27196": "dsrL; reverse-Dsr-associated flavoprotein",

    # Sulfide oxidation
    "K17218": "sqr; sulfide:quinone oxidoreductase",
    "K17229": "fccB; flavocytochrome-c sulfide dehydrogenase flavoprotein",
    "K17230": "fccA; flavocytochrome-c sulfide dehydrogenase cytochrome",

    # Sulfite oxidation
    "K05301": "sorA; sulfite dehydrogenase subunit A",
    "K00386": "sorB; sulfite dehydrogenase subunit B",
    "K21307": "soeA; sulfite dehydrogenase quinone subunit A",
    "K21308": "soeB; sulfite dehydrogenase quinone subunit B",
    "K21309": "soeC; sulfite dehydrogenase quinone subunit C",
    "K00387": "SUOX; sulfite oxidase",

    # Low-oxygen respiration
    "K00425": "cydA; cytochrome-bd ubiquinol oxidase subunit I",
    "K00426": "cydB; cytochrome-bd ubiquinol oxidase subunit II",
    "K00424": "cydX; cytochrome-bd accessory subunit X",
    "K22501": "alternative cytochrome-bd small accessory subunit",
    "K00404": "ccoN; cbb3-type cytochrome-c oxidase subunit I",
    "K00405": "ccoO; cbb3-type cytochrome-c oxidase subunit II",
    "K00406": "ccoP; cbb3-type cytochrome-c oxidase subunit III",
    "K00407": "ccoQ; cbb3-type cytochrome-c oxidase subunit IV",
    "K15862": "ccoNO; fused cbb3-type oxidase subunits I and II",
}


KO_PATTERN = re.compile(r"K\d{5}")


# ---------------------------------------------------------------------
# Rule helper
#
# Each route has:
#
# required:
#     Every listed KO must occur in the same MAG.
#
# forbidden:
#     None of the listed KOs may occur in the MAG.
#
# A rule is positive if at least one alternative route is satisfied.
# ---------------------------------------------------------------------

def route(required, forbidden=()):
    return {
        "required": tuple(required),
        "forbidden": tuple(forbidden),
    }


# ---------------------------------------------------------------------
# Candidate rule definitions
# ---------------------------------------------------------------------

RULES = [
    # ================================================================
    # 1. Nitrate respiration
    # ================================================================
    {
        "capability": "Nitrate respiration",
        "rule_name": "Nitrate_NarGHI_strict",
        "routes": [
            route(("K00370", "K00371", "K00374")),
        ],
        "status": "Retain_candidate_complete_system",
        "evidence_basis": "KEGG M00529 and M00530",
        "rationale": (
            "Requires the alpha, beta and membrane subunits of the "
            "respiratory Nar nitrate-reductase complex."
        ),
    },
    {
        "capability": "Nitrate respiration",
        "rule_name": "Nitrate_NarGH_core",
        "routes": [
            route(("K00370", "K00371")),
        ],
        "status": "MAG_tolerant_diagnostic",
        "evidence_basis": "Nar catalytic core",
        "rationale": (
            "Tests the NarGH catalytic core without requiring the "
            "membrane subunit. This is diagnostic only until compared "
            "with the strict NarGHI result."
        ),
    },
    {
        "capability": "Nitrate respiration",
        "rule_name": "Nitrate_NapAB",
        "routes": [
            route(("K02567", "K02568")),
        ],
        "status": "Retain_candidate_complete_system",
        "evidence_basis": "KEGG M00529 and M00530",
        "rationale": (
            "Requires the catalytic and cytochrome subunits of the "
            "periplasmic nitrate-reductase route."
        ),
    },
    {
        "capability": "Nitrate respiration",
        "rule_name": "Redox_Nitrate_respiration_preliminary",
        "routes": [
            route(("K00370", "K00371", "K00374")),
            route(("K02567", "K02568")),
        ],
        "status": "Preliminary_capability_rule",
        "evidence_basis": "KEGG M00529 and M00530",
        "rationale": (
            "Positive when either complete NarGHI or complete NapAB "
            "is detected."
        ),
    },

    # ================================================================
    # 2. Nitrite and/or nitric-oxide reduction
    # ================================================================
    {
        "capability": "Nitrite and nitric-oxide reduction",
        "rule_name": "Nitrite_to_NO_NirK",
        "routes": [
            route(("K00368",)),
        ],
        "status": "Retain_diagnostic_single_gene",
        "evidence_basis": "KEGG M00529",
        "rationale": (
            "NirK is an accepted diagnostic enzyme for reduction of "
            "nitrite to nitric oxide."
        ),
    },
    {
        "capability": "Nitrite and nitric-oxide reduction",
        "rule_name": "Nitrite_to_NO_NirS",
        "routes": [
            route(("K15864",)),
        ],
        "status": "Retain_diagnostic_single_gene",
        "evidence_basis": "KEGG M00529",
        "rationale": (
            "NirS is an alternative diagnostic enzyme for reduction "
            "of nitrite to nitric oxide."
        ),
    },
    {
        "capability": "Nitrite and nitric-oxide reduction",
        "rule_name": "NO_reduction_NorB_diagnostic",
        "routes": [
            route(("K04561",)),
        ],
        "status": "Diagnostic_sensitivity_rule",
        "evidence_basis": "NorB catalytic subunit",
        "rationale": (
            "Tests the catalytic NorB marker alone, corresponding to "
            "the permissive marker-based approach used in some "
            "functional-module studies."
        ),
    },
    {
        "capability": "Nitrite and nitric-oxide reduction",
        "rule_name": "NO_reduction_NorBC_complete",
        "routes": [
            route(("K04561", "K02305")),
        ],
        "status": "Retain_candidate_complete_system",
        "evidence_basis": "KEGG M00529",
        "rationale": (
            "Requires both nitric-oxide reductase subunits specified "
            "by the KEGG denitrification module."
        ),
    },
    {
        "capability": "Nitrite and nitric-oxide reduction",
        "rule_name": "Redox_Nitrite_NO_reduction_preliminary",
        "routes": [
            route(("K00368",)),
            route(("K15864",)),
            route(("K04561", "K02305")),
        ],
        "status": "Preliminary_capability_rule",
        "evidence_basis": "KEGG M00529",
        "rationale": (
            "Positive when a MAG encodes either accepted nitrite-to-NO "
            "reductase or the complete NorBC nitric-oxide reductase."
        ),
    },

    # ================================================================
    # 3. DNRA
    # ================================================================
    {
        "capability": "DNRA",
        "rule_name": "DNRA_NirBD",
        "routes": [
            route(("K00362", "K00363")),
        ],
        "status": "Retain_candidate_complete_route",
        "evidence_basis": "KEGG M00530",
        "rationale": (
            "Requires both subunits of the NADH-dependent nitrite "
            "reductase route."
        ),
    },
    {
        "capability": "DNRA",
        "rule_name": "DNRA_NrfA_diagnostic",
        "routes": [
            route(("K03385",)),
        ],
        "status": "Diagnostic_sensitivity_rule",
        "evidence_basis": "NrfA catalytic marker",
        "rationale": (
            "Tests NrfA alone as a commonly used diagnostic marker "
            "for respiratory nitrite reduction to ammonium."
        ),
    },
    {
        "capability": "DNRA",
        "rule_name": "DNRA_NrfAH",
        "routes": [
            route(("K03385", "K15876")),
        ],
        "status": "Retain_candidate_complete_route",
        "evidence_basis": "KEGG M00530",
        "rationale": (
            "Requires the NrfA catalytic enzyme and NrfH "
            "electron-transfer component."
        ),
    },
    {
        "capability": "DNRA",
        "rule_name": "Redox_DNRA_preliminary",
        "routes": [
            route(("K00362", "K00363")),
            route(("K03385", "K15876")),
        ],
        "status": "Preliminary_capability_rule",
        "evidence_basis": "KEGG M00530",
        "rationale": (
            "Positive when either the complete NirBD or complete "
            "NrfAH route is detected."
        ),
    },

    # ================================================================
    # 4. Nitrous-oxide reduction
    # ================================================================
    {
        "capability": "Nitrous-oxide reduction",
        "rule_name": "Redox_Nitrous_oxide_reduction_preliminary",
        "routes": [
            route(("K00376",)),
        ],
        "status": "Preliminary_capability_rule",
        "evidence_basis": "KEGG M00529",
        "rationale": (
            "NosZ is the diagnostic catalytic enzyme for reduction of "
            "nitrous oxide to dinitrogen."
        ),
    },

    # ================================================================
    # 5. Dissimilatory sulfate reduction
    # ================================================================
    {
        "capability": "Dissimilatory sulfate reduction",
        "rule_name": "Sulfate_reduction_DsrAB_diagnostic",
        "routes": [
            route(("K11180", "K11181")),
        ],
        "status": "Diagnostic_not_direction_specific",
        "evidence_basis": "KEGG M00596 and sulfur metabolism",
        "rationale": (
            "DsrAB indicates dissimilatory sulfur metabolism but may "
            "operate in the reverse direction in sulfur oxidisers. "
            "It must not alone define sulfate reduction."
        ),
    },
    {
        "capability": "Dissimilatory sulfate reduction",
        "rule_name": "Sulfate_reduction_legacy_core",
        "routes": [
            route(
                (
                    "K00958",
                    "K00394",
                    "K00395",
                    "K11180",
                    "K11181",
                ),
                forbidden=("K27196",),
            ),
        ],
        "status": "Primary_candidate_MAG_compatible",
        "evidence_basis": "Core KEGG M00596 reactions",
        "rationale": (
            "Requires sulfate activation, APS reduction and DsrAB, "
            "while excluding the reverse-Dsr-associated DsrL marker."
        ),
    },
    {
        "capability": "Dissimilatory sulfate reduction",
        "rule_name": "Sulfate_reduction_current_KEGG_strict",
        "routes": [
            route(
                (
                    "K00958",
                    "K00394",
                    "K00395",
                    "K11180",
                    "K11181",
                    "K27187",
                    "K27188",
                    "K27189",
                    "K27190",
                    "K27191",
                ),
                forbidden=("K27196",),
            ),
        ],
        "status": "Strict_sensitivity_rule",
        "evidence_basis": "Current KEGG M00596",
        "rationale": (
            "Implements the current complete KEGG module, including "
            "the DsrMKJOP terminal complex and exclusion of DsrL."
        ),
    },
    {
        "capability": "Dissimilatory sulfate reduction",
        "rule_name": "Reverse_Dsr_indicator_K27196",
        "routes": [
            route(("K27196",)),
        ],
        "status": "Directionality_audit_only",
        "evidence_basis": "KEGG M00596 exclusion term",
        "rationale": (
            "DsrL is audited as an indicator that DsrAB may participate "
            "in oxidative rather than sulfate-reducing metabolism."
        ),
    },
    {
        "capability": "Dissimilatory sulfate reduction",
        "rule_name": "Redox_Sulfate_reduction_preliminary",
        "routes": [
            route(
                (
                    "K00958",
                    "K00394",
                    "K00395",
                    "K11180",
                    "K11181",
                ),
                forbidden=("K27196",),
            ),
        ],
        "status": "Preliminary_capability_rule",
        "evidence_basis": "Core KEGG M00596",
        "rationale": (
            "Uses the MAG-compatible complete biochemical core and "
            "excludes MAGs carrying the reverse-Dsr-associated marker."
        ),
    },

    # ================================================================
    # 6. Sulfide oxidation
    # ================================================================
    {
        "capability": "Sulfide oxidation",
        "rule_name": "Sulfide_oxidation_SQR",
        "routes": [
            route(("K17218",)),
        ],
        "status": "Retain_diagnostic_single_enzyme",
        "evidence_basis": "KEGG sulfur metabolism",
        "rationale": (
            "Sulfide:quinone oxidoreductase independently catalyses "
            "oxidation of sulfide."
        ),
    },
    {
        "capability": "Sulfide oxidation",
        "rule_name": "Sulfide_oxidation_FccAB",
        "routes": [
            route(("K17229", "K17230")),
        ],
        "status": "Retain_candidate_complete_route",
        "evidence_basis": "KEGG M00985",
        "rationale": (
            "Requires both flavoprotein and cytochrome components of "
            "flavocytochrome-c sulfide dehydrogenase."
        ),
    },
    {
        "capability": "Sulfide oxidation",
        "rule_name": "Redox_Sulfide_oxidation_preliminary",
        "routes": [
            route(("K17218",)),
            route(("K17229", "K17230")),
        ],
        "status": "Preliminary_capability_rule",
        "evidence_basis": "KEGG sulfur metabolism and M00985",
        "rationale": (
            "Positive when either SQR or the complete FccAB route is "
            "detected."
        ),
    },

    # ================================================================
    # 7. Sulfite oxidation
    # ================================================================
    {
        "capability": "Sulfite oxidation",
        "rule_name": "Sulfite_oxidation_SorAB",
        "routes": [
            route(("K05301", "K00386")),
        ],
        "status": "Retain_candidate_complete_route",
        "evidence_basis": "KEGG sulfur metabolism",
        "rationale": (
            "Requires both subunits of cytochrome-linked sulfite "
            "dehydrogenase."
        ),
    },
    {
        "capability": "Sulfite oxidation",
        "rule_name": "Sulfite_oxidation_SoeABC",
        "routes": [
            route(("K21307", "K21308", "K21309")),
        ],
        "status": "Retain_candidate_complete_route",
        "evidence_basis": "KEGG sulfur metabolism",
        "rationale": (
            "Requires all three subunits of quinone-linked sulfite "
            "dehydrogenase."
        ),
    },
    {
        "capability": "Sulfite oxidation",
        "rule_name": "Sulfite_oxidation_SUOX",
        "routes": [
            route(("K00387",)),
        ],
        "status": "Retain_diagnostic_single_enzyme",
        "evidence_basis": "KEGG sulfur metabolism",
        "rationale": (
            "Sulfite oxidase can independently establish a sulfite-"
            "oxidation route."
        ),
    },
    {
        "capability": "Sulfite oxidation",
        "rule_name": "Redox_Sulfite_oxidation_preliminary",
        "routes": [
            route(("K05301", "K00386")),
            route(("K21307", "K21308", "K21309")),
            route(("K00387",)),
        ],
        "status": "Preliminary_capability_rule",
        "evidence_basis": "KEGG sulfur metabolism",
        "rationale": (
            "Positive when SorAB, SoeABC or sulfite oxidase is "
            "detected."
        ),
    },

    # ================================================================
    # 8. Low-oxygen respiration
    # ================================================================
    {
        "capability": "Low-oxygen respiration",
        "rule_name": "Low_oxygen_Cytochrome_bd_core",
        "routes": [
            route(("K00425", "K00426")),
        ],
        "status": "MAG_tolerant_diagnostic",
        "evidence_basis": "Cytochrome-bd catalytic core",
        "rationale": (
            "Tests the two principal cytochrome-bd oxidase subunits "
            "without requiring a small accessory subunit."
        ),
    },
    {
        "capability": "Low-oxygen respiration",
        "rule_name": "Low_oxygen_Cytochrome_bd_complete",
        "routes": [
            route(("K00425", "K00426", "K00424")),
            route(("K00425", "K00426", "K22501")),
        ],
        "status": "Retain_candidate_complete_system",
        "evidence_basis": "KEGG M00153",
        "rationale": (
            "Implements the complete KEGG cytochrome-bd module with "
            "either accepted small accessory subunit."
        ),
    },
    {
        "capability": "Low-oxygen respiration",
        "rule_name": "Low_oxygen_cbb3_standard",
        "routes": [
            route(("K00404", "K00405", "K00406", "K00407")),
        ],
        "status": "Retain_candidate_complete_system",
        "evidence_basis": "KEGG M00156",
        "rationale": (
            "Requires all four standard cbb3-type cytochrome-c "
            "oxidase subunits."
        ),
    },
    {
        "capability": "Low-oxygen respiration",
        "rule_name": "Low_oxygen_cbb3_fused",
        "routes": [
            route(("K15862", "K00406", "K00407")),
        ],
        "status": "Retain_candidate_complete_system",
        "evidence_basis": "KEGG M00156",
        "rationale": (
            "Alternative complete cbb3 route with fused CcoNO and "
            "separate CcoP and CcoQ."
        ),
    },
    {
        "capability": "Low-oxygen respiration",
        "rule_name": "Redox_Low_oxygen_respiration_strict",
        "routes": [
            route(("K00425", "K00426", "K00424")),
            route(("K00425", "K00426", "K22501")),
            route(("K00404", "K00405", "K00406", "K00407")),
            route(("K15862", "K00406", "K00407")),
        ],
        "status": "Strict_preliminary_capability_rule",
        "evidence_basis": "KEGG M00153 and M00156",
        "rationale": (
            "Positive when a complete cytochrome-bd or complete "
            "cbb3-type oxidase module is detected."
        ),
    },
    {
        "capability": "Low-oxygen respiration",
        "rule_name": "Redox_Low_oxygen_respiration_MAG_tolerant",
        "routes": [
            route(("K00425", "K00426")),
            route(("K00404", "K00405", "K00406", "K00407")),
            route(("K15862", "K00406", "K00407")),
        ],
        "status": "MAG_tolerant_sensitivity_rule",
        "evidence_basis": "KEGG M00153 and M00156",
        "rationale": (
            "Allows the cytochrome-bd catalytic core while retaining "
            "complete cbb3 requirements."
        ),
    },
]


# ---------------------------------------------------------------------
# The eight candidate capability columns
# ---------------------------------------------------------------------

PRELIMINARY_CAPABILITY_RULES = [
    "Redox_Nitrate_respiration_preliminary",
    "Redox_Nitrite_NO_reduction_preliminary",
    "Redox_DNRA_preliminary",
    "Redox_Nitrous_oxide_reduction_preliminary",
    "Redox_Sulfate_reduction_preliminary",
    "Redox_Sulfide_oxidation_preliminary",
    "Redox_Sulfite_oxidation_preliminary",
    "Redox_Low_oxygen_respiration_strict",
]


# ---------------------------------------------------------------------
# Multigene systems for recovery diagnostics
# ---------------------------------------------------------------------

MULTIGENE_SYSTEMS = [
    {
        "system": "NarGHI",
        "components": ("K00370", "K00371", "K00374"),
    },
    {
        "system": "NapAB",
        "components": ("K02567", "K02568"),
    },
    {
        "system": "NorBC",
        "components": ("K04561", "K02305"),
    },
    {
        "system": "NirBD",
        "components": ("K00362", "K00363"),
    },
    {
        "system": "NrfAH",
        "components": ("K03385", "K15876"),
    },
    {
        "system": "Sulfate_reduction_core",
        "components": (
            "K00958",
            "K00394",
            "K00395",
            "K11180",
            "K11181",
        ),
    },
    {
        "system": "DsrMKJOP",
        "components": (
            "K27187",
            "K27188",
            "K27189",
            "K27190",
            "K27191",
        ),
    },
    {
        "system": "FccAB",
        "components": ("K17229", "K17230"),
    },
    {
        "system": "SorAB",
        "components": ("K05301", "K00386"),
    },
    {
        "system": "SoeABC",
        "components": ("K21307", "K21308", "K21309"),
    },
    {
        "system": "Cytochrome_bd_core",
        "components": ("K00425", "K00426"),
    },
    {
        "system": "Cytochrome_bd_CydABX",
        "components": ("K00425", "K00426", "K00424"),
    },
    {
        "system": "Cytochrome_bd_alternative_accessory",
        "components": ("K00425", "K00426", "K22501"),
    },
    {
        "system": "cbb3_standard",
        "components": ("K00404", "K00405", "K00406", "K00407"),
    },
    {
        "system": "cbb3_fused_CcoNO",
        "components": ("K15862", "K00406", "K00407"),
    },
]


# ---------------------------------------------------------------------
# General utility functions
# ---------------------------------------------------------------------

def clean_text(value):
    if value is None:
        return ""

    return str(value).strip()


def safe_percent(numerator, denominator):
    if denominator == 0:
        return 0.0

    return numerator / denominator * 100.0


def route_is_positive(ko_set, route_definition):
    required = route_definition["required"]
    forbidden = route_definition["forbidden"]

    required_present = all(
        ko in ko_set
        for ko in required
    )

    forbidden_absent = all(
        ko not in ko_set
        for ko in forbidden
    )

    return required_present and forbidden_absent


def rule_is_positive(ko_set, routes):
    return any(
        route_is_positive(
            ko_set,
            route_definition,
        )
        for route_definition in routes
    )


def format_route(route_definition):
    required_text = " + ".join(
        route_definition["required"]
    )

    forbidden = route_definition["forbidden"]

    if forbidden:
        forbidden_text = " + ".join(
            f"NOT {ko}"
            for ko in forbidden
        )

        return (
            f"({required_text}) AND ({forbidden_text})"
        )

    if len(route_definition["required"]) > 1:
        return f"({required_text})"

    return required_text


def format_rule_logic(routes):
    return " OR ".join(
        format_route(route_definition)
        for route_definition in routes
    )


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
# Read the master table
# ---------------------------------------------------------------------

def read_master_table():
    """
    Return:
        mag_kos
        mag_metadata
        metadata_columns
        observed_ko_definitions
        gene_rows_per_mag
    """
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
    observed_ko_definitions = defaultdict(set)
    gene_rows_per_mag = defaultdict(int)

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
                mag_kos[mag].add(ko)

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
        observed_ko_definitions,
        gene_rows_per_mag,
    )


# ---------------------------------------------------------------------
# Calculate all binary rule calls
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

        calls[mag][
            "Redox_preliminary_capability_count_0_8"
        ] = sum(
            calls[mag][rule_name]
            for rule_name in PRELIMINARY_CAPABILITY_RULES
        )

        calls[mag][
            "Redox_preliminary_strategy_score_0_1"
        ] = (
            calls[mag][
                "Redox_preliminary_capability_count_0_8"
            ]
            / 8.0
        )

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
        "Included_in_preliminary_8_capability_score",
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

            writer.writerow(
                {
                    "Capability": (
                        rule_definition["capability"]
                    ),
                    "Rule_name": rule_name,
                    "Rule_logic": format_rule_logic(
                        rule_definition["routes"]
                    ),
                    "Decision_status": (
                        rule_definition["status"]
                    ),
                    "Evidence_basis": (
                        rule_definition["evidence_basis"]
                    ),
                    "Biological_rationale": (
                        rule_definition["rationale"]
                    ),
                    "Included_in_preliminary_8_capability_score": (
                        "Yes"
                        if rule_name
                        in PRELIMINARY_CAPABILITY_RULES
                        else "No"
                    ),
                }
            )


# ---------------------------------------------------------------------
# Write rule-level audit
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
        "Included_in_preliminary_8_capability_score",
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

            writer.writerow(
                {
                    "Capability": (
                        rule_definition["capability"]
                    ),
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
                    "Decision_status": (
                        rule_definition["status"]
                    ),
                    "Evidence_basis": (
                        rule_definition["evidence_basis"]
                    ),
                    "Biological_rationale": (
                        rule_definition["rationale"]
                    ),
                    "Included_in_preliminary_8_capability_score": (
                        "Yes"
                        if rule_name
                        in PRELIMINARY_CAPABILITY_RULES
                        else "No"
                    ),
                }
            )


# ---------------------------------------------------------------------
# Write MAG-by-rule presence matrix
# ---------------------------------------------------------------------

def write_presence_matrix(
    mag_kos,
    mag_metadata,
    metadata_columns,
    calls,
    gene_rows_per_mag,
):
    rule_columns = [
        rule_definition["rule_name"]
        for rule_definition in RULES
    ]

    output_columns = (
        ["Final_MAG"]
        + metadata_columns
        + ["Master_gene_rows"]
        + rule_columns
        + [
            "Redox_preliminary_capability_count_0_8",
            "Redox_preliminary_strategy_score_0_1",
            "Genome_level_completeness_correction",
        ]
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
                "Redox_preliminary_strategy_score_0_1"
            ] = (
                f"{calls[mag]['Redox_preliminary_strategy_score_0_1']:.4f}"
            )

            writer.writerow(
                output_row
            )


# ---------------------------------------------------------------------
# Write individual KO prevalence
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
# Write observed-versus-expected KO definition audit
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
                mag_count = distribution[
                    number_detected
                ]

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

    print("Preliminary redox-flexibility audit completed.")
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
    print(
        f"Rule audit: {RULE_AUDIT_OUTPUT}"
    )
    print(
        f"Presence matrix: {PRESENCE_MATRIX_OUTPUT}"
    )
    print(
        f"KO prevalence: {KO_PREVALENCE_OUTPUT}"
    )
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

    print("Eight preliminary redox capabilities")
    print("------------------------------------")

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
            f"{rule_name:<52}"
            f"{positive:>5} MAGs "
            f"({prevalence:>6.2f}%)"
        )

    score_distribution = defaultdict(int)

    for mag in mag_kos:
        capability_count = calls[mag][
            "Redox_preliminary_capability_count_0_8"
        ]

        score_distribution[
            capability_count
        ] += 1

    print()
    print("Preliminary redox capability-count distribution")
    print("-----------------------------------------------")

    for capability_count in range(0, 9):
        mag_count = score_distribution[
            capability_count
        ]

        prevalence = safe_percent(
            mag_count,
            total_mags,
        )

        print(
            f"{capability_count}/8: "
            f"{mag_count} MAGs "
            f"({prevalence:.2f}%)"
        )

    print()
    print(
        "These scores are preliminary audit outputs and must not "
        "yet be used as the final redox-flexibility strategy."
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
        observed_ko_definitions,
        gene_rows_per_mag,
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
        calls,
        gene_rows_per_mag,
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

