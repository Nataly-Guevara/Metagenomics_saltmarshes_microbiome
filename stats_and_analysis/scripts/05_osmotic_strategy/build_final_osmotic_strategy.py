#!/usr/bin/env python3

"""
Build final osmotic-adaptation ecological strategy tables.

Framework:
    KO annotations
        ->
    (1) module/pathway evidence
        ->
    (2) biological capability calls (0/1)
        ->
    (3) standardized osmotic-strategy breadth (0-1)

The five equally weighted osmotic capabilities are:

    1. Trehalose biosynthesis
    2. Compatible-solute uptake
    3. Ectoine biosynthesis
    4. Mechanosensitive osmotic protection
    5. Ion homeostasis

Important:
- Individual KOs are NOT counted as independent ecological traits.
- Alternative biological routes are collapsed into one capability.
- A genome carrying multiple alternative routes does not receive
  extra weight for the same ecological capability.
- Module-evidence values are retained for audit/QC.
- A stricter trehalose definition is retained as a sensitivity analysis.
- Individual MAG strategy scores are NOT corrected for genome completeness.
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
      "osmotic_adaptation"
)

TABLE_DIR = OUTPUT_DIR / "tables"
ITOL_DIR = OUTPUT_DIR / "itol"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
ITOL_DIR.mkdir(parents=True, exist_ok=True)

FINAL_TABLE = (
    TABLE_DIR
    / "osmotic_strategy_final_373_MAGs.tsv"
)

SUMMARY_TABLE = (
    TABLE_DIR
    / "osmotic_strategy_final_summary.tsv"
)

MODULE_TABLE = (
    TABLE_DIR
    / "osmotic_module_evidence_373_MAGs.tsv"
)

ITOL_MAIN = (
    ITOL_DIR
    / "itol_osmotic_strategy_score_352.txt"
)

ITOL_STRICT = (
    ITOL_DIR
    / "itol_osmotic_strategy_score_strict_trehalose_352.txt"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def parse_ko_field(value):
    """
    Return all KEGG K numbers contained in a KO field.
    Works whether the field contains one KO or several.
    """

    if not value:
        return set()

    return set(
        re.findall(r"K\d{5}", value)
    )


def module_all(kos_present, required):
    """
    Multi-gene module:
    all listed genes are required for a complete module.

    Returns:
        detected marker count
        total required markers
        fractional module evidence
        binary complete/not-complete call
    """

    required = set(required)

    detected = len(
        kos_present & required
    )

    total = len(required)

    evidence = (
        detected / total
        if total > 0
        else 0
    )

    complete = int(
        required <= kos_present
    )

    return detected, total, evidence, complete


def module_any(kos_present, alternatives):
    """
    Alternative diagnostic markers:
    any listed KO is sufficient.

    Returns:
        detected marker count
        total possible markers
        binary any-marker evidence
        binary capability evidence
    """

    alternatives = set(alternatives)

    detected = len(
        kos_present & alternatives
    )

    total = len(alternatives)

    evidence = int(
        detected >= 1
    )

    complete = evidence

    return detected, total, evidence, complete


def singleton(kos_present, ko):
    """
    Single diagnostic marker.
    """

    present = int(
        ko in kos_present
    )

    return present, 1, float(present), present


# ============================================================
# READ MASTER TABLE
# ============================================================

mag_kos = defaultdict(set)
mag_metadata = {}

with INPUT.open("r", newline="") as infile:

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
            + ", ".join(sorted(missing))
        )

    for row in reader:

        mag = row["Final_MAG"].strip()

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
            mag_kos[mag].add(ko)


mags = sorted(
    mag_metadata
)

if len(mags) != 373:

    print(
        "WARNING: expected 373 MAGs, "
        f"but found {len(mags)}."
    )


# ============================================================
# BIOLOGICAL MODULE DEFINITIONS
# ============================================================

#
# ------------------------------------------------------------
# 1. TREHALOSE BIOSYNTHESIS
# ------------------------------------------------------------
#
# Recognised alternative trehalose-producing routes.
#
# Broad definition:
#     OtsAB OR TreYZ OR TreS OR TreT OR K16055
#
# Strict sensitivity definition:
#     OtsAB OR TreYZ OR K16055
#
# TreS and TreT are retained in the broad biosynthetic
# capability but excluded from the stricter sensitivity
# definition because their ecological role may be more
# context-dependent.
#

TREHALOSE_MODULES = {

    "OtsAB": {
        "type": "all",
        "kos": [
            "K00697",
            "K01087",
        ],
    },

    "TreYZ": {
        "type": "all",
        "kos": [
            "K06044",
            "K01236",
        ],
    },

    "TreS": {
        "type": "single",
        "kos": [
            "K05343",
        ],
    },

    "TreT": {
        "type": "single",
        "kos": [
            "K13057",
        ],
    },

    "Bifunctional_K16055": {
        "type": "single",
        "kos": [
            "K16055",
        ],
    },
}


#
# ------------------------------------------------------------
# 2. COMPATIBLE-SOLUTE UPTAKE
# ------------------------------------------------------------
#
# Alternative complete transporter systems.
#
# A MAG does not need both systems.
#

COMPATIBLE_SOLUTE_MODULES = {

    "Osmoprotectant_transporter": {
        "type": "all",
        "kos": [
            "K05845",
            "K05846",
            "K05847",
        ],
    },

    "ProVWX_transporter": {
        "type": "all",
        "kos": [
            "K02000",
            "K02001",
            "K02002",
        ],
    },
}


#
# ------------------------------------------------------------
# 3. ECTOINE BIOSYNTHESIS
# ------------------------------------------------------------
#
# EctABC represents the ectoine-specific biosynthetic core.
#

ECTOINE_MODULES = {

    "EctABC": {
        "type": "all",
        "kos": [
            "K06718",
            "K00836",
            "K06720",
        ],
    },
}


#
# ------------------------------------------------------------
# 4. MECHANOSENSITIVE OSMOTIC PROTECTION
# ------------------------------------------------------------
#
# These markers represent recognised mechanosensitive
# channel families.
#
# Any recognised channel is sufficient to provide evidence
# for mechanosensitive osmotic protection.
#

MECHANOSENSOR_MARKERS = [
    "K03282",
    "K03442",
    "K05802",
    "K16053",
    "K22044",
]


#
# ------------------------------------------------------------
# 5. ION HOMEOSTASIS
# ------------------------------------------------------------
#
# Several alternative Na+/K+ homeostasis mechanisms are
# retained separately.
#
# Complete multicomponent systems require all specified KOs.
# Single diagnostic transporters/antiporters are treated
# independently.
#

ION_MODULES = {

    "Ktr_Trk_system": {
        "type": "all",
        "kos": [
            "K03498",
            "K03499",
        ],
    },

    "NatAB_system": {
        "type": "all",
        "kos": [
            "K09696",
            "K09697",
        ],
    },

    "Kdp_system": {
        "type": "all",
        "kos": [
            "K01545",
            "K01546",
            "K01547",
            "K01548",
        ],
    },

    "Ion_marker_K15987": {
        "type": "single",
        "kos": [
            "K15987",
        ],
    },

    "Ion_marker_K24160": {
        "type": "single",
        "kos": [
            "K24160",
        ],
    },

    "Ion_marker_K03549": {
        "type": "single",
        "kos": [
            "K03549",
        ],
    },

    "Ion_marker_K11105": {
        "type": "single",
        "kos": [
            "K11105",
        ],
    },
}


# ============================================================
# GENERIC MODULE EVALUATION
# ============================================================

def evaluate_module(
    kos_present,
    definition
):

    module_type = definition["type"]
    kos = definition["kos"]

    if module_type == "all":

        return module_all(
            kos_present,
            kos
        )

    elif module_type == "single":

        return singleton(
            kos_present,
            kos[0]
        )

    elif module_type == "any":

        return module_any(
            kos_present,
            kos
        )

    else:

        raise ValueError(
            f"Unknown module type: "
            f"{module_type}"
        )


# ============================================================
# SCORE EACH MAG
# ============================================================

results = []

for mag in mags:

    kos_present = mag_kos[mag]

    row = {
        "Final_MAG":
            mag,
        **mag_metadata[mag],
    }


    # ========================================================
    # 1. TREHALOSE
    # ========================================================

    trehalose_positive_routes = []

    for name, definition in (
        TREHALOSE_MODULES.items()
    ):

        detected, total, evidence, complete = (
            evaluate_module(
                kos_present,
                definition
            )
        )

        row[
            f"Trehalose_{name}_markers_detected"
        ] = detected

        row[
            f"Trehalose_{name}_markers_total"
        ] = total

        row[
            f"Trehalose_{name}_module_evidence"
        ] = round(evidence, 3)

        row[
            f"Trehalose_{name}_complete"
        ] = complete

        if complete:
            trehalose_positive_routes.append(
                name
            )


    row[
        "Trehalose_routes_detected"
    ] = len(
        trehalose_positive_routes
    )

    row[
        "Trehalose_route_names"
    ] = (
        ";".join(
            trehalose_positive_routes
        )
        if trehalose_positive_routes
        else "None"
    )


    # Broad biological capability
    row[
        "Capability_Trehalose_broad"
    ] = int(
        len(
            trehalose_positive_routes
        ) >= 1
    )


    # Strict sensitivity capability
    strict_routes = {
        "OtsAB",
        "TreYZ",
        "Bifunctional_K16055",
    }

    row[
        "Capability_Trehalose_strict"
    ] = int(
        any(
            route in strict_routes
            for route
            in trehalose_positive_routes
        )
    )


    # Primary score currently uses broad trehalose
    row[
        "Capability_Trehalose"
    ] = row[
        "Capability_Trehalose_broad"
    ]


    # ========================================================
    # 2. COMPATIBLE-SOLUTE UPTAKE
    # ========================================================

    solute_positive_modules = []

    for name, definition in (
        COMPATIBLE_SOLUTE_MODULES.items()
    ):

        detected, total, evidence, complete = (
            evaluate_module(
                kos_present,
                definition
            )
        )

        row[
            f"Solute_{name}_markers_detected"
        ] = detected

        row[
            f"Solute_{name}_markers_total"
        ] = total

        row[
            f"Solute_{name}_module_evidence"
        ] = round(evidence, 3)

        row[
            f"Solute_{name}_complete"
        ] = complete

        if complete:
            solute_positive_modules.append(
                name
            )


    row[
        "Compatible_solute_systems_detected"
    ] = len(
        solute_positive_modules
    )

    row[
        "Compatible_solute_system_names"
    ] = (
        ";".join(
            solute_positive_modules
        )
        if solute_positive_modules
        else "None"
    )

    row[
        "Capability_Compatible_solute_uptake"
    ] = int(
        len(
            solute_positive_modules
        ) >= 1
    )


    # ========================================================
    # 3. ECTOINE BIOSYNTHESIS
    # ========================================================

    ectoine_positive = []

    for name, definition in (
        ECTOINE_MODULES.items()
    ):

        detected, total, evidence, complete = (
            evaluate_module(
                kos_present,
                definition
            )
        )

        row[
            f"Ectoine_{name}_markers_detected"
        ] = detected

        row[
            f"Ectoine_{name}_markers_total"
        ] = total

        row[
            f"Ectoine_{name}_module_evidence"
        ] = round(evidence, 3)

        row[
            f"Ectoine_{name}_complete"
        ] = complete

        if complete:
            ectoine_positive.append(
                name
            )


    row[
        "Capability_Ectoine_biosynthesis"
    ] = int(
        len(
            ectoine_positive
        ) >= 1
    )


    # ========================================================
    # 4. MECHANOSENSITIVE PROTECTION
    # ========================================================

    mechanosensor_detected = [
        ko
        for ko in MECHANOSENSOR_MARKERS
        if ko in kos_present
    ]

    row[
        "Mechanosensor_markers_detected"
    ] = len(
        mechanosensor_detected
    )

    row[
        "Mechanosensor_marker_names"
    ] = (
        ";".join(
            mechanosensor_detected
        )
        if mechanosensor_detected
        else "None"
    )

    row[
        "Mechanosensor_module_evidence"
    ] = round(
        len(
            mechanosensor_detected
        )
        /
        len(
            MECHANOSENSOR_MARKERS
        ),
        3,
    )

    row[
        "Capability_Mechanosensitive_protection"
    ] = int(
        len(
            mechanosensor_detected
        ) >= 1
    )


    # ========================================================
    # 5. ION HOMEOSTASIS
    # ========================================================

    ion_positive_modules = []

    for name, definition in (
        ION_MODULES.items()
    ):

        detected, total, evidence, complete = (
            evaluate_module(
                kos_present,
                definition
            )
        )

        row[
            f"Ion_{name}_markers_detected"
        ] = detected

        row[
            f"Ion_{name}_markers_total"
        ] = total

        row[
            f"Ion_{name}_module_evidence"
        ] = round(
            evidence,
            3
        )

        row[
            f"Ion_{name}_complete"
        ] = complete

        if complete:
            ion_positive_modules.append(
                name
            )


    row[
        "Ion_systems_detected"
    ] = len(
        ion_positive_modules
    )

    row[
        "Ion_system_names"
    ] = (
        ";".join(
            ion_positive_modules
        )
        if ion_positive_modules
        else "None"
    )

    row[
        "Capability_Ion_homeostasis"
    ] = int(
        len(
            ion_positive_modules
        ) >= 1
    )


    # ========================================================
    # FINAL OSMOTIC STRATEGY
    # ========================================================

    capability_columns = [
        "Capability_Trehalose",
        "Capability_Compatible_solute_uptake",
        "Capability_Ectoine_biosynthesis",
        "Capability_Mechanosensitive_protection",
        "Capability_Ion_homeostasis",
    ]

    capability_count = sum(
        row[col]
        for col
        in capability_columns
    )

    row[
        "Osmotic_capabilities_detected"
    ] = capability_count

    row[
        "Osmotic_strategy_score"
    ] = round(
        capability_count / 5,
        3
    )


    # ========================================================
    # STRICT TREHALOSE SENSITIVITY SCORE
    # ========================================================

    strict_capability_count = (
        row[
            "Capability_Trehalose_strict"
        ]
        +
        row[
            "Capability_Compatible_solute_uptake"
        ]
        +
        row[
            "Capability_Ectoine_biosynthesis"
        ]
        +
        row[
            "Capability_Mechanosensitive_protection"
        ]
        +
        row[
            "Capability_Ion_homeostasis"
        ]
    )

    row[
        "Osmotic_capabilities_detected_strict_trehalose"
    ] = strict_capability_count

    row[
        "Osmotic_strategy_score_strict_trehalose"
    ] = round(
        strict_capability_count / 5,
        3
    )

    results.append(row)


# ============================================================
# WRITE COMPLETE FINAL TABLE
# ============================================================

if not results:
    raise ValueError(
        "No MAG results were generated."
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
    writer.writerows(results)


# ============================================================
# WRITE MODULE-EVIDENCE TABLE
# ============================================================

module_columns = [
    col
    for col in fieldnames
    if (
        col == "Final_MAG"
        or
        col in {
            "Completeness",
            "Contamination",
            "Site",
            "Treatment",
            "Compartment",
            "GTDB_phylum",
            "GTDB_order",
        }
        or
        "markers_detected" in col
        or
        "markers_total" in col
        or
        "module_evidence" in col
        or
        col.endswith("_complete")
        or
        col.endswith("_route_names")
        or
        col.endswith("_system_names")
        or
        col == "Mechanosensor_marker_names"
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
# SUMMARY STATISTICS
# ============================================================

summary_rows = []

capability_summary = [
    (
        "Trehalose_broad",
        "Capability_Trehalose_broad"
    ),
    (
        "Trehalose_strict",
        "Capability_Trehalose_strict"
    ),
    (
        "Compatible_solute_uptake",
        "Capability_Compatible_solute_uptake"
    ),
    (
        "Ectoine_biosynthesis",
        "Capability_Ectoine_biosynthesis"
    ),
    (
        "Mechanosensitive_protection",
        "Capability_Mechanosensitive_protection"
    ),
    (
        "Ion_homeostasis",
        "Capability_Ion_homeostasis"
    ),
]

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


# Main strategy-score distribution

score_distribution = Counter(
    row[
        "Osmotic_capabilities_detected"
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
            f"{count_value}_of_5",
        "MAGs":
            n,
        "Percent":
            f"{100 * n / len(results):.2f}",
    })


# Strict-trehalose sensitivity distribution

strict_distribution = Counter(
    row[
        "Osmotic_capabilities_detected_strict_trehalose"
    ]
    for row in results
)

for count_value in sorted(
    strict_distribution
):

    n = strict_distribution[
        count_value
    ]

    summary_rows.append({
        "Summary_type":
            "Strict_trehalose_sensitivity_distribution",
        "Category":
            f"{count_value}_of_5",
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
# EXTRACT BACTERIAL TREE TIP NAMES
# ============================================================

tree_tips = set()

if TREE.exists():

    newick = TREE.read_text().strip()

    # Typical Newick tip labels occurring after "(" or ","
    # and before ":".
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
            tree_tips.add(item)

else:

    print(
        f"WARNING: tree not found: {TREE}"
    )


tree_results = [
    row
    for row in results
    if row["Final_MAG"] in tree_tips
]


# ============================================================
# WRITE iTOL DATASET — MAIN SCORE
# ============================================================

with ITOL_MAIN.open("w") as outfile:

    outfile.write(
        "DATASET_GRADIENT\n"
    )

    outfile.write(
        "SEPARATOR TAB\n"
    )

    outfile.write(
        "DATASET_LABEL\t"
        "Osmotic adaptation\n"
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
        "Osmotic strategy breadth\n"
    )

    outfile.write(
        "DATA\n"
    )

    for row in tree_results:

        outfile.write(
            f"{row['Final_MAG']}\t"
            f"{row['Osmotic_strategy_score']}\n"
        )


# ============================================================
# WRITE iTOL DATASET — STRICT TREHALOSE SENSITIVITY
# ============================================================

with ITOL_STRICT.open("w") as outfile:

    outfile.write(
        "DATASET_GRADIENT\n"
    )

    outfile.write(
        "SEPARATOR TAB\n"
    )

    outfile.write(
        "DATASET_LABEL\t"
        "Osmotic adaptation - strict trehalose\n"
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
        "Osmotic strategy breadth\n"
    )

    outfile.write(
        "DATA\n"
    )

    for row in tree_results:

        outfile.write(
            f"{row['Final_MAG']}\t"
            f"{row['Osmotic_strategy_score_strict_trehalose']}\n"
        )


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print("=" * 72)
print("FINAL OSMOTIC-ADAPTATION STRATEGY ANALYSIS")
print("=" * 72)

print()
print(
    f"MAGs assessed: {len(results)}"
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
print("Capability prevalence:")
print()

for label, column in capability_summary:

    count = sum(
        int(row[column])
        for row in results
    )

    percent = (
        100 * count / len(results)
    )

    print(
        f"  {label:<35}"
        f"{count:>4} MAGs "
        f"({percent:6.2f}%)"
    )


print()
print(
    "Main osmotic strategy distribution "
    "(broad trehalose definition):"
)

for capability_count in sorted(
    score_distribution
):

    n = score_distribution[
        capability_count
    ]

    print(
        f"  {capability_count}/5 : "
        f"{n:>4} MAGs "
        f"({100*n/len(results):6.2f}%)"
    )


print()
print(
    "Sensitivity distribution "
    "(strict trehalose definition):"
)

for capability_count in sorted(
    strict_distribution
):

    n = strict_distribution[
        capability_count
    ]

    print(
        f"  {capability_count}/5 : "
        f"{n:>4} MAGs "
        f"({100*n/len(results):6.2f}%)"
    )


print()
print("Outputs:")
print(
    f"  Final MAG table: {FINAL_TABLE}"
)

print(
    f"  Module evidence: {MODULE_TABLE}"
)

print(
    f"  Summary:         {SUMMARY_TABLE}"
)

print(
    f"  iTOL main:       {ITOL_MAIN}"
)

print(
    f"  iTOL sensitivity:{ITOL_STRICT}"
)

print()
print(
    "Interpretation:"
)

print(
    "  Osmotic_strategy_score = "
    "number of detected biological capabilities / 5"
)

print(
    "  The score measures standardized "
    "osmotic capability breadth, not physiological performance."
)

print(
    "  Genome completeness is retained in the table "
    "but is not used to inflate individual MAG scores."
)

print("=" * 72)

