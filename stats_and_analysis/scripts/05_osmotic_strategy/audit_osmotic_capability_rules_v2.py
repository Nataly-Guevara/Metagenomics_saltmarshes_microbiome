#!/usr/bin/env python3

"""
Second audit of osmotic-adaptation capability rules.

This script:

1. Reads the final 373-MAG gene master table.
2. Collects all KOs detected in each MAG.
3. Tests corrected and alternative osmotic-capability rules.
4. Separates sodium and potassium homeostasis mechanisms.
5. Compares strict and MAG-tolerant Kdp rules.
6. Checks complete trehalose/maltose transporter architecture.
7. Reports individual mechanosensitive-channel classes.
8. Writes:
   - a rule-level audit table;
   - a MAG-by-rule presence matrix;
   - individual candidate-KO prevalence;
   - multigene system component diagnostics.

The script uses only the Python standard library.
"""

import csv
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

OUTPUT_DIR = os.path.join(
    PROJECT_ROOT,
    "results/ecological_strategies/osmotic_adaptation/tables",
)

AUDIT_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_capability_rule_audit_v2.tsv",
)

PRESENCE_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_capability_presence_v2.tsv",
)

KO_PREVALENCE_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_candidate_KO_prevalence_v2.tsv",
)

SYSTEM_DIAGNOSTIC_OUTPUT = os.path.join(
    OUTPUT_DIR,
    "osmotic_multigene_system_diagnostics_v2.tsv",
)


# ---------------------------------------------------------------------
# Rule definitions
#
# Each rule contains one or more alternative routes.
# A MAG is positive when it contains every KO in at least one route.
#
# Example:
# alternatives = [
#     ("K00697", "K01087"),
#     ("K16055",),
# ]
#
# means:
# (K00697 AND K01087) OR K16055
# ---------------------------------------------------------------------

RULES = [
    {
        "capability": "Trehalose biosynthesis",
        "rule_name": "Trehalose_OtsAB",
        "alternatives": [
            ("K00697", "K01087"),
        ],
        "status": "Retain",
        "rationale": (
            "Complete OtsAB route requiring trehalose-6-phosphate "
            "synthase and trehalose-6-phosphate phosphatase."
        ),
    },
    {
        "capability": "Trehalose biosynthesis",
        "rule_name": "Trehalose_bifunctional_K16055",
        "alternatives": [
            ("K16055",),
        ],
        "status": "Retain_provisionally",
        "rationale": (
            "Single bifunctional trehalose-6-phosphate "
            "synthase/phosphatase marker; retain subject to "
            "annotation and taxonomy review."
        ),
    },
    {
        "capability": "Trehalose biosynthesis",
        "rule_name": "Trehalose_biosynthesis_any",
        "alternatives": [
            ("K00697", "K01087"),
            ("K16055",),
        ],
        "status": "Provisional_final_rule",
        "rationale": (
            "Positive when either the complete OtsAB route or the "
            "bifunctional K16055 route is detected."
        ),
    },
    {
        "capability": "Trehalose uptake",
        "rule_name": "Trehalose_transporter_core_K10236_K10237_K10238",
        "alternatives": [
            ("K10236", "K10237", "K10238"),
        ],
        "status": "Diagnostic_only",
        "rationale": (
            "Binding and permease components are present, but the "
            "ATP-binding component is not required by this diagnostic rule."
        ),
    },
    {
        "capability": "Trehalose uptake",
        "rule_name": "Trehalose_transporter_complete_with_ATPase",
        "alternatives": [
            ("K10236", "K10237", "K10238", "K10111"),
            ("K10236", "K10237", "K10238", "K10112"),
        ],
        "status": "Candidate_complete_system",
        "rationale": (
            "Requires binding protein, two permease components, and "
            "either recognised multiple-sugar transporter ATPase."
        ),
    },
    {
        "capability": "Compatible-solute uptake",
        "rule_name": "Osmoprotectant_transporter_K05845_K05846_K05847",
        "alternatives": [
            ("K05845", "K05846", "K05847"),
        ],
        "status": "Retain",
        "rationale": (
            "Requires substrate-binding, permease, and ATP-binding "
            "components of the osmoprotectant ABC transporter."
        ),
    },
    {
        "capability": "Compatible-solute uptake",
        "rule_name": "ProVWX_transporter_K02000_K02001_K02002",
        "alternatives": [
            ("K02000", "K02001", "K02002"),
        ],
        "status": "Retain",
        "rationale": (
            "Requires complete ATP-binding, permease, and "
            "substrate-binding components of the ProVWX system."
        ),
    },
    {
        "capability": "Compatible-solute uptake",
        "rule_name": "Compatible_solute_uptake_any",
        "alternatives": [
            ("K05845", "K05846", "K05847"),
            ("K02000", "K02001", "K02002"),
        ],
        "status": "Provisional_final_rule",
        "rationale": (
            "Positive when either complete osmoprotectant uptake "
            "system is detected."
        ),
    },
    {
        "capability": "Ectoine biosynthesis",
        "rule_name": "Ectoine_EctABC",
        "alternatives": [
            ("K06718", "K00836", "K06720"),
        ],
        "status": "Retain",
        "rationale": (
            "Requires all three dedicated ectoine-biosynthesis "
            "enzymes EctA, EctB, and EctC."
        ),
    },

    # Individual mechanosensitive-channel classes
    {
        "capability": "Mechanosensitive channels",
        "rule_name": "Channel_MscL_K03282",
        "alternatives": [
            ("K03282",),
        ],
        "status": "Retain_channel_subtype",
        "rationale": "Large-conductance mechanosensitive channel.",
    },
    {
        "capability": "Mechanosensitive channels",
        "rule_name": "Channel_MscS_K03442",
        "alternatives": [
            ("K03442",),
        ],
        "status": "Retain_channel_subtype",
        "rationale": "Small-conductance mechanosensitive channel.",
    },
    {
        "capability": "Mechanosensitive channels",
        "rule_name": "Channel_MscK_K05802",
        "alternatives": [
            ("K05802",),
        ],
        "status": "Retain_channel_subtype",
        "rationale": (
            "Potassium-dependent mechanosensitive channel."
        ),
    },
    {
        "capability": "Mechanosensitive channels",
        "rule_name": "Channel_MscM_K16053",
        "alternatives": [
            ("K16053",),
        ],
        "status": "Retain_channel_subtype",
        "rationale": "Miniconductance mechanosensitive channel.",
    },
    {
        "capability": "Mechanosensitive channels",
        "rule_name": "Channel_YbiO_K22044",
        "alternatives": [
            ("K22044",),
        ],
        "status": "Retain_channel_subtype",
        "rationale": "Moderate-conductance mechanosensitive channel.",
    },
    {
        "capability": "Mechanosensitive channels",
        "rule_name": "Mechanosensitive_channel_any",
        "alternatives": [
            ("K03282",),
            ("K03442",),
            ("K05802",),
            ("K16053",),
            ("K22044",),
        ],
        "status": "Retain_but_report_saturation",
        "rationale": (
            "Positive when at least one recognised mechanosensitive "
            "channel class is detected."
        ),
    },

    # Potassium systems
    {
        "capability": "Potassium homeostasis",
        "rule_name": "Potassium_TrkKtr_K03498_K03499",
        "alternatives": [
            ("K03498", "K03499"),
        ],
        "status": "Retain_provisionally",
        "rationale": (
            "Requires both membrane and regulatory components of "
            "the Trk/Ktr potassium-uptake system."
        ),
    },
    {
        "capability": "Potassium homeostasis",
        "rule_name": "Potassium_KdpFABC_strict",
        "alternatives": [
            ("K01545", "K01546", "K01547", "K01548"),
        ],
        "status": "Strict_sensitivity_rule",
        "rationale": (
            "Strict rule requiring all four annotated KdpFABC "
            "components."
        ),
    },
    {
        "capability": "Potassium homeostasis",
        "rule_name": "Potassium_KdpABC_without_KdpF",
        "alternatives": [
            ("K01546", "K01547", "K01548"),
        ],
        "status": "Intermediate_sensitivity_rule",
        "rationale": (
            "Requires KdpA, KdpB, and KdpC but does not require "
            "the very small KdpF accessory subunit."
        ),
    },
    {
        "capability": "Potassium homeostasis",
        "rule_name": "Potassium_KdpAB_functional_core",
        "alternatives": [
            ("K01546", "K01547"),
        ],
        "status": "MAG_tolerant_candidate",
        "rationale": (
            "Requires potassium-binding KdpA and ATP-binding KdpB, "
            "the central transport machinery."
        ),
    },
    {
        "capability": "Potassium homeostasis",
        "rule_name": "Potassium_KUP_K03549",
        "alternatives": [
            ("K03549",),
        ],
        "status": "Candidate_single_component_system",
        "rationale": (
            "KUP-family potassium-uptake protein; evaluate whether "
            "it is sufficiently diagnostic for the final capability."
        ),
    },
    {
        "capability": "Potassium homeostasis",
        "rule_name": "Potassium_H_antiporter_K11105",
        "alternatives": [
            ("K11105",),
        ],
        "status": "Candidate_single_component_system",
        "rationale": (
            "Potassium/hydrogen antiporter; analyse separately from "
            "sodium antiport."
        ),
    },

    # Sodium systems
    {
        "capability": "Sodium homeostasis",
        "rule_name": "Sodium_NatAB_K09696_K09697",
        "alternatives": [
            ("K09696", "K09697"),
        ],
        "status": "Retain",
        "rationale": (
            "Requires permease and ATP-binding components of the "
            "NatAB sodium transport system."
        ),
    },
    {
        "capability": "Sodium homeostasis",
        "rule_name": "Sodium_H_antiporter_K24160",
        "alternatives": [
            ("K24160",),
        ],
        "status": "Candidate_single_component_system",
        "rationale": (
            "Recognised sodium/hydrogen antiporter; retain separately "
            "from potassium antiporters."
        ),
    },
    {
        "capability": "Sodium homeostasis",
        "rule_name": "Sodium_pumping_pyrophosphatase_K15987",
        "alternatives": [
            ("K15987",),
        ],
        "status": "Candidate_single_component_system",
        "rationale": (
            "Sodium-translocating pyrophosphatase candidate; evaluate "
            "as an independent sodium-homeostasis mechanism."
        ),
    },

    # Composite diagnostic rules
    {
        "capability": "Sodium homeostasis",
        "rule_name": "Sodium_homeostasis_candidate_any",
        "alternatives": [
            ("K09696", "K09697"),
            ("K24160",),
            ("K15987",),
        ],
        "status": "Diagnostic_not_final",
        "rationale": (
            "Broad sodium-only diagnostic combining complete NatAB, "
            "sodium/hydrogen antiport, or sodium-pumping pyrophosphatase."
        ),
    },
    {
        "capability": "Potassium homeostasis",
        "rule_name": "Potassium_homeostasis_candidate_any",
        "alternatives": [
            ("K03498", "K03499"),
            ("K01546", "K01547"),
            ("K03549",),
            ("K11105",),
        ],
        "status": "Diagnostic_not_final",
        "rationale": (
            "Broad potassium-only diagnostic combining Trk/Ktr, "
            "KdpAB core, KUP, or potassium/hydrogen antiport."
        ),
    },
    {
        "capability": "Ion homeostasis",
        "rule_name": "Ion_homeostasis_candidate_any",
        "alternatives": [
            ("K09696", "K09697"),
            ("K24160",),
            ("K15987",),
            ("K03498", "K03499"),
            ("K01546", "K01547"),
            ("K03549",),
            ("K11105",),
        ],
        "status": "Diagnostic_not_final",
        "rationale": (
            "Broad union of candidate sodium and potassium systems. "
            "This must not yet be used in the final five-capability score."
        ),
    },
]


# ---------------------------------------------------------------------
# Multigene systems for component-recovery diagnostics
# ---------------------------------------------------------------------

MULTIGENE_SYSTEMS = [
    {
        "system": "Trehalose_OtsAB",
        "components": ("K00697", "K01087"),
    },
    {
        "system": "Trehalose_transporter_core",
        "components": ("K10236", "K10237", "K10238"),
    },
    {
        "system": "Osmoprotectant_K05845_K05846_K05847",
        "components": ("K05845", "K05846", "K05847"),
    },
    {
        "system": "ProVWX_K02000_K02001_K02002",
        "components": ("K02000", "K02001", "K02002"),
    },
    {
        "system": "EctABC",
        "components": ("K06718", "K00836", "K06720"),
    },
    {
        "system": "TrkKtr_K03498_K03499",
        "components": ("K03498", "K03499"),
    },
    {
        "system": "NatAB_K09696_K09697",
        "components": ("K09696", "K09697"),
    },
    {
        "system": "KdpFABC",
        "components": ("K01545", "K01546", "K01547", "K01548"),
    },
]


KO_PATTERN = re.compile(r"K\d{5}")


def rule_is_positive(ko_set, alternatives):
    """Return True if every KO in at least one alternative route is present."""
    for route in alternatives:
        if all(ko in ko_set for ko in route):
            return True
    return False


def format_rule_logic(alternatives):
    """Return a readable representation of alternative KO rules."""
    route_strings = []

    for route in alternatives:
        route_strings.append(" + ".join(route))

    return " OR ".join(route_strings)


def safe_percent(numerator, denominator):
    if denominator == 0:
        return 0.0

    return 100.0 * numerator / denominator


def find_metadata_columns(fieldnames):
    """
    Select MAG-level metadata columns when present.

    The script preserves only known MAG-level columns and avoids
    gene-specific columns that vary among rows.
    """
    preferred_columns = [
        "Completeness",
        "Contamination",
        "Quality_class",
        "Site",
        "Treatment",
        "Compartment",
        "GTDB_taxonomy",
        "GTDB taxonomy",
        "GTDB_Taxonomy",
        "Classification",
        "classification",
        "Taxonomy",
        "taxonomy",
    ]

    return [
        column
        for column in preferred_columns
        if column in fieldnames
    ]


def read_master_table(input_file):
    """
    Read the gene master table.

    Returns:
        mag_kos:
            dictionary mapping each MAG to a set of KOs.

        mag_metadata:
            dictionary mapping each MAG to its MAG-level metadata.

        metadata_columns:
            ordered list of metadata fields preserved in output.
    """
    mag_kos = defaultdict(set)
    mag_metadata = {}

    try:
        handle = open(
            input_file,
            "r",
            encoding="utf-8-sig",
            newline="",
        )
    except OSError as exc:
        sys.exit(
            "ERROR: Could not open input file:\n"
            f"{input_file}\n"
            f"{exc}"
        )

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            sys.exit("ERROR: Input table has no header.")

        required_columns = {"Final_MAG", "KO"}
        missing_columns = required_columns.difference(reader.fieldnames)

        if missing_columns:
            sys.exit(
                "ERROR: Input table is missing required column(s): "
                + ", ".join(sorted(missing_columns))
            )

        metadata_columns = find_metadata_columns(reader.fieldnames)

        for row in reader:
            mag = (row.get("Final_MAG") or "").strip()

            if not mag:
                continue

            ko_text = row.get("KO") or ""

            for ko in KO_PATTERN.findall(ko_text):
                mag_kos[mag].add(ko)

            if mag not in mag_metadata:
                mag_metadata[mag] = {
                    column: (row.get(column) or "").strip()
                    for column in metadata_columns
                }

    return mag_kos, mag_metadata, metadata_columns


def collect_candidate_kos():
    """Collect every KO referenced by any rule or diagnostic system."""
    candidate_kos = set()

    for rule in RULES:
        for route in rule["alternatives"]:
            candidate_kos.update(route)

    for system in MULTIGENE_SYSTEMS:
        candidate_kos.update(system["components"])

    return sorted(candidate_kos)


def calculate_rule_calls(mag_kos):
    """
    Calculate binary rule calls for every MAG.

    Also calculates the number of mechanosensitive-channel classes.
    """
    calls = {}

    channel_kos = (
        "K03282",
        "K03442",
        "K05802",
        "K16053",
        "K22044",
    )

    for mag, ko_set in mag_kos.items():
        calls[mag] = {}

        for rule in RULES:
            positive = rule_is_positive(
                ko_set,
                rule["alternatives"],
            )

            calls[mag][rule["rule_name"]] = int(positive)

        calls[mag]["Mechanosensitive_channel_class_count"] = sum(
            1 for ko in channel_kos if ko in ko_set
        )

    return calls


def write_rule_audit(mag_kos, calls):
    total_mags = len(mag_kos)

    columns = [
        "Capability",
        "Rule_name",
        "Rule_logic",
        "Positive_MAGs",
        "Total_MAGs",
        "Prevalence_percent",
        "Decision_status",
        "Biological_rationale",
    ]

    with open(
        AUDIT_OUTPUT,
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

        for rule in RULES:
            rule_name = rule["rule_name"]

            positive_count = sum(
                calls[mag][rule_name]
                for mag in mag_kos
            )

            writer.writerow(
                {
                    "Capability": rule["capability"],
                    "Rule_name": rule_name,
                    "Rule_logic": format_rule_logic(
                        rule["alternatives"]
                    ),
                    "Positive_MAGs": positive_count,
                    "Total_MAGs": total_mags,
                    "Prevalence_percent": (
                        f"{safe_percent(positive_count, total_mags):.2f}"
                    ),
                    "Decision_status": rule["status"],
                    "Biological_rationale": rule["rationale"],
                }
            )


def write_presence_matrix(
    mag_kos,
    mag_metadata,
    metadata_columns,
    calls,
):
    rule_columns = [
        rule["rule_name"]
        for rule in RULES
    ]

    output_columns = (
        ["Final_MAG"]
        + metadata_columns
        + rule_columns
        + ["Mechanosensitive_channel_class_count"]
    )

    with open(
        PRESENCE_OUTPUT,
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
            output_row = {"Final_MAG": mag}

            output_row.update(
                mag_metadata.get(mag, {})
            )

            output_row.update(
                calls[mag]
            )

            writer.writerow(output_row)


def write_ko_prevalence(mag_kos, candidate_kos):
    total_mags = len(mag_kos)

    columns = [
        "KO",
        "MAGs_detected",
        "MAGs_total",
        "Prevalence_percent",
    ]

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

        prevalence_rows = []

        for ko in candidate_kos:
            detected = sum(
                1
                for ko_set in mag_kos.values()
                if ko in ko_set
            )

            prevalence_rows.append(
                {
                    "KO": ko,
                    "MAGs_detected": detected,
                    "MAGs_total": total_mags,
                    "Prevalence_percent": safe_percent(
                        detected,
                        total_mags,
                    ),
                }
            )

        prevalence_rows.sort(
            key=lambda row: (
                -row["MAGs_detected"],
                row["KO"],
            )
        )

        for row in prevalence_rows:
            writer.writerow(
                {
                    "KO": row["KO"],
                    "MAGs_detected": row["MAGs_detected"],
                    "MAGs_total": row["MAGs_total"],
                    "Prevalence_percent": (
                        f"{row['Prevalence_percent']:.2f}"
                    ),
                }
            )


def write_system_diagnostics(mag_kos):
    """
    For every multigene system, count how many MAGs contain:

    - zero components;
    - one component;
    - two components;
    - etc.;
    - all required components.

    This helps distinguish complete systems from isolated components.
    """
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
        SYSTEM_DIAGNOSTIC_OUTPUT,
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

        for system in MULTIGENE_SYSTEMS:
            system_name = system["system"]
            components = system["components"]
            number_required = len(components)

            count_distribution = defaultdict(int)

            for ko_set in mag_kos.values():
                components_detected = sum(
                    1
                    for ko in components
                    if ko in ko_set
                )

                count_distribution[components_detected] += 1

            for components_detected in range(
                number_required + 1
            ):
                mag_count = count_distribution[
                    components_detected
                ]

                if components_detected == 0:
                    recovery_class = "No_components"

                elif components_detected == number_required:
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
                        "Components_detected": components_detected,
                        "MAGs": mag_count,
                        "MAGs_total": total_mags,
                        "Prevalence_percent": (
                            f"{safe_percent(mag_count, total_mags):.2f}"
                        ),
                        "Recovery_class": recovery_class,
                    }
                )


def print_summary(mag_kos, calls, candidate_kos):
    total_mags = len(mag_kos)

    print("Second osmotic capability audit completed.")
    print()
    print(f"MAGs assessed: {total_mags}")
    print(f"Candidate KOs assessed: {len(candidate_kos)}")
    print()
    print(f"Rule audit: {AUDIT_OUTPUT}")
    print(f"Presence matrix: {PRESENCE_OUTPUT}")
    print(f"KO prevalence: {KO_PREVALENCE_OUTPUT}")
    print(
        "System diagnostics: "
        f"{SYSTEM_DIAGNOSTIC_OUTPUT}"
    )
    print()

    header_rule = "Capability rule"
    header_positive = "Positive MAGs"
    header_prevalence = "Prevalence %"

    print(
        f"{header_rule:<56}"
        f"{header_positive:>15}"
        f"{header_prevalence:>15}"
    )

    print("-" * 86)

    for rule in RULES:
        rule_name = rule["rule_name"]

        positive_count = sum(
            calls[mag][rule_name]
            for mag in mag_kos
        )

        prevalence = safe_percent(
            positive_count,
            total_mags,
        )

        print(
            f"{rule_name:<56}"
            f"{positive_count:>15}"
            f"{prevalence:>15.2f}"
        )

    channel_count_distribution = defaultdict(int)

    for mag in mag_kos:
        class_count = calls[mag][
            "Mechanosensitive_channel_class_count"
        ]

        channel_count_distribution[class_count] += 1

    print()
    print("Mechanosensitive-channel class diversity")
    print("-----------------------------------------")

    for class_count in sorted(channel_count_distribution):
        mag_count = channel_count_distribution[class_count]

        print(
            f"{class_count} channel classes: "
            f"{mag_count} MAGs "
            f"({safe_percent(mag_count, total_mags):.2f}%)"
        )


def main():
    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    mag_kos, mag_metadata, metadata_columns = (
        read_master_table(INPUT_FILE)
    )

    if not mag_kos:
        sys.exit(
            "ERROR: No MAGs were recovered from the input table."
        )

    candidate_kos = collect_candidate_kos()

    calls = calculate_rule_calls(mag_kos)

    write_rule_audit(
        mag_kos,
        calls,
    )

    write_presence_matrix(
        mag_kos,
        mag_metadata,
        metadata_columns,
        calls,
    )

    write_ko_prevalence(
        mag_kos,
        candidate_kos,
    )

    write_system_diagnostics(
        mag_kos,
    )

    print_summary(
        mag_kos,
        calls,
        candidate_kos,
    )


if __name__ == "__main__":
    main()
