#!/usr/bin/env python3

"""
Audit Dsr pathway directionality in MAGs currently called positive
for the Sat-AprAB-DsrAB sulfate-reduction core.

PURPOSE
-------
This script does NOT change the final redox score.

It performs a targeted genomic-context audit for the Dsr-positive MAGs
and asks whether their Dsr system is more consistent with:

    1. reductive sulfate/sulfite metabolism,
    2. oxidative/reverse-Dsr sulfur metabolism, or
    3. unresolved directionality.

The output is intended to determine which MAGs require DsrAB
phylogenetic analysis before final interpretation.

Evidence layers
---------------
A. Core genes already used in the redox pipeline:
   sat, aprA, aprB, dsrA, dsrB

B. Dsr-associated context:
   dsrC
   dsrD
   dsrE
   dsrF
   dsrH
   dsrL
   DsrMKJOP

C. Independent sulfur-oxidation context:
   Sox system
   SQR
   SoxCD

Important:
- DsrMKJOP is NOT directional by itself.
- DsrL is recorded but NOT used alone to determine direction.
- DsrC annotations can overlap with TusE-like proteins and are
  therefore retained as supporting evidence only.
- The automated interpretation is deliberately conservative.
- Final directionality can be resolved using DsrAB phylogeny.
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

MASTER = (
    BASE
    / "results/functional_annotation/kofamscan/"
      "final_373_MAGs_gene_master_table.tsv"
)

REDOX_TABLE = (
    BASE
    / "results/ecological_strategies/"
      "redox_flexibility/tables/"
      "redox_strategy_final_v2_373_MAGs.tsv"
)

OUTPUT_DIR = (
    BASE
    / "results/ecological_strategies/"
      "redox_flexibility/dsr_directionality"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_TABLE = (
    OUTPUT_DIR
    / "dsr_directionality_audit.tsv"
)

SUMMARY_TABLE = (
    OUTPUT_DIR
    / "dsr_directionality_summary.tsv"
)


# ============================================================
# KNOWN KO MARKERS
# ============================================================

# Core pathway used in the redox score
CORE_KOS = {
    "sat": "K00958",
    "aprA": "K00394",
    "aprB": "K00395",
    "dsrA": "K11180",
    "dsrB": "K11181",
}

# DsrC-related KOs.
# These can overlap with TusE/DsrC-like proteins, so they are
# supporting evidence only and are NOT used for directionality.
DSRC_KOS = {
    "K11179",
    "K23077",
}

# DsrL
DSRL_KOS = {
    "K27196",
}

# DsrMKJOP
DSRMKJOP_KOS = {
    "dsrM": "K27187",
    "dsrK": "K27188",
    "dsrJ": "K27189",
    "dsrO": "K27190",
    "dsrP": "K27191",
}

# Sox core already used in redox pipeline
SOX_CORE_KOS = {
    "soxA": "K17222",
    "soxX": "K17223",
    "soxB": "K17224",
    "soxY": "K17226",
    "soxZ": "K17227",
}

# Optional SoxCD oxidation components
SOXCD_KOS = {
    "soxC": "K17225",
    "soxD": "K22622",
}

# Sulfide:quinone oxidoreductase
SQR_KO = "K17218"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def extract_kos(value):
    if not value:
        return set()

    return set(
        re.findall(
            r"K\d{5}",
            value
        )
    )


def text_has_gene_symbol(text, symbol):
    """
    Conservative gene-symbol search across the complete annotation row.

    Examples:
      dsrD
      dsrE
      dsrF
      dsrH
      dsrL

    Word boundaries prevent matches inside unrelated longer strings.
    """

    pattern = (
        r"(?i)(?<![A-Za-z0-9])"
        + re.escape(symbol)
        + r"(?![A-Za-z0-9])"
    )

    return bool(
        re.search(
            pattern,
            text
        )
    )


# ============================================================
# IDENTIFY THE DSR-POSITIVE MAGs FROM FINAL REDOX TABLE
# ============================================================

candidate_mags = set()
redox_metadata = {}

with REDOX_TABLE.open(
    "r",
    newline=""
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t",
    )

    required = {
        "Final_MAG",
        "Capability_Dissimilatory_sulfate_reduction_core",
        "Completeness",
        "Contamination",
        "Site",
        "Treatment",
        "Compartment",
        "GTDB_phylum",
        "GTDB_order",
    }

    if reader.fieldnames is None:
        raise ValueError(
            "Redox table has no header."
        )

    missing = (
        required
        - set(reader.fieldnames)
    )

    if missing:
        raise ValueError(
            "Missing columns in redox table: "
            + ", ".join(
                sorted(missing)
            )
        )

    for row in reader:

        if (
            row[
                "Capability_Dissimilatory_sulfate_reduction_core"
            ].strip()
            == "1"
        ):

            mag = row[
                "Final_MAG"
            ].strip()

            candidate_mags.add(
                mag
            )

            redox_metadata[mag] = {
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


# ============================================================
# READ ALL ANNOTATION ROWS FOR THOSE MAGs
# ============================================================

mag_kos = defaultdict(set)
mag_annotation_text = defaultdict(list)

with MASTER.open(
    "r",
    newline=""
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t",
    )

    if reader.fieldnames is None:
        raise ValueError(
            "Master table has no header."
        )

    if "Final_MAG" not in reader.fieldnames:
        raise ValueError(
            "Master table does not contain Final_MAG."
        )

    if "KO" not in reader.fieldnames:
        raise ValueError(
            "Master table does not contain KO."
        )

    for row in reader:

        mag = row[
            "Final_MAG"
        ].strip()

        if mag not in candidate_mags:
            continue

        for ko in extract_kos(
            row.get(
                "KO",
                ""
            )
        ):
            mag_kos[mag].add(
                ko
            )

        # Search all annotation columns because dsrD/dsrEFH may
        # not have reliable KO assignments in the current table.
        complete_row_text = " | ".join(
            str(value)
            for value in row.values()
            if value is not None
        )

        mag_annotation_text[mag].append(
            complete_row_text
        )


# ============================================================
# AUDIT EACH CANDIDATE MAG
# ============================================================

results = []

for mag in sorted(
    candidate_mags
):

    kos = mag_kos[
        mag
    ]

    annotation_text = "\n".join(
        mag_annotation_text[
            mag
        ]
    )

    row = {
        "Final_MAG":
            mag,
        **redox_metadata[mag],
    }


    # --------------------------------------------------------
    # Core pathway
    # --------------------------------------------------------

    for gene, ko in CORE_KOS.items():

        row[
            f"{gene}_{ko}"
        ] = int(
            ko in kos
        )


    # --------------------------------------------------------
    # DsrC supporting evidence
    # --------------------------------------------------------

    dsrC_hits = sorted(
        kos
        & DSRC_KOS
    )

    row[
        "DsrC_related_KOs"
    ] = (
        ";".join(
            dsrC_hits
        )
        if dsrC_hits
        else "None"
    )

    row[
        "DsrC_related_present"
    ] = int(
        bool(
            dsrC_hits
        )
    )


    # --------------------------------------------------------
    # dsrD / dsrEFH annotation-text evidence
    #
    # We deliberately search annotation text rather than
    # pretending that these genes have simple universal KO
    # identifiers in our current annotation framework.
    # --------------------------------------------------------

    dsrD = int(
        text_has_gene_symbol(
            annotation_text,
            "dsrD"
        )
    )

    dsrE = int(
        text_has_gene_symbol(
            annotation_text,
            "dsrE"
        )
    )

    dsrF = int(
        text_has_gene_symbol(
            annotation_text,
            "dsrF"
        )
    )

    dsrH = int(
        text_has_gene_symbol(
            annotation_text,
            "dsrH"
        )
    )

    row[
        "dsrD_annotation_evidence"
    ] = dsrD

    row[
        "dsrE_annotation_evidence"
    ] = dsrE

    row[
        "dsrF_annotation_evidence"
    ] = dsrF

    row[
        "dsrH_annotation_evidence"
    ] = dsrH

    dsrEFH_count = (
        dsrE
        + dsrF
        + dsrH
    )

    row[
        "dsrEFH_members_detected"
    ] = dsrEFH_count

    row[
        "dsrEFH_any"
    ] = int(
        dsrEFH_count > 0
    )


    # --------------------------------------------------------
    # DsrL
    # --------------------------------------------------------

    dsrL_ko = int(
        bool(
            kos
            & DSRL_KOS
        )
    )

    dsrL_text = int(
        text_has_gene_symbol(
            annotation_text,
            "dsrL"
        )
    )

    row[
        "dsrL_K27196"
    ] = dsrL_ko

    row[
        "dsrL_annotation_evidence"
    ] = dsrL_text

    row[
        "dsrL_any_evidence"
    ] = int(
        dsrL_ko
        or dsrL_text
    )


    # --------------------------------------------------------
    # DsrMKJOP
    #
    # Important: this complex is NOT interpreted as directional.
    # --------------------------------------------------------

    dsrmkjop_count = 0

    for gene, ko in DSRMKJOP_KOS.items():

        present = int(
            ko in kos
        )

        row[
            f"{gene}_{ko}"
        ] = present

        dsrmkjop_count += (
            present
        )

    row[
        "DsrMKJOP_members_detected"
    ] = dsrmkjop_count

    row[
        "DsrMKJOP_complete"
    ] = int(
        dsrmkjop_count
        == len(
            DSRMKJOP_KOS
        )
    )


    # --------------------------------------------------------
    # Sox sulfur-oxidation context
    # --------------------------------------------------------

    sox_core_count = 0

    for gene, ko in SOX_CORE_KOS.items():

        present = int(
            ko in kos
        )

        row[
            f"{gene}_{ko}"
        ] = present

        sox_core_count += (
            present
        )

    row[
        "Sox_core_members_detected"
    ] = sox_core_count

    row[
        "Sox_core_complete"
    ] = int(
        sox_core_count
        == len(
            SOX_CORE_KOS
        )
    )


    # --------------------------------------------------------
    # SoxCD
    # --------------------------------------------------------

    soxcd_count = 0

    for gene, ko in SOXCD_KOS.items():

        present = int(
            ko in kos
        )

        row[
            f"{gene}_{ko}"
        ] = present

        soxcd_count += (
            present
        )

    row[
        "SoxCD_members_detected"
    ] = soxcd_count

    row[
        "SoxCD_complete"
    ] = int(
        soxcd_count
        == len(
            SOXCD_KOS
        )
    )


    # --------------------------------------------------------
    # SQR
    # --------------------------------------------------------

    sqr = int(
        SQR_KO
        in kos
    )

    row[
        "sqr_K17218"
    ] = sqr


    # --------------------------------------------------------
    # Independent oxidative-context summary
    # --------------------------------------------------------

    oxidative_context = int(
        row[
            "Sox_core_complete"
        ]
        or row[
            "SoxCD_complete"
        ]
        or sqr
    )

    row[
        "Independent_sulfur_oxidation_context"
    ] = oxidative_context


    # ========================================================
    # CONSERVATIVE DIRECTIONALITY INTERPRETATION
    # ========================================================
    #
    # These are deliberately called "support" rather than
    # definitive classifications.
    #
    # Strong reductive-context pattern:
    #
    #   dsrD detected
    #   AND no dsrEFH annotation evidence
    #
    # Oxidative-context pattern:
    #
    #   dsrEFH evidence
    #   OR independent sulfur-oxidation systems
    #
    # Mixed evidence remains unresolved.
    #
    # DsrAB phylogeny is recommended for unresolved cases and
    # is the final confirmation for directionality.
    # ========================================================

    reductive_context = int(
        dsrD
        and not row[
            "dsrEFH_any"
        ]
    )

    oxidative_dsr_context = int(
        row[
            "dsrEFH_any"
        ]
    )

    row[
        "Reductive_context_support"
    ] = reductive_context

    row[
        "Oxidative_Dsr_context_support"
    ] = oxidative_dsr_context


    # --------------------------------------------------------
    # Provisional interpretation
    # --------------------------------------------------------

    if (
        reductive_context
        and not oxidative_dsr_context
        and not oxidative_context
    ):

        interpretation = (
            "reductive_context_supported"
        )

        phylogeny_priority = (
            "confirm_with_DsrAB_phylogeny"
        )

    elif (
        oxidative_dsr_context
        and not reductive_context
    ):

        interpretation = (
            "oxidative_reverse_Dsr_context_supported"
        )

        phylogeny_priority = (
            "confirm_with_DsrAB_phylogeny"
        )

    elif (
        reductive_context
        and (
            oxidative_dsr_context
            or oxidative_context
        )
    ):

        interpretation = (
            "mixed_context_unresolved"
        )

        phylogeny_priority = (
            "high_priority_DsrAB_phylogeny"
        )

    elif (
        oxidative_context
        and not reductive_context
    ):

        interpretation = (
            "oxidative_context_present_direction_unresolved"
        )

        phylogeny_priority = (
            "high_priority_DsrAB_phylogeny"
        )

    else:

        interpretation = (
            "insufficient_context_direction_unresolved"
        )

        phylogeny_priority = (
            "high_priority_DsrAB_phylogeny"
        )


    row[
        "Provisional_directionality"
    ] = interpretation

    row[
        "DsrAB_phylogeny_recommendation"
    ] = phylogeny_priority

    results.append(
        row
    )


# ============================================================
# WRITE AUDIT TABLE
# ============================================================

if not results:
    raise ValueError(
        "No Dsr-positive MAGs were found."
    )

fieldnames = list(
    results[0].keys()
)

with OUTPUT_TABLE.open(
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
# SUMMARY
# ============================================================

categories = defaultdict(
    int
)

for row in results:

    categories[
        row[
            "Provisional_directionality"
        ]
    ] += 1


with SUMMARY_TABLE.open(
    "w",
    newline=""
) as outfile:

    writer = csv.writer(
        outfile,
        delimiter="\t",
    )

    writer.writerow([
        "Category",
        "MAGs",
        "Percent",
    ])

    for category in sorted(
        categories
    ):

        count = categories[
            category
        ]

        writer.writerow([
            category,
            count,
            f"{100 * count / len(results):.2f}",
        ])


# ============================================================
# PRINT RESULTS
# ============================================================

print()
print(
    "=" * 72
)

print(
    "DSR DIRECTIONALITY GENOMIC-CONTEXT AUDIT"
)

print(
    "=" * 72
)

print()

print(
    f"Dsr-core candidate MAGs assessed: "
    f"{len(results)}"
)

print()

print(
    "Provisional genomic-context interpretation:"
)

for category in sorted(
    categories
):

    count = categories[
        category
    ]

    print(
        f"  {category:<48}"
        f"{count:>3} MAGs "
        f"({100 * count / len(results):6.2f}%)"
    )


print()
print(
    "Per-MAG summary:"
)

print()

for row in results:

    print(
        f"{row['Final_MAG']:<30} "
        f"{row['Provisional_directionality']:<45} "
        f"dsrD={row['dsrD_annotation_evidence']} "
        f"dsrEFH={row['dsrEFH_members_detected']} "
        f"dsrL={row['dsrL_any_evidence']} "
        f"Sox={row['Sox_core_complete']} "
        f"SQR={row['sqr_K17218']}"
    )


print()
print(
    "Outputs:"
)

print(
    f"  Audit table:   "
    f"{OUTPUT_TABLE}"
)

print(
    f"  Summary table: "
    f"{SUMMARY_TABLE}"
)

print()

print(
    "Important:"
)

print(
    "  This audit provides genomic-context evidence only."
)

print(
    "  DsrMKJOP is not treated as directional because it can "
    "participate in both reductive and oxidative Dsr systems."
)

print(
    "  DsrD and DsrEFH annotations are searched across all "
    "available annotation fields rather than inferred from "
    "unverified KO assignments."
)

print(
    "  DsrAB phylogeny remains the recommended final test "
    "for pathway directionality."
)

print(
    "=" * 72
)
