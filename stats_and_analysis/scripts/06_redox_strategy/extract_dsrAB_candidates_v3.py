#!/usr/bin/env python3

"""
Prepare the FINAL saltmarsh DsrA/DsrB pairs for phylogenetic analysis.

VERSION 3
---------

This script uses the successfully extracted V2 DsrA/DsrB gene manifest.

PAIRING LOGIC
-------------

For each MAG:

1. Prefer an adjacent DsrA-DsrB pair.
2. Otherwise prefer a nearby same-contig pair.
3. Otherwise prefer another unique same-contig pair.
4. If the MAG contains exactly one DsrA and one DsrB on
   different contigs, retain them as a lower-confidence
   cross-contig unique pair.
5. Never make an arbitrary choice when multiple equally
   plausible pairs exist.

The final FASTA files use the same Pair_ID for DsrA and DsrB,
allowing them to be aligned separately and concatenated later.

This script does NOT infer reductive versus oxidative direction.
That will be determined from DsrAB phylogeny plus genomic context.
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

INPUT_MANIFEST = (
    BASE
    / "results/ecological_strategies/"
      "redox_flexibility/dsr_directionality/"
      "dsrAB_phylogeny/"
      "dsrAB_gene_manifest_v2.tsv"
)

DSR_AUDIT = (
    BASE
    / "results/ecological_strategies/"
      "redox_flexibility/dsr_directionality/"
      "dsr_directionality_audit.tsv"
)

OUTPUT_DIR = (
    BASE
    / "results/ecological_strategies/"
      "redox_flexibility/dsr_directionality/"
      "dsrAB_phylogeny"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

FINAL_DSRA_FASTA = (
    OUTPUT_DIR
    / "dsrA_final_MAG_pairs_v3.faa"
)

FINAL_DSRB_FASTA = (
    OUTPUT_DIR
    / "dsrB_final_MAG_pairs_v3.faa"
)

PAIR_MANIFEST = (
    OUTPUT_DIR
    / "dsrAB_pair_manifest_v3.tsv"
)

MAG_SUMMARY = (
    OUTPUT_DIR
    / "dsrAB_MAG_summary_v3.tsv"
)


# ============================================================
# HELPERS
# ============================================================

def read_fasta_record(path, target_id):
    """
    Retrieve one protein sequence from a FASTA file.
    """

    current_id = None
    sequence = []

    with Path(path).open("r") as handle:

        for line in handle:

            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):

                if (
                    current_id == target_id
                ):

                    return "".join(
                        sequence
                    )

                current_id = (
                    line[1:]
                    .split()[0]
                )

                sequence = []

            else:

                sequence.append(
                    line.strip()
                )

        if current_id == target_id:

            return "".join(
                sequence
            )

    return ""


def parse_position(protein_id):
    """
    Example:

        k127_382766_24

    becomes:

        contig = k127_382766
        gene_number = 24
    """

    match = re.match(
        r"(.+)_([0-9]+)$",
        protein_id
    )

    if not match:

        return (
            protein_id,
            None,
        )

    return (
        match.group(1),
        int(
            match.group(2)
        ),
    )


def classify_pair(a, b):
    """
    Classify DsrA-DsrB physical relationship.

    Returns:
        rank
        relationship
        gene_distance

    Lower rank = stronger pairing evidence.
    """

    if (
        a["contig"]
        == b["contig"]
        and
        a["gene_number"]
        is not None
        and
        b["gene_number"]
        is not None
    ):

        distance = abs(
            a["gene_number"]
            - b["gene_number"]
        )

        if distance == 1:

            return (
                0,
                "adjacent_same_contig",
                distance,
            )

        elif distance <= 3:

            return (
                1,
                "nearby_same_contig",
                distance,
            )

        else:

            return (
                2,
                "same_contig_distant",
                distance,
            )

    return (
        3,
        "cross_contig_unique_candidate",
        "",
    )


# ============================================================
# READ DIRECTIONALITY METADATA
# ============================================================

audit_metadata = {}

with DSR_AUDIT.open(
    "r",
    newline=""
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t",
    )

    for row in reader:

        mag = row[
            "Final_MAG"
        ].strip()

        audit_metadata[
            mag
        ] = row


# ============================================================
# READ V2 GENE MANIFEST
# ============================================================

genes_by_mag = defaultdict(
    list
)

with INPUT_MANIFEST.open(
    "r",
    newline=""
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t",
    )

    required = {
        "Final_MAG",
        "Dsr_subunit",
        "protein_id",
        "protein_fasta",
        "sequence_found",
    }

    if reader.fieldnames is None:

        raise ValueError(
            "Input gene manifest has no header."
        )

    missing = (
        required
        - set(reader.fieldnames)
    )

    if missing:

        raise ValueError(
            "Missing required columns: "
            + ", ".join(
                sorted(missing)
            )
        )

    for row in reader:

        if row[
            "sequence_found"
        ].strip() != "1":

            continue

        protein_id = row[
            "protein_id"
        ].strip()

        contig, gene_number = (
            parse_position(
                protein_id
            )
        )

        sequence = read_fasta_record(
            row[
                "protein_fasta"
            ],
            protein_id,
        )

        if not sequence:

            raise ValueError(
                "Could not recover sequence for "
                f"{row['Final_MAG']} | "
                f"{protein_id}"
            )

        gene = dict(
            row
        )

        gene[
            "contig"
        ] = contig

        gene[
            "gene_number"
        ] = gene_number

        gene[
            "sequence"
        ] = sequence

        genes_by_mag[
            row[
                "Final_MAG"
            ]
        ].append(
            gene
        )


# ============================================================
# PAIR DsrA AND DsrB WITHIN EACH MAG
# ============================================================

selected_pairs = []
summary_rows = []


for mag in sorted(
    genes_by_mag
):

    genes = genes_by_mag[
        mag
    ]

    dsra = [
        gene
        for gene in genes
        if gene[
            "Dsr_subunit"
        ] == "DsrA"
    ]

    dsrb = [
        gene
        for gene in genes
        if gene[
            "Dsr_subunit"
        ] == "DsrB"
    ]

    candidate_pairs = []

    for a in dsra:

        for b in dsrb:

            (
                rank,
                relationship,
                distance,
            ) = classify_pair(
                a,
                b,
            )

            candidate_pairs.append({
                "Final_MAG":
                    mag,

                "DsrA":
                    a,

                "DsrB":
                    b,

                "rank":
                    rank,

                "relationship":
                    relationship,

                "gene_distance":
                    distance,
            })


    selected = None
    decision = ""


    if not dsra:

        decision = (
            "no_DsrA"
        )


    elif not dsrb:

        decision = (
            "no_DsrB"
        )


    elif candidate_pairs:

        best_rank = min(
            pair[
                "rank"
            ]
            for pair in candidate_pairs
        )

        best_pairs = [
            pair
            for pair in candidate_pairs
            if pair[
                "rank"
            ] == best_rank
        ]


        # ----------------------------------------------
        # Unique best pair
        # ----------------------------------------------

        if len(
            best_pairs
        ) == 1:

            selected = (
                best_pairs[0]
            )

            if (
                selected[
                    "relationship"
                ]
                == "cross_contig_unique_candidate"
            ):

                # Cross-contig pairing is accepted only
                # when exactly one DsrA and one DsrB exist.
                if (
                    len(dsra) == 1
                    and
                    len(dsrb) == 1
                ):

                    decision = (
                        "selected_cross_contig_unique"
                    )

                else:

                    selected = None

                    decision = (
                        "cross_contig_not_unique"
                    )

            else:

                decision = (
                    "selected_"
                    + selected[
                        "relationship"
                    ]
                )


        # ----------------------------------------------
        # More than one equally good pair
        # ----------------------------------------------

        else:

            selected = None

            decision = (
                "ambiguous_multiple_best_pairs"
            )


    if selected is not None:

        selected[
            "pairing_decision"
        ] = decision

        selected_pairs.append(
            selected
        )


    metadata = audit_metadata.get(
        mag,
        {}
    )

    summary_rows.append({
        "Final_MAG":
            mag,

        "GTDB_phylum":
            metadata.get(
                "GTDB_phylum",
                ""
            ),

        "GTDB_order":
            metadata.get(
                "GTDB_order",
                ""
            ),

        "Provisional_directionality":
            metadata.get(
                "Provisional_directionality",
                ""
            ),

        "DsrA_hits":
            len(
                dsra
            ),

        "DsrB_hits":
            len(
                dsrb
            ),

        "Candidate_pair_combinations":
            len(
                candidate_pairs
            ),

        "Pair_selected":
            int(
                selected
                is not None
            ),

        "Pairing_decision":
            decision,

        "Pair_relationship":
            (
                selected[
                    "relationship"
                ]
                if selected
                is not None
                else ""
            ),

        "Ready_for_phylogeny":
            int(
                selected
                is not None
            ),
    })


# ============================================================
# ASSIGN FINAL PAIR IDs
# ============================================================

for index, pair in enumerate(
    selected_pairs,
    start=1,
):

    pair[
        "Pair_ID"
    ] = (
        f"Saltmarsh_DsrAB_"
        f"{index:03d}"
    )


# ============================================================
# WRITE FINAL FASTA FILES
# ============================================================

with FINAL_DSRA_FASTA.open(
    "w"
) as dsra_out, FINAL_DSRB_FASTA.open(
    "w"
) as dsrb_out:

    for pair in selected_pairs:

        pair_id = pair[
            "Pair_ID"
        ]

        mag = pair[
            "Final_MAG"
        ]

        relationship = pair[
            "relationship"
        ]


        dsra_out.write(
            f">{pair_id} "
            f"{mag}"
            f"|relationship={relationship}\n"
            f"{pair['DsrA']['sequence']}\n"
        )


        dsrb_out.write(
            f">{pair_id} "
            f"{mag}"
            f"|relationship={relationship}\n"
            f"{pair['DsrB']['sequence']}\n"
        )


# ============================================================
# WRITE PAIR MANIFEST
# ============================================================

pair_rows = []

for pair in selected_pairs:

    pair_rows.append({
        "Pair_ID":
            pair[
                "Pair_ID"
            ],

        "Final_MAG":
            pair[
                "Final_MAG"
            ],

        "Pair_relationship":
            pair[
                "relationship"
            ],

        "Gene_distance":
            pair[
                "gene_distance"
            ],

        "Pairing_decision":
            pair[
                "pairing_decision"
            ],

        "DsrA_protein_id":
            pair[
                "DsrA"
            ][
                "protein_id"
            ],

        "DsrB_protein_id":
            pair[
                "DsrB"
            ][
                "protein_id"
            ],

        "DsrA_contig":
            pair[
                "DsrA"
            ][
                "contig"
            ],

        "DsrB_contig":
            pair[
                "DsrB"
            ][
                "contig"
            ],

        "DsrA_KOfam_score":
            pair[
                "DsrA"
            ].get(
                "score",
                ""
            ),

        "DsrB_KOfam_score":
            pair[
                "DsrB"
            ].get(
                "score",
                ""
            ),

        "DsrA_E_value":
            pair[
                "DsrA"
            ].get(
                "E_value",
                ""
            ),

        "DsrB_E_value":
            pair[
                "DsrB"
            ].get(
                "E_value",
                ""
            ),
    })


pair_fields = [
    "Pair_ID",
    "Final_MAG",
    "Pair_relationship",
    "Gene_distance",
    "Pairing_decision",
    "DsrA_protein_id",
    "DsrB_protein_id",
    "DsrA_contig",
    "DsrB_contig",
    "DsrA_KOfam_score",
    "DsrB_KOfam_score",
    "DsrA_E_value",
    "DsrB_E_value",
]


with PAIR_MANIFEST.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=pair_fields,
        delimiter="\t",
    )

    writer.writeheader()

    writer.writerows(
        pair_rows
    )


# ============================================================
# WRITE MAG SUMMARY
# ============================================================

summary_fields = [
    "Final_MAG",
    "GTDB_phylum",
    "GTDB_order",
    "Provisional_directionality",
    "DsrA_hits",
    "DsrB_hits",
    "Candidate_pair_combinations",
    "Pair_selected",
    "Pairing_decision",
    "Pair_relationship",
    "Ready_for_phylogeny",
]


with MAG_SUMMARY.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=summary_fields,
        delimiter="\t",
    )

    writer.writeheader()

    writer.writerows(
        summary_rows
    )


# ============================================================
# PRINT SUMMARY
# ============================================================

relationship_counts = defaultdict(
    int
)

for pair in selected_pairs:

    relationship_counts[
        pair[
            "relationship"
        ]
    ] += 1


print()
print(
    "=" * 72
)

print(
    "FINAL DsrA / DsrB PHYLOGENY PAIR PREPARATION - V3"
)

print(
    "=" * 72
)

print()

print(
    f"MAGs containing extracted Dsr genes: "
    f"{len(genes_by_mag)}"
)

print(
    f"MAGs with a final DsrAB pair: "
    f"{len(selected_pairs)}"
)

print()

print(
    "Pair relationships:"
)

for relationship in sorted(
    relationship_counts
):

    print(
        f"  {relationship:<35}"
        f"{relationship_counts[relationship]:>3}"
    )


print()
print(
    "MAG-level pairing:"
)

for row in summary_rows:

    print(
        f"  {row['Final_MAG']:<75} "
        f"A={row['DsrA_hits']} "
        f"B={row['DsrB_hits']} "
        f"{row['Pairing_decision']}"
    )


print()
print(
    "Outputs:"
)

print(
    f"  Final DsrA FASTA: "
    f"{FINAL_DSRA_FASTA}"
)

print(
    f"  Final DsrB FASTA: "
    f"{FINAL_DSRB_FASTA}"
)

print(
    f"  Pair manifest:    "
    f"{PAIR_MANIFEST}"
)

print(
    f"  MAG summary:      "
    f"{MAG_SUMMARY}"
)

print()

print(
    "Interpretation:"
)

print(
    "  Adjacent same-contig DsrAB pairs have the strongest "
    "genomic-context support."
)

print(
    "  Cross-contig pairs are retained only when the MAG "
    "contains exactly one DsrA and exactly one DsrB."
)

print(
    "  Cross-contig status must be retained as lower-confidence "
    "metadata during phylogenetic interpretation."
)

print()

print(
    "Next step:"
)

print(
    "  Add validated reductive and oxidative DsrAB reference "
    "sequences, align DsrA and DsrB separately, concatenate "
    "the aligned pairs, and reconstruct the DsrAB phylogeny."
)

print(
    "=" * 72
)
