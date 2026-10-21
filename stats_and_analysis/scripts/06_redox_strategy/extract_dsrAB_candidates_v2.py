#!/usr/bin/env python3

"""
Extract DsrA and DsrB proteins from Dsr-core candidate MAGs.

VERSION 2
---------
Correctly searches separate protein directories for:

    MetaBAT2:
    results/checkm2_all_bins/protein_files/

    VAMB:
    results/checkm2_vamb_bins/protein_files/

PURPOSE
-------
The redox analysis identified MAGs containing the
Sat-AprAB-DsrAB core, but DsrAB can occur in both:

    reductive sulfate/sulfite metabolism

and

    oxidative/reverse-Dsr sulfur metabolism.

This script prepares DsrA/DsrB protein sequences for
phylogenetic analysis.

Workflow:
    Dsr-core candidate MAGs
        ->
    identify all K11180 (DsrA) and K11181 (DsrB)
        ->
    identify correct protein FASTA by binner
        ->
    extract amino-acid sequences
        ->
    identify neighbouring DsrA/DsrB pairs
        ->
    prepare DsrA and DsrB FASTA files for phylogeny

Important:
- Multiple DsrA/DsrB copies are retained.
- No arbitrary "first hit" selection is made.
- Adjacent genes are considered strongest pairing evidence.
- Nearby genes on the same contig are retained.
- Ambiguous pairing is explicitly reported.
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

DSR_AUDIT = (
    BASE
    / "results/ecological_strategies/"
      "redox_flexibility/dsr_directionality/"
      "dsr_directionality_audit.tsv"
)

METABAT_PROTEIN_DIR = (
    BASE
    / "results/checkm2_all_bins/protein_files"
)

VAMB_PROTEIN_DIR = (
    BASE
    / "results/checkm2_vamb_bins/protein_files"
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

ALL_DSRA_FASTA = (
    OUTPUT_DIR
    / "dsrA_all_MAG_candidates_v2.faa"
)

ALL_DSRB_FASTA = (
    OUTPUT_DIR
    / "dsrB_all_MAG_candidates_v2.faa"
)

PAIRED_DSRA_FASTA = (
    OUTPUT_DIR
    / "dsrA_paired_MAG_candidates_v2.faa"
)

PAIRED_DSRB_FASTA = (
    OUTPUT_DIR
    / "dsrB_paired_MAG_candidates_v2.faa"
)

GENE_MANIFEST = (
    OUTPUT_DIR
    / "dsrAB_gene_manifest_v2.tsv"
)

PAIR_MANIFEST = (
    OUTPUT_DIR
    / "dsrAB_pair_manifest_v2.tsv"
)

MAG_SUMMARY = (
    OUTPUT_DIR
    / "dsrAB_MAG_summary_v2.tsv"
)


# ============================================================
# CONSTANTS
# ============================================================

DSRA_KO = "K11180"
DSRB_KO = "K11181"

# DsrA and DsrB commonly occur next to one another.
# Retain pairs separated by up to three predicted genes.
MAX_PAIR_DISTANCE = 3


# ============================================================
# FASTA READER
# ============================================================

def read_fasta(path):
    """
    Read protein FASTA.

    Returns:
        dict:
            sequence_id -> amino-acid sequence

    Only the first whitespace-delimited part of the FASTA
    header is used as the sequence ID.
    """

    sequences = {}

    current_id = None
    current_sequence = []

    with path.open("r") as handle:

        for line in handle:

            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:

                    sequences[current_id] = "".join(
                        current_sequence
                    )

                current_id = (
                    line[1:]
                    .split()[0]
                )

                current_sequence = []

            else:

                current_sequence.append(
                    line.strip()
                )

        if current_id is not None:

            sequences[current_id] = "".join(
                current_sequence
            )

    return sequences


# ============================================================
# PARSE PRODIGAL-LIKE PROTEIN ID
# ============================================================

def parse_protein_position(protein_id):
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

        return protein_id, None

    return (
        match.group(1),
        int(match.group(2)),
    )


# ============================================================
# FIND CORRECT PROTEIN FASTA
# ============================================================

def find_protein_fasta(original_bin):
    """
    Locate protein FASTA according to binner.

    VAMB bins contain '__vamb_bin.' in Original_bin.
    Other bins are treated as MetaBAT2 bins.
    """

    filename = (
        f"{original_bin}.faa"
    )

    if "__vamb_bin." in original_bin:

        preferred = (
            VAMB_PROTEIN_DIR
            / filename
        )

        alternative = (
            METABAT_PROTEIN_DIR
            / filename
        )

    else:

        preferred = (
            METABAT_PROTEIN_DIR
            / filename
        )

        alternative = (
            VAMB_PROTEIN_DIR
            / filename
        )

    if preferred.exists():
        return preferred

    # Safety fallback in case a bin ended up in the other directory.
    if alternative.exists():
        return alternative

    return None


# ============================================================
# 1. READ DSR-CORE CANDIDATE MAGs
# ============================================================

candidate_mags = set()
audit_metadata = {}

with DSR_AUDIT.open(
    "r",
    newline=""
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t"
    )

    if reader.fieldnames is None:

        raise ValueError(
            "Dsr directionality audit has no header."
        )

    for row in reader:

        mag = row[
            "Final_MAG"
        ].strip()

        if not mag:
            continue

        candidate_mags.add(
            mag
        )

        audit_metadata[
            mag
        ] = row


if not candidate_mags:

    raise ValueError(
        "No candidate MAGs found in Dsr audit."
    )


# ============================================================
# 2. FIND ALL DsrA / DsrB GENES
# ============================================================

genes = []

with MASTER.open(
    "r",
    newline=""
) as infile:

    reader = csv.DictReader(
        infile,
        delimiter="\t"
    )

    if reader.fieldnames is None:

        raise ValueError(
            "Master table has no header."
        )

    required_columns = {
        "gene_id",
        "Final_MAG",
        "Original_bin",
        "KO",
        "score",
        "E_value",
        "KO_definition",
    }

    missing = (
        required_columns
        - set(reader.fieldnames)
    )

    if missing:

        raise ValueError(
            "Missing columns in master table: "
            + ", ".join(
                sorted(missing)
            )
        )

    for row in reader:

        mag = row[
            "Final_MAG"
        ].strip()

        if mag not in candidate_mags:
            continue

        ko = row[
            "KO"
        ].strip()

        if ko not in {
            DSRA_KO,
            DSRB_KO,
        }:
            continue

        gene_id = row[
            "gene_id"
        ].strip()

        original_bin = row[
            "Original_bin"
        ].strip()

        # gene_id format:
        #
        # Final_MAG__protein_id
        #
        # Example:
        # metabat2__...__bin.78__k127_382766_24
        #
        # Protein FASTA:
        # >k127_382766_24

        prefix = (
            mag
            + "__"
        )

        if gene_id.startswith(
            prefix
        ):

            protein_id = gene_id[
                len(prefix):
            ]

        else:

            # Conservative fallback.
            match = re.search(
                r"(k\d+_.+_[0-9]+)$",
                gene_id
            )

            if match:

                protein_id = (
                    match.group(1)
                )

            else:

                protein_id = (
                    gene_id
                )

        contig, gene_number = (
            parse_protein_position(
                protein_id
            )
        )

        genes.append({
            "gene_id":
                gene_id,

            "protein_id":
                protein_id,

            "Final_MAG":
                mag,

            "Original_bin":
                original_bin,

            "KO":
                ko,

            "Dsr_subunit":
                (
                    "DsrA"
                    if ko == DSRA_KO
                    else "DsrB"
                ),

            "contig":
                contig,

            "gene_number":
                gene_number,

            "score":
                row.get(
                    "score",
                    ""
                ),

            "E_value":
                row.get(
                    "E_value",
                    ""
                ),

            "KO_definition":
                row.get(
                    "KO_definition",
                    ""
                ),
        })


if not genes:

    raise ValueError(
        "No DsrA/DsrB genes found."
    )


# ============================================================
# 3. EXTRACT PROTEIN SEQUENCES
# ============================================================

genes_by_bin = defaultdict(
    list
)

for gene in genes:

    genes_by_bin[
        gene[
            "Original_bin"
        ]
    ].append(
        gene
    )


missing_fasta_files = []
missing_sequences = []

used_protein_files = set()


for original_bin, bin_genes in (
    genes_by_bin.items()
):

    faa = find_protein_fasta(
        original_bin
    )

    if faa is None:

        missing_fasta_files.append(
            original_bin
        )

        for gene in bin_genes:

            gene[
                "sequence"
            ] = ""

            gene[
                "sequence_found"
            ] = 0

            gene[
                "protein_fasta"
            ] = "NOT_FOUND"

        continue


    used_protein_files.add(
        str(faa)
    )

    sequences = read_fasta(
        faa
    )


    for gene in bin_genes:

        protein_id = gene[
            "protein_id"
        ]

        sequence = sequences.get(
            protein_id,
            ""
        )

        gene[
            "sequence"
        ] = sequence

        gene[
            "sequence_found"
        ] = int(
            bool(sequence)
        )

        gene[
            "protein_fasta"
        ] = str(
            faa
        )

        if not sequence:

            missing_sequences.append(
                (
                    gene[
                        "Final_MAG"
                    ],
                    protein_id,
                    str(faa),
                )
            )


# ============================================================
# 4. WRITE ALL DsrA AND DsrB FASTA FILES
# ============================================================

with ALL_DSRA_FASTA.open(
    "w"
) as dsra_out, ALL_DSRB_FASTA.open(
    "w"
) as dsrb_out:

    for gene in genes:

        if not gene[
            "sequence_found"
        ]:
            continue

        header = (
            f"{gene['Final_MAG']}"
            f"|{gene['Dsr_subunit']}"
            f"|{gene['protein_id']}"
        )

        fasta_record = (
            f">{header}\n"
            f"{gene['sequence']}\n"
        )

        if gene[
            "Dsr_subunit"
        ] == "DsrA":

            dsra_out.write(
                fasta_record
            )

        else:

            dsrb_out.write(
                fasta_record
            )


# ============================================================
# 5. IDENTIFY CANDIDATE DsrA-DsrB PAIRS
# ============================================================

genes_by_mag = defaultdict(
    list
)

for gene in genes:

    genes_by_mag[
        gene[
            "Final_MAG"
        ]
    ].append(
        gene
    )


all_candidate_pairs = []

for mag in sorted(
    candidate_mags
):

    mag_genes = genes_by_mag.get(
        mag,
        []
    )

    dsra_genes = [
        gene
        for gene in mag_genes
        if gene[
            "Dsr_subunit"
        ] == "DsrA"
    ]

    dsrb_genes = [
        gene
        for gene in mag_genes
        if gene[
            "Dsr_subunit"
        ] == "DsrB"
    ]


    for a in dsra_genes:

        for b in dsrb_genes:

            # Must occur on the same contig.
            if (
                a[
                    "contig"
                ]
                != b[
                    "contig"
                ]
            ):
                continue

            if (
                a[
                    "gene_number"
                ]
                is None
                or
                b[
                    "gene_number"
                ]
                is None
            ):
                continue

            distance = abs(
                a[
                    "gene_number"
                ]
                - b[
                    "gene_number"
                ]
            )

            if (
                distance
                <= MAX_PAIR_DISTANCE
            ):

                all_candidate_pairs.append({
                    "Final_MAG":
                        mag,

                    "contig":
                        a[
                            "contig"
                        ],

                    "DsrA_protein_id":
                        a[
                            "protein_id"
                        ],

                    "DsrB_protein_id":
                        b[
                            "protein_id"
                        ],

                    "DsrA_gene_number":
                        a[
                            "gene_number"
                        ],

                    "DsrB_gene_number":
                        b[
                            "gene_number"
                        ],

                    "gene_distance":
                        distance,

                    "adjacent":
                        int(
                            distance == 1
                        ),

                    "DsrA_sequence_found":
                        a[
                            "sequence_found"
                        ],

                    "DsrB_sequence_found":
                        b[
                            "sequence_found"
                        ],

                    "DsrA_sequence":
                        a[
                            "sequence"
                        ],

                    "DsrB_sequence":
                        b[
                            "sequence"
                        ],
                })


# ============================================================
# 6. IDENTIFY AMBIGUOUS PAIRING
# ============================================================

partners_for_a = defaultdict(
    int
)

partners_for_b = defaultdict(
    int
)


for pair in all_candidate_pairs:

    key_a = (
        pair[
            "Final_MAG"
        ],
        pair[
            "DsrA_protein_id"
        ],
    )

    key_b = (
        pair[
            "Final_MAG"
        ],
        pair[
            "DsrB_protein_id"
        ],
    )

    partners_for_a[
        key_a
    ] += 1

    partners_for_b[
        key_b
    ] += 1


for pair in all_candidate_pairs:

    key_a = (
        pair[
            "Final_MAG"
        ],
        pair[
            "DsrA_protein_id"
        ],
    )

    key_b = (
        pair[
            "Final_MAG"
        ],
        pair[
            "DsrB_protein_id"
        ],
    )

    pair[
        "pairing_ambiguous"
    ] = int(
        partners_for_a[
            key_a
        ] > 1
        or
        partners_for_b[
            key_b
        ] > 1
    )


    if (
        pair[
            "adjacent"
        ]
        and not pair[
            "pairing_ambiguous"
        ]
    ):

        confidence = (
            "high_adjacent_unique"
        )

    elif (
        pair[
            "adjacent"
        ]
        and pair[
            "pairing_ambiguous"
        ]
    ):

        confidence = (
            "adjacent_but_ambiguous"
        )

    elif (
        pair[
            "gene_distance"
        ]
        <= MAX_PAIR_DISTANCE
        and not pair[
            "pairing_ambiguous"
        ]
    ):

        confidence = (
            "nearby_unique"
        )

    else:

        confidence = (
            "nearby_ambiguous"
        )

    pair[
        "pair_confidence"
    ] = confidence


# ============================================================
# 7. SELECT PAIRS FOR INITIAL PHYLOGENY
# ============================================================

selected_pairs = [
    pair
    for pair in all_candidate_pairs
    if (
        not pair[
            "pairing_ambiguous"
        ]
        and pair[
            "DsrA_sequence_found"
        ]
        and pair[
            "DsrB_sequence_found"
        ]
    )
]


# ============================================================
# 8. WRITE PAIRED FASTA FILES
# ============================================================
#
# DsrA and DsrB FASTA files use the SAME pair ID.
#
# We will later:
#
#   align DsrA separately
#   align DsrB separately
#   concatenate aligned sequences by pair ID
#
# ============================================================

with PAIRED_DSRA_FASTA.open(
    "w"
) as dsra_out, PAIRED_DSRB_FASTA.open(
    "w"
) as dsrb_out:

    for number, pair in enumerate(
        selected_pairs,
        start=1
    ):

        pair_id = (
            f"Saltmarsh_DsrAB_"
            f"{number:03d}"
        )

        pair[
            "pair_id"
        ] = pair_id

        metadata = (
            f"{pair['Final_MAG']}"
            f"|contig={pair['contig']}"
            f"|distance={pair['gene_distance']}"
            f"|confidence={pair['pair_confidence']}"
        )

        dsra_out.write(
            f">{pair_id} {metadata}\n"
            f"{pair['DsrA_sequence']}\n"
        )

        dsrb_out.write(
            f">{pair_id} {metadata}\n"
            f"{pair['DsrB_sequence']}\n"
        )


for pair in all_candidate_pairs:

    if "pair_id" not in pair:

        pair[
            "pair_id"
        ] = ""


# ============================================================
# 9. WRITE GENE MANIFEST
# ============================================================

gene_fields = [
    "Final_MAG",
    "Original_bin",
    "Dsr_subunit",
    "KO",
    "gene_id",
    "protein_id",
    "contig",
    "gene_number",
    "score",
    "E_value",
    "KO_definition",
    "sequence_found",
    "protein_fasta",
]


with GENE_MANIFEST.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=gene_fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()

    for gene in genes:

        writer.writerow(
            gene
        )


# ============================================================
# 10. WRITE PAIR MANIFEST
# ============================================================

pair_fields = [
    "pair_id",
    "Final_MAG",
    "contig",
    "DsrA_protein_id",
    "DsrB_protein_id",
    "DsrA_gene_number",
    "DsrB_gene_number",
    "gene_distance",
    "adjacent",
    "pairing_ambiguous",
    "pair_confidence",
    "DsrA_sequence_found",
    "DsrB_sequence_found",
]


with PAIR_MANIFEST.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=pair_fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()

    for pair in all_candidate_pairs:

        writer.writerow(
            pair
        )


# ============================================================
# 11. MAG-LEVEL SUMMARY
# ============================================================

summary_rows = []

selected_mags = {
    pair[
        "Final_MAG"
    ]
    for pair in selected_pairs
}


for mag in sorted(
    candidate_mags
):

    mag_genes = genes_by_mag.get(
        mag,
        []
    )

    dsra_count = sum(
        1
        for gene in mag_genes
        if gene[
            "Dsr_subunit"
        ] == "DsrA"
    )

    dsrb_count = sum(
        1
        for gene in mag_genes
        if gene[
            "Dsr_subunit"
        ] == "DsrB"
    )

    mag_pairs = [
        pair
        for pair in all_candidate_pairs
        if pair[
            "Final_MAG"
        ] == mag
    ]

    mag_selected = [
        pair
        for pair in selected_pairs
        if pair[
            "Final_MAG"
        ] == mag
    ]

    adjacent_pairs = sum(
        pair[
            "adjacent"
        ]
        for pair in mag_pairs
    )

    ambiguous_pairs = sum(
        pair[
            "pairing_ambiguous"
        ]
        for pair in mag_pairs
    )

    summary_rows.append({
        "Final_MAG":
            mag,

        "GTDB_phylum":
            audit_metadata[
                mag
            ].get(
                "GTDB_phylum",
                ""
            ),

        "GTDB_order":
            audit_metadata[
                mag
            ].get(
                "GTDB_order",
                ""
            ),

        "Provisional_directionality":
            audit_metadata[
                mag
            ].get(
                "Provisional_directionality",
                ""
            ),

        "DsrA_hits":
            dsra_count,

        "DsrB_hits":
            dsrb_count,

        "same_contig_candidate_pairs":
            len(
                mag_pairs
            ),

        "adjacent_pairs":
            adjacent_pairs,

        "ambiguous_candidate_pairs":
            ambiguous_pairs,

        "selected_phylogeny_pairs":
            len(
                mag_selected
            ),

        "ready_for_initial_phylogeny":
            int(
                mag in selected_mags
            ),
    })


summary_fields = [
    "Final_MAG",
    "GTDB_phylum",
    "GTDB_order",
    "Provisional_directionality",
    "DsrA_hits",
    "DsrB_hits",
    "same_contig_candidate_pairs",
    "adjacent_pairs",
    "ambiguous_candidate_pairs",
    "selected_phylogeny_pairs",
    "ready_for_initial_phylogeny",
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
# 12. PRINT SUMMARY
# ============================================================

dsra_total = sum(
    1
    for gene in genes
    if gene[
        "Dsr_subunit"
    ] == "DsrA"
)

dsrb_total = sum(
    1
    for gene in genes
    if gene[
        "Dsr_subunit"
    ] == "DsrB"
)

sequence_success = sum(
    gene[
        "sequence_found"
    ]
    for gene in genes
)

adjacent_total = sum(
    pair[
        "adjacent"
    ]
    for pair in all_candidate_pairs
)

ambiguous_total = sum(
    pair[
        "pairing_ambiguous"
    ]
    for pair in all_candidate_pairs
)


print()
print(
    "=" * 72
)

print(
    "DsrA / DsrB EXTRACTION AND PAIRING - V2"
)

print(
    "=" * 72
)

print()

print(
    f"Candidate MAGs: "
    f"{len(candidate_mags)}"
)

print(
    f"DsrA hits (K11180): "
    f"{dsra_total}"
)

print(
    f"DsrB hits (K11181): "
    f"{dsrb_total}"
)

print(
    f"Protein sequences successfully extracted: "
    f"{sequence_success}/{len(genes)}"
)

print(
    f"Protein FASTA files used: "
    f"{len(used_protein_files)}"
)

print()

print(
    f"Same-contig candidate DsrAB pairs "
    f"(distance <= {MAX_PAIR_DISTANCE} genes): "
    f"{len(all_candidate_pairs)}"
)

print(
    f"Adjacent DsrA-DsrB pairs: "
    f"{adjacent_total}"
)

print(
    f"Ambiguous candidate pairs: "
    f"{ambiguous_total}"
)

print(
    f"Non-ambiguous pairs selected for initial phylogeny: "
    f"{len(selected_pairs)}"
)

print(
    f"MAGs represented by selected pairs: "
    f"{len(selected_mags)}/{len(candidate_mags)}"
)


if missing_fasta_files:

    print()
    print(
        "WARNING - protein FASTA not found for bins:"
    )

    for original_bin in missing_fasta_files:

        print(
            f"  {original_bin}"
        )


if missing_sequences:

    print()
    print(
        "WARNING - protein IDs not found in located FASTA:"
    )

    for mag, protein_id, faa in missing_sequences:

        print(
            f"  {mag} | "
            f"{protein_id} | "
            f"{faa}"
        )


print()
print(
    "Outputs:"
)

print(
    f"  All DsrA FASTA:       "
    f"{ALL_DSRA_FASTA}"
)

print(
    f"  All DsrB FASTA:       "
    f"{ALL_DSRB_FASTA}"
)

print(
    f"  Paired DsrA FASTA:    "
    f"{PAIRED_DSRA_FASTA}"
)

print(
    f"  Paired DsrB FASTA:    "
    f"{PAIRED_DSRB_FASTA}"
)

print(
    f"  Gene manifest:        "
    f"{GENE_MANIFEST}"
)

print(
    f"  Pair manifest:        "
    f"{PAIR_MANIFEST}"
)

print(
    f"  MAG summary:          "
    f"{MAG_SUMMARY}"
)

print()

print(
    "Next:"
)

print(
    "  If sequence recovery is complete, add validated "
    "reductive and oxidative reference DsrAB sequences, "
    "align DsrA and DsrB separately, concatenate aligned "
    "pairs, and reconstruct the DsrAB phylogeny."
)

print(
    "=" * 72
)
