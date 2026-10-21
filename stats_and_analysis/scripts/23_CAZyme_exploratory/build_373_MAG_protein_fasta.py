#!/usr/bin/env python3

"""
Build a combined protein FASTA + matching GFF3 for the
373 final dereplicated MAGs.

Purpose
-------
Prepare dbCAN inputs while preserving exact MAG identity.

Outputs
-------
1. final_373_MAGs_proteins.faa
   Protein headers:
       >Final_MAG|protein_id

2. final_373_MAGs_genes.gff3
   Gene IDs exactly match the protein FASTA IDs.

3. final_373_MAGs_protein_manifest.tsv
   Mapping:
       Final_MAG
       Original_bin
       protein_id
       dbCAN_gene_id
       contig
       coordinates
       strand
       source FASTA

This allows us to run:
    - CAZyme annotation
    - substrate-oriented CAZyme analysis
    - later CGC analysis

The script only uses MAGs present in the final 373-MAG master table.
"""

import csv
import re
from collections import OrderedDict
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

# Preferred source because filenames should correspond directly
# to Final_MAG names.
PRODIGAL_PROTEIN_DIR = (
    BASE
    / "results/functional_annotation/prodigal/proteins"
)

# Fallback locations
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
      "salicornia_cazyme/input"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_FASTA = (
    OUTPUT_DIR
    / "final_373_MAGs_proteins.faa"
)

OUTPUT_GFF = (
    OUTPUT_DIR
    / "final_373_MAGs_genes.gff3"
)

OUTPUT_MANIFEST = (
    OUTPUT_DIR
    / "final_373_MAGs_protein_manifest.tsv"
)

OUTPUT_MAG_MANIFEST = (
    OUTPUT_DIR
    / "final_373_MAGs_source_files.tsv"
)


# ============================================================
# READ FINAL MAG LIST
# ============================================================

mag_info = OrderedDict()

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

    required = {
        "Final_MAG",
        "Original_bin",
        "Completeness",
        "Contamination",
        "Site",
        "Treatment",
        "Compartment",
        "GTDB_phylum",
        "GTDB_order",
    }

    missing = (
        required
        - set(reader.fieldnames)
    )

    if missing:
        raise ValueError(
            "Missing master-table columns: "
            + ", ".join(
                sorted(missing)
            )
        )

    for row in reader:

        mag = row[
            "Final_MAG"
        ].strip()

        if not mag:
            continue

        if mag not in mag_info:

            mag_info[mag] = {
                "Original_bin":
                    row[
                        "Original_bin"
                    ].strip(),

                "Completeness":
                    row[
                        "Completeness"
                    ],

                "Contamination":
                    row[
                        "Contamination"
                    ],

                "Site":
                    row[
                        "Site"
                    ],

                "Treatment":
                    row[
                        "Treatment"
                    ],

                "Compartment":
                    row[
                        "Compartment"
                    ],

                "GTDB_phylum":
                    row[
                        "GTDB_phylum"
                    ],

                "GTDB_order":
                    row[
                        "GTDB_order"
                    ],
            }


if len(mag_info) != 373:

    print(
        "WARNING: expected 373 final MAGs, "
        f"found {len(mag_info)}."
    )


# ============================================================
# LOCATE PROTEIN FASTA
# ============================================================

def find_protein_fasta(
    final_mag,
    original_bin
):
    """
    Try the final Prodigal protein directory first,
    then the CheckM2 binner-specific directories.
    """

    candidates = [

        # Preferred final-MAG protein set
        PRODIGAL_PROTEIN_DIR
        / f"{final_mag}.faa",

        # MetaBAT2 / generic CheckM2
        METABAT_PROTEIN_DIR
        / f"{original_bin}.faa",

        # VAMB CheckM2
        VAMB_PROTEIN_DIR
        / f"{original_bin}.faa",
    ]

    for path in candidates:

        if path.exists():
            return path

    return None


# ============================================================
# FASTA PARSER
# ============================================================

def read_prodigal_fasta(path):
    """
    Parse Prodigal protein FASTA.

    Expected header style:

    >k127_1296662_1 # 3 # 920 # 1 # ID=...

    Returns dictionaries containing:
      protein_id
      start
      end
      strand
      original_header
      sequence
    """

    records = []

    current_header = None
    sequence_parts = []

    def finish_record(
        header,
        sequence
    ):

        if header is None:
            return

        parts = [
            item.strip()
            for item in header.split("#")
        ]

        protein_id = (
            parts[0]
            .lstrip(">")
            .split()[0]
        )

        start = None
        end = None
        strand_raw = None

        if len(parts) >= 4:

            try:
                start = int(
                    parts[1]
                )

                end = int(
                    parts[2]
                )

                strand_raw = int(
                    parts[3]
                )

            except ValueError:
                pass

        if strand_raw == 1:
            strand = "+"

        elif strand_raw == -1:
            strand = "-"

        else:
            strand = "."

        # Contig = protein ID without final gene-number suffix
        match = re.match(
            r"(.+)_([0-9]+)$",
            protein_id
        )

        if match:
            contig = match.group(1)
            gene_number = match.group(2)

        else:
            contig = protein_id
            gene_number = ""

        records.append({
            "protein_id":
                protein_id,

            "contig":
                contig,

            "gene_number":
                gene_number,

            "start":
                start,

            "end":
                end,

            "strand":
                strand,

            "original_header":
                header,

            "sequence":
                "".join(
                    sequence
                ),
        })


    with path.open("r") as infile:

        for line in infile:

            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):

                finish_record(
                    current_header,
                    sequence_parts
                )

                current_header = line

                sequence_parts = []

            else:

                sequence_parts.append(
                    line.strip()
                )

        finish_record(
            current_header,
            sequence_parts
        )

    return records


# ============================================================
# BUILD COMBINED DATASET
# ============================================================

protein_manifest_rows = []
mag_manifest_rows = []

missing_fasta_mags = []
duplicate_dbcan_ids = set()
seen_dbcan_ids = set()

total_proteins = 0
mags_written = 0
genes_with_coordinates = 0


with OUTPUT_FASTA.open(
    "w"
) as fasta_out, OUTPUT_GFF.open(
    "w"
) as gff_out:

    gff_out.write(
        "##gff-version 3\n"
    )

    for final_mag, metadata in (
        mag_info.items()
    ):

        original_bin = metadata[
            "Original_bin"
        ]

        protein_fasta = (
            find_protein_fasta(
                final_mag,
                original_bin
            )
        )

        if protein_fasta is None:

            missing_fasta_mags.append(
                final_mag
            )

            mag_manifest_rows.append({
                "Final_MAG":
                    final_mag,

                "Original_bin":
                    original_bin,

                "Protein_FASTA":
                    "NOT_FOUND",

                "Proteins":
                    0,

                "Status":
                    "missing_protein_fasta",
            })

            continue


        records = read_prodigal_fasta(
            protein_fasta
        )


        if not records:

            mag_manifest_rows.append({
                "Final_MAG":
                    final_mag,

                "Original_bin":
                    original_bin,

                "Protein_FASTA":
                    str(
                        protein_fasta
                    ),

                "Proteins":
                    0,

                "Status":
                    "empty_protein_fasta",
            })

            continue


        mags_written += 1


        for record in records:

            protein_id = record[
                "protein_id"
            ]

            dbcan_gene_id = (
                f"{final_mag}|"
                f"{protein_id}"
            )

            if dbcan_gene_id in seen_dbcan_ids:

                duplicate_dbcan_ids.add(
                    dbcan_gene_id
                )

                continue

            seen_dbcan_ids.add(
                dbcan_gene_id
            )


            # ----------------------------------------------
            # Protein FASTA
            # ----------------------------------------------

            fasta_out.write(
                f">{dbcan_gene_id}\n"
                f"{record['sequence']}\n"
            )

            total_proteins += 1


            # ----------------------------------------------
            # GFF3
            # ----------------------------------------------

            contig_id = (
                f"{final_mag}|"
                f"{record['contig']}"
            )

            start = record[
                "start"
            ]

            end = record[
                "end"
            ]

            strand = record[
                "strand"
            ]


            if (
                start is not None
                and
                end is not None
            ):

                genes_with_coordinates += 1

                # Make sure GFF coordinates are ascending.
                gff_start = min(
                    start,
                    end
                )

                gff_end = max(
                    start,
                    end
                )

                attributes = (
                    f"ID={dbcan_gene_id};"
                    f"Name={dbcan_gene_id};"
                    f"Final_MAG={final_mag};"
                    f"Original_protein={protein_id}"
                )

                gff_out.write(
                    f"{contig_id}\t"
                    f"Prodigal\t"
                    f"CDS\t"
                    f"{gff_start}\t"
                    f"{gff_end}\t"
                    f".\t"
                    f"{strand}\t"
                    f".\t"
                    f"{attributes}\n"
                )


            # ----------------------------------------------
            # Protein manifest
            # ----------------------------------------------

            protein_manifest_rows.append({
                "Final_MAG":
                    final_mag,

                "Original_bin":
                    original_bin,

                "protein_id":
                    protein_id,

                "dbCAN_gene_id":
                    dbcan_gene_id,

                "contig":
                    record[
                        "contig"
                    ],

                "dbCAN_contig_id":
                    contig_id,

                "gene_number":
                    record[
                        "gene_number"
                    ],

                "start":
                    (
                        start
                        if start is not None
                        else ""
                    ),

                "end":
                    (
                        end
                        if end is not None
                        else ""
                    ),

                "strand":
                    strand,

                "protein_length_aa":
                    len(
                        record[
                            "sequence"
                        ]
                    ),

                "protein_fasta":
                    str(
                        protein_fasta
                    ),

                "Site":
                    metadata[
                        "Site"
                    ],

                "Treatment":
                    metadata[
                        "Treatment"
                    ],

                "Compartment":
                    metadata[
                        "Compartment"
                    ],

                "GTDB_phylum":
                    metadata[
                        "GTDB_phylum"
                    ],

                "GTDB_order":
                    metadata[
                        "GTDB_order"
                    ],

                "Completeness":
                    metadata[
                        "Completeness"
                    ],

                "Contamination":
                    metadata[
                        "Contamination"
                    ],
            })


        mag_manifest_rows.append({
            "Final_MAG":
                final_mag,

            "Original_bin":
                original_bin,

            "Protein_FASTA":
                str(
                    protein_fasta
                ),

            "Proteins":
                len(
                    records
                ),

            "Status":
                "ok",
        })


# ============================================================
# WRITE PROTEIN MANIFEST
# ============================================================

protein_fields = [
    "Final_MAG",
    "Original_bin",
    "protein_id",
    "dbCAN_gene_id",
    "contig",
    "dbCAN_contig_id",
    "gene_number",
    "start",
    "end",
    "strand",
    "protein_length_aa",
    "protein_fasta",
    "Site",
    "Treatment",
    "Compartment",
    "GTDB_phylum",
    "GTDB_order",
    "Completeness",
    "Contamination",
]


with OUTPUT_MANIFEST.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=protein_fields,
        delimiter="\t",
    )

    writer.writeheader()

    writer.writerows(
        protein_manifest_rows
    )


# ============================================================
# WRITE MAG SOURCE MANIFEST
# ============================================================

mag_fields = [
    "Final_MAG",
    "Original_bin",
    "Protein_FASTA",
    "Proteins",
    "Status",
]


with OUTPUT_MAG_MANIFEST.open(
    "w",
    newline=""
) as outfile:

    writer = csv.DictWriter(
        outfile,
        fieldnames=mag_fields,
        delimiter="\t",
    )

    writer.writeheader()

    writer.writerows(
        mag_manifest_rows
    )


# ============================================================
# PRINT QC
# ============================================================

print()
print(
    "=" * 72
)

print(
    "FINAL 373-MAG PROTEIN + GFF PREPARATION"
)

print(
    "=" * 72
)

print()

print(
    f"Final MAGs expected: "
    f"{len(mag_info)}"
)

print(
    f"MAGs with proteins written: "
    f"{mags_written}"
)

print(
    f"Total proteins written: "
    f"{total_proteins}"
)

print(
    f"Genes with GFF coordinates: "
    f"{genes_with_coordinates}"
)

print(
    f"Unique dbCAN gene IDs: "
    f"{len(seen_dbcan_ids)}"
)

print(
    f"Duplicate dbCAN IDs skipped: "
    f"{len(duplicate_dbcan_ids)}"
)

print(
    f"MAGs missing protein FASTA: "
    f"{len(missing_fasta_mags)}"
)


if missing_fasta_mags:

    print()

    print(
        "Missing protein FASTAs:"
    )

    for mag in missing_fasta_mags:

        print(
            f"  {mag}"
        )


if duplicate_dbcan_ids:

    print()

    print(
        "Duplicate IDs:"
    )

    for gene_id in sorted(
        duplicate_dbcan_ids
    ):

        print(
            f"  {gene_id}"
        )


print()
print(
    "Outputs:"
)

print(
    f"  Combined proteins: "
    f"{OUTPUT_FASTA}"
)

print(
    f"  Combined GFF3:     "
    f"{OUTPUT_GFF}"
)

print(
    f"  Protein manifest:  "
    f"{OUTPUT_MANIFEST}"
)

print(
    f"  MAG source table:  "
    f"{OUTPUT_MAG_MANIFEST}"
)


print()
print(
    "Next:"
)

print(
    "  Run dbCAN CAZyme annotation on the combined protein FASTA."
)

print(
    "  The matching GFF3 will later allow CGC detection while "
    "preserving MAG and contig boundaries."
)

print(
    "=" * 72
)
