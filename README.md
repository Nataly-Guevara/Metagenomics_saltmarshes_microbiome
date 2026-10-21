# Metagenomics Saltmarsh Microbiome

This repository contains the reproducible bioinformatics workflow for shotgun metagenomic analysis of saltmarsh sediment microbial communities.

The repository is designed to store **scripts, documentation, software notes, metadata templates, and small summary outputs**. Raw sequencing reads and large intermediate files are not stored here.

## Repository structure

| Folder | Purpose |
| --- | --- |
| [`metagenomics_project/`](metagenomics_project/) | Existing preprocessing and bioinformatics workflow structure. |
| [`stats_and_analysis/`](stats_and_analysis/) | Chapter 2 statistics and ecological analyses, with scripts, inputs, and outputs organised by analysis. |

The statistical analysis collection is being assembled. See its [README](stats_and_analysis/README.md) for the layout and current status.

## Data availability

Raw metagenomic sequencing data are not stored in this GitHub repository because of file size and reproducibility best practices. Raw reads should be deposited in NCBI and referenced using BioProject, BioSample, and SRA accession numbers in `metadata/NCBI_accessions.md`.

Large intermediate files, including trimmed reads, host-filtered reads, assemblies, BAM files, MAG FASTA files, GTDB-Tk databases, CheckM2 databases, and temporary job outputs, should remain on HPC or institutional storage.

## Workflow overview

The current workflow includes:

1. Quality control and read trimming
2. Host and contaminant read removal
3. Assembly with MEGAHIT
4. Read mapping and depth generation
5. Genome binning with MetaBAT2
6. MAG quality assessment with CheckM2
7. MAG filtering and dereplication with dRep
8. Taxonomic classification with GTDB-Tk

## Notes on reproducibility

PBS scripts are written for Gadi/NCI and may require adjustment for other HPC systems. Software installation notes and database locations should be documented in `envs/software_notes.md`.

## What is not tracked

The repository intentionally excludes:

- raw FASTQ files
- trimmed FASTQ files
- host-filtered FASTQ files
- BAM/SAM/CRAM alignment files
- large assembly FASTA files
- MAG/bin FASTA files
- GTDB-Tk, CheckM2, or other reference databases
- failed run folders and large intermediate outputs
- PBS `.o`, `.e`, `.log`, `.out`, and `.err` files

## Project status

This repository is under active development. Scripts and documentation will be added as the metagenomics workflow progresses.
