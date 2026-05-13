
# Metagenomics Saltmarsh Microbiome

This repository contains the computational workflow for shotgun metagenomic analysis of saltmarsh sediment microbial communities.

The project includes scripts and documentation for quality control, read trimming, lane merging, host or contaminant read removal, taxonomic profiling, functional annotation, assembly, binning, and summary reporting.

Raw sequencing data and large intermediate files are not stored in this GitHub repository.

## Repository structure

```text
metagenomics_project/
├── 00_raw_data/
├── 01_metadata/
├── 02_qc/
│   ├── fastqc_raw/
│   ├── multiqc_raw/
│   ├── trimmed_reads/
│   ├── fastqc_trimmed/
│   └── multiqc_trimmed/
├── 03_host_removal/
├── 04_taxonomy/
├── 05_function/
├── 06_assembly/
├── 07_bins/
├── 08_results/
│   ├── tables/
│   ├── figures/
│   └── summary_reports/
├── 09_logs/
├── 10_scripts/
│   ├── bash/
│   ├── pbs/
│   ├── R/
│   └── python/
└── README.md

**Data policy**

The repository excludes large sequencing and intermediate files, including FASTQ files, BAM/SAM files, large compressed outputs, and large R objects. These files remain on the HPC system or institutional storage.

**Project status**

This repository is under active development. Scripts and documentation will be added as the metagenomics workflow progresses.


