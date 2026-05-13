cat > README.md << 'EOF'
# Metagenomics Saltmarsh Microbiome

This repository contains the computational workflow for shotgun metagenomic analysis of saltmarsh sediment microbial communities.

Raw sequencing data and large intermediate files are not stored in this GitHub repository. This repository is intended to store scripts, small metadata files, selected QC reports, summary tables, figures, and documentation.

## Repository structure

```text
metagenomics_project/
├── 00_raw_data/
├── 01_metadata/
├── 02_qc/
│   ├── fastqc_raw/
│   ├── multiqc_raw/
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
```

## Folder descriptions

- `00_raw_data/`: placeholder only. Raw FASTQ files remain on HPC storage and are not uploaded to GitHub.
- `01_metadata/`: sample metadata, experimental design files, and soil parameter tables, if small and non-sensitive.
- `02_qc/`: selected quality-control reports from FastQC and MultiQC.
- `03_host_removal/`: scripts and small summary outputs related to host or contaminant read removal.
- `04_taxonomy/`: scripts and small taxonomic profiling summaries.
- `05_function/`: scripts and small functional annotation summaries.
- `06_assembly/`: scripts and small assembly statistics. Large contig files are not uploaded.
- `07_bins/`: scripts and small binning/MAG quality summaries. MAG FASTA files are not uploaded.
- `08_results/`: selected final tables, figures, and summary reports.
- `09_logs/`: selected small job logs or run summaries.
- `10_scripts/`: bash, PBS, R, and Python scripts used in the workflow.

## Data policy

The following files are excluded from this repository:

- raw sequencing reads
- trimmed reads
- host-filtered reads
- large assembly files
- large bin/MAG files
- large compressed outputs
- large R objects

These files should remain on Gadi, Gandalf, or institutional storage.

## Project status

This repository is under active development. Scripts and documentation will be added as the metagenomics workflow progresses.
EOF
