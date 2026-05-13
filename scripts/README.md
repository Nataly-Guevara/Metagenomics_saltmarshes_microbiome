# Metagenomics analysis scripts

This folder contains scripts used for shotgun metagenomic preprocessing and downstream analysis of saltmarsh sediment samples.

## Folder structure

- `qc/`: FastQC and MultiQC quality-control scripts.
- `trimming/`: read trimming and filtering scripts using fastp.
- `merging/`: scripts for merging cleaned lane-level reads by biological sample.
- `host_removal/`: host or contaminant read removal scripts.
- `taxonomy/`: taxonomic profiling scripts.
- `functional_annotation/`: functional profiling and annotation scripts.
- `assembly_binning/`: assembly, binning, and MAG-related scripts.

Raw sequencing data are not included in this repository.
