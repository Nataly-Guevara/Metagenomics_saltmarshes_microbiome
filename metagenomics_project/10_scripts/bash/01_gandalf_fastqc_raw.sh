#!/bin/bash

module load fastqc/0.11.9

INPUT_DIR="/path/to/metagenomics_project/data/raw_metagenomic_reads"
OUTPUT_DIR="/path/to/metagenomics_project/qc/fastqc_raw"

mkdir -p "$OUTPUT_DIR"

fastqc "$INPUT_DIR"/*.fastq.gz -o "$OUTPUT_DIR"
