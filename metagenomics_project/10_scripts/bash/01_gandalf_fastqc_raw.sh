#!/bin/bash

module load fastqc/0.11.9

INPUT_DIR="/home/gnat/metagenomics-saltmarshes/data/raw_metagenomic_reads"
OUTPUT_DIR="/home/gnat/metagenomics-saltmarshes/qc/fastqc_raw"

mkdir -p "$OUTPUT_DIR"

fastqc "$INPUT_DIR"/*.fastq.gz -o "$OUTPUT_DIR"
