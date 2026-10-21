#!/bin/bash

# Run FastQC and MultiQC on fastp-trimmed lane-level FASTQ files

TRIMMED_DIR="data/trimmed_lanes"
FASTQC_OUT="qc/fastqc_trimmed"
MULTIQC_OUT="qc/multiqc_trimmed"

mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

fastqc "$TRIMMED_DIR"/*.fastq.gz \
  --outdir "$FASTQC_OUT" \
  --threads 8

multiqc "$FASTQC_OUT" qc/fastp_reports \
  --outdir "$MULTIQC_OUT" \
  --filename multiqc_trimmed_report.html
