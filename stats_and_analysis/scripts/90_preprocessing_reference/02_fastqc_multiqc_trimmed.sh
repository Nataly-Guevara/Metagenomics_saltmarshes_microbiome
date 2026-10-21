#!/bin/bash

# Run FastQC and MultiQC on fastp-trimmed lane-level FASTQ files

TRIMMED_DIR="data/trimmed_lanes_v2"
FASTQC_OUT="qc/fastqc_trimmed_v2"
MULTIQC_OUT="qc/multiqc_trimmed_v2"

mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

fastqc "$TRIMMED_DIR"/*.fastq.gz \
  --outdir "$FASTQC_OUT" \
  --threads 8

multiqc "$FASTQC_OUT" qc/fastp_reports_v2 \
  --outdir "$MULTIQC_OUT" \
  --filename multiqc_trimmed_v2_report.html
