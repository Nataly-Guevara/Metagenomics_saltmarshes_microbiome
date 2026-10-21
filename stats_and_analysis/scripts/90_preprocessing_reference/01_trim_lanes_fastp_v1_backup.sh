#!/bin/bash

RAW_DIR="data/raw_metagenomic_reads"
OUT_DIR="data/trimmed_lanes"
REPORT_DIR="qc/fastp_reports"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$REPORT_DIR" "$LOG_DIR"

for R1 in "$RAW_DIR"/*_R1.fastq.gz
do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

    BASENAME=$(basename "$R1" _R1.fastq.gz)

    OUT_R1="$OUT_DIR/${BASENAME}_trimmed_R1.fastq.gz"
    OUT_R2="$OUT_DIR/${BASENAME}_trimmed_R2.fastq.gz"

    HTML="$REPORT_DIR/${BASENAME}_fastp.html"
    JSON="$REPORT_DIR/${BASENAME}_fastp.json"

    echo "Processing $BASENAME"

    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "$OUT_R1" \
      -O "$OUT_R2" \
      --detect_adapter_for_pe \
      --trim_poly_g \
      --cut_front \
      --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40 \
      --length_required 50 \
      --thread 8 \
      --html "$HTML" \
      --json "$JSON" \
      2> "$LOG_DIR/${BASENAME}_fastp.log"

done
