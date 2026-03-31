#!/bin/bash
# 1_Fastp_BWA_opt.sh
# Usage: bash 1_Fastp_BWA_opt.sh -r ./rawdata -g /path/to/hg19.fa -t 16

set -euo pipefail
shopt -s nullglob

# --- 1. Default settings ---
RAW_DIR="./rawdata"
CLEAN_DIR="./clean_fastq"
BAM_DIR="./bam_files"
GENOME_INDEX="/public/users/pghzeng/Database/AA_data_repo/hg19/hg19full.fa"
THREADS=16

TOOLS_ENV=""

run_cmd() {
    # Usage: run_cmd <envNameOrEmpty> <command> [args...]
    local env_name="$1"
    shift
    local exe="$1"

    if command -v "$exe" >/dev/null 2>&1; then
        "$@"
        return
    fi

    if [ -z "$env_name" ]; then
        echo "Error: '$exe' not found in PATH and no env specified to run it." >&2
        exit 127
    fi

    if command -v mamba >/dev/null 2>&1; then
        mamba run -n "$env_name" "$@"
        return
    fi
    if command -v conda >/dev/null 2>&1; then
        conda run -n "$env_name" "$@"
        return
    fi

    echo "Error: '$exe' not found in PATH; and neither mamba nor conda is available to run env '$env_name'." >&2
    exit 127
}

# --- 2. Help function ---
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -r <dir>    Raw data directory [Default: $RAW_DIR]"
    echo "  -c <dir>    Clean FASTQ output directory [Default: $CLEAN_DIR]"
    echo "  -b <dir>    BAM output directory [Default: $BAM_DIR]"
    echo "  -g <file>   BWA genome index path [Default: $GENOME_INDEX]"
    echo "  -t <int>    Number of threads [Default: $THREADS]"
    echo "  -e <name>   conda/mamba environment name to use when fastp/bwa/samtools are not in PATH [Default: empty]"
    echo "  -h          Show this help message"
    exit 1
}

# --- 3. Parse command-line arguments ---
while getopts "r:c:b:g:t:e:h" opt; do
    case $opt in
        r) RAW_DIR="$OPTARG" ;;
        c) CLEAN_DIR="$OPTARG" ;;
        b) BAM_DIR="$OPTARG" ;;
        g) GENOME_INDEX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        e) TOOLS_ENV="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# --- 4. Validate required inputs ---
if [ ! -d "$RAW_DIR" ]; then
    echo "Error: Raw directory '$RAW_DIR' not found."
    exit 1
fi

# Check the BWA index path.
# This script currently prints the configured path instead of validating all index sidecar files.
echo "========================================"
echo "Configuration:"
echo "Raw Dir:      $RAW_DIR"
echo "Clean Dir:    $CLEAN_DIR"
echo "BAM Dir:      $BAM_DIR"
echo "Genome Index: $GENOME_INDEX"
echo "Threads:      $THREADS"
echo "Tools Env:    ${TOOLS_ENV:-<use PATH>}"
echo "========================================"

mkdir -p "$CLEAN_DIR" "$BAM_DIR"

# --- 5. Main loop ---
# This assumes input files use the exact suffix .1.fq.gz.
# Update the glob if your FASTQ naming convention differs.
r1_files=("${RAW_DIR}"/*.1.fq.gz)
if [ "${#r1_files[@]}" -eq 0 ]; then
    echo "Error: No *.1.fq.gz files found in $RAW_DIR"
    exit 1
fi

for r1 in "${r1_files[@]}"; do
    base=$(basename "$r1" .1.fq.gz)
    r2="${RAW_DIR}/${base}.2.fq.gz"

    # Check whether the matching R2 file exists.
    if [ ! -f "$r2" ]; then
        echo "Warning: Paired file $r2 not found, skipping $base."
        continue
    fi

    echo ">>> [$(date '+%T')] Processing ${base}..."

    clean_r1="${CLEAN_DIR}/${base}_R1.clean.fq.gz"
    clean_r2="${CLEAN_DIR}/${base}_R2.clean.fq.gz"
    fastp_html="${CLEAN_DIR}/${base}.html"
    fastp_json="${CLEAN_DIR}/${base}.json"
    fastp_log="${CLEAN_DIR}/${base}.fastp.log"

    out_bam="${BAM_DIR}/${base}.sorted.bam"
    out_bai="${BAM_DIR}/${base}.sorted.bam.bai"

    # fastp with resume-friendly behavior: skip if clean FASTQs already exist.
    if [ -s "$clean_r1" ] && [ -s "$clean_r2" ]; then
        echo "    -> fastp outputs exist, skipping fastp."
    else
        echo "    -> Running fastp..."
        run_cmd "$TOOLS_ENV" fastp -i "$r1" -I "$r2" \
            -o "$clean_r1" -O "$clean_r2" \
            --detect_adapter_for_pe --fix_mgi_id --umi --umi_loc=per_read \
            --umi_len=8 --umi_prefix=UMI --trim_front1 8 --trim_front2 8 \
            -w "$THREADS" -h "$fastp_html" -j "$fastp_json" &> "$fastp_log"
    fi

    # BWA + samtools with resume-friendly behavior: skip if BAM and BAI already exist.
    if [ -s "$out_bam" ] && [ -s "$out_bai" ]; then
        echo "    -> BAM exists, skipping BWA-MEM & sorting."
    else
        if [ ! -s "$clean_r1" ] || [ ! -s "$clean_r2" ]; then
            echo "Error: clean fastq not found for ${base}. Expected: $clean_r1 and $clean_r2"
            exit 1
        fi

        echo "    -> Running BWA-MEM & Sorting..."
        run_cmd "$TOOLS_ENV" bwa mem -5SP -t "$THREADS" "$GENOME_INDEX" \
            "$clean_r1" \
            "$clean_r2" | \
            run_cmd "$TOOLS_ENV" samtools view -bhS - | run_cmd "$TOOLS_ENV" samtools sort -@ "$THREADS" -o "$out_bam"
    fi

    # If BAM exists but the index is missing, generate only the index.
    if [ -s "$out_bam" ] && [ ! -s "$out_bai" ]; then
        echo "    -> Indexing BAM..."
        run_cmd "$TOOLS_ENV" samtools index "$out_bam"
    fi
    echo ">>> ${base} Done."
done