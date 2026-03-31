#!/bin/bash
# 2_Add_Tags.sh

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

BAM_DIR="./bam_files"
TAGGED_DIR="./tagged_bam"
ADD_UMI_PY="${SCRIPT_DIR}/add_umi_tag.py"
THREADS=4

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

usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -i <dir>   Input BAM directory containing *.sorted.bam from 1_Fastp_BWA.sh [Default: $BAM_DIR]"
    echo "  -o <dir>   Output directory for tagged BAM files [Default: $TAGGED_DIR]"
    echo "  -p <file>  Path to add_umi_tag.py [Default: $ADD_UMI_PY]"
    echo "  -t <int>   Number of threads for samtools index [Default: $THREADS]"
    echo "  -e <name>  conda/mamba environment name to use when python/samtools are not in PATH [Default: empty]"
    echo "  -h         Show this help message"
    exit 1
}

while getopts "i:o:p:t:e:h" opt; do
    case $opt in
        i) BAM_DIR="$OPTARG" ;;
        o) TAGGED_DIR="$OPTARG" ;;
        p) ADD_UMI_PY="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        e) TOOLS_ENV="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

if [ ! -d "$BAM_DIR" ]; then
    echo "Error: BAM directory '$BAM_DIR' not found."
    exit 1
fi
if [ ! -f "$ADD_UMI_PY" ]; then
    echo "Error: add_umi_tag.py not found at '$ADD_UMI_PY'"
    exit 1
fi

mkdir -p "$TAGGED_DIR"

bams=("${BAM_DIR}"/*.sorted.bam)
if [ "${#bams[@]}" -eq 0 ]; then
    echo "Error: No *.sorted.bam files found in $BAM_DIR"
    exit 1
fi

for bam in "${bams[@]}"; do
    base=$(basename "$bam" .sorted.bam)
    out_bam="${TAGGED_DIR}/${base}.tagged.bam"
    out_bai="${TAGGED_DIR}/${base}.tagged.bam.bai"

    # Resume-friendly behavior: skip if both tagged BAM and BAI already exist.
    if [ -s "$out_bam" ] && [ -s "$out_bai" ]; then
        echo ">>> ${base}: tagged BAM exists, skipping."
        continue
    fi

    if [ -s "$out_bam" ] && [ ! -s "$out_bai" ]; then
        echo ">>> ${base}: BAM exists but index missing, indexing..."
        run_cmd "$TOOLS_ENV" samtools index -@ "$THREADS" "$out_bam"
        continue
    fi

    echo ">>> Injecting RX tags for ${base}..."
    run_cmd "$TOOLS_ENV" python "$ADD_UMI_PY" "$bam" "$out_bam"
    run_cmd "$TOOLS_ENV" samtools index -@ "$THREADS" "$out_bam"
done
