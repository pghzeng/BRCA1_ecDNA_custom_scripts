#!/bin/bash
# 3_Pairtools_Process.sh

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

TAGGED_DIR="./tagged_bam"
PAIRS_DIR="./pairs_output"
CHROMS_SIZE="/data0/database/hg19_zpgh/Genome/Genome.sizes"

THREADS=12
MIN_MAPQ=30
WALKS_POLICY="5unique"
CHROM_SUBSET=""

# Force a working compressor so that *.gz outputs are truly gzip/BGZF.
# Some installations may not have bgzip/pbgzip available, causing pairtools to
# silently write plain text even when the filename ends with .gz.
PAIRTOOLS_CMD_OUT="gzip -c"
PAIRTOOLS_CMD_IN=""

PAIRTOOLS_ENV=""
SAMTOOLS_ENV=""

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

gzip_ok() {
    # Return 0 if a .gz file exists and passes gzip integrity test.
    # For BGZF files, `gzip -t` is still expected to work.
    local path="$1"
    if [ ! -s "$path" ]; then
        return 1
    fi
    if [[ "$path" != *.gz ]]; then
        return 0
    fi
    if command -v gzip >/dev/null 2>&1; then
        gzip -t "$path" >/dev/null 2>&1
        return $?
    fi
    # Fallback: best-effort read test of the first few lines
    zcat "$path" >/dev/null 2>&1
}

usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -i <dir>   Input tagged BAM directory containing *.tagged.bam from 2_Add_Tags.sh [Default: $TAGGED_DIR]"
    echo "  -o <dir>   Output directory for pairs files and final BAM files [Default: $PAIRS_DIR]"
    echo "  -s <file>  Path to chrom.sizes file passed to pairtools --chroms-path [Default: $CHROMS_SIZE]"
    echo "  -t <int>   Number of threads for pairtools sort and samtools [Default: $THREADS]"
    echo "  -q <int>   Minimum MAPQ for pairtools parse [Default: $MIN_MAPQ]"
    echo "  -w <str>   Value for pairtools parse --walks-policy [Default: $WALKS_POLICY]"
    echo "  -u <str>   Optional chromosome subset passed to pairtools parse --chrom-subset [Default: empty]"
    echo "  -Z <str>   Output compression command passed to pairtools --cmd-out, e.g. 'bgzip -c' or 'gzip -c' [Default: $PAIRTOOLS_CMD_OUT]"
    echo "  -z <str>   Input decompression command passed to pairtools --cmd-in, e.g. 'bgzip -dc' or 'gzip -dc' [Default: auto]"
    echo "  -e <name>  conda/mamba environment name for pairtools when pairtools is not in PATH [Default: empty]"
    echo "  -E <name>  conda/mamba environment name for samtools when samtools is not in PATH [Default: empty]"
    echo "  -h         Show this help message"
    exit 1
}

while getopts "i:o:s:t:q:w:u:Z:z:e:E:h" opt; do
    case $opt in
        i) TAGGED_DIR="$OPTARG" ;;
        o) PAIRS_DIR="$OPTARG" ;;
        s) CHROMS_SIZE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MIN_MAPQ="$OPTARG" ;;
        w) WALKS_POLICY="$OPTARG" ;;
        u) CHROM_SUBSET="$OPTARG" ;;
        Z) PAIRTOOLS_CMD_OUT="$OPTARG" ;;
        z) PAIRTOOLS_CMD_IN="$OPTARG" ;;
        e) PAIRTOOLS_ENV="$OPTARG" ;;
        E) SAMTOOLS_ENV="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Auto-select cmd-in based on cmd-out when it is not provided.
if [ -z "$PAIRTOOLS_CMD_IN" ]; then
    case "$PAIRTOOLS_CMD_OUT" in
        *pbgzip*) PAIRTOOLS_CMD_IN="pbgzip -dc" ;;
        *bgzip*)  PAIRTOOLS_CMD_IN="bgzip -dc" ;;
        *gzip*)   PAIRTOOLS_CMD_IN="gzip -dc" ;;
        *)        PAIRTOOLS_CMD_IN="" ;;
    esac
fi

if [ ! -d "$TAGGED_DIR" ]; then
    echo "Error: tagged BAM directory '$TAGGED_DIR' not found."
    exit 1
fi
if [ -z "$CHROMS_SIZE" ] || [ ! -f "$CHROMS_SIZE" ]; then
    echo "Error: chrom.sizes file not found: '$CHROMS_SIZE' (use -s to override)."
    exit 1
fi

mkdir -p "$PAIRS_DIR"

bams=("${TAGGED_DIR}"/*.tagged.bam)
if [ "${#bams[@]}" -eq 0 ]; then
    echo "Error: No *.tagged.bam files found in $TAGGED_DIR"
    exit 1
fi

for bam in "${bams[@]}"; do
    base=$(basename "$bam" .tagged.bam)
    echo ">>> Processing pairs for ${base}..."

    qname_bam="${PAIRS_DIR}/${base}.qname_sorted.bam"
    raw_pairs="${PAIRS_DIR}/${base}.raw.pairs.gz"
    sorted_pairs="${PAIRS_DIR}/${base}.sorted.pairs.gz"
    dedup_pairs="${PAIRS_DIR}/${base}.deduped.pairs.gz"
    dedup_stats="${PAIRS_DIR}/${base}.dedup_stats.txt"
    tmp_sam="${PAIRS_DIR}/${base}.tmp.sam"
    final_bam="${PAIRS_DIR}/${base}.final.bam"
    final_bai="${PAIRS_DIR}/${base}.final.bam.bai"

    # Resume-friendly behavior: skip the whole sample if the final BAM already exists.
    if [ -s "$final_bam" ] && [ -s "$final_bai" ]; then
        echo "    -> final BAM exists, skipping ${base}."
        continue
    fi

    # Step 1: parse. Skip if raw pairs already exist and pass integrity checks.
    if gzip_ok "$raw_pairs"; then
        echo "    -> raw pairs exist, skipping parse."
    else
        if [ -e "$raw_pairs" ]; then
            echo "    -> raw pairs exists but seems incomplete/corrupt, re-generating..."
            rm -f "$raw_pairs"
        fi
        if [ ! -s "$qname_bam" ] || [ "$bam" -nt "$qname_bam" ]; then
            echo "    -> Name-sorting tagged BAM for pairtools parse..."
            qname_tmp="${qname_bam}.tmp.$$"
            rm -f "$qname_tmp"
            run_cmd "$SAMTOOLS_ENV" samtools sort -n -@ "$THREADS" -o "$qname_tmp" "$bam"
            mv -f "$qname_tmp" "$qname_bam"
        else
            echo "    -> qname-sorted BAM exists, skipping name-sort."
        fi
        echo "    -> Running pairtools parse..."
        raw_tmp="${raw_pairs}.tmp.$$"
        parse_args=(
            parse
            --walks-policy "$WALKS_POLICY"
            --min-mapq "$MIN_MAPQ"
            --add-columns RX
            --chroms-path "$CHROMS_SIZE"
            --cmd-out "$PAIRTOOLS_CMD_OUT"
        )
        if [ -n "$CHROM_SUBSET" ]; then
            parse_args+=(--chrom-subset "$CHROM_SUBSET")
        fi
        parse_args+=(-o "$raw_tmp" "$qname_bam")
        rm -f "$raw_tmp"
        run_cmd "$PAIRTOOLS_ENV" pairtools "${parse_args[@]}"
        if ! gzip_ok "$raw_tmp"; then
            echo "Error: pairtools parse produced an invalid output: $raw_tmp" >&2
            rm -f "$raw_tmp"
            exit 1
        fi
        mv -f "$raw_tmp" "$raw_pairs"
        rm -f "$qname_bam"
    fi

    # Step 2: sort. Skip if sorted pairs already exist and pass integrity checks.
    if gzip_ok "$sorted_pairs"; then
        echo "    -> sorted pairs exist, skipping sort."
    else
        if [ -e "$sorted_pairs" ]; then
            echo "    -> sorted pairs exists but seems incomplete/corrupt, re-generating..."
            rm -f "$sorted_pairs"
        fi
        if [ ! -s "$raw_pairs" ]; then
            echo "Error: raw pairs not found for ${base}: $raw_pairs"
            exit 1
        fi
        echo "    -> Running pairtools sort..."
        sorted_tmp="${sorted_pairs}.tmp.$$"
        rm -f "$sorted_tmp"
        sort_args=(sort --nproc "$THREADS" "$raw_pairs" -o "$sorted_tmp" --cmd-out "$PAIRTOOLS_CMD_OUT")
        if [[ "$raw_pairs" == *.gz ]] && [ -n "$PAIRTOOLS_CMD_IN" ]; then
            sort_args+=(--cmd-in "$PAIRTOOLS_CMD_IN")
        fi
        run_cmd "$PAIRTOOLS_ENV" pairtools "${sort_args[@]}"
        if ! gzip_ok "$sorted_tmp"; then
            echo "Error: pairtools sort produced an invalid output: $sorted_tmp" >&2
            rm -f "$sorted_tmp"
            exit 1
        fi
        mv -f "$sorted_tmp" "$sorted_pairs"
    fi

    # Step 3: dedup. Skip if dedup pairs and stats already exist and are valid.
    if gzip_ok "$dedup_pairs" && [ -s "$dedup_stats" ]; then
        echo "    -> dedup outputs exist, skipping dedup."
    else
        if [ -e "$dedup_pairs" ] || [ -e "$dedup_stats" ]; then
            echo "    -> dedup outputs exist but seem incomplete/corrupt, re-generating..."
            rm -f "$dedup_pairs" "$dedup_stats"
        fi
        if [ ! -s "$sorted_pairs" ]; then
            echo "Error: sorted pairs not found for ${base}: $sorted_pairs"
            exit 1
        fi
        echo "    -> Running pairtools dedup (UMI-aware)..."
        # pairtools parse --add-columns RX creates RX1/RX2, so dedup must use that column pair.
        dedup_tmp="${dedup_pairs}.tmp.$$"
        stats_tmp="${dedup_stats}.tmp.$$"
        rm -f "$dedup_tmp" "$stats_tmp"
        dedup_args=(dedup --extra-col-pair RX1 RX2 --mark-dups --cmd-out "$PAIRTOOLS_CMD_OUT" --output-stats "$stats_tmp" -o "$dedup_tmp" "$sorted_pairs")
        if [[ "$sorted_pairs" == *.gz ]] && [ -n "$PAIRTOOLS_CMD_IN" ]; then
            dedup_args+=(--cmd-in "$PAIRTOOLS_CMD_IN")
        fi
        run_cmd "$PAIRTOOLS_ENV" pairtools "${dedup_args[@]}"
        if ! gzip_ok "$dedup_tmp" || [ ! -s "$stats_tmp" ]; then
            echo "Error: pairtools dedup produced invalid output for ${base}." >&2
            rm -f "$dedup_tmp" "$stats_tmp"
            exit 1
        fi
        mv -f "$dedup_tmp" "$dedup_pairs"
        mv -f "$stats_tmp" "$dedup_stats"
    fi

    # Step 4: split -> BAM. If final BAM already exists, skip; otherwise generate it from dedup pairs.
    if [ -s "$final_bam" ] && [ ! -s "$final_bai" ]; then
        echo "    -> final BAM exists but index missing, indexing..."
        run_cmd "$SAMTOOLS_ENV" samtools index -@ "$THREADS" "$final_bam"
        continue
    fi

    if [ -s "$final_bam" ] && [ -s "$final_bai" ]; then
        echo "    -> final BAM exists, skipping split/BAM generation."
        continue
    fi

    if [ ! -s "$dedup_pairs" ]; then
        echo "Error: dedup pairs not found for ${base}: $dedup_pairs"
        exit 1
    fi

    echo "    -> Converting deduped pairs back to BAM..."
    rm -f "$tmp_sam"
    split_args=(split --output-sam "$tmp_sam" "$dedup_pairs")
    if [[ "$dedup_pairs" == *.gz ]] && [ -n "$PAIRTOOLS_CMD_IN" ]; then
        split_args+=(--cmd-in "$PAIRTOOLS_CMD_IN")
    fi
    run_cmd "$PAIRTOOLS_ENV" pairtools "${split_args[@]}"
    run_cmd "$SAMTOOLS_ENV" samtools view -bhS "$tmp_sam" | run_cmd "$SAMTOOLS_ENV" samtools sort -@ "$THREADS" -o "$final_bam"
    run_cmd "$SAMTOOLS_ENV" samtools index -@ "$THREADS" "$final_bam"
    rm -f "$tmp_sam"
done