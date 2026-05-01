#!/usr/bin/env bash
# joint_coverage_bed.sh
#
# Compute regions where coverage in BOTH BAMs exceeds a threshold.
#
# Strategy:
#   1. Run mosdepth on each BAM with --quantize so it emits a BED of
#      "high coverage" intervals (>= MIN_COV) directly. This is fast
#      because mosdepth streams the BAM once and never writes per-base depth.
#   2. bedtools intersect the two high-coverage BEDs to get the joint set.
#
# Requires: mosdepth, bedtools, bgzip, tabix.

set -euo pipefail

usage() {
    cat <<EOF
Usage: $0 -a BAM1 -b BAM2 -o OUT_PREFIX [-c MIN_COV] [-t THREADS] [-r REGIONS_BED]

  -a BAM1         First BAM file (indexed)
  -b BAM2         Second BAM file (indexed)
  -o OUT_PREFIX   Output prefix (e.g., results/joint)
  -c MIN_COV      Minimum coverage threshold (default: 10)
  -t THREADS      Threads for mosdepth (default: 4)
  -r REGIONS_BED  Optional BED to restrict analysis (e.g., autosomes only)
  -m MAPQ         Minimum mapping quality (default: 20)

Output:
  OUT_PREFIX.bam1.highcov.bed.gz   intervals where BAM1 coverage >= MIN_COV
  OUT_PREFIX.bam2.highcov.bed.gz   intervals where BAM2 coverage >= MIN_COV
  OUT_PREFIX.joint.highcov.bed     intersection (intervals high-cov in BOTH)
  OUT_PREFIX.joint.highcov.bed.gz  bgzipped + tabix-indexed
EOF
    exit 1
}

MIN_COV=10
THREADS=4
MAPQ=20
REGIONS=""
BAM1=""
BAM2=""
OUT=""

while getopts "a:b:o:c:t:r:m:h" opt; do
    case $opt in
        a) BAM1="$OPTARG" ;;
        b) BAM2="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        c) MIN_COV="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        r) REGIONS="$OPTARG" ;;
        m) MAPQ="$OPTARG" ;;
        h|*) usage ;;
    esac
done

[[ -z "$BAM1" || -z "$BAM2" || -z "$OUT" ]] && usage

# Sanity checks
for f in "$BAM1" "$BAM2"; do
    [[ -f "$f" ]] || { echo "ERROR: BAM not found: $f" >&2; exit 1; }
    [[ -f "$f.bai" || -f "${f%.bam}.bai" ]] || {
        echo "ERROR: BAM index missing for $f. Run: samtools index $f" >&2
        exit 1
    }
done

mkdir -p "$(dirname "$OUT")"

# mosdepth --quantize cuts coverage into bins. We use two bins:
#   "LOW"  for [0, MIN_COV)
#   "HIGH" for [MIN_COV, infinity)
# The output .quantized.bed.gz has every interval labeled, and we just
# grep for HIGH.
QUANTIZE="0:${MIN_COV}:LOW ${MIN_COV}:HIGH"
# mosdepth quantize syntax: thresholds separated by colons; we set
# breakpoints at 0 and MIN_COV with labels via MOSDEPTH_Q0/Q1 env vars.
export MOSDEPTH_Q0=NO_COVERAGE
export MOSDEPTH_Q1=LOW_COVERAGE
export MOSDEPTH_Q2=HIGH_COVERAGE

REGIONS_ARG=()
if [[ -n "$REGIONS" ]]; then
    [[ -f "$REGIONS" ]] || { echo "ERROR: regions BED not found: $REGIONS" >&2; exit 1; }
    REGIONS_ARG=(--by "$REGIONS")
fi

run_mosdepth() {
    local bam="$1"
    local prefix="$2"
    echo "[$(date +%H:%M:%S)] mosdepth on $bam -> $prefix" >&2
    mosdepth \
        --threads "$THREADS" \
        --no-per-base \
        --mapq "$MAPQ" \
        --quantize "0:1:${MIN_COV}:" \
        "${REGIONS_ARG[@]}" \
        "$prefix" \
        "$bam"
}

# Run mosdepth on each BAM
run_mosdepth "$BAM1" "${OUT}.bam1"
run_mosdepth "$BAM2" "${OUT}.bam2"

# Extract HIGH_COVERAGE intervals from each .quantized.bed.gz
extract_high() {
    local prefix="$1"
    local out_bed="$2"
    zcat "${prefix}.quantized.bed.gz" \
        | awk -v OFS='\t' '$4 == "HIGH_COVERAGE" {print $1, $2, $3}' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        | bgzip -c > "$out_bed"
    tabix -f -p bed "$out_bed"
}

echo "[$(date +%H:%M:%S)] Extracting HIGH_COVERAGE intervals" >&2
extract_high "${OUT}.bam1" "${OUT}.bam1.highcov.bed.gz"
extract_high "${OUT}.bam2" "${OUT}.bam2.highcov.bed.gz"

# Intersect
echo "[$(date +%H:%M:%S)] Computing joint high-coverage intervals" >&2
bedtools intersect \
    -a "${OUT}.bam1.highcov.bed.gz" \
    -b "${OUT}.bam2.highcov.bed.gz" \
    -sorted \
    > "${OUT}.joint.highcov.bed"

bgzip -kf "${OUT}.joint.highcov.bed"
tabix -f -p bed "${OUT}.joint.highcov.bed.gz"

# Summary
N_INTERVALS=$(wc -l < "${OUT}.joint.highcov.bed")
TOTAL_BP=$(awk '{s+=$3-$2} END {print s+0}' "${OUT}.joint.highcov.bed")

echo ""
echo "=== Joint high-coverage summary (>= ${MIN_COV}x in BOTH BAMs) ==="
echo "Intervals: ${N_INTERVALS}"
echo "Total bp:  ${TOTAL_BP}"
echo "Output:    ${OUT}.joint.highcov.bed[.gz]"
