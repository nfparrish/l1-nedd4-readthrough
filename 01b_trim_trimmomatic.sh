#!/bin/bash
#SBATCH --job-name=trim_tm
#SBATCH --output=/work/nfp8/2026_03_02/GSE226189/logs/trim_%A_%a.out
#SBATCH --error=/work/nfp8/2026_03_02/GSE226189/logs/trim_%A_%a.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=48G
#SBATCH -c 8
#SBATCH --time=04:00:00
#SBATCH --array=1-82
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"
module load Trimmomatic/0.39

# ---- Resolve sample ----
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -z "$SRR" ]] && { echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: trimming ${SRR}"

mkdir -p "$FASTQ_DIR"

TRIMMOMATIC_JAR="$TRIMMOMATIC"
ADAPTERS_DIR="$(dirname "$TRIMMOMATIC_JAR")/adapters"
THREADS=8

if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"

    # Trimmomatic PE output naming (with .P suffix for paired, .U for unpaired):
    # We keep only paired-end reads and rename to match expected pattern
    TRIM_R1_PAIRED="${FASTQ_DIR}/${SRR}_1_trimmed_P.fq.gz"
    TRIM_R2_PAIRED="${FASTQ_DIR}/${SRR}_2_trimmed_P.fq.gz"
    TRIM_R1_UNPAIRED="${FASTQ_DIR}/${SRR}_1_trimmed_U.fq.gz"
    TRIM_R2_UNPAIRED="${FASTQ_DIR}/${SRR}_2_trimmed_U.fq.gz"

    # Check if already trimmed
    if [[ -f "$TRIM_R1_PAIRED" && -f "$TRIM_R2_PAIRED" ]]; then
        echo "  Already trimmed — skipping"; exit 0
    fi

    [[ -f "$R1" ]] || { echo "ERROR: $R1 not found" >&2; exit 1; }
    [[ -f "$R2" ]] || { echo "ERROR: $R2 not found" >&2; exit 1; }

    echo "  Trimming paired-end reads with Trimmomatic..."
    java -jar "$TRIMMOMATIC_JAR" PE -threads "$THREADS" \
        "$R1" "$R2" \
        "$TRIM_R1_PAIRED" "$TRIM_R1_UNPAIRED" \
        "$TRIM_R2_PAIRED" "$TRIM_R2_UNPAIRED" \
        ILLUMINACLIP:"${ADAPTERS_DIR}/TruSeq3-PE.fa":2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:${TRIM_QUALITY} \
        MINLEN:${TRIM_MIN_LENGTH}

    # Clean up unpaired (we only use paired-end reads for quantification)
    rm -f "$TRIM_R1_UNPAIRED" "$TRIM_R2_UNPAIRED"

    echo "[$(date -Is)] Trimmed R1: $(du -sh "$TRIM_R1_PAIRED" 2>/dev/null | cut -f1 || echo 'N/A')"
    echo "[$(date -Is)] Trimmed R2: $(du -sh "$TRIM_R2_PAIRED" 2>/dev/null | cut -f1 || echo 'N/A')"

else
    # Single-end
    RAW="${FASTQ_DIR}/${SRR}.fastq.gz"
    TRIM_SE="${FASTQ_DIR}/${SRR}_trimmed.fq.gz"

    if [[ -f "$TRIM_SE" ]]; then
        echo "  Already trimmed — skipping"; exit 0
    fi

    [[ -f "$RAW" ]] || { echo "ERROR: $RAW not found" >&2; exit 1; }

    echo "  Trimming single-end reads with Trimmomatic..."
    java -jar "$TRIMMOMATIC_JAR" SE -threads "$THREADS" \
        "$RAW" "$TRIM_SE" \
        ILLUMINACLIP:"${ADAPTERS_DIR}/TruSeq3-SE.fa":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:${TRIM_QUALITY} \
        MINLEN:${TRIM_MIN_LENGTH}

    echo "[$(date -Is)] Trimmed: $(du -sh "$TRIM_SE" 2>/dev/null | cut -f1 || echo 'N/A')"
fi

echo "[$(date -Is)] Trimmomatic complete: ${SRR}"
