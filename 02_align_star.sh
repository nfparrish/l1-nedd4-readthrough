#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --output=/work/nfp8/2026_03_02/GSE226189/logs/align_%A_%a.out
#SBATCH --error=/work/nfp8/2026_03_02/GSE226189/logs/align_%A_%a.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=06:00:00
# Array range set dynamically by launcher; fallback:
#SBATCH --array=1-82
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"
module load "$MOD_STAR"
module load "$MOD_SAMTOOLS"

# ---- Resolve sample ----
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
if [[ -z "$SRR" ]]; then
    echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1
fi
echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: aligning ${SRR}"

OUTPREFIX="${BAM_DIR}/${SRR}_"
FINAL_BAM="${OUTPREFIX}Aligned.sortedByCoord.out.bam"

# ---- Skip if already done ----
if [[ -f "${FINAL_BAM}.bai" ]]; then
    echo "Already aligned & indexed — skipping"; exit 0
fi

mkdir -p "${BAM_DIR}"

# ---- Determine input files — prefer Trimmomatic output, fall back to raw ----
if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
    # Trimmomatic PE naming: SRR_1_trimmed_P.fq.gz / SRR_2_trimmed_P.fq.gz
    TRIM_R1="${FASTQ_DIR}/${SRR}_1_trimmed_P.fq.gz"
    TRIM_R2="${FASTQ_DIR}/${SRR}_2_trimmed_P.fq.gz"
    RAW_R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    RAW_R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
    if [[ -f "$TRIM_R1" && -f "$TRIM_R2" ]]; then
        READ_FILES="${TRIM_R1} ${TRIM_R2}"
        echo "  Using Trimmomatic outputs"
    elif [[ -f "$RAW_R1" && -f "$RAW_R2" ]]; then
        READ_FILES="${RAW_R1} ${RAW_R2}"
        echo "  WARNING: trimmed FASTQs not found — using raw"
    else
        echo "ERROR: no input FASTQs found for ${SRR}" >&2; exit 1
    fi
else
    TRIM_SE="${FASTQ_DIR}/${SRR}_trimmed.fq.gz"
    RAW_SE="${FASTQ_DIR}/${SRR}.fastq.gz"
    if [[ -f "$TRIM_SE" ]]; then
        READ_FILES="$TRIM_SE"
    elif [[ -f "$RAW_SE" ]]; then
        READ_FILES="$RAW_SE"
        echo "  WARNING: trimmed FASTQ not found — using raw"
    else
        echo "ERROR: no input FASTQ found for ${SRR}" >&2; exit 1
    fi
fi

# ---- STAR alignment ----
STAR \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn $READ_FILES \
    --readFilesCommand zcat \
    --runThreadN "$STAR_THREADS" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax "$STAR_MULTIMAP_MAX" \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --outFileNamePrefix "$OUTPREFIX" \
    --outTmpDir "${BAM_DIR}/tmp_STAR_${SRR}_$$"

# ---- Index BAM ----
samtools index "$FINAL_BAM"

echo "[$(date -Is)] Alignment complete: ${SRR}"
echo "  BAM: $(du -sh "$FINAL_BAM" | cut -f1)"
echo "  Reads: $(samtools flagstat "$FINAL_BAM" | head -1)"
