#!/bin/bash
#SBATCH --job-name=tophat_align
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=48G
#SBATCH -c 8
#SBATCH --time=48:00:00
#SBATCH --array=1-4
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config_test_mcf7_tophat.sh}"

source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
module load "$MOD_PYTHON2"
module load "$MOD_TOPHAT"
module load "$MOD_BOWTIE2"
module load "$MOD_SAMTOOLS"

[[ -f "$TOPHAT_INSERT_SIZE_ENV" ]] || { echo "ERROR: missing insert-size env ${TOPHAT_INSERT_SIZE_ENV}" >&2; exit 1; }
source "$TOPHAT_INSERT_SIZE_ENV"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -n "$SRR" ]] || { echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

echo "[$(date -Is)] TopHat2 aligning ${SRR}"

if [[ "$TOPHAT_USE_TRIMMED" == "1" ]]; then
    R1="${FASTQ_DIR}/${SRR}_1_trimmed_P.fq.gz"
    R2="${FASTQ_DIR}/${SRR}_2_trimmed_P.fq.gz"
else
    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
fi

[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: missing FASTQ for ${SRR}" >&2; exit 1; }

mkdir -p "$BAM_DIR" "$LOG_DIR"

SAMPLE_OUT_DIR="${BAM_DIR}/${SRR}_tophat"
FINAL_BAM="${BAM_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"

if [[ -f "${FINAL_BAM}.bai" ]]; then
    echo "[$(date -Is)] ${SRR} already aligned in TopHat2 branch — skipping"
    exit 0
fi

rm -rf "$SAMPLE_OUT_DIR"

tophat2 \
    --num-threads "${ALIGN_CPUS:-8}" \
    --read-realign-edit-dist "$TOPHAT_READ_REALIGN_EDIT_DIST" \
    --mate-inner-dist "$TOPHAT_MATE_INNER_DIST" \
    --mate-std-dev "$TOPHAT_MATE_STD_DEV" \
    --output-dir "$SAMPLE_OUT_DIR" \
    $TOPHAT_EXTRA_ARGS \
    "$BOWTIE2_INDEX_PREFIX" \
    "$R1" "$R2"

samtools sort -@ "${ALIGN_CPUS:-8}" -o "$FINAL_BAM" "${SAMPLE_OUT_DIR}/accepted_hits.bam"
samtools index "$FINAL_BAM"

cp "${SAMPLE_OUT_DIR}/align_summary.txt" "${BAM_DIR}/${SRR}_align_summary.txt"

echo "[$(date -Is)] TopHat2 alignment complete: ${SRR}"
echo "  Final BAM: ${FINAL_BAM}"
echo "  Using mate inner distance ${TOPHAT_MATE_INNER_DIST} and std dev ${TOPHAT_MATE_STD_DEV}"
