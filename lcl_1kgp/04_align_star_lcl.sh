#!/bin/bash
#SBATCH --job-name=star_lcl
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --time=06:00:00
#SBATCH --array=1-60
#SBATCH --output=/work/nfp8/LCL_1kgp/logs/star_%A_%a.out
#SBATCH --error=/work/nfp8/LCL_1kgp/logs/star_%A_%a.err
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
REF_DIR="/hpc/group/parrishlab/refs"
STAR_INDEX="${REF_DIR}/star_GRCh38_v45_sjdb149"
SRR_LIST="/work/nfp8/LCL_1kgp/metadata/srr_list.txt"
FASTQ_DIR="/work/nfp8/LCL_1kgp/fastq"
BAM_DIR="/work/nfp8/LCL_1kgp/bam_star"
LOG_DIR="/work/nfp8/LCL_1kgp/logs"

mkdir -p "$BAM_DIR" "$LOG_DIR"

source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
module load STAR/2.7.11a
module load samtools/1.21

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -n "$SRR" ]] || { echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: STAR aligning ${SRR}"

R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: missing FASTQs for ${SRR}" >&2; exit 1; }

OUTPREFIX="${BAM_DIR}/${SRR}_"
FINAL_BAM="${OUTPREFIX}Aligned.sortedByCoord.out.bam"

if [[ -f "${FINAL_BAM}.bai" ]]; then
    echo "[$(date -Is)] ${SRR} already aligned — skipping"
    exit 0
fi

STAR \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --runThreadN 8 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 1 \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --outFileNamePrefix "$OUTPREFIX" \
    --outTmpDir "${BAM_DIR}/tmp_STAR_${SRR}_$$"

samtools index "$FINAL_BAM"

echo "[$(date -Is)] STAR complete: ${SRR}"
echo "  BAM: $(du -sh "$FINAL_BAM" | cut -f1)"
echo "  Mapped: $(samtools flagstat "$FINAL_BAM" | head -1)"
