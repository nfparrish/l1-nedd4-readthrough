#!/bin/bash
#SBATCH --job-name=tophat_lcl
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=48G
#SBATCH -c 8
#SBATCH --time=12:00:00
#SBATCH --array=1-60
#SBATCH --output=/work/nfp8/LCL_1kgp/logs/tophat_%A_%a.out
#SBATCH --error=/work/nfp8/LCL_1kgp/logs/tophat_%A_%a.err
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
REF_DIR="/hpc/group/parrishlab/refs"
BOWTIE2_INDEX="${REF_DIR}/bowtie2_GRCh38_primary_assembly"
GENCODE_GTF="${REF_DIR}/gencode.v45.primary_assembly.annotation.gtf"
SRR_LIST="/work/nfp8/LCL_1kgp/metadata/srr_list.txt"
FASTQ_DIR="/work/nfp8/LCL_1kgp/fastq"
BAM_DIR="/work/nfp8/LCL_1kgp/bam_tophat"
LOG_DIR="/work/nfp8/LCL_1kgp/logs"

# Insert size: typical for NEBNext Ultra II polyA 2x150 PE
# (250 bp mean insert, fragment length ~400 bp; inner dist = fragment - 2*150 = 100)
# Adjust if picard CollectInsertSizeMetrics gives a different value.
MATE_INNER_DIST=100
MATE_STD_DEV=50

# TopHat2 reports MAPQ=50 for unique alignments
UNIQUE_MAPQ=50

# Strandedness: NEBNext dUTP is fr-firststrand; empirically may be fr-secondstrand
# Verify with infer_experiment.py after first alignment. Use "fr-firststrand" for now.
LIBRARY_TYPE="fr-firststrand"

mkdir -p "$BAM_DIR" "$LOG_DIR"

source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
module load Python/2.7.11
module load TopHat/2.1.1
module load Bowtie2/2.4.4-rhel8
module load samtools/1.21

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -n "$SRR" ]] || { echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: TopHat2 aligning ${SRR}"

R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: missing FASTQs for ${SRR}" >&2; exit 1; }

FINAL_BAM="${BAM_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"
if [[ -f "${FINAL_BAM}.bai" ]]; then
    echo "[$(date -Is)] ${SRR} already aligned — skipping"
    exit 0
fi

SAMPLE_OUT="${BAM_DIR}/${SRR}_tophat_tmp"
rm -rf "$SAMPLE_OUT"

tophat2 \
    --num-threads 8 \
    --library-type "$LIBRARY_TYPE" \
    --mate-inner-dist "$MATE_INNER_DIST" \
    --mate-std-dev "$MATE_STD_DEV" \
    --read-realign-edit-dist 0 \
    --output-dir "$SAMPLE_OUT" \
    "$BOWTIE2_INDEX" \
    "$R1" "$R2"

samtools sort -@ 8 -o "$FINAL_BAM" "${SAMPLE_OUT}/accepted_hits.bam"
samtools index "$FINAL_BAM"

cp "${SAMPLE_OUT}/align_summary.txt" "${BAM_DIR}/${SRR}_align_summary.txt"
rm -rf "$SAMPLE_OUT"

echo "[$(date -Is)] TopHat2 complete: ${SRR}"
echo "  BAM: $(du -sh "$FINAL_BAM" | cut -f1)"
echo "  Summary:"
cat "${BAM_DIR}/${SRR}_align_summary.txt"
