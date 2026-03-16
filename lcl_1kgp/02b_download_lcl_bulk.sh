#!/bin/bash
#SBATCH --job-name=dl_lcl_bulk
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH --time=06:00:00
#SBATCH --array=1-719
#SBATCH --output=/work/nfp8/LCL_1kgp/logs/download_bulk_%A_%a.out
#SBATCH --error=/work/nfp8/LCL_1kgp/logs/download_bulk_%A_%a.err
set -euo pipefail

SRR_LIST="/work/nfp8/LCL_1kgp/metadata/srr_remaining.txt"
FASTQ_DIR="/work/nfp8/LCL_1kgp/fastq"
SRA_CACHE="/work/nfp8/LCL_1kgp/sra_cache"

mkdir -p "$FASTQ_DIR" "$SRA_CACHE" /work/nfp8/LCL_1kgp/logs

source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
module load SRA-Toolkit/3.0.0-rhel8

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -n "$SRR" ]] || { echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: downloading ${SRR}"

R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"
if [[ -f "$R1" && -f "$R2" ]]; then
    echo "[$(date -Is)] ${SRR} already downloaded — skipping"
    exit 0
fi

fasterq-dump \
    --split-files \
    --threads 8 \
    --outdir "$FASTQ_DIR" \
    --temp "$SRA_CACHE" \
    "$SRR"

gzip -f "${FASTQ_DIR}/${SRR}_1.fastq"
gzip -f "${FASTQ_DIR}/${SRR}_2.fastq"

rm -rf "${SRA_CACHE}/${SRR}" "${SRA_CACHE}/${SRR}.sra"

echo "[$(date -Is)] Done: ${SRR}  R1=$(du -sh ${R1} | cut -f1)  R2=$(du -sh ${R2} | cut -f1)"
