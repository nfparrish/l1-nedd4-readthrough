#!/bin/bash
#SBATCH --job-name=igv_report
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=01:00:00
set -euo pipefail

source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
module load samtools/1.21

PYTHON=/hpc/group/parrishlab/envs/l1nedd4_py39/bin/python3
SCRIPT=/hpc/home/nfp8/copilot/2026_03_03/scripts/make_igv_report.py
BAM_DIR=/work/nfp8/2026_03_02/GSE226189/bam
REGION="chr15:55958955-55959955"
OUT=/work/nfp8/2026_03_02/GSE226189/results/igv_report.html
TITLE="GSE226189 — IGV-style coverage at ${REGION}"

echo "[$(date -Is)] Starting IGV report for ${REGION}"
echo "  BAMs: ${BAM_DIR}"
echo "  Output: ${OUT}"

"$PYTHON" "$SCRIPT" \
    --bam-dir "$BAM_DIR" \
    --region "$REGION" \
    --out "$OUT" \
    --title "$TITLE"

echo "[$(date -Is)] Done: $OUT"
