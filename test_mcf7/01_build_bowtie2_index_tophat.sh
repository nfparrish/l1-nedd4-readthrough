#!/bin/bash
#SBATCH --job-name=bt2_index_mcf7
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=48G
#SBATCH -c 12
#SBATCH --time=08:00:00
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config_test_mcf7_tophat.sh}"

source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
module load "$MOD_BOWTIE2"

mkdir -p "$(dirname "$BOWTIE2_INDEX_PREFIX")"

if [[ -f "${BOWTIE2_INDEX_PREFIX}.1.bt2" || -f "${BOWTIE2_INDEX_PREFIX}.1.bt2l" ]]; then
    echo "[$(date -Is)] Bowtie2 index already present at ${BOWTIE2_INDEX_PREFIX}"
    exit 0
fi

echo "[$(date -Is)] Building Bowtie2 index from ${GENOME_FA}"
bowtie2-build --threads "${INDEX_CPUS:-12}" "$GENOME_FA" "$BOWTIE2_INDEX_PREFIX"
echo "[$(date -Is)] Bowtie2 index complete: ${BOWTIE2_INDEX_PREFIX}"
