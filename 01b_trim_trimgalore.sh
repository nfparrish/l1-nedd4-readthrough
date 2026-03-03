#!/bin/bash
#SBATCH --job-name=trim_tg
#SBATCH --output=/work/nfp8/2026_03_02/GSE226189/logs/trim_%A_%a.out
#SBATCH --error=/work/nfp8/2026_03_02/GSE226189/logs/trim_%A_%a.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH --time=02:00:00
#SBATCH --array=1-82
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"
module load "$MOD_TRIMGALORE"

# ---- Resolve sample ----
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -z "$SRR" ]] && { echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: trimming ${SRR}"

mkdir -p "$FASTQ_DIR"

if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
    R1="${FASTQ_DIR}/${SRR}_1.fastq.gz"
    R2="${FASTQ_DIR}/${SRR}_2.fastq.gz"

    # Trim Galore PE output naming: SRR_1_val_1.fq.gz, SRR_2_val_2.fq.gz
    TRIM_R1="${FASTQ_DIR}/${SRR}_1_val_1.fq.gz"
    TRIM_R2="${FASTQ_DIR}/${SRR}_2_val_2.fq.gz"

    if [[ -f "$TRIM_R1" && -f "$TRIM_R2" ]]; then
        echo "  Already trimmed — skipping"; exit 0
    fi

    [[ -f "$R1" ]] || { echo "ERROR: $R1 not found" >&2; exit 1; }
    [[ -f "$R2" ]] || { echo "ERROR: $R2 not found" >&2; exit 1; }

    trim_galore \
        --path_to_cutadapt "${VENV_DIR}/bin/cutadapt" \
        --quality "${TRIM_QUALITY}" \
        --length  "${TRIM_MIN_LENGTH}" \
        --paired \
        --fastqc \
        --output_dir "$FASTQ_DIR" \
        "$R1" "$R2"

    echo "[$(date -Is)] Trimmed R1: $(du -sh "$TRIM_R1" | cut -f1)"
    echo "[$(date -Is)] Trimmed R2: $(du -sh "$TRIM_R2" | cut -f1)"

else
    # Single-end — output: ${SRR}_trimmed.fq.gz
    RAW="${FASTQ_DIR}/${SRR}.fastq.gz"
    TRIM_SE="${FASTQ_DIR}/${SRR}_trimmed.fq.gz"

    if [[ -f "$TRIM_SE" ]]; then
        echo "  Already trimmed — skipping"; exit 0
    fi

    [[ -f "$RAW" ]] || { echo "ERROR: $RAW not found" >&2; exit 1; }

    trim_galore \
        --path_to_cutadapt "${VENV_DIR}/bin/cutadapt" \
        --quality "${TRIM_QUALITY}" \
        --length  "${TRIM_MIN_LENGTH}" \
        --fastqc \
        --output_dir "$FASTQ_DIR" \
        "$RAW"
fi

echo "[$(date -Is)] Trim Galore complete: ${SRR}"
