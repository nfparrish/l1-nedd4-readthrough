#!/bin/bash
#SBATCH --job-name=dl_mcf7_ena
#SBATCH --output=/work/nfp8/MCF7_ctrl/logs/download_%j.out
#SBATCH --error=/work/nfp8/MCF7_ctrl/logs/download_%j.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=4G
#SBATCH -c 4
#SBATCH --time=04:00:00
#SBATCH --array=1-2
set -euo pipefail

# Source test config (injected via PIPELINE_CONFIG env var by launcher)
source "${PIPELINE_CONFIG}"

mkdir -p "$FASTQ_DIR" "$LOG_DIR"

# ---- Resolve accession for this array task ----
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -z "$ACC" ]] && { echo "ERROR: no accession at task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "[$(date -Is)] Downloading ${ACC} from ENA FTP"

# ---- Build ENA FTP base URL ----
# Convention: /vol1/fastq/{first6}/{accession}/ for 9-char ERR+6-digit accessions
# For 10-char accessions (ERR+7), path is /vol1/fastq/{first6}/{last3}/{accession}/
ENA_BASE="https://ftp.sra.ebi.ac.uk/vol1/fastq"
ACC_LEN="${#ACC}"
if [[ $ACC_LEN -le 9 ]]; then
    PREFIX="${ACC:0:6}"
    URL_DIR="${ENA_BASE}/${PREFIX}/${ACC}"
else
    PREFIX="${ACC:0:6}"
    SUBDIR="${ACC: -3}"
    URL_DIR="${ENA_BASE}/${PREFIX}/${SUBDIR}/${ACC}"
fi

# ---- Download ----
download_file() {
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        echo "  Already exists: $dest"
        return
    fi
    echo "  Fetching: $url"
    wget -q --show-progress --continue -O "${dest}.part" "$url" && mv "${dest}.part" "$dest"
    echo "  Done: $(du -sh "$dest" | cut -f1)"
}

if [[ "$LIBRARY_LAYOUT" == "PE" ]]; then
    download_file "${URL_DIR}/${ACC}_1.fastq.gz" "${FASTQ_DIR}/${ACC}_1.fastq.gz"
    download_file "${URL_DIR}/${ACC}_2.fastq.gz" "${FASTQ_DIR}/${ACC}_2.fastq.gz"
else
    download_file "${URL_DIR}/${ACC}.fastq.gz" "${FASTQ_DIR}/${ACC}.fastq.gz"
fi

echo "[$(date -Is)] Download complete: ${ACC}"
