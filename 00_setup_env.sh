#!/bin/bash
#SBATCH --job-name=setup_env
#SBATCH --output=/hpc/home/nfp8/copilot/2026_03_03/logs/setup_env_%j.out
#SBATCH --error=/hpc/home/nfp8/copilot/2026_03_03/logs/setup_env_%j.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=4G
#SBATCH -c 2
#SBATCH --time=01:00:00
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"

echo "[$(date -Is)] Setting up Python virtual environment and directory tree"

# ---- Create all directories ----
for d in "$BAM_DIR" "$COVERAGE_DIR" "$COUNTS_DIR" "$RESULTS_DIR" "$LOG_DIR" \
         "$PERSISTENT_BAM" "$PERSISTENT_RESULTS" \
         "$REF_DIR" "${STAR_INDEX}" \
         "${SCRIPT_DIR}/logs"; do
    mkdir -p "$d"
done

# ---- Python venv with pandas, matplotlib, seaborn, multiqc ----
if [[ ! -f "${VENV_DIR}/bin/python" ]]; then
    echo "Creating Python venv at ${VENV_DIR}"
    python3 -m venv "${VENV_DIR}"
fi

# Use venv Python directly (no activate script sourcing)
"${VENV_DIR}/bin/pip" install --upgrade pip
# Note: cutadapt pinned to 4.x for compatibility with Trim Galore 0.4.3
"${VENV_DIR}/bin/pip" install pandas matplotlib seaborn multiqc 'cutadapt>=3.7,<5'

echo "[$(date -Is)] Python venv ready at ${VENV_DIR}"
echo "  Packages: $(${VENV_DIR}/bin/pip list --format=columns 2>/dev/null | grep -iE 'pandas|matplotlib|seaborn|multiqc' | tr '\n' ' ')"

echo "[$(date -Is)] Setup complete"
