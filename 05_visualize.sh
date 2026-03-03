#!/bin/bash
#SBATCH --job-name=visualize_rt
#SBATCH --output=/work/nfp8/2026_03_02/GSE226189/logs/visualize_%j.out
#SBATCH --error=/work/nfp8/2026_03_02/GSE226189/logs/visualize_%j.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=8G
#SBATCH -c 2
#SBATCH --time=01:00:00
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"

FIG_DIR="${PERSISTENT_RESULTS}/figures"
VENV_PYTHON="${VENV_DIR}/bin/python3"
mkdir -p "$FIG_DIR"

${VENV_PYTHON} "${SCRIPT_DIR}/05_visualize.py" \
    --results-dir  "$RESULTS_DIR" \
    --coverage-dir "$COVERAGE_DIR" \
    --rt-window    "$RT_WINDOW" \
    --output-dir   "$FIG_DIR" \
    --gse          "$GSE"

echo "[$(date -Is)] Visualization job complete — figures in ${FIG_DIR}"
