#!/usr/bin/env bash
# ============================================================================
# run_test_mcf7.sh — MCF7 positive-control test run
# ============================================================================
# Downloads ERR973734 + ERR973735 (MCF-7 untreated, positive control) and
# ERR973728 + ERR973729 (MCF-7 ORF1 shRNA sh1085, negative control) from ENA,
# then runs the full L1-NEDD4 pipeline on these 4 samples.
#
# Expected result:
#   Untreated (ERR973734/5): detectable RT CPM with decay gradient
#   ORF1 shRNA (ERR973728/9): reduced RT CPM (L1 transcript degraded by RISC)
# ============================================================================
set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(cd "${TEST_DIR}/.." && pwd)"
export PIPELINE_CONFIG="${TEST_DIR}/config_test_mcf7.sh"

# Source the test config so we can reference its vars here
source "$PIPELINE_CONFIG"

# ---- Log directory in home (always writable from login node) ----
LOG_DIR_HOME="${TEST_DIR}/logs"
mkdir -p "$LOG_DIR_HOME"

echo "=================================================================="
echo " L1-NEDD4 MCF7 Positive-Control Test"
echo " $(date)"
echo " Config:  ${PIPELINE_CONFIG}"
echo " Samples: $(cat "$SRR_LIST" | tr '\n' ' ')"
echo " Work:    ${WORK_BASE}"
echo " Logs:    ${LOG_DIR_HOME}"
echo "=================================================================="

# Verify STAR index and venv exist before submitting
if [[ ! -f "${STAR_INDEX}/SA" ]]; then
    echo "WARNING: STAR index not found at ${STAR_INDEX}"
    echo "         Submitting anyway — alignment will fail if index is missing."
fi
if [[ ! -f "${VENV_DIR}/bin/activate" ]]; then
    echo "WARNING: Python venv not found at ${VENV_DIR}"
    echo "         Run 00_setup_env.sh first."
fi

# Helper: submit with optional afterok dependency, log to home dir
DEP_ON=""
submit() {
    local label="$1" script="$2"
    shift 2
    local dep_flag=""
    [[ -n "$DEP_ON" ]] && dep_flag="--dependency=afterok:${DEP_ON}"
    local jid
    jid=$(sbatch \
        --export=ALL,PIPELINE_CONFIG="${PIPELINE_CONFIG}" \
        --output="${LOG_DIR_HOME}/${label}_%A_%a.out" \
        --error="${LOG_DIR_HOME}/${label}_%A_%a.err" \
        $dep_flag "$@" "$script" | awk '{print $NF}')
    echo "$jid"
}

# Helper: submit with an explicit dual dependency string
submit_dep2() {
    local label="$1" dep_str="$2" script="$3"
    shift 3
    local jid
    jid=$(sbatch \
        --export=ALL,PIPELINE_CONFIG="${PIPELINE_CONFIG}" \
        --output="${LOG_DIR_HOME}/${label}_%A_%a.out" \
        --error="${LOG_DIR_HOME}/${label}_%A_%a.err" \
        --dependency="afterok:${dep_str}" \
        "$@" "$script" | awk '{print $NF}')
    echo "$jid"
}

# ------ Step 0: Download from ENA FTP ----------------------
JOB_DL=$(submit "dl" "${TEST_DIR}/00_download_mcf7.sh" \
    --array=1-4 \
    --job-name=dl_mcf7)
echo "[Step 0] Download (array 1-4): job ${JOB_DL}"
DEP_ON="$JOB_DL"

# ------ Step 0.5: Setup env (venv + dirs) if needed -------
JOB_SETUP=$(submit "setup" "${PIPELINE_DIR}/00_setup_env.sh" \
    --job-name=setup_mcf7 \
    --mem=4G -c 2 --time=01:00:00)
echo "[Step 0.5] Setup env: job ${JOB_SETUP}"
DEP_ON="$JOB_SETUP"

# ------ Step 1: Build STAR index --------------------------
# Runs after setup; trim runs in parallel
JOB_IDX=$(submit "star_index" "${PIPELINE_DIR}/01_build_star_index.sh" \
    --job-name=index_mcf7 \
    --mem=48G -c 12 --time=04:00:00)
echo "[Step 1] Build STAR index: job ${JOB_IDX}"
# Do NOT advance DEP_ON here — trim also depends on setup (DEP_ON still = JOB_SETUP)

# ------ Step 1.5: Adapter trimming (parallel with index) --
JOB_TRIM=$(submit "trim" "${PIPELINE_DIR}/01b_trim_trimgalore.sh" \
    --array=1-4 \
    --job-name=trim_mcf7 \
    --mem=8G -c 4 --time=02:00:00)
echo "[Step 1.5] Trim Galore (array 1-4): job ${JOB_TRIM}"

# ------ Step 2: STAR alignment — depends on BOTH index AND trim
JOB_ALIGN=$(submit_dep2 "align" "${JOB_IDX}:${JOB_TRIM}" \
    "${PIPELINE_DIR}/02_align_star.sh" \
    --array=1-4 \
    --job-name=align_mcf7 \
    --mem=36G -c 8 --time=06:00:00)
echo "[Step 2] Align (array 1-4): job ${JOB_ALIGN}"
DEP_ON="$JOB_ALIGN"

# ------ Step 3: Post-alignment ----------------------------
JOB_POST=$(submit "post" "${PIPELINE_DIR}/03_post_align.sh" \
    --array=1-4 \
    --job-name=post_mcf7 \
    --mem=16G -c 4 --time=02:00:00)
echo "[Step 3] Post-align (array 1-4): job ${JOB_POST}"
DEP_ON="$JOB_POST"

# ------ Step 4: Collect results + TMM --------------------------
JOB_COLLECT=$(submit "collect" "${PIPELINE_DIR}/04_collect_results.sh" \
    --job-name=collect_mcf7 \
    --mem=16G -c 4 --time=01:00:00)
echo "[Step 4] Collect: job ${JOB_COLLECT}"
DEP_ON="$JOB_COLLECT"

# ------ Step 5: Visualize ----------------------------
JOB_VIS=$(submit "visualize" "${PIPELINE_DIR}/05_visualize.sh" \
    --job-name=vis_mcf7 \
    --mem=8G -c 2 --time=00:30:00)
echo "[Step 5] Visualize: job ${JOB_VIS}"
DEP_ON="$JOB_VIS"

# ------ Step 6: Sanity check — assert MCF7 RT CPM > 0 ----------
JOB_CHECK=$(submit "check" "${TEST_DIR}/06_sanity_check.sh" \
    --job-name=check_mcf7 \
    --mem=2G -c 1 --time=00:10:00)
echo "[Step 6] Sanity check: job ${JOB_CHECK}"

echo "=================================================================="
echo " All jobs submitted. Track with:"
echo "   squeue -u nfp8 --sort=i"
echo " Final results will be in:"
echo "   ${PERSISTENT_RESULTS}/"
echo "=================================================================="
