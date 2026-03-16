#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(cd "${TEST_DIR}/.." && pwd)"
export PIPELINE_CONFIG="${TEST_DIR}/config_test_mcf7_tophat.sh"

source "$PIPELINE_CONFIG"

LOG_DIR_HOME="${TEST_DIR}/logs"
mkdir -p "$LOG_DIR_HOME"

echo "=================================================================="
echo " MCF7 TopHat2 Troubleshooting Remap"
echo " $(date)"
echo " Config:  ${PIPELINE_CONFIG}"
echo " Samples: $(tr '\n' ' ' < "$SRR_LIST")"
echo " Work:    ${WORK_BASE}"
echo " Logs:    ${LOG_DIR_HOME}"
echo "=================================================================="

DEP_ON=""

submit() {
    local label="$1" script="$2"
    shift 2
    local dep_flag=""
    [[ -n "$DEP_ON" ]] && dep_flag="--dependency=afterok:${DEP_ON}"
    sbatch \
        --export=ALL,PIPELINE_CONFIG="${PIPELINE_CONFIG}" \
        --output="${LOG_DIR_HOME}/${label}_%A_%a.out" \
        --error="${LOG_DIR_HOME}/${label}_%A_%a.err" \
        $dep_flag "$@" "$script" | awk '{print $NF}'
}

submit_dep2() {
    local label="$1" dep_str="$2" script="$3"
    shift 3
    sbatch \
        --export=ALL,PIPELINE_CONFIG="${PIPELINE_CONFIG}" \
        --output="${LOG_DIR_HOME}/${label}_%A_%a.out" \
        --error="${LOG_DIR_HOME}/${label}_%A_%a.err" \
        --dependency="afterok:${dep_str}" \
        "$@" "$script" | awk '{print $NF}'
}

JOB_SETUP=$(submit "setup_tophat" "${PIPELINE_DIR}/00_setup_env.sh" \
    --job-name=setup_mcf7_tophat \
    --mem=4G -c 2 --time=01:00:00)
echo "[Step 0] Setup env: job ${JOB_SETUP}"
DEP_ON="$JOB_SETUP"

JOB_INDEX=$(submit "bt2_index" "${TEST_DIR}/01_build_bowtie2_index_tophat.sh" \
    --job-name=bt2_index_mcf7 \
    --mem=48G -c 12 --time=08:00:00)
echo "[Step 1] Bowtie2 index: job ${JOB_INDEX}"

JOB_INSERT=$(submit "insert" "${TEST_DIR}/01c_estimate_insert_size_tophat.sh" \
    --job-name=insert_mcf7_tophat \
    --mem=8G -c 2 --time=01:00:00)
echo "[Step 2] Insert-size estimation: job ${JOB_INSERT}"

JOB_ALIGN=$(submit_dep2 "align_tophat" "${JOB_INDEX}:${JOB_INSERT}" \
    "${TEST_DIR}/02_align_tophat.sh" \
    --array=1-4 \
    --job-name=align_mcf7_tophat \
    --mem=48G -c 8 --time=24:00:00)
echo "[Step 3] TopHat2 align (array 1-4): job ${JOB_ALIGN}"
DEP_ON="$JOB_ALIGN"

JOB_POST=$(submit "post_tophat" "${PIPELINE_DIR}/03_post_align.sh" \
    --array=1-4 \
    --job-name=post_mcf7_tophat \
    --mem=64G -c 4 --time=08:00:00)
echo "[Step 4] Post-align (array 1-4): job ${JOB_POST}"
DEP_ON="$JOB_POST"

JOB_COLLECT=$(submit "collect_tophat" "${PIPELINE_DIR}/04_collect_results.sh" \
    --job-name=collect_mcf7_tophat \
    --mem=16G -c 4 --time=02:00:00)
echo "[Step 5] Collect results: job ${JOB_COLLECT}"
DEP_ON="$JOB_COLLECT"

JOB_VIS=$(submit "vis_tophat" "${PIPELINE_DIR}/05_visualize.sh" \
    --job-name=vis_mcf7_tophat \
    --mem=8G -c 2 --time=01:00:00)
echo "[Step 6] Visualize: job ${JOB_VIS}"

echo "=================================================================="
echo " All TopHat2 troubleshooting jobs submitted. Track with:"
echo "   squeue -u nfp8 --sort=i"
echo " Results will land under:"
echo "   ${RESULTS_DIR}"
echo "   ${PERSISTENT_RESULTS}"
echo "=================================================================="
