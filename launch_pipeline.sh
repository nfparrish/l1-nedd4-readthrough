#!/usr/bin/env bash
# ============================================================================
# launch_pipeline.sh — L1-NEDD4 Read-Through Detection Pipeline Orchestrator
# ============================================================================
# Submits all SLURM jobs with dependency chains.
#
# Usage:
#   bash launch_pipeline.sh                  # start from scratch
#   bash launch_pipeline.sh --after JOBID    # chain after an existing job (e.g., download)
#   bash launch_pipeline.sh --skip-index     # skip genome build if index exists
#   bash launch_pipeline.sh --skip-setup     # skip Python env setup if venv exists
#   bash launch_pipeline.sh --from STEP      # resume from step number (0-5)
#
# Steps:
#   0   Setup Python env + directories
#   1   Build STAR genome index
#   1.5 Adapter trimming (Trim Galore, array) [runs parallel to index]
#   2   STAR alignment (array) [depends on both index AND trim]
#   3   Post-alignment QC (array)
#   4   Collect results + TMM normalization
#   5   Visualization
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# ── Parse arguments ────────────────────────────────────────────────────────
AFTER_JOB=""
SKIP_INDEX=false
SKIP_SETUP=false
FROM_STEP=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --after)     AFTER_JOB="$2"; shift 2 ;;
        --skip-index) SKIP_INDEX=true; shift ;;
        --skip-setup) SKIP_SETUP=true; shift ;;
        --from)      FROM_STEP="$2"; shift 2 ;;
        -h|--help)
            head -20 "$0" | grep "^#" | sed 's/^# *//'
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Helper: submit with optional dependency (pass as second arg or use DEP_ON)
submit() {
    local script="$1"
    shift
    local dep_flag=""
    if [[ -n "${DEP_ON:-}" ]]; then
        dep_flag="--dependency=afterok:${DEP_ON}"
    fi
    local jid
    jid=$(sbatch $dep_flag "$@" "$script" | awk '{print $NF}')
    echo "$jid"
}

# Helper: submit with an explicit dependency string (for dual-dep cases)
submit_dep() {
    local dep_str="$1" script="$2"
    shift 2
    local jid
    jid=$(sbatch --dependency="afterok:${dep_str}" "$@" "$script" | awk '{print $NF}')
    echo "$jid"
}

echo "=================================================================="
echo " L1-NEDD4 Pipeline — ${GSE}"
echo " $(date)"
echo " Script dir:  ${SCRIPT_DIR}"
echo " FASTQ dir:   ${FASTQ_DIR}"
echo " Work base:   ${WORK_BASE}"
echo " Persistent:  ${PERSISTENT_BASE}"
echo "=================================================================="

# Track the dependency chain
DEP_ON="${AFTER_JOB}"
JOB_INDEX=""
JOB_TRIM=""

# ── Step 0: Setup environment ─────────────────────────────────────────────
if [[ $FROM_STEP -le 0 ]]; then
    if $SKIP_SETUP && [[ -d "${VENV_DIR}" ]]; then
        echo "[Step 0] SKIP — venv already exists at ${VENV_DIR}"
    else
        JOB_SETUP=$(submit "${SCRIPT_DIR}/00_setup_env.sh")
        echo "[Step 0] Setup env: job ${JOB_SETUP}"
        DEP_ON="$JOB_SETUP"
    fi
fi

# ── Step 1: Build STAR index ──────────────────────────────────────────────
if [[ $FROM_STEP -le 1 ]]; then
    if $SKIP_INDEX && [[ -f "${STAR_INDEX}/SA" ]]; then
        echo "[Step 1] SKIP — STAR index exists at ${STAR_INDEX}"
    else
        JOB_INDEX=$(submit "${SCRIPT_DIR}/01_build_star_index.sh")
        echo "[Step 1] Build STAR index: job ${JOB_INDEX}"
        # Note: do NOT advance DEP_ON here — trim runs in parallel with index
    fi
fi

# ── Step 1.5: Adapter trimming (parallel with index build) ───────────────
if [[ $FROM_STEP -le 1 ]]; then
    JOB_TRIM=$(submit "${SCRIPT_DIR}/01b_trim_trimgalore.sh" \
        --array="1-${NUM_SAMPLES}%16")
    echo "[Step 1.5] Trim Galore (array 1-${NUM_SAMPLES}): job ${JOB_TRIM}"
fi

# ── Step 2: STAR alignment (array) — requires BOTH index AND trim ────────
if [[ $FROM_STEP -le 2 ]]; then
    # Build dual dependency: index + trim (either may be empty if skipped)
    ALIGN_DEP="${DEP_ON}"
    [[ -n "$JOB_INDEX" && -n "$JOB_TRIM" ]] && ALIGN_DEP="${JOB_INDEX}:${JOB_TRIM}"
    [[ -n "$JOB_INDEX" && -z "$JOB_TRIM" ]] && ALIGN_DEP="$JOB_INDEX"
    [[ -z "$JOB_INDEX" && -n "$JOB_TRIM" ]] && ALIGN_DEP="$JOB_TRIM"

    if [[ -n "$ALIGN_DEP" ]]; then
        JOB_ALIGN=$(submit_dep "$ALIGN_DEP" "${SCRIPT_DIR}/02_align_star.sh" \
            --array="1-${NUM_SAMPLES}%16")
    else
        JOB_ALIGN=$(sbatch --array="1-${NUM_SAMPLES}%16" \
            "${SCRIPT_DIR}/02_align_star.sh" | awk '{print $NF}')
    fi
    echo "[Step 2] STAR align (array 1-${NUM_SAMPLES}): job ${JOB_ALIGN}"
    DEP_ON="$JOB_ALIGN"
fi

# ── Step 3: Post-alignment QC (array) ─────────────────────────────────────
if [[ $FROM_STEP -le 3 ]]; then
    JOB_POST=$(submit "${SCRIPT_DIR}/03_post_align.sh" \
        --array="1-${NUM_SAMPLES}%32")
    echo "[Step 3] Post-align QC (array 1-${NUM_SAMPLES}): job ${JOB_POST}"
    DEP_ON="$JOB_POST"
fi

# ── Step 4: Collect results + normalization ────────────────────────────────
if [[ $FROM_STEP -le 4 ]]; then
    JOB_COLLECT=$(submit "${SCRIPT_DIR}/04_collect_results.sh")
    echo "[Step 4] Collect results: job ${JOB_COLLECT}"
    DEP_ON="$JOB_COLLECT"
fi

# ── Step 5: Visualization ─────────────────────────────────────────────────
if [[ $FROM_STEP -le 5 ]]; then
    JOB_VIS=$(submit "${SCRIPT_DIR}/05_visualize.sh")
    echo "[Step 5] Visualize: job ${JOB_VIS}"
fi

echo "=================================================================="
echo " All jobs submitted. Monitor with:"
echo "   squeue -u nfp8 --sort=i"
echo "   sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS"
echo "=================================================================="
