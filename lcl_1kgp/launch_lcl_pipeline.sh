#!/usr/bin/env bash
# launch_lcl_pipeline.sh
# Submit: (1) download, then (2a) TopHat2 + (2b) STAR in parallel after download.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRR_LIST="/work/nfp8/LCL_1kgp/metadata/srr_list.txt"
N=$(wc -l < "$SRR_LIST")
echo "Submitting pipeline for ${N} samples"

# 1 — Download
DL_JOB=$(sbatch \
    --parsable \
    --array="1-${N}" \
    "${SCRIPT_DIR}/02_download_lcl.sh")
echo "Download job: ${DL_JOB}"

# 2a — TopHat2 (after all downloads complete)
TH_JOB=$(sbatch \
    --parsable \
    --array="1-${N}" \
    --dependency="afterok:${DL_JOB}" \
    "${SCRIPT_DIR}/03_align_tophat_lcl.sh")
echo "TopHat2 job: ${TH_JOB}"

# 2b — STAR (after all downloads complete, runs in parallel with TopHat2)
ST_JOB=$(sbatch \
    --parsable \
    --array="1-${N}" \
    --dependency="afterok:${DL_JOB}" \
    "${SCRIPT_DIR}/04_align_star_lcl.sh")
echo "STAR job: ${ST_JOB}"

echo ""
echo "Pipeline submitted. Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f /work/nfp8/LCL_1kgp/logs/download_${DL_JOB}_1.out"
echo ""
echo "After alignments complete, run infer_experiment.py on a 1/1 sample BAM to"
echo "confirm strandedness, then generate IGV reports."
echo ""
echo "TopHat2 BAMs:  /work/nfp8/LCL_1kgp/bam_tophat/"
echo "STAR BAMs:     /work/nfp8/LCL_1kgp/bam_star/"
echo "Metadata:      /work/nfp8/LCL_1kgp/metadata/"
