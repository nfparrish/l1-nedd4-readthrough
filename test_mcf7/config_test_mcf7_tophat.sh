#!/usr/bin/env bash
# ============================================================================
# config_test_mcf7_tophat.sh — MCF7 TopHat2 troubleshooting remap
# ============================================================================

_PARENT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config_test_mcf7.sh"
source "$_PARENT"

export GSE="MCF7_ctrl_tophat"

export WORK_BASE="/work/nfp8/${GSE}"
export FASTQ_DIR="/work/nfp8/MCF7_ctrl/fastq"
export BAM_DIR="${WORK_BASE}/bam"
export COVERAGE_DIR="${WORK_BASE}/coverage"
export COUNTS_DIR="${WORK_BASE}/counts"
export RESULTS_DIR="${WORK_BASE}/results"
export LOG_DIR="${WORK_BASE}/logs"

export PERSISTENT_BAM="/hpc/group/parrishlab/L1_NEDD4/bam_nedd4_locus/${GSE}"
export PERSISTENT_RESULTS="/hpc/group/parrishlab/L1_NEDD4/results/${GSE}"

export MOD_TOPHAT="TopHat/2.1.1"
export MOD_BOWTIE2="Bowtie2/2.4.4-rhel8"
export MOD_PICARD="Picard/2.18.2"
export MOD_PYTHON2="Python/2.7.11"
export MOD_SAMTOOLS="samtools/1.21"

export BOWTIE2_INDEX_PREFIX="${REF_DIR}/bowtie2_GRCh38_primary_assembly"

# TopHat2 uses mate inner distance, not total insert size.
export TOPHAT_INSERT_METRICS_DIR="${RESULTS_DIR}/insert_size"
export TOPHAT_INSERT_SIZE_ENV="${TOPHAT_INSERT_METRICS_DIR}/tophat_insert_size.env"
export TOPHAT_INSERT_SIZE_SOURCE_BAM_DIR="/work/nfp8/MCF7_ctrl/bam"
export TOPHAT_INSERT_SIZE_SOURCE_SAMPLES="ERR973734 ERR973735"

# TopHat2 reports MAPQ=50 for unique alignments (STAR uses 255)
export UNIQUE_MAPQ=50

# Use trimmed FASTQ for troubleshooting; adapter trimming was likely performed
# in the paper even though it was not described explicitly.
export TOPHAT_USE_TRIMMED="${TOPHAT_USE_TRIMMED:-1}"
export TOPHAT_READ_REALIGN_EDIT_DIST=0
export TOPHAT_EXTRA_ARGS="${TOPHAT_EXTRA_ARGS:-}"
