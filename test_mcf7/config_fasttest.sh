#!/bin/bash
# Fast test config using 1 trimmed MCF7 sample (~1.7GB total)
# This allows testing downstream pipeline steps in parallel while full trim runs

export PIPELINE_DIR="/hpc/home/nfp8/copilot/2026_03_03"
export WORK_DIR="/work/nfp8/2026_03_03_fasttest"
export RESULTS_DIR="${WORK_DIR}/persistent_results"
export LOG_DIR="${WORK_DIR}/logs"

# Module versions
export MOD_STAR="STAR/2.7.11a"
export MOD_SAMTOOLS="samtools/1.21"
export MOD_SUBREAD="Subread/2.0.3"
export MOD_BEDTOOLS="Bedtools/2.30.0"
export MOD_PICARD="Picard/2.18.2"
export MOD_R="R/4.4.3"
export MOD_FASTQC="FastQC/0.11.7"

# Library metadata
export LIBRARY_LAYOUT="PE"
export STRANDEDNESS="reverse"

# Single sample for fast testing
export SAMPLES=("ERR973728")
export SRR_LIST="${WORK_DIR}/ERR_list.txt"
export NUM_SAMPLES=1

# Use trimmed Paired-end reads (already created by Trimmomatic, ~850MB each)
export FASTQ_DIR="/work/nfp8/MCF7_ctrl/fastq"
export TRIM_SUFFIX="_trimmed_P"  # Trimmomatic PE output
export BAM_DIR="${WORK_DIR}/bam"
export COVERAGE_DIR="${WORK_DIR}/coverage"
export COUNTS_DIR="${WORK_DIR}/counts"
export LOG_DIR="${WORK_DIR}/logs"

# Use existing STAR index (built during MCF7 test, stored in persistent ref dir)
export STAR_INDEX="/hpc/group/parrishlab/refs/star_GRCh38_v45_sjdb149"
export GTF="/hpc/group/parrishlab/refs/gencode.v45.primary_assembly.annotation.gtf"

# NEDD4 / L1 insertion locus (hg38)
# L1 breakpoint: chr15:55,958,934-55,958,936 (plus strand)
# RT_WINDOW: 1 kb downstream of L1 3' breakpoint
export RT_WINDOW="chr15:55958936-55959936"
export BAM_SUBSET_REGION="chr15:54958936-56959936"
export BACKGROUND_WINDOWS_BED="/hpc/group/parrishlab/refs/nedd4_background_windows.bed"
export NEDD4_EXONS_BED="/hpc/group/parrishlab/refs/nedd4_exons.bed"

# Trimming parameters (not used in this test, but kept for completeness)
export TRIM_QUALITY=20
export TRIM_MIN_LENGTH=36

