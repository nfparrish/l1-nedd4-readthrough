#!/usr/bin/env bash
# ============================================================================
# config_test_mcf7.sh — Override config for MCF7 positive-control test
# ============================================================================
# Sources the parent config.sh for shared settings (modules, genomic coords,
# VENV_DIR, REF_DIR, STAR_INDEX) then overrides the dataset-specific vars.
# Used by: run_test_mcf7.sh  (export PIPELINE_CONFIG before sbatch)
# ============================================================================

# ---- Load shared / infrastructure settings first ----
_PARENT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
source "$_PARENT"

# ---- Dataset identity ----
export GSE="MCF7_ctrl"

# ---- ENA accessions (E-MTAB-3788 / PRJEB8975) ----
# Positive controls — MCF-7 untreated (replicates 1-2):
#   ERR973734, ERR973735
# Negative controls — MCF-7 ORF1 shRNA sh1085 (replicates 1-2):
#   ERR973728, ERR973729
#   The shRNA targets L1 ORF1. Transcripts initiated from the L1 sense promoter
#   (including read-through) contain ORF1 sequence and are degraded by RISC,
#   so RT CPM should be reduced relative to untreated cells.
# (Scramble shRNA replicates ERR973731-733 are available as a vehicle control if needed.)
#
# Library: TruSeq Stranded Total RNA, first-strand cDNA (SDRF: LIBRARY_STRAND=first strand)
# Read length: 150 bp PE (SPOT_LENGTH=300), confirming STAR_SJDB_OVERHANG=149 is correct
# STRANDEDNESS: "reverse" (fr-firststrand / RF). In STAR ReadsPerGene.out.tab,
#   correct column is "antisense" (column 4).
export SRR_LIST="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/ERR_list.txt"
export NUM_SAMPLES=4
export LIBRARY_LAYOUT="PE"
export STRANDEDNESS="reverse"

# ---- Directory layout (separate work area from GSE226189) ----
export WORK_BASE="/work/nfp8/MCF7_ctrl"
export FASTQ_DIR="${WORK_BASE}/fastq"
export BAM_DIR="${WORK_BASE}/bam"
export COVERAGE_DIR="${WORK_BASE}/coverage"
export COUNTS_DIR="${WORK_BASE}/counts"
export RESULTS_DIR="${WORK_BASE}/results"
export LOG_DIR="${WORK_BASE}/logs"

# ---- Persistent storage ----
export PERSISTENT_BAM="/hpc/group/parrishlab/L1_NEDD4/bam_nedd4_locus/MCF7_ctrl"
export PERSISTENT_RESULTS="/hpc/group/parrishlab/L1_NEDD4/results/MCF7_ctrl"

# ---- STAR params ----
# PRJEB8975 reads are 100 bp. The shared STAR index was built with
# sjdbOverhang=149 (for 150 bp NovaSeq). This is fine: STAR 2-pass mode
# recomputes junctions, and a larger overhang does not degrade shorter reads.
# Do NOT override STAR_SJDB_OVERHANG here — it is baked into the genome index.
