#!/usr/bin/env bash
# ============================================================================
# L1-NEDD4 Read-Through Detection Pipeline — Configuration
# ============================================================================
# Edit this file to target a new dataset. All pipeline scripts source this.
# ============================================================================

# --- Dataset identity ---
export GSE="GSE226189"
export BIOPROJECT="PRJNA939148"
export LIBRARY_LAYOUT="PE"       # PE (paired-end) or SE (single-end)
export STRANDEDNESS="unstranded"  # unstranded | forward | reverse (verify w/ infer_experiment.py)

# --- SRR list (one accession per line) ---
export SRR_LIST="/work/nfp8/2026_03_02/${GSE}/SRR_list.txt"
export NUM_SAMPLES=$(wc -l < "$SRR_LIST" 2>/dev/null || echo 82)

# --- Directory layout ---
# Transient work (auto-purged after 75 days on /work)
export WORK_BASE="/work/nfp8/2026_03_02/${GSE}"
export FASTQ_DIR="${WORK_BASE}/fastq"
export BAM_DIR="${WORK_BASE}/bam"
export COVERAGE_DIR="${WORK_BASE}/coverage"
export COUNTS_DIR="${WORK_BASE}/counts"
export RESULTS_DIR="${WORK_BASE}/results"
export LOG_DIR="${WORK_BASE}/logs"

# Persistent storage (no auto-purge)
export PERSISTENT_BASE="/hpc/group/parrishlab/L1_NEDD4"
export PERSISTENT_BAM="${PERSISTENT_BASE}/bam_nedd4_locus/${GSE}"
export PERSISTENT_RESULTS="${PERSISTENT_BASE}/results/${GSE}"

# Scripts (in home directory for persistence)
export SCRIPT_DIR="/hpc/home/nfp8/copilot/2026_03_03"

# Reference genome & STAR index (persistent, shared across datasets)
export REF_DIR="/hpc/group/parrishlab/refs"
export GENOME_FA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
export GENCODE_GTF="${REF_DIR}/gencode.v45.primary_assembly.annotation.gtf"
export STAR_INDEX="${REF_DIR}/star_GRCh38_v45_sjdb149"

# Python virtual environment (persistent)
export VENV_DIR="/hpc/group/parrishlab/envs/l1nedd4_py39"

# --- SLURM defaults ---
export SLURM_PARTITION="common"
export SLURM_ACCOUNT="parrishlab"

# --- DCC modules ---
export MOD_STAR="STAR/2.7.11a"
export MOD_SAMTOOLS="samtools/1.21"
export MOD_R="R/4.4.3"
export MOD_FASTQC="FastQC/0.11.7"
# Note: Trimmomatic (now used instead of legacy Trim Galore 0.4.3) loads inline in trim script
export MOD_SUBREAD="Subread/2.0.3"       # featureCounts
export MOD_BEDTOOLS="Bedtools/2.30.0"

# --- L1-NEDD4 genomic coordinates (GRCh38) ---
# NEDD4 gene: chr15:55,826,917–55,993,660 (minus strand)
# L1-NEDD4 insertion breakpoint: chr15:55,958,934–55,958,936 (plus strand)
# Read-through window: 1 kb downstream of L1 3' breakpoint on + strand
export RT_WINDOW="chr15:55958936-55959936"
export L1_BREAKPOINT_CHR="chr15"
export L1_BREAKPOINT_3P=55958936
export PAPER_DOWNSTREAM_WINDOW="${RT_WINDOW}"
export PAPER_UPSTREAM_WINDOW="chr15:55957935-55958935"
export BAM_SUBSET_REGION="chr15:54958936-56959936"   # NEDD4 locus ±1 Mb

# Proxy SNP (for annotation only — NOT for window placement)
export PROXY_SNP_POS="chr15:55862437"    # rs16976600

# --- Adapter trimming parameters (applies to Trimmomatic) ---
export TRIM_QUALITY=20          # Phred quality cutoff
export TRIM_MIN_LENGTH=36       # discard reads shorter than this post-trim

# --- featureCounts strand code (derived from STRANDEDNESS at runtime) ---
# 0=unstranded  1=forward (fr-secondstrand)  2=reverse (fr-firststrand/TruSeq)
# Set automatically in 03_post_align.sh from $STRANDEDNESS

# --- featureCounts scope ---
# SKIP   = skip featureCounts entirely
# PANEL  = ~100 housekeeping genes only (fast; useful for QC / normalization)
# CHR15  = all chr15 genes (good balance for NEDD4-focused analysis)
# GENOME = full annotation, all chromosomes (default; most complete)
export FC_SCOPE="${FC_SCOPE:-GENOME}"  # respect override from --export
export FC_PANEL_GTF="${REF_DIR}/housekeeping_panel.gtf"
export FC_CHR15_GTF="${REF_DIR}/gencode.v45.chr15.gtf"

# --- Background windows for RT CPM flanking normalization ---
# 10 × 1 kb intronic windows: 5 upstream and 5 downstream of L1, verified within NEDD4 intron 19
# Exon positions are subtracted at runtime via bedtools using NEDD4_EXONS_BED.
export BACKGROUND_WINDOWS_BED="${REF_DIR}/nedd4_background_windows.bed"
export NEDD4_EXONS_BED="${REF_DIR}/nedd4_exons.bed"
export NEDD4_ENSEMBL_ID="ENSG00000069869"

# --- STAR alignment parameters ---
export STAR_SJDB_OVERHANG=149    # max read length - 1 (NovaSeq 150bp)
export STAR_THREADS=8
export STAR_MULTIMAP_MAX=1       # uniquely mapped only for read-through

# --- Resource requests ---
export ALIGN_MEM="36G"
export ALIGN_CPUS=8
export ALIGN_TIME="06:00:00"
export INDEX_MEM="48G"
export INDEX_CPUS=12
export INDEX_TIME="04:00:00"
export QUANT_MEM="8G"
export QUANT_CPUS=2
export QUANT_TIME="01:00:00"
