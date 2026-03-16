# Session Handoff — 2026-03-16

## Project Overview

This project measures **L1-NEDD4 readthrough transcription** as a function of
L1 insertion genotype at chr15:55,958,935 (hg38, plus strand).  The L1 element
produces sense-strand read-through into the 3′ flanking region.  We quantify
this via **RT_CPM**: the number of sense-strand (R1-forward) reads in a 1 kb
window downstream of the L1 breakpoint, normalized per million uniquely mapped
reads.

---

## 1. What Was Done This Session

### 1a. MCF7 Positive Control — RT_CPM Calculation
- Computed RT_CPM for **4 MCF7 samples × 2 aligners (STAR + TopHat2) = 8 values**.
- Script: `lcl_1kgp/compute_rt_cpm.py`
- Result: untreated MCF7 shows ~0.76 CPM (STAR mean), shRNA knockdown ~0.19 CPM → **4× reduction**, confirming the metric works as expected.

### 1b. GSE226189 Fibroblasts — RT_CPM (82 samples)
- Computed RT_CPM for all **82 skin fibroblast samples** in GSE226189.
- Script: `lcl_1kgp/compute_rt_cpm_GSE226189.py`
- Output TSV: `lcl_1kgp/GSE226189_rt_cpm.tsv`
- **Distribution**: min=0.07, mean=0.94, max=12.31 (SRR23630232 is a major outlier)

### 1c. 3-Genotype Gaussian Mixture Model Fit
- Fit log(RT_CPM) to a 3-component log-normal mixture with **Hardy-Weinberg equilibrium weights** (q=0.302 from 1000 Genomes: p²=0.487, 2pq=0.422, q²=0.091).
- Script: `lcl_1kgp/fit_rt_cpm_model.py`
- Figure: `lcl_1kgp/GSE226189_rt_cpm_histogram.png`
- **Key findings**:
  - **Free 3G model** (3 free means + shared σ) wins by AIC/BIC over Null and Dosage
  - **Linear dosage constraint is REJECTED** (LRT χ²(2)=18.40, p=0.0001)
  - Fitted means: 0/0→0.36 CPM, 0/1→0.47 CPM, 1/1→3.29 CPM
  - Fold-changes: 0/1 is 1.3× (not 2×), 1/1 is 9.2× (not 4×) — super-linear at homozygous
  - MAP assignments: 0/0=41, 0/1=30, 1/1=11 (close to HWE expectation 40/35/7.5)

---

## 2. SLURM Jobs — Current State (as of 2026-03-16T~15:00)

| Job ID | Name | Status | Notes |
|--------|------|--------|-------|
| **44236980** | `dl_lcl` (60 selected) | ~55/60 COMPLETED, 5 FAILED | Failed tasks: 3,8,33,44,53 |
| **44236981** | `tophat_lcl` (60 align) | PENDING (Dependency on 44236980) | Will NOT auto-start — see §3 |
| **44236982** | `star_lcl` (60 align) | PENDING (Dependency on 44236980) | Same dependency issue |
| **44237285** | `dl_lcl_bulk` (719 remaining) | RUNNING ~66/719 | Background accumulation |

### ⚠️ CRITICAL: Alignment jobs will NOT start automatically

Because 5 tasks of 44236980 **FAILED**, the SLURM dependency (`afterok:44236980`)
will **never be satisfied**.  The alignment jobs will stay PENDING forever.

**You must fix this** — see §3 below.

---

## 3. Immediate Next Steps (When You Reconnect)

### Step 1: Resubmit 5 failed downloads

The 5 missing SRR accessions (SLURM array tasks 3, 8, 33, 44, 53):
- `SRR19762770`  (task 3)
- `SRR19762740`  (task 8)
- `SRR19762244`  (task 33)
- `SRR19762790`  (task 44)
- `SRR19762677`  (task 53)

```bash
cd /hpc/home/nfp8/copilot/2026_03_03/lcl_1kgp
sbatch --array=3,8,33,44,53 02_download_lcl.sh
# Note the new job ID: <RETRY_JOBID>
```

Verify completion:
```bash
cd /work/nfp8/LCL_1kgp/fastq
for srr in SRR19762770 SRR19762740 SRR19762244 SRR19762790 SRR19762677; do
  [[ -f "${srr}_1.fastq.gz" && -f "${srr}_2.fastq.gz" ]] && echo "OK: $srr" || echo "MISSING: $srr"
done
```

### Step 2: Cancel stale alignment jobs and resubmit with new dependency

```bash
scancel 44236981 44236982

# Resubmit with dependency on the retry job
sbatch --dependency=afterok:<RETRY_JOBID> 03_align_tophat_lcl.sh
sbatch --dependency=afterok:<RETRY_JOBID> 04_align_star_lcl.sh
```

### Step 3: Verify LCL strandedness after first alignment completes

Once a STAR BAM exists for a 1/1 genotype sample, run:
```bash
module load RSeQC
infer_experiment.py -r /hpc/group/parrishlab/refs/hg38_RefSeq.bed \
  -i /work/nfp8/LCL_1kgp/bam_star/<SAMPLE>_Aligned.sortedByCoord.out.bam
```

MCF7 empirically showed SENSE = R1-forward (opposite of expected NEBNext dUTP).
LCL data is from a different lab (PRJNA851328); strandedness must be verified.

### Step 4: Compute RT_CPM for LCL cohort

Adapt `compute_rt_cpm.py` for the 60 LCL BAMs:
- BAM dir: wherever the alignment scripts output (check `03_align_tophat_lcl.sh` and `04_align_star_lcl.sh`)
- Log dir: same (STAR `_Log.final.out`)
- Genotype metadata: `/work/nfp8/LCL_1kgp/metadata/selected_samples.tsv`
- Group by genotype (0/0, 0/1, 1/1) and compare distributions
- Run the same 3-genotype GMM model from `fit_rt_cpm_model.py`

### Step 5: Generate LCL IGV reports by genotype

Use `scripts/make_igv_report.py` to make grouped strand-separated coverage plots:
- One report for STAR BAMs, one for TopHat2
- Group/label by genotype using `--labels-file`
- Expect 1/1 samples to show the "blue curve" readthrough signal

---

## 4. Key Technical Reference

### Readthrough Metric Definition
- **RT_WINDOW:** `chr15:55958936-55959936` (1 kb downstream of L1 3′ breakpoint)
- **SENSE flags:** `-f 64 -F 16` (R1 required, reverse NOT set = R1 forward = plus-strand RNA)
- **MAPQ filter:** `-q 20` (captures STAR MAPQ=255 and TopHat2 MAPQ=50 uniques)
- **Denominator:** Uniquely mapped reads from `*_Log.final.out` (STAR) or `align_summary.txt` (TopHat2)
- **RT_CPM = window_reads / total_mapped × 1,000,000**

### Hardy-Weinberg Frequencies (1000 Genomes)
- Source: PRJNA851328 genotypes: 1232 (0/0), 1028 (0/1), 242 (1/1), N=2502
- Allele freq: q = 0.3022, p = 0.6978
- HWE weights: p²=0.4870, 2pq=0.4217, q²=0.0913

### Genome References
- Genome: `/hpc/group/parrishlab/refs/GRCh38.primary_assembly.genome.fa`
- Annotation: `/hpc/group/parrishlab/refs/gencode.v45.primary_assembly.annotation.gtf`
- STAR index: `/hpc/group/parrishlab/refs/star_GRCh38_v45_sjdb149/`
- Bowtie2 index: `/hpc/group/parrishlab/refs/bowtie2_GRCh38_primary_assembly`

---

## 5. Data Locations

| Path | Contents | Size |
|------|----------|------|
| `/work/nfp8/MCF7_ctrl/bam/` | STAR BAMs for 4 MCF7 samples | — |
| `/work/nfp8/MCF7_ctrl_tophat/bam/` | TopHat2 BAMs for 4 MCF7 samples | — |
| `/work/nfp8/2026_03_02/GSE226189/bam/` | STAR BAMs for 82 GSE226189 fibroblasts | — |
| `/work/nfp8/LCL_1kgp/fastq/` | LCL downloads (selected + bulk) | ~1 TB |
| `/work/nfp8/LCL_1kgp/metadata/` | Genotypes, selected samples, SRR lists | — |

### ⚠️ /work auto-purge: 75 days — files created 2026-03-03 purge ~2026-05-17

---

## 6. File Inventory (this commit)

```
lcl_1kgp/
├── 01_select_samples.py           # Select 20 per genotype from PRJNA851328
├── 02_download_lcl.sh             # SLURM: download 60 selected (job 44236980)
├── 02b_download_lcl_bulk.sh       # SLURM: download 719 remaining (job 44237285)
├── 03_align_tophat_lcl.sh         # SLURM: TopHat2 alignment (job 44236981)
├── 04_align_star_lcl.sh           # SLURM: STAR alignment (job 44236982)
├── launch_lcl_pipeline.sh         # Orchestration script (submits all above)
├── compute_rt_cpm.py              # RT_CPM for 4 MCF7 ctrl samples
├── compute_rt_cpm_GSE226189.py    # RT_CPM for 82 GSE226189 fibroblasts
├── fit_rt_cpm_model.py            # 3-genotype GMM fit + histogram + LRT
├── GSE226189_rt_cpm.tsv           # 82-sample RT_CPM results
├── GSE226189_rt_cpm_histogram.png # 3-panel figure
├── GSE226189_rt_cpm_histogram.html# HTML wrapper for browser viewing
└── HANDOFF.md                     # This file
```

---

## 7. Conda Environment

All scripts run in: `l1nedd4_py39` (Python 3.9)
- Key packages: numpy, scipy (1.13.1), matplotlib (3.9.4), pysam
- samtools/1.21 available via `module load samtools/1.21`

---

## 8. GitHub Repo

- Repo: `nfparrish/l1-nedd4-readthrough`
- Branch: `main`
- Remote: `https://github.com/nfparrish/l1-nedd4-readthrough.git`
