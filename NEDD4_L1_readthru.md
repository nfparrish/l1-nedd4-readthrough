# L1-NEDD4 Read-Through Transcription Detection

**Last updated:** 2026-03-02  
**Author:** nfp8 (Parrish Lab, Duke)  
**Platform:** Duke DCC (SLURM HPC)

---

## 1. Scientific Background

### 1.1 NEDD4 and LINE-1

NEDD4 (Neural precursor cell Expressed Developmentally Down-regulated 4) is an E3
ubiquitin ligase on chromosome 15 (GRCh38: chr15:55,826,917–55,993,660, **minus
strand**). It regulates membrane protein turnover via ubiquitin-mediated endocytosis and
is implicated in neurodevelopment, immune signaling, and cancer.

A full-length L1HS LINE-1 retrotransposon (~6 kb) is inserted **antisense** into NEDD4
intron 19. The insertion breakpoint spans **chr15:55,958,934–55,958,936** on the plus
strand. This L1 element retains an intact internal sense promoter that can drive
transcription through (read past) its 3′ end into flanking genomic sequence — a
phenomenon called **L1 read-through transcription**.

### 1.2 Read-Through Window

When L1's sense promoter fires, the resulting transcript extends beyond the L1 3′ poly-A
signal. Since the L1 is on the plus strand and NEDD4 is on the minus strand, L1
read-through transcripts travel in the **antisense** direction relative to NEDD4 mRNA.

The **read-through detection window** is defined as the 1 kb region immediately
downstream (plus-strand sense) of the L1 3′ breakpoint:

- **Window:** chr15:55,958,936–55,959,936 (GRCh38, plus strand)
- **Rationale:** Read-through signal decays with distance from L1. A 1 kb window
  captures the strongest signal while minimizing intronic noise.

### 1.3 Proxy SNP

rs16976600 (chr15:55,862,437) is a proxy SNP in the NEDD4 locus associated with somatic
L1 retrotransposition. It is used for **annotation only** and does **not** define the
read-through window placement.

---

## 2. What We Are Measuring

### 2.1 The Core Question

**Is the L1 element at NEDD4 intron 19 producing read-through transcripts in a given
sample, and how strong is the signal?**

The 1 kb window sits in intronic sequence within NEDD4. Because NEDD4 is on the minus
strand, legitimate NEDD4 transcription does not produce plus-strand reads here. The only
expected source of plus-strand coverage in this window is L1 sense-promoter-driven
read-through. (For unstranded libraries both strands are counted, but true read-through
signal dominates the background.)

### 2.2 Read-Through CPM (RT CPM)

For each sample:

$$
\text{RT CPM} = \frac{\bar{d}_{\text{window}}}{N_{\text{primary mapped}}} \times 10^6
$$

Where:
- $\bar{d}_{\text{window}}$ = mean per-base read depth across the 1 kb read-through
  window (from `samtools depth -a -d 0`)
- $N_{\text{primary mapped}}$ = total primary mapped reads (from `samtools flagstat`,
  the "primary mapped" line)

**Key design choices:**
- Normalization uses **total primary mapped reads**, not NEDD4-specific exon counts. This
  avoids circular dependence (NEDD4 expression is itself modulated by the L1 insertion).
- `samtools depth -a` forces output at all positions in the window, including zero-depth
  positions, so the mean reflects true signal density.
- `-d 0` removes the default depth cap of 8000, relevant for deeply sequenced libraries.

### 2.3 TMM Normalization

A secondary normalization using Trimmed Mean of M-values (edgeR) is computed across all
expressed genes. TMM factors correct for composition bias between libraries. This is
applied to gene-level count matrices but is **not** folded back into RT CPM — the two
metrics are independent:

| Metric | Purpose | Denominator |
|--------|---------|-------------|
| RT CPM | Read-through signal intensity | Primary mapped reads |
| TMM-CPM | Gene expression levels | TMM-adjusted library size |

---

## 3. Expected Results: Active vs. Silent

### 3.1 Sample with active L1-NEDD4 read-through

This individual carries the L1HS insertion and its promoter is unmethylated/active.

**Typical metrics:**
| Field | Value |
|-------|-------|
| rt_mean_depth | ~8–15 |
| rt_max_depth | ~30–50 |
| total_mapped | ~40,000,000 |
| rt_cpm | ~0.20–0.35 |

**Coverage profile across the 1 kb window:**
```
Depth
  40 │██
     │████
  25 │██████
     │█████████
  10 │█████████████████
     │██████████████████████████████▄▂▁▁▁▁▁▁
   0 └──────────────────────────────────────────→
     55,958,936                              55,959,936
     ↑ L1 3′ break                    (plus strand →)
```

The **signature of genuine read-through**: a left-to-right decay gradient. RNA polymerase
initiating from the L1 sense promoter reads through the 3′ poly-A, with processivity
dropping as distance increases. Some transcripts terminate at 200 bp, others at 500 bp,
a few reach 1 kb. This produces the characteristic "waterfall" profile — peak depth at
the breakpoint, monotonically decreasing.

### 3.2 Sample without read-through

Either this individual lacks the L1 insertion at this locus (it's polymorphic — absent in
many people) or carries it but the promoter is fully methylated.

**Typical metrics:**
| Field | Value |
|-------|-------|
| rt_mean_depth | 0.0–0.1 |
| rt_max_depth | 0–2 |
| total_mapped | ~50,000,000 |
| rt_cpm | 0.000000–0.002 |

**Coverage profile:**
```
Depth
  40 │
  25 │
  10 │
     │
   0 │▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁
     └──────────────────────────────────────────→
     55,958,936                              55,959,936
```

Flat at zero. An occasional depth-1 position from a mismapped read, but no gradient, no
spatial structure.

### 3.3 Cohort-Level Interpretation

- **Strip plot:** Most fibroblast samples cluster near zero (non-carriers or silenced).
  A handful of elevated points indicate active carriers. If the L1 insertion frequency is
  ~5–15% in the population and penetrance is partial, expect 3–10 positive samples in an
  82-sample cohort.
- **Aggregated coverage profile:** The mean line is low (diluted by negatives) but should
  still show a slight decay gradient. The upper IQR bound reveals the positive tail.
- **Ranked bar plot:** A sharp knee separating a few high bars from baseline noise. The
  inflection point tells you how many samples have active L1-NEDD4 transcription.

### 3.4 MCF7 as Positive Control

MCF7 is a breast adenocarcinoma cell line with global L1 hypomethylation. RT CPM should
be clearly positive — likely well above typical fibroblast levels — with a strong,
unambiguous decay gradient. If the pipeline fails to detect signal in MCF7, something is
wrong. If it succeeds, we have a calibrated true positive before interpreting the subtler
primary fibroblast data.

---

## 4. Datasets

### 4.1 GSE226189 (Primary Analysis)

- **Title:** Aging effects on primary human skin fibroblasts
- **BioProject:** PRJNA939148
- **Samples:** 82 paired-end (NovaSeq 6000, 150 bp)
- **Library:** Unstranded (verify with `infer_experiment.py`)
- **SRR list:** `/work/nfp8/2026_03_02/GSE226189/SRR_list.txt`
- **Scientific question:** Is L1-NEDD4 read-through detectable in primary fibroblasts,
  and does signal vary with donor age?

### 4.2 MCF7 Control (Pipeline Validation)

- **Accessions:** ERR973734, ERR973735 (ENA, Project PRJEB8975)
- **Cell line:** MCF7 (breast adenocarcinoma)
- **Samples:** 2 paired-end (100 bp)
- **Library:** TruSeq Stranded Total RNA → `STRANDEDNESS="reverse"` (fr-firststrand / RF)
- **Purpose:** Positive control. Validates pipeline end-to-end.
- **Expected result:** Both samples should show measurable RT CPM > 0 with decay profile.

### 4.3 Future: GSE51518

- **Description:** 10 donors, single-end old/young skin fibroblasts
- **Status:** Not yet on this cluster; planned follow-up
- **Config change needed:** `LIBRARY_LAYOUT="SE"`, new SRR list, verify strandedness

---

## 5. Pipeline Architecture

All scripts are in `/hpc/home/nfp8/copilot/2026_03_03/`. Every pipeline script sources
its config via `${PIPELINE_CONFIG:-config.sh}`, making the pipeline reusable across
datasets by exporting a different config file.

### 5.1 Script Inventory

| Step | Script | SLURM Type | Description |
|------|--------|------------|-------------|
| — | `config.sh` | — | Central config (paths, coords, modules, resources) |
| 0 | `00_setup_env.sh` | Single job | Creates directories + Python venv (pandas, matplotlib, seaborn) |
| 1 | `01_build_star_index.sh` | Single job (48G, 12 CPU) | Downloads GRCh38 + GENCODE v45 GTF, builds STAR index |
| 2 | `02_align_star.sh` | Array job (36G, 8 CPU) | STAR 2-pass alignment, BAM sort + index, GeneCounts |
| 3 | `03_post_align.sh` | Array job (8G) | flagstat, RT window depth (-d 0), BAM subset (±1 Mb), count copy |
| 4 | `04_collect_results.sh` | Single job (16G) | Assemble count matrix, compute RT CPM, TMM normalization (R/edgeR) |
| 5 | `05_visualize.sh` + `.py` | Single job | Strip plot, coverage profile, ranked bar plot, summary stats |
| — | `launch_pipeline.sh` | Orchestrator | Submits all steps with `--dependency=afterok` chains |

### 5.2 Test Harness (MCF7)

| Script | Purpose |
|--------|---------|
| `test_mcf7/config_test_mcf7.sh` | Config override for MCF7 (inherits parent config) |
| `test_mcf7/ERR_list.txt` | 2 ENA accessions |
| `test_mcf7/00_download_mcf7.sh` | Downloads from ENA FTP (array 1-2) |
| `test_mcf7/run_test_mcf7.sh` | End-to-end launcher with SLURM dependency chain |
| `test_mcf7/06_sanity_check.sh` | Automated PASS/FAIL assertions on RT CPM |

### 5.3 Dependency Chain

```
download → [setup_env → build_index →] align (array) → post_align (array) → collect → visualize → [sanity_check]
```

Steps in brackets are skippable if resources already exist (venv, STAR index).

---

## 6. Storage Layout

### 6.1 Transient (auto-purged after 75 days)

```
/work/nfp8/<dataset>/
  ├── fastq/          # Raw FASTQ downloads
  ├── bam/            # Full-genome BAMs + STAR outputs
  ├── coverage/       # flagstat, _readthrough.depth files
  ├── counts/         # ReadsPerGene.out.tab per sample
  ├── results/        # readthrough_signal.csv, count_matrix.csv, TMM outputs
  └── logs/           # SLURM stdout/stderr
```

### 6.2 Persistent (no purge)

```
/hpc/group/parrishlab/
  ├── refs/                            # GRCh38 genome, GENCODE GTF, STAR index (~30 GB)
  ├── envs/l1nedd4_py39/              # Python 3.9 venv
  └── L1_NEDD4/
      ├── bam_nedd4_locus/<dataset>/   # Subset BAMs (NEDD4 ±1 Mb per sample)
      └── results/<dataset>/           # Final CSVs + figures/
```

### 6.3 Scripts

```
/hpc/home/nfp8/copilot/2026_03_03/    # All pipeline scripts (persistent home dir)
```

**Important:** `/hpc/group/parrishlab` and `/work` are writable **only from compute
nodes**. Do NOT use `/hpc/group/kolab` (storage quota full).

---

## 7. Key Genomic Coordinates (GRCh38)

| Feature | Coordinates | Strand |
|---------|-------------|--------|
| NEDD4 gene | chr15:55,826,917–55,993,660 | Minus |
| L1HS insertion breakpoint | chr15:55,958,934–55,958,936 | Plus |
| Read-through window (1 kb) | chr15:55,958,936–55,959,936 | Plus |
| BAM subset region (±1 Mb) | chr15:54,958,936–56,959,936 | — |
| Proxy SNP rs16976600 | chr15:55,862,437 | — |
| NEDD4 Ensembl ID | ENSG00000069869 | — |

---

## 8. Software Environment

| Tool | DCC Module / Source | Version |
|------|---------------------|---------|
| STAR | `STAR/2.7.11a` | 2.7.11a |
| samtools | `samtools/1.21` | 1.21 |
| R | `R/4.4.3` | 4.4.3 |
| edgeR | BiocManager install | ≥4.0 |
| FastQC | `FastQC/0.11.7` | 0.11.7 |
| Python | System `/usr/bin/python3` | 3.9.25 |
| pandas/matplotlib/seaborn | venv pip install | latest |

---

## 9. Outputs

| File | Description |
|------|-------------|
| `readthrough_signal.csv` | Per-sample: SRR, rt_mean_depth, rt_max_depth, total_mapped, rt_cpm |
| `count_matrix.csv` | Gene × sample raw counts from STAR GeneCounts |
| `tmm_norm_factors.csv` | Per-sample TMM normalization factors (edgeR) |
| `tmm_cpm_matrix.csv` | TMM-normalized CPM for all expressed genes |
| `nedd4_counts.csv` | NEDD4 (ENSG00000069869) counts extracted separately |
| `*_rt_cpm_strip.png` | Strip plot of RT CPM across all samples |
| `*_rt_coverage_profile.png` | Mean ± IQR per-base depth across RT window |
| `*_rt_cpm_ranked.png` | Ranked bar plot colored by median threshold |
| `*_summary_stats.txt` | Descriptive statistics for RT CPM |

---

## 10. Known Issues and Future Work

### 10.1 Open Questions
- **GSE226189 strandedness:** Configured as "unstranded" based on NovaSeq default; verify
  empirically with `infer_experiment.py` on a completed BAM.
- **Strand-specific depth:** `samtools depth` counts reads on both strands. For stranded
  libraries, consider filtering by strand flag to isolate plus-strand reads only, matching
  the L1 transcription direction.
- **MCF7 genomic rearrangements:** MCF7 has massive structural variation. Confirm the
  L1-NEDD4 read-through signal by manual IGV inspection of subset BAMs.

### 10.2 Planned Extensions
- `infer_experiment.py` (RSeQC) strandedness check as a pipeline pre-step
- Strand-specific depth extraction for stranded libraries
- Integration of GSE51518 (old vs. young fibroblasts, single-end)
- Differential RT CPM analysis between age groups (Wilcoxon / permutation test)
- Chimeric junction analysis from STAR `Chimeric.out.junction` files

---

## 11. How to Run

```bash
# === MCF7 positive-control test (2 samples) ===
bash /hpc/home/nfp8/copilot/2026_03_03/test_mcf7/run_test_mcf7.sh

# === Full GSE226189 pipeline (82 samples) ===
bash /hpc/home/nfp8/copilot/2026_03_03/launch_pipeline.sh --after <DOWNLOAD_JOBID>

# Skip env + index if already built:
bash /hpc/home/nfp8/copilot/2026_03_03/launch_pipeline.sh --skip-setup --skip-index --after <JOBID>

# Resume from a specific step:
bash /hpc/home/nfp8/copilot/2026_03_03/launch_pipeline.sh --from 3

# Monitor:
squeue -u nfp8 --sort=i
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS
```
