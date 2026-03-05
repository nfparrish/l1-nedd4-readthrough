# L1-NEDD4 Readthrough Pipeline — Handoff, 2026-03-04

**Session closed:** 2026-03-04 ~20:15 EST  
**Operator:** Copilot (Claude Sonnet 4.6)  
**Repo:** `git@github.com:nfparrish/l1-nedd4-readthrough.git`  
**Branch:** `main` — latest commit `ac227d8`

---

## 1. Overview

Pipeline detects L1-element-driven NEDD4 readthrough transcription from bulk RNA-seq.  
Two cohorts in flight:

| Cohort | Samples | Status |
|--------|---------|--------|
| **MCF7** (GSE47006) | 4 ERR samples — positive controls | Post-align RUNNING (tasks 1-3); task 4 pending |
| **GSE226189** | 82 fibroblast FASTQs | Trim RUNNING; full chain queued |

---

## 2. Environment

```
HPC:            Duke DCC, SLURM, account=parrishlab, partition=common
Pipeline dir:   /hpc/home/nfp8/copilot/2026_03_03/
MCF7 config:    /hpc/home/nfp8/copilot/2026_03_03/test_mcf7/config_test_mcf7.sh
Default config: /hpc/home/nfp8/copilot/2026_03_03/config.sh
MCF7 work dir:  /work/nfp8/MCF7_ctrl/
GSE226189 dir:  /work/nfp8/2026_03_02/GSE226189/
Persistent out: /hpc/group/parrishlab/L1_NEDD4/
STAR index:     /hpc/group/parrishlab/refs/star_GRCh38_v45_sjdb149
GENCODE GTF:    /hpc/group/parrishlab/refs/gencode.v45.primary_assembly.annotation.gtf
Panel GTF:      /hpc/group/parrishlab/refs/housekeeping_panel.gtf  ← BUILT (15 MB, 100 genes)
Chr15 GTF:      /hpc/group/parrishlab/refs/gencode.v45.chr15.gtf   ← BUILT
Python venv:    /hpc/group/parrishlab/envs/l1nedd4_py39
RT window:      chr15:55958936-55959936  (1 kb downstream of NEDD4 stop)
```

---

## 3. Pipeline Scripts (current state)

### `00_setup_env.sh`
Creates venv, directories, chr15 GTF subset, housekeeping panel GTF.  
**Fix in this session:** housekeeping panel `PATTERN` used `\|` (invalid ERE) — replaced with `paste -sd'|'`.

### `01b_trim_trimmomatic.sh`
Adapter trimming (Trimmomatic). Sources config via `PIPELINE_CONFIG` env var.

### `02_align_star.sh`
STAR 2-pass alignment (150 bp SJDB). Sources config.

### `03_post_align.sh`  ← most complex, most recent changes
Steps performed per sample:
1. **MarkDup** (Picard) — produces `{SRR}_dedup.bam`
2. **Flagstat** — alignment QC
3. **Plus-strand RT depth** — samtools depth at `RT_WINDOW` (forward coordinates for `reverse`-stranded)
3b. **Minus-strand RT depth** *(NEW)* — sanity control; also captures NEDD4 intron retention from minus-strand reads → `{SRR}_readthrough_minus.depth`
4. **Background depth** — samtools depth over `BACKGROUND_WINDOWS_BED`
5. **Subset BAM** — extract reads overlapping RT_WINDOW for IGV
6. **featureCounts** — strand-aware, controlled by `FC_SCOPE` (see below)
7. **Copy STAR GeneCounts** — retains both quantification outputs

**Fixes in this session:**
- `samtools index` added before `samtools depth -r` (stranded libs only)  
- `|| true` on intermediate temp-BAM `rm -f` lines (NFS `.nfs*` lock handles)  
- **Missing `fi`** at end of FC_SCOPE outer if-block — caused `unexpected end of file` at line 261 for all jobs reaching step 6 (fixed commit `ac227d8`)

### `04_collect_results.sh`
Aggregates all sample metrics into `readthrough_signal.csv`.  
**New columns:** `rt_minus_mean_depth`, `rt_minus_cpm`

### `launch_pipeline.sh`
Orchestrates full pipeline as SLURM dependency chain.  
**Fix in this session:** `submit()` and `submit_dep()` now pass `--export=ALL,PIPELINE_CONFIG="${PIPELINE_CONFIG}"` — previously scripts couldn't find `config.sh` from SLURM's temp dir.

---

## 4. FC_SCOPE Feature (new)

`FC_SCOPE` is set in `config.sh` (default: `GENOME`) or overridden via `--export`:

| Value | GTF used | Speed | Use case |
|-------|----------|-------|----------|
| `SKIP` | none | instant | skip featureCounts |
| `PANEL` | `housekeeping_panel.gtf` (100 genes) | ~2 min | QC normalization |
| `CHR15` | chr15 subset | ~15 min | NEDD4-focused analysis |
| `GENOME` | full GENCODE v45 | ~45 min | complete quantification |

Panel genes cover: ACTB/GAPDH, ribosomes (RPL/RPS), glycolysis, chaperones (HSP90/HSPA), chromatin, HIF, MYC, EGFR/ERBB2, AKT/mTOR/PIK3CA, STAT3/JAK, Wnt, Notch, TFRC, UPR.

---

## 5. Current SLURM Queue (at session close)

```
MCF7 cohort:
  43895927_[1-3]   post_align    RUNNING  ~2h elapsed (steps 1-4 of 7 likely)
  43896864         mcf7_panel_resubmit  PENDING afterany:43895927
                   → deletes GENOME featureCounts for ERR973734/35/28
                   → resubmits tasks 1-3 with FC_SCOPE=PANEL
  43897191_[4]     post_align    PENDING (Priority)
                   → ERR973729, FC_SCOPE=PANEL, no dependency (GTF already built)

GSE226189 cohort:
  43897100_[1-82]  trim_tm       RUNNING (13 concurrent, ~82 tasks total)
  43897101_[1-82]  star_align    PENDING afterok:43897100
  43897102_[1-82]  post_align    PENDING afterok:43897101
  43897103         collect_results  PENDING afterok:43897102
  43897104         visualize_rt  PENDING afterok:43897103
```

**Monitor with:**
```bash
squeue -u nfp8 --format="%.10i %.22j %.8T %.10M %R" --sort=i
```

**If any task fails, check logs:**
- MCF7: `/work/nfp8/MCF7_ctrl/logs/`
- GSE226189: `/work/nfp8/2026_03_02/GSE226189/logs/`
- Pipeline logs: `/hpc/home/nfp8/copilot/2026_03_03/logs/`

---

## 6. MCF7 Biology Notes

ERR973729 is the **non-readthrough negative within the positive-control set**. All other 3 MCF7 samples (ERR973734, ERR973735, ERR973728) show plus-strand readthrough depth; ERR973729 shows zero plus-strand depth in `RT_WINDOW`. This is biologically real — the 112 reads appearing in `samtools view` over the window are split-junction NEDD4 exonic reads whose junctions span the window through long introns (no M bases inside the window). The new step 3b minus-strand depth file confirms NEDD4 intronic reads on the minus strand (as expected for a reverse-stranded library).

---

## 7. Commits This Session

| Commit | Description |
|--------|-------------|
| `d24c6b9` | Throttle post-align array to %16 |
| `0ad9529` | Fix missing samtools index before depth for stranded libs |
| `3df18b4` | Add `\|\| true` to stranded temp-BAM `rm -f` (NFS locks) |
| `901c3e8` | Add step 3b (minus-strand RT depth) + FC_SCOPE feature |
| `0b68037` | Fix ERE pattern `\|` → `paste \|` and pass PIPELINE_CONFIG in launch_pipeline.sh |
| `ac227d8` | Fix missing `fi` closing outer FC_SCOPE if-block in step 6 |

---

## 8. Next Steps After Jobs Complete

1. **MCF7 PANEL featureCounts** — once 43895927 + 43896864 + 43897191 all complete, run `04_collect_results.sh` with `PIPELINE_CONFIG=config_test_mcf7.sh` to generate the updated `readthrough_signal.csv` with minus-strand columns.

2. **GSE226189 results** — jobs are fully chained; `collect_results` and `visualize` will auto-run. Results land in `/hpc/group/parrishlab/L1_NEDD4/GSE226189/`.

3. **Strandedness validation for GSE226189** — config uses `unstranded`; confirm with `infer_experiment.py` from RSeQC on a completed BAM before interpreting featureCounts direction. If stranded, re-run post_align with corrected `STRANDEDNESS`.

4. **MultiQC** — run `multiqc` across the GSE226189 logs once alignment completes (or add as a step after `04_collect_results.sh`).

5. **FC_SCOPE for GSE226189** — currently set to `GENOME` in `config.sh`. Consider switching to `PANEL` for faster iteration, then re-running with `GENOME` for publication-quality counts.

---

## 9. Key File Paths Quick Reference

```bash
# Pipeline
/hpc/home/nfp8/copilot/2026_03_03/03_post_align.sh    # main workhorse
/hpc/home/nfp8/copilot/2026_03_03/config.sh           # GSE226189 defaults
/hpc/home/nfp8/copilot/2026_03_03/test_mcf7/config_test_mcf7.sh  # MCF7

# Outputs
/work/nfp8/MCF7_ctrl/coverage/*_readthrough.depth      # plus-strand RT
/work/nfp8/MCF7_ctrl/coverage/*_readthrough_minus.depth # minus-strand RT (new)
/work/nfp8/MCF7_ctrl/counts/*_featureCounts.txt        # gene counts
/work/nfp8/MCF7_ctrl/results/readthrough_signal.csv    # summary table

# Reference files
/hpc/group/parrishlab/refs/housekeeping_panel.gtf      # 100-gene PANEL
/hpc/group/parrishlab/refs/gencode.v45.chr15.gtf       # chr15 subset
```
