# Session Handoff — 2026-03-03

> **Last update:** ~2026-03-03T03:00 UTC-5  
> **Next session:** Check progress in ~2-3 hours

---

## Current Pipeline Status

### Active Jobs (as of last check)

```
43870250       index_mcf7  RUNNING  ~1h 5m elapsed / 4h limit  (STAR genome index)
43870677_[1-4] trim_mcf7   RUNNING  All 4 tasks on cores       (Trim Galore with cutadapt 4.9)
43870678_[1-4] align_mcf7  PENDING  Depends on 43870250+43870677
43870679_[1-4] post_mcf7   PENDING  Depends on align
43870680      collect_mcf7 PENDING  Depends on post
43870681      vis_mcf7     PENDING  Depends on collect
43870682      check_mcf7   PENDING  Depends on vis
```

### Expected Timeline

| Task | Expected Duration | Status |
|---|---|---|
| Index (43870250) | 4 hours | ~1h done, ~3h remaining |
| Trim (43870677_[1-4]) | 30 min each | Currently running |
| Align (43870678_[1-4]) | 6 hours | Will start after index+trim |
| Post-align (43870679_[1-4]) | 2 hours | After align |
| Collect (43870680) | 1 hour | After post |
| Vis (43870681) | 30 min | After collect |
| Check (43870682) | 10 min | After vis |

**Full pipeline ETA:** ~13-14 hours from when index started.

---

## What Was Done This Session

### Peer-Review Implementation (✅ COMPLETE)

Applied all 9 scientific recommendations from external review:

1. **Critical:** PCR duplicate marking → `03_post_align.sh` samtools markdup pipeline
2. **Critical:** Strand-specific depth (MAPQ≥255) → handles reverse/forward/unstranded libraries
3. **Critical:** MAPQ filter on depth calculation → prevents multimapping artifacts
4. **Moderate:** Adapter trimming → new `01b_trim_trimgalore.sh` (PE/SE)
5. **Moderate:** featureCounts alongside STAR GeneCounts → both retained for comparison
6. **Moderate:** STAR alignment QC summary → `04_collect_results.sh` Part 0 parses Log.final.out
7. **Moderate:** Background normalization → flanking window `rt_bg_fold` + `rt_bg_log2fold` metrics
8. **Moderate:** Removed stale STAR flags → `--outSAMstrandField intronMotif`, chimeric params
9. Global CPM retained, with background fold-enrichment added

### Infrastructure Setup

1. **Prompt history logging:** `save_prompt_history.sh` + updated `agents.md` with mandatory pre-summary logging
2. **Trim dependency fix:** `cutadapt` version pinned to 4.x (5.x incompatible with Trim Galore 0.4.3)
3. **Dual dependencies:** `launch_pipeline.sh` and `run_test_mcf7.sh` handle both index+trim deps on align

### Updated/Created Files

```
/hpc/home/nfp8/copilot/2026_03_03/
├── 01b_trim_trimgalore.sh           [NEW] Trim Galore adapter trimming
├── 03_post_align.sh                 [UPDATED] markdup, strand depth, featureCounts
├── 04_collect_results.sh            [UPDATED] STAR QC, featureCounts matrix, bg fold
├── 02_align_star.sh                 [UPDATED] Trimmed input preference, stale flags removed
├── 01_build_star_index.sh           [UPDATED] NEDD4 exon + background BED generation
├── launch_pipeline.sh               [UPDATED] Step 1.5 (trim), dual align dep
├── config.sh                        [UPDATED] New modules, trim params, BED paths
├── 00_setup_env.sh                  [UPDATED] cutadapt 4.x pinned
├── agents.md                        [NEW] Storage policy + prompt logging mandate
├── save_prompt_history.sh           [NEW] Session prompt archiver
├── 2026_03_03_AI_prompt_hx.md       [NEW] Prompt history log (append-only)
└── test_mcf7/
    └── run_test_mcf7.sh             [UPDATED] Trim step, dual align dep, post mem 16G
```

---

## Known Issues & Fixes Applied

### ✅ Issue: Trim Galore 0.4.3 requires cutadapt CLI compatible version

**Symptom:** All trim tasks exit code 2: `cutadapt: error: unrecognized arguments: -f /path/to/fastq`

**Root Cause:** cutadapt 5.2 has completely different argument structure than 4.x

**Fix Applied:**
- Downgraded to cutadapt 4.9 in venv: `/hpc/group/parrishlab/envs/l1nedd4_py39/`
- Updated `00_setup_env.sh` to pin `'cutadapt>=3.7,<5'`
- Updated `01b_trim_trimgalore.sh` to pass `--path_to_cutadapt "${VENV_DIR}/bin/cutadapt"`

**Status:** ✅ FIXED — new trim tasks (43870677) running successfully

---

## Monitoring & Next Steps

### When Returning (in ~2-3 hours)

1. **Check job status:**
   ```bash
   squeue -u nfp8 --sort=i --format="%.10i %.20j %.8T %.12M %R"
   ```

2. **Expected progression:**
   - Trim (43870677) should be COMPLETED
   - Align (43870678) should be RUNNING (waiting on index 43870250)
   - If index still RUNNING, alignment will wait
   - If index COMPLETED before trim, alignment will start as soon as trim finishes

3. **If anything FAILED:**
   ```bash
   # Check recent errors
   cat /hpc/home/nfp8/copilot/2026_03_03/test_mcf7/logs/POST_JOB_123.err
   # or for array jobs
   sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed
   ```

4. **Expected final output files** (when complete):
   - `/work/nfp8/2026_03_02/persistent_results/`
     - `readthrough_signal.csv` (rt_mean, rt_cpm, rt_bg_fold, rt_bg_log2fold)
     - `star_alignment_qc.csv` (unique %, multi %, mismatch rate, etc.)
     - `fc_count_matrix.csv` (featureCounts gene counts)
     - `count_matrix.csv` (STAR GeneCounts for comparison)
     - `nedd4_fc_counts.csv` + `nedd4_counts.csv` (NEDD4-specific counts)
     - `tmm_cpm_matrix.csv` (TMM-normalized expression)

### Critical Config Values

| Variable | Value | File |
|---|---|---|
| VENV_DIR | `/hpc/group/parrishlab/envs/l1nedd4_py39` | config.sh |
| cutadapt version | 4.9 (pinned 3.7-<5) | 00_setup_env.sh |
| STAR modules | Trim_galore/0.4.3, Subread/2.0.3, STAR/2.7.11a | config.sh |
| Strandedness | reverse (fr-firststrand, TruSeq) | config.sh |
| Library layout | PE (paired-end) | config.sh |

### Git Status

All changes pushed to: `https://github.com/nfparrish/l1-nedd4-readthrough`

Latest commits:
1. `1378cd7` — Pin cutadapt to 4.x for Trim Galore 0.4.3 compatibility
2. `cbd7dbc` — Add prompt history logging
3. `8a8f3fc` — Fix: add --path_to_cutadapt to Trim Galore
4. `30ff0f8` — Peer-review fixes: markdup, strand depth, trim, featureCounts, STAR QC, background fold norm

---

## Resumption Checklist

- [ ] SSH back in
- [ ] Run `squeue -u nfp8` to check job status
- [ ] If all COMPLETED: check `/work/nfp8/2026_03_02/persistent_results/` for output files
- [ ] If jobs still PENDING or RUNNING: wait and recheck in 30 min
- [ ] If any FAILED: check error logs in `test_mcf7/logs/` 
- [ ] If errors: escalate to agent with job ID and error output
- [ ] Update `2026_03_03_AI_prompt_hx.md` with session continuation prompts

---

## Glossary

- **rt_cpm:** Read-through CPM (original metric, retained)
- **rt_bg_fold:** Read-through / background fold enrichment (raw ratio)
- **rt_bg_log2fold:** log₂ of fold enrichment (for statistical testing)
- **MAPQ≥255:** Uniquely mapped reads (default STAR setting)
- **markdup -s:** Mark (not remove) duplicates; output metrics
- **featureCounts:** HTSeq-like gene quantification; more sensitive than STAR GeneCounts

---

**Session prepared by:** GitHub Copilot (Claude Opus 4)  
**Safe to disconnect:** ✅ YES
