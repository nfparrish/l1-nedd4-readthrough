#!/bin/bash
#SBATCH --job-name=collect_results
#SBATCH --output=/work/nfp8/2026_03_02/GSE226189/logs/collect_%j.out
#SBATCH --error=/work/nfp8/2026_03_02/GSE226189/logs/collect_%j.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=02:00:00
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"

mkdir -p "$RESULTS_DIR" "$PERSISTENT_RESULTS"

# Use venv python directly (no activate script needed)
VENV_PYTHON="${VENV_DIR}/bin/python3"

echo "[$(date -Is)] Starting results collection"

# ===========================================================================
# Part 0: STAR alignment QC summary (parse Log.final.out for each sample)
# ===========================================================================
${VENV_PYTHON} - <<'PYEOF'
import os, sys, re
import pandas as pd

BAM_DIR     = os.environ["BAM_DIR"]
SRR_LIST    = os.environ["SRR_LIST"]
RESULTS_DIR = os.environ["RESULTS_DIR"]

with open(SRR_LIST) as f:
    srrs = [line.strip() for line in f if line.strip()]

fields = {
    "total_reads":          r"Number of input reads \|\s+(\d+)",
    "pct_unique_mapped":    r"Uniquely mapped reads % \|\s+([\d.]+)%",
    "pct_multi_mapped":     r"% of reads mapped to multiple loci \|\s+([\d.]+)%",
    "pct_unmapped_short":   r"% of reads unmapped: too short \|\s+([\d.]+)%",
    "num_splices":          r"Number of splices: Total \|\s+(\d+)",
    "mismatch_rate":        r"Mismatch rate per base, % \|\s+([\d.]+)%",
}

rows = []
for srr in srrs:
    log_path = os.path.join(BAM_DIR, f"{srr}_Log.final.out")
    if not os.path.exists(log_path):
        print(f"  SKIP {srr}: no STAR Log.final.out", file=sys.stderr)
        continue
    with open(log_path) as fh:
        text = fh.read()
    row = {"srr": srr}
    for key, pat in fields.items():
        m = re.search(pat, text)
        row[key] = float(m.group(1)) if m else None
    rows.append(row)

df = pd.DataFrame(rows)
out = os.path.join(RESULTS_DIR, "star_alignment_qc.csv")
df.to_csv(out, index=False)
print(f"STAR QC: {len(df)} samples -> {out}")
if not df.empty:
    print(df[["srr","total_reads","pct_unique_mapped","pct_multi_mapped"]].to_string(index=False))
PYEOF

echo "[$(date -Is)] STAR QC summary done"

# ===========================================================================
# Part A: Read-through signal (CPM) + background-normalized fold enrichment
# ===========================================================================
${VENV_PYTHON} - <<'PYEOF'
import os, sys, glob
import pandas as pd
import numpy as np

COVERAGE_DIR = os.environ["COVERAGE_DIR"]
SRR_LIST     = os.environ["SRR_LIST"]
RESULTS_DIR  = os.environ["RESULTS_DIR"]
RT_WINDOW    = os.environ["RT_WINDOW"]

with open(SRR_LIST) as f:
    srrs = [line.strip() for line in f if line.strip()]


def mean_depth(path: str):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", header=None, names=["chr", "pos", "depth"])
    return df["depth"].mean() if len(df) > 0 else 0.0

results = []
paper_results = []
for srr in srrs:
    rt_file       = os.path.join(COVERAGE_DIR, f"{srr}_readthrough.depth")
    rt_minus_file = os.path.join(COVERAGE_DIR, f"{srr}_readthrough_minus.depth")
    paper_sense_file = os.path.join(COVERAGE_DIR, f"{srr}_paper_downstream_sense.depth")
    paper_antisense_file = os.path.join(COVERAGE_DIR, f"{srr}_paper_upstream_antisense.depth")
    fs_file = os.path.join(COVERAGE_DIR, f"{srr}_flagstat.txt")
    bg_file = os.path.join(COVERAGE_DIR, f"{srr}_background.depth")

    if not os.path.exists(rt_file) or not os.path.exists(fs_file):
        print(f"  SKIP {srr}: missing depth or flagstat", file=sys.stderr)
        continue

    # Read-through mean depth (plus strand for stranded; all reads for unstranded)
    rt = pd.read_csv(rt_file, sep="\t", header=None, names=["chr","pos","depth"])
    rt_mean = rt["depth"].mean() if len(rt) > 0 else 0.0
    rt_max  = int(rt["depth"].max()) if len(rt) > 0 else 0

    # Minus-strand depth (sanity control / NEDD4 intron retention proxy)
    same_window_minus_mean = None
    same_window_minus_cpm  = None
    if os.path.exists(rt_minus_file):
        rtm = pd.read_csv(rt_minus_file, sep="\t", header=None, names=["chr","pos","depth"])
        same_window_minus_mean = rtm["depth"].mean() if len(rtm) > 0 else 0.0

    # Total primary mapped reads from flagstat (CPM denominator for RT windows)
    total_mapped = 0
    with open(fs_file) as fh:
        lines = fh.readlines()
    for line in lines:
        if "primary mapped (" in line:
            total_mapped = int(line.split()[0])
            break
    if total_mapped == 0:   # fallback for older samtools
        for line in lines:
            if "mapped (" in line and "primary" not in line:
                total_mapped = int(line.split()[0])
                break

    # Paper-style CPM denominator: total R1 reads with MAPQ>=20 (matches Philippe et al.)
    r1_count_file = os.path.join(COVERAGE_DIR, f"{srr}_r1_mapq20_count.txt")
    r1_mapped = 0
    if os.path.exists(r1_count_file):
        try:
            r1_mapped = int(open(r1_count_file).read().strip())
        except ValueError:
            pass
    # Fall back to total_mapped if R1 count not available
    paper_denom = r1_mapped if r1_mapped > 0 else total_mapped

    rt_cpm = (rt_mean / total_mapped) * 1e6 if total_mapped > 0 else 0.0
    if same_window_minus_mean is not None and total_mapped > 0:
        same_window_minus_cpm = (same_window_minus_mean / total_mapped) * 1e6

    paper_downstream_mean = mean_depth(paper_sense_file)
    paper_upstream_antisense_mean = mean_depth(paper_antisense_file)
    paper_downstream_cpm = None
    paper_upstream_antisense_cpm = None
    if paper_downstream_mean is not None and paper_denom > 0:
        paper_downstream_cpm = (paper_downstream_mean / paper_denom) * 1e6
    if paper_upstream_antisense_mean is not None and paper_denom > 0:
        paper_upstream_antisense_cpm = (paper_upstream_antisense_mean / paper_denom) * 1e6

    # Background flanking window mean depth (NEDD4 exons excluded)
    bg_mean = None
    rt_bg_fold = None
    if os.path.exists(bg_file):
        bg = pd.read_csv(bg_file, sep="\t", header=None, names=["chr","pos","depth"])
        if len(bg) > 0:
            bg_mean = bg["depth"].mean()
            # Fold enrichment of read-through over local background
            # log2 scale friendlier for diff testing; use raw fold in CSV
            rt_bg_fold = rt_mean / bg_mean if bg_mean > 0 else None

    results.append({
        "srr":                 srr,
        "rt_window":           RT_WINDOW,
        "rt_mean_depth":       round(rt_mean, 4),
        "rt_max_depth":        rt_max,
        "total_mapped":        total_mapped,
        "rt_cpm":              round(rt_cpm, 6),
        "same_window_minus_mean_depth": round(same_window_minus_mean, 4) if same_window_minus_mean is not None else None,
        "same_window_minus_cpm":        round(same_window_minus_cpm, 6) if same_window_minus_cpm is not None else None,
        "bg_mean_depth":       round(bg_mean, 4) if bg_mean is not None else None,
        "rt_bg_fold":          round(rt_bg_fold, 4) if rt_bg_fold is not None else None,
        "rt_bg_log2fold":      round(np.log2(rt_bg_fold), 4) if rt_bg_fold is not None and rt_bg_fold > 0 else None,
    })

    paper_results.append({
        "srr": srr,
        "paper_downstream_window": os.environ.get("PAPER_DOWNSTREAM_WINDOW", RT_WINDOW),
        "paper_upstream_window": os.environ.get("PAPER_UPSTREAM_WINDOW", ""),
        "paper_downstream_sense_mean_depth": round(paper_downstream_mean, 4) if paper_downstream_mean is not None else None,
        "paper_downstream_sense_cpm": round(paper_downstream_cpm, 6) if paper_downstream_cpm is not None else None,
        "paper_upstream_antisense_mean_depth": round(paper_upstream_antisense_mean, 4) if paper_upstream_antisense_mean is not None else None,
        "paper_upstream_antisense_cpm": round(paper_upstream_antisense_cpm, 6) if paper_upstream_antisense_cpm is not None else None,
    })

df = pd.DataFrame(results)
out_csv = os.path.join(RESULTS_DIR, "readthrough_signal.csv")
df.to_csv(out_csv, index=False)
print(f"Wrote {len(df)} samples to {out_csv}")
print(df[["srr","rt_mean_depth","rt_cpm","same_window_minus_mean_depth","same_window_minus_cpm","bg_mean_depth","rt_bg_fold"]].to_string(index=False))

paper_df = pd.DataFrame(paper_results)
paper_out_csv = os.path.join(RESULTS_DIR, "paper_style_signal.csv")
paper_df.to_csv(paper_out_csv, index=False)
print(f"Wrote {len(paper_df)} samples to {paper_out_csv}")
if not paper_df.empty:
    print(paper_df[["srr","paper_downstream_sense_mean_depth","paper_upstream_antisense_mean_depth"]].to_string(index=False))
PYEOF

echo "[$(date -Is)] Read-through CSV done"

# ===========================================================================
# Part B1: featureCounts count matrix
# ===========================================================================
${VENV_PYTHON} - <<'PYEOF'
import os, sys
import pandas as pd

COUNTS_DIR  = os.environ["COUNTS_DIR"]
SRR_LIST    = os.environ["SRR_LIST"]
RESULTS_DIR = os.environ["RESULTS_DIR"]
FC_SCOPE    = os.environ.get("FC_SCOPE", "GENOME")

with open(SRR_LIST) as f:
    srrs = [line.strip() for line in f if line.strip()]

counts = {}
for srr in srrs:
    fc_file = os.path.join(COUNTS_DIR, f"{srr}_featureCounts_{FC_SCOPE}.txt")
    if not os.path.exists(fc_file):
        print(f"  SKIP {srr}: no featureCounts file", file=sys.stderr)
        continue
    # featureCounts output: header=2 lines, cols: Geneid,Chr,Start,End,Strand,Length,<bam>
    df = pd.read_csv(fc_file, sep="\t", comment="#", header=0, index_col=0)
    # Last column is the count column (BAM path as header)
    count_col = df.columns[-1]
    counts[srr] = df[count_col]

if not counts:
    print("WARNING: no featureCounts files found — skipping featureCounts matrix", file=sys.stderr)
else:
    mat = pd.DataFrame(counts)
    out_csv = os.path.join(RESULTS_DIR, "fc_count_matrix.csv")
    mat.to_csv(out_csv)
    print(f"featureCounts matrix: {mat.shape[0]} genes x {mat.shape[1]} samples -> {out_csv}")
    nedd4_ids = [g for g in mat.index if "ENSG00000069869" in str(g)]
    if nedd4_ids:
        nedd4 = mat.loc[nedd4_ids].T
        nedd4.to_csv(os.path.join(RESULTS_DIR, "nedd4_fc_counts.csv"))
        print(f"NEDD4 featureCounts saved ({len(nedd4_ids)} gene IDs)")
PYEOF

echo "[$(date -Is)] featureCounts matrix done"

# ===========================================================================
# Part B2: STAR GeneCounts count matrix (retained for cross-validation)
# ===========================================================================
${VENV_PYTHON} - <<'PYEOF'
import os, sys
import pandas as pd

COUNTS_DIR   = os.environ["COUNTS_DIR"]
SRR_LIST     = os.environ["SRR_LIST"]
RESULTS_DIR  = os.environ["RESULTS_DIR"]
STRANDEDNESS = os.environ.get("STRANDEDNESS", "unstranded")

with open(SRR_LIST) as f:
    srrs = [line.strip() for line in f if line.strip()]

col_map  = {"unstranded": "unstranded", "forward": "sense", "reverse": "antisense"}
use_col  = col_map.get(STRANDEDNESS, "unstranded")

counts = {}
for srr in srrs:
    tab = os.path.join(COUNTS_DIR, f"{srr}_ReadsPerGene.out.tab")
    if not os.path.exists(tab):
        print(f"  SKIP {srr}: no ReadsPerGene file", file=sys.stderr)
        continue
    df = pd.read_csv(tab, sep="\t", header=None,
                     names=["gene_id","unstranded","sense","antisense"],
                     skiprows=4)
    counts[srr] = df.set_index("gene_id")[use_col]

if not counts:
    print("WARNING: no STAR GeneCounts files found", file=sys.stderr)
else:
    mat = pd.DataFrame(counts)
    out_csv = os.path.join(RESULTS_DIR, "count_matrix.csv")
    mat.to_csv(out_csv)
    print(f"STAR GeneCounts matrix: {mat.shape[0]} genes x {mat.shape[1]} samples -> {out_csv}")
    nedd4_ids = [g for g in mat.index if "ENSG00000069869" in str(g)]
    if nedd4_ids:
        mat.loc[nedd4_ids].T.to_csv(os.path.join(RESULTS_DIR, "nedd4_counts.csv"))
        print(f"NEDD4 GeneCounts saved ({len(nedd4_ids)} gene IDs)")
PYEOF

echo "[$(date -Is)] STAR GeneCounts matrix done"

# ===========================================================================
# Part C: TMM normalization (R / edgeR) — uses featureCounts matrix when available
# ===========================================================================
module load "$MOD_R"
export R_LIBS_USER="${HOME}/R/libs"
mkdir -p "${R_LIBS_USER}"

Rscript - <<'REOF'
r_libs <- Sys.getenv("R_LIBS_USER", unset=file.path(Sys.getenv("HOME"), "R", "libs"))
dir.create(r_libs, showWarnings=FALSE, recursive=TRUE)
# Include system R library so module-installed Bioconductor packages are visible
sys_lib <- "/admin/apps/rhel9/R-4.4.3/lib64/R/library"
.libPaths(c(r_libs, sys_lib, .libPaths()))

suppressPackageStartupMessages({
    if (!requireNamespace("edgeR", quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
        BiocManager::install("edgeR", ask=FALSE, quiet=TRUE)
    }
    library(edgeR)
})

results_dir <- Sys.getenv("RESULTS_DIR")

# Prefer featureCounts matrix; fall back to STAR GeneCounts
fc_file     <- file.path(results_dir, "fc_count_matrix.csv")
star_file   <- file.path(results_dir, "count_matrix.csv")
counts_file <- if (file.exists(fc_file)) fc_file else star_file

if (!file.exists(counts_file)) stop("No count matrix found for TMM normalization")
cat(sprintf("TMM input: %s\n", counts_file))

counts <- read.csv(counts_file, row.names=1, check.names=FALSE)
cat(sprintf("Loaded: %d genes x %d samples\n", nrow(counts), ncol(counts)))

min_samples <- min(3, ncol(counts))
keep <- rowSums(counts > 0) >= min_samples
counts_filt <- counts[keep, ]
cat(sprintf("After low-expression filter: %d genes\n", nrow(counts_filt)))

dge <- DGEList(counts=counts_filt)
dge <- calcNormFactors(dge, method="TMM")

write.csv(dge$samples, file.path(results_dir, "tmm_norm_factors.csv"))
cpm_mat <- cpm(dge, normalized.lib.sizes=TRUE)
write.csv(cpm_mat, file.path(results_dir, "tmm_cpm_matrix.csv"))
cat(sprintf("TMM-CPM matrix: %d genes x %d samples saved\n", nrow(cpm_mat), ncol(cpm_mat)))
REOF

echo "[$(date -Is)] TMM normalization done"

# ===========================================================================
# Part D: Copy key results to persistent storage
# ===========================================================================
for f in readthrough_signal.csv star_alignment_qc.csv \
          count_matrix.csv fc_count_matrix.csv \
          nedd4_counts.csv nedd4_fc_counts.csv \
          tmm_norm_factors.csv tmm_cpm_matrix.csv; do
    src="${RESULTS_DIR}/${f}"
    [[ -f "$src" ]] && cp "$src" "${PERSISTENT_RESULTS}/"
done

echo "[$(date -Is)] Results copied to persistent storage: ${PERSISTENT_RESULTS}"
echo "[$(date -Is)] Collection complete"
