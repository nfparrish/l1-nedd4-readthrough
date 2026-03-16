#!/usr/bin/env python3
"""
05_visualize.py — L1-NEDD4 Read-Through Visualization
======================================================
Run within the venv (pandas, matplotlib, seaborn required).
Generates figures:
  1. Strip-plot of read-through CPM across samples
  2. Per-base coverage across the 1 kb read-through window
    3. Ranked bar plot of read-through CPM
    4. Paper-style downstream/upstream summaries when available
    5. Summary statistics table

Usage (executed by SLURM or interactively):
  source /hpc/group/parrishlab/envs/l1nedd4_py39/bin/activate
  python3 05_visualize.py --results-dir /work/nfp8/2026_03_02/GSE226189/results \
                          --coverage-dir /work/nfp8/2026_03_02/GSE226189/coverage \
                          --rt-window chr15:55958936-55959936 \
                          --output-dir   /hpc/group/parrishlab/L1_NEDD4/results/GSE226189/figures
"""
import argparse, os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ── CLI arguments ──────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="L1-NEDD4 read-through figures")
    p.add_argument("--results-dir", required=True, help="Dir with readthrough_signal.csv + matrices")
    p.add_argument("--coverage-dir", required=True, help="Dir with per-sample _readthrough.depth files")
    p.add_argument("--rt-window", default="chr15:55958936-55959936", help="Read-through window coord")
    p.add_argument("--output-dir", required=True, help="Where to write figures")
    p.add_argument("--gse", default=None, help="Dataset label for plot titles (auto-detected if omitted)")
    return p.parse_args()

# ── Figure 1: Strip-plot of read-through CPM ──────────────────────────────
def fig_strip_cpm(df, out_dir, gse):
    fig, ax = plt.subplots(figsize=(4, 6))
    sns.stripplot(y="rt_cpm", data=df, jitter=0.3, size=5, alpha=0.7, ax=ax, color="steelblue")
    ax.axhline(df["rt_cpm"].median(), color="red", linestyle="--", linewidth=0.8, label="median")
    ax.set_ylabel("Read-Through CPM")
    ax.set_title(f"L1-NEDD4 Read-Through Signal\n{gse} (n={len(df)})")
    ax.legend(fontsize=8)
    sns.despine()
    fig.tight_layout()
    path = os.path.join(out_dir, f"{gse}_rt_cpm_strip.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path}")

# ── Figure 2: Aggregated per-base coverage across read-through window ─────
def fig_coverage_profile(coverage_dir, srrs, rt_window, out_dir, gse):
    chrom, coords = rt_window.split(":")
    start, end = int(coords.split("-")[0]), int(coords.split("-")[1])
    positions = np.arange(start, end + 1)
    all_depths = np.zeros((len(srrs), len(positions)), dtype=float)

    loaded_rows = []  # track which row indices were populated
    for i, srr in enumerate(srrs):
        fp = os.path.join(coverage_dir, f"{srr}_readthrough.depth")
        if not os.path.exists(fp):
            continue
        d = pd.read_csv(fp, sep="\t", header=None, names=["chr","pos","depth"])
        if len(d) == 0:
            continue
        # Map positions to array indices
        idx_map = {pos: j for j, pos in enumerate(positions)}
        for _, row in d.iterrows():
            j = idx_map.get(int(row["pos"]))
            if j is not None:
                all_depths[i, j] = row["depth"]
        loaded_rows.append(i)

    if not loaded_rows:
        print("  WARN: no depth files loaded, skipping coverage profile")
        return

    loaded_data = all_depths[loaded_rows, :]  # only rows that were actually populated
    mean_depth = np.nanmean(loaded_data, axis=0)
    q25 = np.nanpercentile(loaded_data, 25, axis=0)
    q75 = np.nanpercentile(loaded_data, 75, axis=0)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.fill_between(positions, q25, q75, alpha=0.25, color="steelblue", label="IQR")
    ax.plot(positions, mean_depth, linewidth=1.0, color="steelblue", label="Mean")

    # Mark the L1 breakpoint (start of window)
    ax.axvline(start, color="red", linestyle="--", linewidth=0.8, label="L1 3′ break")

    ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel("Read Depth")
    ax.set_title(f"L1-NEDD4 Read-Through Window Coverage\n{gse} (n={len(loaded_rows)})")
    ax.legend(fontsize=8)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    sns.despine()
    fig.tight_layout()
    path = os.path.join(out_dir, f"{gse}_rt_coverage_profile.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path}")

# ── Figure 3: RT CPM ranked bar plot ──────────────────────────────────────
def fig_ranked_bar(df, out_dir, gse):
    df_sorted = df.sort_values("rt_cpm", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(max(8, len(df_sorted)*0.15), 4))
    colors = ["firebrick" if x > df["rt_cpm"].median() else "steelblue" for x in df_sorted["rt_cpm"]]
    ax.bar(range(len(df_sorted)), df_sorted["rt_cpm"], color=colors, width=0.8)
    ax.set_xlabel("Sample rank")
    ax.set_ylabel("Read-Through CPM")
    ax.set_title(f"L1-NEDD4 RT CPM — Ranked Samples\n{gse}")
    sns.despine()
    fig.tight_layout()
    path = os.path.join(out_dir, f"{gse}_rt_cpm_ranked.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path}")


def fig_paper_style(df, out_dir, gse):
    value_cols = [
        ("paper_downstream_sense_cpm", "Downstream sense CPM", "#bd4b32"),
        ("paper_upstream_antisense_cpm", "Upstream antisense CPM", "#2a7f62"),
    ]
    available = [col for col, _, _ in value_cols if col in df.columns and df[col].notna().any()]
    if not available:
        note_path = os.path.join(out_dir, f"{gse}_paper_style_note.txt")
        with open(note_path, "w") as fh:
            fh.write(
                "Paper-style outputs unavailable for this dataset.\n"
                "Reason: no strand-resolved mate1-only paper-style metrics were produced, which is expected for unstranded libraries.\n"
            )
        print(f"  Saved: {note_path}")
        return False

    plot_df = df[["srr"] + available].melt(id_vars="srr", var_name="metric", value_name="value")
    plot_df = plot_df.dropna(subset=["value"])
    if plot_df.empty:
        note_path = os.path.join(out_dir, f"{gse}_paper_style_note.txt")
        with open(note_path, "w") as fh:
            fh.write(
                "Paper-style outputs unavailable for this dataset.\n"
                "Reason: all paper-style summary values were empty after collection.\n"
            )
        print(f"  Saved: {note_path}")
        return False

    label_map = {col: label for col, label, _ in value_cols}
    color_map = {col: color for col, _, color in value_cols}
    plot_df["metric_label"] = plot_df["metric"].map(label_map)

    fig, ax = plt.subplots(figsize=(5.5, 6))
    sns.stripplot(
        data=plot_df,
        x="metric_label",
        y="value",
        jitter=0.25,
        size=5,
        alpha=0.7,
        palette=[color_map[m] for m in available],
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("CPM")
    ax.set_title(f"Paper-Style Stranded Proxies\n{gse}")
    sns.despine()
    fig.tight_layout()
    path = os.path.join(out_dir, f"{gse}_paper_style_strip.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved: {path}")
    return True

# ── Summary stats ──────────────────────────────────────────────────────────
def write_summary(df, paper_df, out_dir, gse):
    stats = df["rt_cpm"].describe()
    path = os.path.join(out_dir, f"{gse}_summary_stats.txt")
    with open(path, "w") as fh:
        fh.write(f"L1-NEDD4 Read-Through Summary — {gse}\n")
        fh.write("=" * 50 + "\n")
        fh.write(f"Samples analysed: {len(df)}\n")
        fh.write(f"Median RT CPM:    {df['rt_cpm'].median():.6f}\n")
        fh.write(f"Mean RT CPM:      {df['rt_cpm'].mean():.6f}\n")
        fh.write(f"SD RT CPM:        {df['rt_cpm'].std():.6f}\n")
        fh.write(f"Min:              {df['rt_cpm'].min():.6f}\n")
        fh.write(f"Max:              {df['rt_cpm'].max():.6f}\n")
        fh.write(f"IQR:              {df['rt_cpm'].quantile(0.25):.6f} – {df['rt_cpm'].quantile(0.75):.6f}\n\n")
        fh.write("Non-zero RT signal:\n")
        nz = (df["rt_cpm"] > 0).sum()
        fh.write(f"  {nz}/{len(df)} samples ({nz/len(df)*100:.1f}%)\n")
        fh.write(f"\nTotal mapped reads — median: {df['total_mapped'].median():.0f}\n")
        if paper_df is None or paper_df.empty:
            fh.write("\nPaper-style summary: not available\n")
        else:
            downstream_col = "paper_downstream_sense_cpm"
            upstream_col = "paper_upstream_antisense_cpm"
            downstream_present = downstream_col in paper_df.columns and paper_df[downstream_col].notna().any()
            upstream_present = upstream_col in paper_df.columns and paper_df[upstream_col].notna().any()
            fh.write("\nPaper-style summary:\n")
            if downstream_present:
                fh.write(f"  Median downstream sense CPM: {paper_df[downstream_col].fillna(0).median():.6f}\n")
            else:
                fh.write("  Downstream sense CPM: not available\n")
            if upstream_present:
                fh.write(f"  Median upstream antisense CPM: {paper_df[upstream_col].fillna(0).median():.6f}\n")
            else:
                fh.write("  Upstream antisense CPM: not available\n")
    print(f"  Saved: {path}")

# ── Main ───────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    if args.gse:
        gse = args.gse
    else:
        gse = os.path.basename(os.path.dirname(args.results_dir.rstrip("/")))
        if gse in ("results", "figures", "", "."):
            gse = "unknown_dataset"

    rt_csv = os.path.join(args.results_dir, "readthrough_signal.csv")
    paper_csv = os.path.join(args.results_dir, "paper_style_signal.csv")
    if not os.path.exists(rt_csv):
        sys.exit(f"ERROR: {rt_csv} not found. Run 04_collect_results.sh first.")

    df = pd.read_csv(rt_csv)
    paper_df = pd.read_csv(paper_csv) if os.path.exists(paper_csv) else None
    print(f"Loaded {len(df)} samples from {rt_csv}")
    if paper_df is not None:
        print(f"Loaded {len(paper_df)} samples from {paper_csv}")

    fig_strip_cpm(df, args.output_dir, gse)
    fig_coverage_profile(args.coverage_dir, df["srr"].tolist(), args.rt_window, args.output_dir, gse)
    fig_ranked_bar(df, args.output_dir, gse)
    if paper_df is not None:
        fig_paper_style(paper_df, args.output_dir, gse)
    write_summary(df, paper_df, args.output_dir, gse)
    print("Visualization complete.")

if __name__ == "__main__":
    main()
