#!/usr/bin/env python3
"""
make_igv_report.py — IGV-style per-base coverage HTML for a genomic window.

Usage:
    python3 make_igv_report.py \
        --bam-dir /work/nfp8/2026_03_02/GSE226189/bam \
        --region chr15:55958955-55959955 \
        --out /work/nfp8/2026_03_02/GSE226189/results/igv_report.html \
        --title "GSE226189 chr15:55958955-55959955"
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bam-dir", required=True)
    p.add_argument("--region", required=True, help="e.g. chr15:55958955-55959955")
    p.add_argument("--out", required=True)
    p.add_argument("--title", default="Coverage profile")
    p.add_argument("--bam-suffix", default="_dedup.bam")
    p.add_argument("--min-mapq", type=int, default=20,
                   help="Minimum MAPQ. Default 20 matches Philippe et al. igvtools count threshold.")
    p.add_argument("--labels-file",
                   help="Optional 2-column TSV (no header): sample_name<TAB>display_label")
    return p.parse_args()


def get_strand_coverage(bam_path: Path, region: str, flag_pairs: list, min_mapq: int) -> dict:
    """
    Get per-base depth by summing separate samtools view | samtools depth pipelines,
    one per (include_flags, exclude_flags) pair.  No shell process substitution needed.
    """
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    combined: dict = {}

    for inc, exc in flag_pairs:
        view_cmd = ["samtools", "view", "-b"]
        if inc:
            view_cmd += ["-f", str(inc)]
        if exc:
            view_cmd += ["-F", str(exc)]
        if min_mapq:
            view_cmd += ["-q", str(min_mapq)]
        view_cmd += [str(bam_path), region]

        p_view = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p_depth = subprocess.Popen(
            ["samtools", "depth", "-"],
            stdin=p_view.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        p_view.stdout.close()  # allow p_view to receive SIGPIPE if p_depth exits

        for line in p_depth.stdout:
            cols = line.split("\t")
            if len(cols) >= 3:
                pos = int(cols[1])
                if start <= pos <= end:
                    combined[pos] = combined.get(pos, 0) + int(cols[2])

        p_depth.wait()
        p_view.wait()

    return combined


def main():
    args = parse_args()

    # Load optional sample label map
    label_map = {}
    if args.labels_file:
        with open(args.labels_file) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2:
                    label_map[parts[0].strip()] = parts[1].strip()

    bam_dir = Path(args.bam_dir)
    bam_files = sorted(bam_dir.glob(f"*{args.bam_suffix}"))
    if not bam_files:
        sys.exit(f"ERROR: no BAMs matching *{args.bam_suffix} in {bam_dir}")
    print(f"Found {len(bam_files)} BAM files", flush=True)

    # Parse region
    chrom, coords = args.region.split(":")
    start, end = map(int, coords.split("-"))
    positions = list(range(start, end + 1))

    # Fragment-level coverage: R1 only (one count per fragment), matching Philippe et al.
    # Empirically verified from the MCF7 data (Philippe et al.):
    #   R1 is the SENSE read (same orientation as the mRNA) — consistent with
    #   fr-secondstrand library chemistry or R1/R2 labelling swapped in GEO submission.
    #   Sense (plus-strand RNA):      R1 forward  (-f 64 -F 16 = read1+NOT reverse)
    #   Antisense (minus-strand RNA): R1 reverse  (-f 80      = read1+reverse)
    SENSE_FLAGS     = [(64, 16)]   # R1-forward → plus-strand RNA (e.g. L1 readthrough)
    ANTISENSE_FLAGS = [(80, 0)]    # R1-reverse → minus-strand RNA (e.g. NEDD4 exons)

    samples = []
    all_sense = []
    all_antisense = []

    for i, bam in enumerate(bam_files):
        sample = bam.name.replace(args.bam_suffix, "")
        print(f"  [{i+1}/{len(bam_files)}] {sample}", flush=True)
        sample = label_map.get(sample, sample)
        sense_cov = get_strand_coverage(bam, args.region, SENSE_FLAGS, args.min_mapq)
        anti_cov  = get_strand_coverage(bam, args.region, ANTISENSE_FLAGS, args.min_mapq)
        sense_depths = [sense_cov.get(p, 0) for p in positions]
        anti_depths  = [anti_cov.get(p, 0)  for p in positions]
        samples.append(sample)
        all_sense.append(sense_depths)
        all_antisense.append(anti_depths)

    # Sort descending by combined max signal (highest signal on top = control/scramble samples)
    order = sorted(
        range(len(samples)),
        key=lambda i: max(all_sense[i]) + max(all_antisense[i]),
        reverse=True,
    )
    samples      = [samples[i]      for i in order]
    all_sense    = [all_sense[i]    for i in order]
    all_antisense = [all_antisense[i] for i in order]

    n = len(samples)
    print(f"Building Plotly figure with {n} tracks...", flush=True)

    fig = make_subplots(
        rows=n,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.002,
    )

    zeros = [0] * len(positions)

    for i, (sample, sense, antisense) in enumerate(zip(samples, all_sense, all_antisense)):
        row = i + 1
        neg_anti = [-d for d in antisense]  # mirror below zero

        # --- Sense (blue, above zero) ---
        # Hidden zero-baseline first, then fill="tonexty" on coverage trace
        fig.add_trace(
            go.Scatter(
                x=positions, y=zeros,
                mode="lines",
                line=dict(width=0, color="rgba(0,0,0,0)"),
                hoverinfo="skip", showlegend=False,
            ),
            row=row, col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=positions, y=sense,
                mode="lines",
                fill="tonexty",
                line=dict(width=0.8, color="rgba(33,102,172,1)"),
                fillcolor="rgba(33,102,172,0.5)",
                hoverinfo="skip", showlegend=False,
            ),
            row=row, col=1,
        )

        # --- Antisense (red, below zero) ---
        fig.add_trace(
            go.Scatter(
                x=positions, y=zeros,
                mode="lines",
                line=dict(width=0, color="rgba(0,0,0,0)"),
                hoverinfo="skip", showlegend=False,
            ),
            row=row, col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=positions, y=neg_anti,
                mode="lines",
                fill="tonexty",
                line=dict(width=0.8, color="rgba(178,24,43,1)"),
                fillcolor="rgba(178,24,43,0.5)",
                hoverinfo="skip", showlegend=False,
            ),
            row=row, col=1,
        )
        ymax = 40
        fig.update_yaxes(
            range=[-ymax, ymax],
            showticklabels=False,
            zeroline=True, zerolinewidth=0.5, zerolinecolor="#888",
            title_text="<span style='font-size:7px'>±40</span>",
            title_standoff=1,
            row=row, col=1,
        )

    # Add horizontal sample labels using each subplot's actual y-domain
    for i, sample in enumerate(samples):
        row = i + 1
        axis_key = "yaxis" if row == 1 else f"yaxis{row}"
        domain = fig.layout[axis_key].domain
        y_center = (domain[0] + domain[1]) / 2
        fig.add_annotation(
            x=-0.01, xref="paper",
            y=y_center, yref="paper",
            text=sample,
            showarrow=False,
            xanchor="right",
            yanchor="middle",
            font=dict(size=9),
        )

    # Global layout
    fig.update_layout(
        title=dict(text=args.title, font=dict(size=14)),
        height=max(200, n * 30),  # ~30px per sample
        margin=dict(l=260, r=20, t=50, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
        hovermode="x unified",
    )
    fig.update_xaxes(
        title_text=args.region,
        tickformat=",d",
        row=n, col=1,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(out_path), full_html=True, include_plotlyjs="cdn")
    print(f"Wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
