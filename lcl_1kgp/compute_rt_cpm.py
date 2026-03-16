#!/usr/bin/env python3
"""
Compute normalized readthrough (RT) signal for MCF7 samples.

Metric: reads in RT_WINDOW with SENSE flags (-f 64 -F 16 -q 20)
        divided by total uniquely-mapped reads × 1,000,000 = RT_CPM

RT_WINDOW: chr15:55958936-55959936  (1 kb downstream of L1 3' breakpoint)
SENSE = R1 forward = plus-strand RNA = -f 64 (R1) -F 16 (not reverse) -q 20
"""

import subprocess
import sys

REGION = "chr15:55958936-55959936"

SAMPLES = [
    # (srr, condition)
    ("ERR973734", "untreated_R1"),
    ("ERR973735", "untreated_R2"),
    ("ERR973728", "shRNA_R1"),
    ("ERR973729", "shRNA_R2"),
]

ALIGNERS = {
    "STAR":    "/work/nfp8/MCF7_ctrl/bam/{srr}_Aligned.sortedByCoord.out.bam",
    "TopHat2": "/work/nfp8/MCF7_ctrl_tophat/bam/{srr}_Aligned.sortedByCoord.out.bam",
}

# Total uniquely-mapped reads from alignment logs (pre-computed)
# STAR: "Uniquely mapped reads number" from *_Log.final.out
# TopHat2: Left reads - multiple alignments (MAPQ≥20 ≈ unique Left reads)
TOTAL_MAPPED = {
    "STAR": {
        "ERR973734": 88657457,
        "ERR973735": 112376114,
        "ERR973728": 63742872,
        "ERR973729": 24738303,
    },
    "TopHat2": {
        # Left_mapped - multiple_alignments  from align_summary.txt
        "ERR973734": 94697431 - 19554996,   # = 75142435
        "ERR973735": 120752438 - 25635097,  # = 95117341
        "ERR973728": 66258413 - 13398703,   # = 52859710
        "ERR973729": 25124565 - 3922004,    # = 21202561
    },
}


def count_sense_reads(bam, region):
    """Count R1 forward-strand reads overlapping region (MAPQ >= 20)."""
    cmd = ["samtools", "view", "-c", "-f", "64", "-F", "16", "-q", "20", bam, region]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return int(result.stdout.strip())


def main():
    header = f"{'sample':<15} {'condition':<15} {'aligner':<10} {'window_reads':>13} {'total_mapped':>14} {'RT_CPM':>10}"
    print(header)
    print("-" * len(header))

    rows = []
    for aligner, bam_template in ALIGNERS.items():
        for srr, condition in SAMPLES:
            bam = bam_template.format(srr=srr)
            window_reads = count_sense_reads(bam, REGION)
            total = TOTAL_MAPPED[aligner][srr]
            rt_cpm = window_reads / total * 1_000_000
            rows.append((srr, condition, aligner, window_reads, total, rt_cpm))
            print(f"{srr:<15} {condition:<15} {aligner:<10} {window_reads:>13,} {total:>14,} {rt_cpm:>10.4f}")

    print()
    print("Note: SENSE = R1 forward (-f 64 -F 16 -q 20), region =", REGION)
    print("      TopHat2 denominator = unique Left reads (Left_mapped - multi_alignments)")
    print("      STAR denominator    = Uniquely mapped reads number (Log.final.out)")


if __name__ == "__main__":
    main()
