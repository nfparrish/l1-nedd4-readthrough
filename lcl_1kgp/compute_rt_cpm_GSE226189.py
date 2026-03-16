#!/usr/bin/env python3
"""
Compute normalized readthrough (RT) CPM for all 82 samples in GSE226189
(Transcriptomes of human primary skin fibroblasts, GEO: GSE226189).

Metric: RT_CPM = (sense reads in RT_WINDOW) / (total uniquely mapped reads) × 1e6

RT_WINDOW:  chr15:55958936-55959936  (1 kb downstream of L1 3' breakpoint)
SENSE:      R1 forward = plus-strand RNA = -f 64 -F 16 -q 20
TOTAL:      Uniquely mapped reads number from STAR *_Log.final.out
"""

import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REGION = "chr15:55958936-55959936"
BAM_DIR = Path("/work/nfp8/2026_03_02/GSE226189/bam")
SRR_LIST = Path("/hpc/home/nfp8/copilot/2026_03_02/SRR_list.txt")
OUT_TSV = Path("/hpc/home/nfp8/copilot/2026_03_03/lcl_1kgp/GSE226189_rt_cpm.tsv")


def get_total_mapped(srr: str) -> int:
    log = BAM_DIR / f"{srr}_Log.final.out"
    for line in log.read_text().splitlines():
        if "Uniquely mapped reads number" in line:
            return int(line.split("|")[1].strip())
    raise ValueError(f"Could not parse total mapped from {log}")


def count_sense_reads(srr: str) -> int:
    bam = BAM_DIR / f"{srr}_Aligned.sortedByCoord.out.bam"
    cmd = ["samtools", "view", "-c", "-f", "64", "-F", "16", "-q", "20", str(bam), REGION]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return int(result.stdout.strip())


def process_sample(srr: str):
    total = get_total_mapped(srr)
    window = count_sense_reads(srr)
    rt_cpm = window / total * 1_000_000
    return srr, window, total, rt_cpm


def main():
    srr_list = [s.strip() for s in SRR_LIST.read_text().splitlines() if s.strip()]
    print(f"Processing {len(srr_list)} samples...", flush=True)

    results = []
    # Using 8 threads — samtools window counts are fast I/O
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(process_sample, srr): srr for srr in srr_list}
        for i, fut in enumerate(as_completed(futures), 1):
            srr = futures[fut]
            try:
                res = fut.result()
                results.append(res)
                print(f"  [{i:2d}/82] {srr}  window={res[1]}  RT_CPM={res[3]:.4f}", flush=True)
            except Exception as e:
                print(f"  [{i:2d}/82] {srr}  ERROR: {e}", file=sys.stderr, flush=True)

    # Sort by SRR accession
    results.sort(key=lambda r: r[0])

    # Write TSV
    with OUT_TSV.open("w") as fh:
        fh.write("sample\twindow_reads\ttotal_mapped\tRT_CPM\n")
        for srr, window, total, rt_cpm in results:
            fh.write(f"{srr}\t{window}\t{total}\t{rt_cpm:.6f}\n")

    print(f"\nWrote {len(results)} rows to {OUT_TSV}")
    print(f"\nSummary:")
    cpms = [r[3] for r in results]
    print(f"  min RT_CPM : {min(cpms):.4f}")
    print(f"  max RT_CPM : {max(cpms):.4f}")
    print(f"  mean RT_CPM: {sum(cpms)/len(cpms):.4f}")
    print(f"\nNote: SENSE = R1 forward (-f 64 -F 16 -q 20), region = {REGION}")


if __name__ == "__main__":
    main()
