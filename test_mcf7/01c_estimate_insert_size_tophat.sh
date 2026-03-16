#!/bin/bash
#SBATCH --job-name=insize_mcf7
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=16G
#SBATCH -c 2
#SBATCH --time=01:00:00
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config_test_mcf7_tophat.sh}"

source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
module load samtools/1.21

mkdir -p "$TOPHAT_INSERT_METRICS_DIR"

read_length_from_fastq() {
    local sample="$1"
    local fastq
    if [[ "$TOPHAT_USE_TRIMMED" == "1" ]]; then
        fastq="${FASTQ_DIR}/${sample}_1_trimmed_P.fq.gz"
    else
        fastq="${FASTQ_DIR}/${sample}_1.fastq.gz"
    fi
    zcat "$fastq" | sed -n '2p' | awk '{print length($0); exit}'
}

first_sample="$(awk 'NF{print; exit}' "$SRR_LIST")"
read_len="$(read_length_from_fastq "$first_sample")"
[[ -n "$read_len" ]] || { echo "ERROR: could not infer read length" >&2; exit 1; }

echo "[$(date -Is)] Estimating insert size from samples: ${TOPHAT_INSERT_SIZE_SOURCE_SAMPLES}"
echo "  Inferred read length: ${read_len}"

export METRICS_DIR="$TOPHAT_INSERT_METRICS_DIR"
export SOURCE_BAM_DIR="$TOPHAT_INSERT_SIZE_SOURCE_BAM_DIR"
export SOURCE_SAMPLES="$TOPHAT_INSERT_SIZE_SOURCE_SAMPLES"
export READ_LEN="$read_len"
export INSERT_ENV="$TOPHAT_INSERT_SIZE_ENV"

"${VENV_DIR}/bin/python3" - <<'PYEOF'
import os
import statistics
import subprocess
from pathlib import Path

metrics_dir = Path(os.environ["METRICS_DIR"])
source_bam_dir = Path(os.environ["SOURCE_BAM_DIR"])
samples = os.environ["SOURCE_SAMPLES"].split()
read_len = int(os.environ["READ_LEN"])
insert_env = Path(os.environ["INSERT_ENV"])

means = []
stdevs = []

for sample in samples:
    bam = source_bam_dir / f"{sample}_Aligned.sortedByCoord.out.bam"
    txt = metrics_dir / f"{sample}_insert_metrics.txt"
    if not bam.exists():
        raise SystemExit(f"Missing BAM for insert-size estimation: {bam}")

    # Extract insert sizes from first 500k paired reads to avoid OOM
    result = subprocess.run(
        f"samtools view -f 0x0001 '{bam}' | awk '{{if($9!=0)print $9}}' | head -500000",
        shell=True,
        capture_output=True,
        text=True,
        check=False,
    )
    
    insert_sizes = [
        abs(int(line.strip()))
        for line in result.stdout.split('\n')
        if line.strip() and 0 < abs(int(line.strip())) <= 2000
    ]
    
    if not insert_sizes:
        raise SystemExit(f"No paired reads found in {bam}. Samtools output: {result.stdout[:500] if result.stdout else '(empty)'}")
    
    mean_insert = statistics.mean(insert_sizes)
    stdev_insert = statistics.stdev(insert_sizes)
    
    means.append(mean_insert)
    stdevs.append(stdev_insert)
    
    # Write metrics file in Picard-like format for reference
    with txt.open("w") as handle:
        handle.write("## samtools-based insert size metrics\n")
        handle.write(f"SAMPLE\t{sample}\n")
        handle.write(f"MEDIAN_INSERT_SIZE\tMEAN_INSERT_SIZE\tSTANDARD_DEVIATION\n")
        handle.write(f"{int(statistics.median(insert_sizes))}\t{mean_insert:.2f}\t{stdev_insert:.2f}\n")

mean_insert = statistics.mean(means)
mean_stdev = statistics.mean(stdevs)
mate_inner = round(mean_insert - 2 * read_len)
mate_stdev = round(mean_stdev)

with insert_env.open("w") as handle:
    handle.write(f"export TOPHAT_ESTIMATED_READ_LENGTH={read_len}\n")
    handle.write(f"export TOPHAT_MEAN_INSERT_SIZE={mean_insert:.2f}\n")
    handle.write(f"export TOPHAT_STDDEV_INSERT_SIZE={mean_stdev:.2f}\n")
    handle.write(f"export TOPHAT_MATE_INNER_DIST={mate_inner}\n")
    handle.write(f"export TOPHAT_MATE_STD_DEV={mate_stdev}\n")

print(f"Wrote {insert_env}")
print(f"Mean insert size: {mean_insert:.2f}")
print(f"Insert size stdev: {mean_stdev:.2f}")
print(f"TopHat mate inner distance (-r): {mate_inner}")
print(f"TopHat mate std dev: {mate_stdev}")
PYEOF

echo "[$(date -Is)] Insert-size estimation complete: ${TOPHAT_INSERT_SIZE_ENV}"
