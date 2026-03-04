#!/bin/bash
#SBATCH --job-name=post_align
#SBATCH --output=/work/nfp8/2026_03_02/GSE226189/logs/post_%A_%a.out
#SBATCH --error=/work/nfp8/2026_03_02/GSE226189/logs/post_%A_%a.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH --time=08:00:00
#SBATCH --array=1-82
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"
module load "$MOD_SAMTOOLS"
module load "$MOD_SUBREAD"
module load "$MOD_BEDTOOLS"

# ---- Resolve sample ----
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")
[[ -z "$SRR" ]] && { echo "ERROR: no SRR at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
echo "[$(date -Is)] Task ${SLURM_ARRAY_TASK_ID}: post-processing ${SRR}"

FINAL_BAM="${BAM_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"
[[ -f "$FINAL_BAM" ]] || { echo "ERROR: BAM not found: $FINAL_BAM" >&2; exit 1; }

DEDUP_BAM="${BAM_DIR}/${SRR}_dedup.bam"
TMP_DIR="${BAM_DIR}/tmp_postAlign_${SRR}_$$"
mkdir -p "$COVERAGE_DIR" "$PERSISTENT_BAM" "$COUNTS_DIR" "$TMP_DIR"

# ===========================================================================
# 1. PCR duplicate marking (samtools markdup via collate -> fixmate -> sort pipeline)
#    Required for coverage-based analyses; % duplication reported per sample.
# ===========================================================================
if [[ ! -f "${DEDUP_BAM}.bai" ]]; then
    echo "  [1] Marking PCR duplicates..."
    # collate: rearrange into name-order for fixmate; fixmate: adds MC/ms tags;
    # sort: back to coordinate order; markdup: mark (not remove) duplicates + stats
    samtools collate -O "$FINAL_BAM" "${TMP_DIR}/collate_${SRR}" \
      | samtools fixmate -m -r - - \
      | samtools sort -T "${TMP_DIR}/sort_${SRR}" - \
      | samtools markdup -s - "$DEDUP_BAM" \
          2> "${COVERAGE_DIR}/${SRR}_markdup.stats"
    samtools index "$DEDUP_BAM"
    echo "    Markdup stats: $(head -5 "${COVERAGE_DIR}/${SRR}_markdup.stats" | tr '\n' '|')"
else
    echo "  [1] Dedup BAM already exists — skipping"
fi

# ===========================================================================
# 2. Flagstat on dedup BAM (primary mapped reads = CPM denominator)
# ===========================================================================
FLAGSTAT="${COVERAGE_DIR}/${SRR}_flagstat.txt"
if [[ ! -f "$FLAGSTAT" ]]; then
    echo "  [2] Running flagstat..."
    samtools flagstat "$DEDUP_BAM" > "$FLAGSTAT"
fi

# ===========================================================================
# 3. Read-through window depth — strand-aware, MAPQ >= 255 (uniquely mapped)
#
#   For stranded libraries the L1 read-through travels on the PLUS strand.
#   We extract only plus-strand reads to avoid NEDD4 minus-strand intronic coverage:
#
#   STRANDEDNESS="reverse" (fr-firststrand / TruSeq Stranded):
#     - R2 forward  (flag: read2 [0x80] + NOT reverse [NOT 0x10]) -> plus strand
#     - R1 reverse  (flag: read1 [0x40] + reverse [0x10])         -> plus strand
#
#   STRANDEDNESS="forward" (fr-secondstrand):
#     - R1 forward  (flag: read1 [0x40] + NOT reverse [NOT 0x10]) -> plus strand
#     - R2 reverse  (flag: read2 [0x80] + reverse [0x10])         -> plus strand
#
#   STRANDEDNESS="unstranded": all reads, MAPQ >= 255 filter only.
#   -Q in samtools depth = minimum mapping quality
# ===========================================================================
RT_DEPTH="${COVERAGE_DIR}/${SRR}_readthrough.depth"
if [[ ! -f "$RT_DEPTH" ]]; then
    echo "  [3] Computing read-through depth (${STRANDEDNESS} library)..."

    if [[ "$STRANDEDNESS" == "unstranded" ]]; then
        samtools depth -a -d 0 -Q 255 -r "$RT_WINDOW" "$DEDUP_BAM" > "$RT_DEPTH"

    else
        R2FWD="${TMP_DIR}/${SRR}_r2fwd.bam"
        R1REV="${TMP_DIR}/${SRR}_r1rev.bam"
        MERGED="${TMP_DIR}/${SRR}_plus_strand.bam"

        if [[ "$STRANDEDNESS" == "reverse" ]]; then
            # TruSeq Stranded (fr-firststrand): plus-strand reads
            samtools view -bh -q 255 -f 128 -F 16 "$DEDUP_BAM" "$RT_WINDOW" > "$R2FWD"
            samtools view -bh -q 255 -f  80       "$DEDUP_BAM" "$RT_WINDOW" > "$R1REV"
        else
            # fr-secondstrand: plus-strand reads
            samtools view -bh -q 255 -f  64 -F 16 "$DEDUP_BAM" "$RT_WINDOW" > "$R2FWD"
            samtools view -bh -q 255 -f 144       "$DEDUP_BAM" "$RT_WINDOW" > "$R1REV"
        fi

        samtools merge -f "$MERGED" "$R2FWD" "$R1REV"
        samtools sort "$MERGED" -o "${MERGED%.bam}_sorted.bam"
        samtools index "${MERGED%.bam}_sorted.bam"
        samtools depth -a -d 0 -r "$RT_WINDOW" "${MERGED%.bam}_sorted.bam" > "$RT_DEPTH"
        rm -f "$R2FWD" "$R1REV" "$MERGED" "${MERGED%.bam}_sorted.bam" "${MERGED%.bam}_sorted.bam.bai" || true
    fi

    echo "    RT window: $(wc -l < "$RT_DEPTH") positions | $(awk '{s+=$3} END{printf "mean depth=%.3f", s/NR}' "$RT_DEPTH")"
fi

# ===========================================================================
# 4. Background window depth (same strand logic; flanking normalization in Step 4)
#    BACKGROUND_WINDOWS_BED built by 01_build_star_index.sh (NEDD4 exons subtracted)
# ===========================================================================
BG_DEPTH="${COVERAGE_DIR}/${SRR}_background.depth"
if [[ ! -f "$BG_DEPTH" && -f "${BACKGROUND_WINDOWS_BED}" ]]; then
    echo "  [4] Computing background window depth..."

    if [[ "$STRANDEDNESS" == "unstranded" ]]; then
        samtools depth -a -d 0 -Q 255 -b "$BACKGROUND_WINDOWS_BED" "$DEDUP_BAM" > "$BG_DEPTH"

    else
        R2FWD="${TMP_DIR}/${SRR}_bg_r2fwd.bam"
        R1REV="${TMP_DIR}/${SRR}_bg_r1rev.bam"
        MERGED="${TMP_DIR}/${SRR}_bg_plus.bam"

        if [[ "$STRANDEDNESS" == "reverse" ]]; then
            samtools view -bh -q 255 -f 128 -F 16 "$DEDUP_BAM" > "$R2FWD"
            samtools view -bh -q 255 -f  80       "$DEDUP_BAM" > "$R1REV"
        else
            samtools view -bh -q 255 -f  64 -F 16 "$DEDUP_BAM" > "$R2FWD"
            samtools view -bh -q 255 -f 144       "$DEDUP_BAM" > "$R1REV"
        fi

        samtools merge -f "$MERGED" "$R2FWD" "$R1REV"
        samtools sort "$MERGED" -o "${MERGED%.bam}_sorted.bam"
        samtools index "${MERGED%.bam}_sorted.bam"
        samtools depth -a -d 0 -b "$BACKGROUND_WINDOWS_BED" "${MERGED%.bam}_sorted.bam" > "$BG_DEPTH"
        rm -f "$R2FWD" "$R1REV" "$MERGED" "${MERGED%.bam}_sorted.bam" "${MERGED%.bam}_sorted.bam.bai" || true
    fi

    echo "    Background: $(wc -l < "$BG_DEPTH") positions | $(awk '{s+=$3} END{printf "mean depth=%.3f", s/NR}' "$BG_DEPTH")"
elif [[ ! -f "${BACKGROUND_WINDOWS_BED}" ]]; then
    echo "  [4] WARNING: ${BACKGROUND_WINDOWS_BED} not found — skipping background depth"
fi

# ===========================================================================
# 5. Subset BAM to NEDD4 locus +/-1 Mb -> persistent storage (for IGV inspection)
# ===========================================================================
SUBSET_BAM="${PERSISTENT_BAM}/${SRR}_nedd4locus.bam"
if [[ ! -f "${SUBSET_BAM}.bai" ]]; then
    echo "  [5] Subsetting BAM to NEDD4 locus..."
    samtools view -b -h -q 255 "$DEDUP_BAM" "$BAM_SUBSET_REGION" > "$SUBSET_BAM"
    samtools index "$SUBSET_BAM"
    echo "    Subset BAM: $(du -sh "$SUBSET_BAM" | cut -f1)"
fi

# ===========================================================================
# 6. featureCounts gene quantification (strand-aware, PE-aware)
#    Retained alongside STAR GeneCounts; use featureCounts for final gene expression.
# ===========================================================================
FC_OUT="${COUNTS_DIR}/${SRR}_featureCounts.txt"
if [[ ! -f "$FC_OUT" ]]; then
    echo "  [6] Running featureCounts..."

    FC_STRAND=0
    [[ "$STRANDEDNESS" == "forward" ]] && FC_STRAND=1
    [[ "$STRANDEDNESS" == "reverse" ]] && FC_STRAND=2

    FC_FLAGS="-T ${STAR_THREADS} -a ${GENCODE_GTF} -s ${FC_STRAND}"
    [[ "$LIBRARY_LAYOUT" == "PE" ]] && FC_FLAGS="${FC_FLAGS} -p --countReadPairs -B"

    featureCounts \
        ${FC_FLAGS} \
        -o "$FC_OUT" \
        "$DEDUP_BAM" \
        2> "${COUNTS_DIR}/${SRR}_featureCounts.log"

    echo "    featureCounts: $(grep 'Assigned' "${COUNTS_DIR}/${SRR}_featureCounts.log" | tail -1 || echo 'see log')"
fi

# ===========================================================================
# 7. Copy STAR GeneCounts to counts dir (retain both quantification methods)
# ===========================================================================
STAR_COUNTS="${BAM_DIR}/${SRR}_ReadsPerGene.out.tab"
[[ -f "$STAR_COUNTS" ]] && cp "$STAR_COUNTS" "${COUNTS_DIR}/${SRR}_ReadsPerGene.out.tab"

# ===========================================================================
# Cleanup
# ===========================================================================
rm -rf "$TMP_DIR"
echo "[$(date -Is)] Post-processing complete: ${SRR}"
