#!/bin/bash
#SBATCH --job-name=post_fasttest
#SBATCH --output=/work/nfp8/2026_03_03_fasttest/logs/post_fasttest.out
#SBATCH --error=/work/nfp8/2026_03_03_fasttest/logs/post_fasttest.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=64G
#SBATCH -c 4
#SBATCH --time=02:00:00

set -euo pipefail

# Source config
CONFIG="/hpc/home/nfp8/copilot/2026_03_03/test_mcf7/config_fasttest.sh"
source "$CONFIG"

module load "$MOD_SAMTOOLS"
module load "$MOD_SUBREAD"
module load "$MOD_BEDTOOLS"

# Single sample
SRR="${SAMPLES[0]}"
echo "[$(date -Is)] Post-processing ${SRR}"

FINAL_BAM="${BAM_DIR}/${SRR}_Aligned.sortedByCoord.out.bam"
if [[ ! -f "$FINAL_BAM" ]]; then
    echo "ERROR: BAM not found: $FINAL_BAM"
    exit 1
fi

DEDUP_BAM="${BAM_DIR}/${SRR}_dedup.bam"
TMP_DIR="${BAM_DIR}/tmp_postAlign_${SRR}_$$"
mkdir -p "$COVERAGE_DIR" "$COUNTS_DIR" "$TMP_DIR"

echo "[$(date -Is)] Step 1: MarkDuplicates..."
samtools collate -Ouf "$FINAL_BAM" "$TMP_DIR/collate" \
  | samtools fixmate -m -r -u - - \
  | samtools sort -u -T "$TMP_DIR/sort" - \
  | samtools markdup -s -d 100 - "$DEDUP_BAM"

samtools index "$DEDUP_BAM"
echo "[$(date -Is)]   Dedup BAM: $DEDUP_BAM"

# Read counts before/after dedup
RAW_COUNT=$(samtools view -c -F 0x904 "$FINAL_BAM")
DEDUP_COUNT=$(samtools view -c -F 0x904 "$DEDUP_BAM")
echo "[$(date -Is)]   Reads: raw=$RAW_COUNT, dedup=$DEDUP_COUNT"

# ===========================================================================
# 2. Flagstat on dedup BAM (primary mapped reads = CPM denominator)
# ===========================================================================
FLAGSTAT="${COVERAGE_DIR}/${SRR}_flagstat.txt"
echo "[$(date -Is)] Step 2: Running flagstat..."
samtools flagstat "$DEDUP_BAM" > "$FLAGSTAT"

# ===========================================================================
# 3. Read-through window depth — strand-aware, MAPQ >= 255 (uniquely mapped)
#
#   The non-reference L1 sits in the 5' region of NEDD4 on the PLUS strand.
#   Its sense promoter drives transcription downstream (increasing coords).
#   RT_WINDOW = 1 kb immediately 3' of the L1 insertion breakpoint.
#
#   STRANDEDNESS="reverse" (fr-firststrand / TruSeq Stranded):
#     Plus-strand reads = R2-forward (flag 128, NOT reverse) + R1-reverse (flag 80)
# ===========================================================================
RT_DEPTH="${COVERAGE_DIR}/${SRR}_readthrough.depth"
echo "[$(date -Is)] Step 3: Computing read-through depth in ${RT_WINDOW}..."

R2FWD="${TMP_DIR}/${SRR}_r2fwd.bam"
R1REV="${TMP_DIR}/${SRR}_r1rev.bam"
MERGED="${TMP_DIR}/${SRR}_plus_strand.bam"

samtools view -bh -q 255 -f 128 -F 16 "$DEDUP_BAM" "$RT_WINDOW" > "$R2FWD"
samtools view -bh -q 255 -f  80       "$DEDUP_BAM" "$RT_WINDOW" > "$R1REV"

samtools merge -f "$MERGED" "$R2FWD" "$R1REV"
samtools sort "$MERGED" -o "${MERGED%.bam}_sorted.bam"
samtools depth -a -d 0 -r "$RT_WINDOW" "${MERGED%.bam}_sorted.bam" > "$RT_DEPTH"
rm -f "$R2FWD" "$R1REV" "$MERGED" "${MERGED%.bam}_sorted.bam"

echo "    RT window: $(wc -l < "$RT_DEPTH") positions | $(awk '{s+=$3} END{printf "mean depth=%.3f", s/NR}' "$RT_DEPTH")"

# ===========================================================================
# 4. Background window depth (same strand logic; flanking normalization)
#    BACKGROUND_WINDOWS_BED: 10x1kb intronic windows flanking the L1 insertion
# ===========================================================================
BG_DEPTH="${COVERAGE_DIR}/${SRR}_background.depth"
echo "[$(date -Is)] Step 4: Computing background window depth..."

R2FWD_BG="${TMP_DIR}/${SRR}_bg_r2fwd.bam"
R1REV_BG="${TMP_DIR}/${SRR}_bg_r1rev.bam"
MERGED_BG="${TMP_DIR}/${SRR}_bg_plus.bam"

samtools view -bh -q 255 -f 128 -F 16 "$DEDUP_BAM" > "$R2FWD_BG"
samtools view -bh -q 255 -f  80       "$DEDUP_BAM" > "$R1REV_BG"

samtools merge -f "$MERGED_BG" "$R2FWD_BG" "$R1REV_BG"
samtools sort "$MERGED_BG" -o "${MERGED_BG%.bam}_sorted.bam"
samtools depth -a -d 0 -b "$BACKGROUND_WINDOWS_BED" "${MERGED_BG%.bam}_sorted.bam" > "$BG_DEPTH"
rm -f "$R2FWD_BG" "$R1REV_BG" "$MERGED_BG" "${MERGED_BG%.bam}_sorted.bam"

echo "    Background: $(wc -l < "$BG_DEPTH") positions | $(awk '{s+=$3} END{printf "mean depth=%.3f", s/NR}' "$BG_DEPTH")"

# ===========================================================================
# 5. Subset BAM to NEDD4 locus +/-1 Mb -> for IGV inspection
# ===========================================================================
SUBSET_BAM="${COVERAGE_DIR}/${SRR}_nedd4locus.bam"
echo "[$(date -Is)] Step 5: Subsetting BAM to NEDD4 locus (${BAM_SUBSET_REGION})..."
samtools view -b -h -q 255 "$DEDUP_BAM" "$BAM_SUBSET_REGION" > "$SUBSET_BAM"
samtools index "$SUBSET_BAM"
echo "    Subset BAM: $(du -sh "$SUBSET_BAM" | cut -f1)"

# Compute mean depths
RT_MEAN=$(awk '{s+=$3} END{if(NR>0) printf "%.3f", s/NR; else print "0"}' "$RT_DEPTH")
BG_MEAN=$(awk '{s+=$3} END{if(NR>0) printf "%.3f", s/NR; else print "0"}' "$BG_DEPTH")
echo "[$(date -Is)]   RT mean depth: $RT_MEAN"
echo "[$(date -Is)]   BG mean depth: $BG_MEAN"

echo "[$(date -Is)] Step 6: featureCounts..."
FC_OUT="${COUNTS_DIR}/${SRR}_featureCounts.txt"
featureCounts \
  -p --countReadPairs \
  -B \
  -a "$GTF" \
  -o "$FC_OUT" \
  -s 2 \
  -T 4 \
  "$DEDUP_BAM" 2>&1 | tee "${COUNTS_DIR}/${SRR}_featureCounts.log"

echo "[$(date -Is)] Step 7: Copying STAR GeneCounts..."
STAR_GENECOUNTS="${BAM_DIR}/${SRR}_ReadsPerGene.out.tab"
if [[ -f "$STAR_GENECOUNTS" ]]; then
    cp "$STAR_GENECOUNTS" "$COUNTS_DIR/"
fi

# Cleanup
rm -rf "$TMP_DIR"

echo "[$(date -Is)] Post-processing complete for ${SRR}"
