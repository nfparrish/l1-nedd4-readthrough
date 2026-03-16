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
# Uniquely-mapped MAPQ threshold: STAR=255, TopHat2=50. Configurable via UNIQUE_MAPQ.
UNIQUE_MAPQ="${UNIQUE_MAPQ:-255}"
# Paper-style MAPQ threshold: Philippe et al. use MAPQ>=20 for igvtools count.
# This recovers junction-flanking read pairs where one mate maps into the L1 body.
PAPER_MAPQ="${PAPER_MAPQ:-20}"

# 3. Read-through window depth — strand-aware, MAPQ >= UNIQUE_MAPQ (uniquely mapped)
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
        samtools depth -a -d 0 -Q "$UNIQUE_MAPQ" -r "$RT_WINDOW" "$DEDUP_BAM" > "$RT_DEPTH"

    else
        R2FWD="${TMP_DIR}/${SRR}_r2fwd.bam"
        R1REV="${TMP_DIR}/${SRR}_r1rev.bam"
        MERGED="${TMP_DIR}/${SRR}_plus_strand.bam"

        if [[ "$STRANDEDNESS" == "reverse" ]]; then
            # TruSeq Stranded (fr-firststrand): plus-strand reads
            samtools view -bh -q "$UNIQUE_MAPQ" -f 128 -F 16 "$DEDUP_BAM" "$RT_WINDOW" > "$R2FWD"
            samtools view -bh -q "$UNIQUE_MAPQ" -f  80       "$DEDUP_BAM" "$RT_WINDOW" > "$R1REV"
        else
            # fr-secondstrand: plus-strand reads
            samtools view -bh -q "$UNIQUE_MAPQ" -f  64 -F 16 "$DEDUP_BAM" "$RT_WINDOW" > "$R2FWD"
            samtools view -bh -q "$UNIQUE_MAPQ" -f 144       "$DEDUP_BAM" "$RT_WINDOW" > "$R1REV"
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
# 3b. Same-window minus-strand depth — opposite-strand background view for the
#     same genomic interval. Mirrors step 3 flag logic but selects minus-strand
#     RNA, which at this locus is largely host-gene background.
#
#   STRANDEDNESS="reverse" (fr-firststrand): minus-strand RNA reads are
#     - R1 forward (flag: read1 [0x40] + NOT reverse [NOT 0x10]) -> minus-strand RNA
#     - R2 reverse (flag: read2 [0x80] + reverse [0x10] = 0x90)  -> minus-strand RNA
#
#   STRANDEDNESS="forward" (fr-secondstrand): minus-strand RNA reads are
#     - R1 reverse (flag: read1 [0x40] + reverse [0x10] = 0x50)  -> minus-strand RNA
#     - R2 forward (flag: read2 [0x80] + NOT reverse [NOT 0x10]) -> minus-strand RNA
#
#   STRANDEDNESS="unstranded": reads aligning to minus genomic orientation
#     (flag 0x10 set); not interpretable as RNA strand-of-origin but useful
#     as a proxy for minus-strand coverage at this locus.
# ===========================================================================
RT_MINUS_DEPTH="${COVERAGE_DIR}/${SRR}_readthrough_minus.depth"
if [[ ! -f "$RT_MINUS_DEPTH" ]]; then
    echo "  [3b] Computing same-window minus-strand background depth (${STRANDEDNESS} library)..."

    if [[ "$STRANDEDNESS" == "unstranded" ]]; then
        MMINUS="${TMP_DIR}/${SRR}_minus_all.bam"
        samtools view -bh -q "$UNIQUE_MAPQ" -f 16 "$DEDUP_BAM" "$RT_WINDOW" \
            | samtools sort -T "${TMP_DIR}/sort_minus_${SRR}" - -o "$MMINUS"
        samtools index "$MMINUS"
        samtools depth -a -d 0 -r "$RT_WINDOW" "$MMINUS" > "$RT_MINUS_DEPTH"
        rm -f "$MMINUS" "${MMINUS}.bai" || true

    else
        MR1="${TMP_DIR}/${SRR}_minus_r1.bam"
        MR2="${TMP_DIR}/${SRR}_minus_r2.bam"
        MMERGED="${TMP_DIR}/${SRR}_minus_strand.bam"

        if [[ "$STRANDEDNESS" == "reverse" ]]; then
            # fr-firststrand: minus-strand RNA -> R1 forward + R2 reverse
            samtools view -bh -q "$UNIQUE_MAPQ" -f  64 -F 16 "$DEDUP_BAM" "$RT_WINDOW" > "$MR1"
            samtools view -bh -q "$UNIQUE_MAPQ" -f 144       "$DEDUP_BAM" "$RT_WINDOW" > "$MR2"
        else
            # fr-secondstrand: minus-strand RNA -> R1 reverse + R2 forward
            samtools view -bh -q "$UNIQUE_MAPQ" -f  80       "$DEDUP_BAM" "$RT_WINDOW" > "$MR1"
            samtools view -bh -q "$UNIQUE_MAPQ" -f 128 -F 16 "$DEDUP_BAM" "$RT_WINDOW" > "$MR2"
        fi

        samtools merge -f "$MMERGED" "$MR1" "$MR2"
        samtools sort "$MMERGED" -o "${MMERGED%.bam}_sorted.bam"
        samtools index "${MMERGED%.bam}_sorted.bam"
        samtools depth -a -d 0 -r "$RT_WINDOW" "${MMERGED%.bam}_sorted.bam" > "$RT_MINUS_DEPTH"
        rm -f "$MR1" "$MR2" "$MMERGED" "${MMERGED%.bam}_sorted.bam" "${MMERGED%.bam}_sorted.bam.bai" || true
    fi

    echo "    RT minus: $(wc -l < "$RT_MINUS_DEPTH") positions | $(awk '{s+=$3} END{printf "mean depth=%.3f", s/NR}' "$RT_MINUS_DEPTH")"
fi

# ===========================================================================
# 3c. Paper-style stranded windows using mate 1 only (Philippe et al. method).
#     - Uses PAPER_MAPQ=20 (igvtools default) not UNIQUE_MAPQ, so junction-
#       flanking reads where a mate lands in the L1 body are retained.
#     - Uses raw FINAL_BAM (no dedup filter) matching paper's igvtools workflow.
#     - Counts R1 only (one tag per fragment) on the correct RNA strand.
#     - Also writes a companion R1_MAPQ20 count for CPM normalization.
#
#   STRANDEDNESS="reverse" (fr-firststrand / TruSeq Stranded):
#     - downstream sense RNA (plus strand):   R1 reverse  (-f 80)
#     - upstream antisense RNA (minus strand): R1 forward  (-f 64 -F 16)
#
#   STRANDEDNESS="forward" (fr-secondstrand):
#     - downstream sense RNA (plus strand):   R1 forward  (-f 64 -F 16)
#     - upstream antisense RNA (minus strand): R1 reverse  (-f 80)
#
#   STRANDEDNESS="unstranded": skipped.
# ===========================================================================
PAPER_SENSE_DEPTH="${COVERAGE_DIR}/${SRR}_paper_downstream_sense.depth"
PAPER_ANTISENSE_DEPTH="${COVERAGE_DIR}/${SRR}_paper_upstream_antisense.depth"
PAPER_R1_COUNT_FILE="${COVERAGE_DIR}/${SRR}_r1_mapq20_count.txt"

# Write R1 MAPQ>=20 count for CPM denominator (once, regardless of strandedness)
if [[ ! -f "$PAPER_R1_COUNT_FILE" ]]; then
    echo "  [3c-pre] Counting R1 reads (MAPQ>=${PAPER_MAPQ}) for CPM denominator..."
    samtools view -c -q "$PAPER_MAPQ" -f 64 "$FINAL_BAM" 2>/dev/null > "$PAPER_R1_COUNT_FILE"
    echo "    R1 count: $(cat $PAPER_R1_COUNT_FILE)"
fi

if [[ ! -f "$PAPER_SENSE_DEPTH" || ! -f "$PAPER_ANTISENSE_DEPTH" ]]; then
    echo "  [3c] Computing paper-style mate1 stranded windows (${STRANDEDNESS}, MAPQ>=${PAPER_MAPQ}, pre-dedup)..."

    if [[ "$STRANDEDNESS" == "unstranded" ]]; then
        echo "    Skipping paper-style windows for unstranded library"
    else
        PAPER_SENSE_BAM="${TMP_DIR}/${SRR}_paper_sense_mate1.bam"
        PAPER_ANTISENSE_BAM="${TMP_DIR}/${SRR}_paper_antisense_mate1.bam"

        if [[ "$STRANDEDNESS" == "reverse" ]]; then
            # fr-firststrand: R1-reverse = plus-strand RNA (sense downstream)
            samtools view -bh -q "$PAPER_MAPQ" -f 80 "$FINAL_BAM" "$PAPER_DOWNSTREAM_WINDOW" > "$PAPER_SENSE_BAM"
            # fr-firststrand: R1-forward = minus-strand RNA (antisense upstream)
            samtools view -bh -q "$PAPER_MAPQ" -f 64 -F 16 "$FINAL_BAM" "$PAPER_UPSTREAM_WINDOW" > "$PAPER_ANTISENSE_BAM"
        else
            # fr-secondstrand: R1-forward = plus-strand RNA (sense downstream)
            samtools view -bh -q "$PAPER_MAPQ" -f 64 -F 16 "$FINAL_BAM" "$PAPER_DOWNSTREAM_WINDOW" > "$PAPER_SENSE_BAM"
            # fr-secondstrand: R1-reverse = minus-strand RNA (antisense upstream)
            samtools view -bh -q "$PAPER_MAPQ" -f 80 "$FINAL_BAM" "$PAPER_UPSTREAM_WINDOW" > "$PAPER_ANTISENSE_BAM"
        fi

        samtools sort "$PAPER_SENSE_BAM" -o "${PAPER_SENSE_BAM%.bam}_sorted.bam"
        samtools index "${PAPER_SENSE_BAM%.bam}_sorted.bam"
        samtools depth -a -d 0 -r "$PAPER_DOWNSTREAM_WINDOW" "${PAPER_SENSE_BAM%.bam}_sorted.bam" > "$PAPER_SENSE_DEPTH"

        samtools sort "$PAPER_ANTISENSE_BAM" -o "${PAPER_ANTISENSE_BAM%.bam}_sorted.bam"
        samtools index "${PAPER_ANTISENSE_BAM%.bam}_sorted.bam"
        samtools depth -a -d 0 -r "$PAPER_UPSTREAM_WINDOW" "${PAPER_ANTISENSE_BAM%.bam}_sorted.bam" > "$PAPER_ANTISENSE_DEPTH"

        rm -f \
            "$PAPER_SENSE_BAM" "${PAPER_SENSE_BAM%.bam}_sorted.bam" "${PAPER_SENSE_BAM%.bam}_sorted.bam.bai" \
            "$PAPER_ANTISENSE_BAM" "${PAPER_ANTISENSE_BAM%.bam}_sorted.bam" "${PAPER_ANTISENSE_BAM%.bam}_sorted.bam.bai" || true

        echo "    Paper downstream sense: $(wc -l < "$PAPER_SENSE_DEPTH") positions | $(awk '{s+=$3} END{printf "mean depth=%.3f", s/NR}' "$PAPER_SENSE_DEPTH")"
        echo "    Paper upstream antisense: $(wc -l < "$PAPER_ANTISENSE_DEPTH") positions | $(awk '{s+=$3} END{printf "mean depth=%.3f", s/NR}' "$PAPER_ANTISENSE_DEPTH")"
    fi
fi

# ===========================================================================
# 4. Background window depth (same strand logic; flanking normalization in Step 4)
#    BACKGROUND_WINDOWS_BED built by 01_build_star_index.sh (NEDD4 exons subtracted)
# ===========================================================================
BG_DEPTH="${COVERAGE_DIR}/${SRR}_background.depth"
if [[ ! -f "$BG_DEPTH" && -f "${BACKGROUND_WINDOWS_BED}" ]]; then
    echo "  [4] Computing background window depth..."

    if [[ "$STRANDEDNESS" == "unstranded" ]]; then
        samtools depth -a -d 0 -Q "$UNIQUE_MAPQ" -b "$BACKGROUND_WINDOWS_BED" "$DEDUP_BAM" > "$BG_DEPTH"

    else
        R2FWD="${TMP_DIR}/${SRR}_bg_r2fwd.bam"
        R1REV="${TMP_DIR}/${SRR}_bg_r1rev.bam"
        MERGED="${TMP_DIR}/${SRR}_bg_plus.bam"

        if [[ "$STRANDEDNESS" == "reverse" ]]; then
            samtools view -bh -q "$UNIQUE_MAPQ" -f 128 -F 16 "$DEDUP_BAM" > "$R2FWD"
            samtools view -bh -q "$UNIQUE_MAPQ" -f  80       "$DEDUP_BAM" > "$R1REV"
        else
            samtools view -bh -q "$UNIQUE_MAPQ" -f  64 -F 16 "$DEDUP_BAM" > "$R2FWD"
            samtools view -bh -q "$UNIQUE_MAPQ" -f 144       "$DEDUP_BAM" > "$R1REV"
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
    samtools view -b -h -q "$UNIQUE_MAPQ" "$DEDUP_BAM" "$BAM_SUBSET_REGION" > "$SUBSET_BAM"
    samtools index "$SUBSET_BAM"
    echo "    Subset BAM: $(du -sh "$SUBSET_BAM" | cut -f1)"
fi

# ===========================================================================
# 6. featureCounts gene quantification (strand-aware, PE-aware)
#    FC_SCOPE (set in config.sh) controls what annotation is used:
#      SKIP   — skip featureCounts entirely
#      PANEL  — ~100 housekeeping genes (FC_PANEL_GTF); fast, good for QC/normalization
#      CHR15  — all chr15 genes (FC_CHR15_GTF); good balance for NEDD4-focused runs
#      GENOME — full GENCODE annotation (GENCODE_GTF); default, most complete
# ===========================================================================
FC_SCOPE="${FC_SCOPE:-GENOME}"
# Scope-tagged filename lets GENOME and PANEL outputs coexist on disk
FC_OUT="${COUNTS_DIR}/${SRR}_featureCounts_${FC_SCOPE}.txt"

if [[ "$FC_SCOPE" == "SKIP" ]]; then
    echo "  [6] featureCounts SKIP (FC_SCOPE=SKIP)"
elif [[ ! -f "$FC_OUT" ]]; then
    echo "  [6] Running featureCounts (FC_SCOPE=${FC_SCOPE})..."

    case "$FC_SCOPE" in
        PANEL)  FC_GTF="${FC_PANEL_GTF:-${REF_DIR}/housekeeping_panel.gtf}" ;;
        CHR15)  FC_GTF="${FC_CHR15_GTF:-${REF_DIR}/gencode.v45.chr15.gtf}" ;;
        GENOME) FC_GTF="${GENCODE_GTF}" ;;
        *)      echo "  WARNING: unknown FC_SCOPE=${FC_SCOPE}, defaulting to GENOME" >&2
                FC_GTF="${GENCODE_GTF}" ;;
    esac

    if [[ ! -f "$FC_GTF" ]]; then
        echo "  WARNING: GTF not found: $FC_GTF — skipping featureCounts" >&2
    else
        FC_STRAND=0
        [[ "$STRANDEDNESS" == "forward" ]] && FC_STRAND=1
        [[ "$STRANDEDNESS" == "reverse" ]] && FC_STRAND=2

        FC_FLAGS="-T ${STAR_THREADS} -a ${FC_GTF} -s ${FC_STRAND}"
        [[ "$LIBRARY_LAYOUT" == "PE" ]] && FC_FLAGS="${FC_FLAGS} -p --countReadPairs -B"

        featureCounts \
            ${FC_FLAGS} \
            -o "$FC_OUT" \
            "$DEDUP_BAM" \
            2> "${COUNTS_DIR}/${SRR}_featureCounts.log"

        echo "    featureCounts: $(grep 'Assigned' "${COUNTS_DIR}/${SRR}_featureCounts.log" | tail -1 || echo 'see log')"
    fi
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
