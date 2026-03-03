#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --output=/hpc/home/nfp8/copilot/2026_03_03/logs/star_index_%j.out
#SBATCH --error=/hpc/home/nfp8/copilot/2026_03_03/logs/star_index_%j.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=48G
#SBATCH -c 12
#SBATCH --time=04:00:00
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"
module load "$MOD_STAR"

echo "[$(date -Is)] Building STAR index"
echo "  Genome:  ${GENOME_FA}"
echo "  GTF:     ${GENCODE_GTF}"
echo "  Index:   ${STAR_INDEX}"
echo "  Overhang: ${STAR_SJDB_OVERHANG}"

mkdir -p "${REF_DIR}" "${STAR_INDEX}"

# ---- Download reference genome if missing ----
if [[ ! -f "${GENOME_FA}" ]]; then
    echo "[$(date -Is)] Downloading GRCh38 primary assembly..."
    wget -q -O "${GENOME_FA}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
    gunzip "${GENOME_FA}.gz"
    echo "[$(date -Is)] Genome downloaded: $(du -sh "${GENOME_FA}" | cut -f1)"
fi

# ---- Download GENCODE annotation if missing ----
if [[ ! -f "${GENCODE_GTF}" ]]; then
    echo "[$(date -Is)] Downloading GENCODE v45 GTF..."
    wget -q -O "${GENCODE_GTF}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz"
    gunzip "${GENCODE_GTF}.gz"
    echo "[$(date -Is)] GTF downloaded: $(du -sh "${GENCODE_GTF}" | cut -f1)"
fi

# ---- Generate STAR index ----
if [[ ! -f "${STAR_INDEX}/Genome" ]]; then
    echo "[$(date -Is)] Running STAR --runMode genomeGenerate..."
    STAR \
        --runMode genomeGenerate \
        --runThreadN "${INDEX_CPUS}" \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GENCODE_GTF}" \
        --sjdbOverhang "${STAR_SJDB_OVERHANG}"
    echo "[$(date -Is)] STAR index built: $(du -sh "${STAR_INDEX}" | cut -f1)"
else
    echo "[$(date -Is)] STAR index already exists — skipping"
fi

# ---- Generate NEDD4 exon BED from GENCODE GTF (for background normalization) ----
if [[ -f "${GENCODE_GTF}" && ! -f "${NEDD4_EXONS_BED}" ]]; then
    echo "[$(date -Is)] Extracting NEDD4 exon coordinates from GTF..."
    module load "${MOD_BEDTOOLS}" 2>/dev/null || true
    awk -v gene="${NEDD4_ENSEMBL_ID}" '
        $3=="exon" && $0 ~ gene {
            split($0, a, "gene_id \""); split(a[2], b, "\"");
            if (b[1] ~ gene) print $1"\t"$4-1"\t"$5
        }
    ' "${GENCODE_GTF}" | sort -k1,1 -k2,2n | uniq > "${NEDD4_EXONS_BED}"
    echo "  $(wc -l < "${NEDD4_EXONS_BED}") NEDD4 exon intervals written to ${NEDD4_EXONS_BED}"
fi

# ---- Generate flanking background windows BED (10 × 1 kb, exons subtracted) ----
if [[ -f "${NEDD4_EXONS_BED}" && ! -f "${BACKGROUND_WINDOWS_BED}" ]]; then
    echo "[$(date -Is)] Building background windows BED..."
    module load "${MOD_BEDTOOLS}" 2>/dev/null || true
    # Candidate windows: 5 upstream (pre-L1 body) + 5 downstream (post-RT window)
    # L1 body ≈ chr15:55,952,934–55,958,934; RT window: chr15:55,958,936–55,959,936
    cat > "${BACKGROUND_WINDOWS_BED}.candidates" << 'BEDEOF'
chr15	55920000	55921000	upstream_1
chr15	55922000	55923000	upstream_2
chr15	55924000	55925000	upstream_3
chr15	55926000	55927000	upstream_4
chr15	55928000	55929000	upstream_5
chr15	55960000	55961000	downstream_1
chr15	55962000	55963000	downstream_2
chr15	55964000	55965000	downstream_3
chr15	55966000	55967000	downstream_4
chr15	55968000	55969000	downstream_5
BEDEOF
    # Subtract NEDD4 exon positions; keep windows with ≥ 500 bp remaining
    bedtools subtract -a "${BACKGROUND_WINDOWS_BED}.candidates" \
                      -b "${NEDD4_EXONS_BED}" \
        | awk '($3-$2) >= 500' \
        > "${BACKGROUND_WINDOWS_BED}"
    rm -f "${BACKGROUND_WINDOWS_BED}.candidates"
    echo "  $(wc -l < "${BACKGROUND_WINDOWS_BED}") background windows retained after exon subtraction"
fi

echo "[$(date -Is)] Done"
