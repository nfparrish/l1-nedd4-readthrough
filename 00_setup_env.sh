#!/bin/bash
#SBATCH --job-name=setup_env
#SBATCH --output=/hpc/home/nfp8/copilot/2026_03_03/logs/setup_env_%j.out
#SBATCH --error=/hpc/home/nfp8/copilot/2026_03_03/logs/setup_env_%j.err
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=4G
#SBATCH -c 2
#SBATCH --time=01:00:00
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${PIPELINE_CONFIG:-${SCRIPT_DIR}/config.sh}"

echo "[$(date -Is)] Setting up Python virtual environment and directory tree"

# ---- Create all directories ----
for d in "$BAM_DIR" "$COVERAGE_DIR" "$COUNTS_DIR" "$RESULTS_DIR" "$LOG_DIR" \
         "$PERSISTENT_BAM" "$PERSISTENT_RESULTS" \
         "$REF_DIR" "${STAR_INDEX}" \
         "${SCRIPT_DIR}/logs"; do
    mkdir -p "$d"
done

# ---- Python venv with pandas, matplotlib, seaborn, multiqc ----
if [[ ! -f "${VENV_DIR}/bin/python" ]]; then
    echo "Creating Python venv at ${VENV_DIR}"
    python3 -m venv "${VENV_DIR}"
fi

# Use venv Python directly (no activate script sourcing)
"${VENV_DIR}/bin/pip" install --upgrade pip
# Note: cutadapt pinned to 3.5 for compatibility with Trim Galore 0.4.3 (4.x+ incompatible)
"${VENV_DIR}/bin/pip" install pandas matplotlib seaborn multiqc 'cutadapt==3.5'

echo "[$(date -Is)] Python venv ready at ${VENV_DIR}"
echo "  Packages: $(${VENV_DIR}/bin/pip list --format=columns 2>/dev/null | grep -iE 'pandas|matplotlib|seaborn|multiqc' | tr '\n' ' ')"

# ---- Create chr15 GTF subset (for FC_SCOPE=CHR15) ----
FC_CHR15_GTF="${FC_CHR15_GTF:-${REF_DIR}/gencode.v45.chr15.gtf}"
if [[ -f "${GENCODE_GTF}" && ! -f "${FC_CHR15_GTF}" ]]; then
    echo "[$(date -Is)] Creating chr15 GTF subset..."
    grep -E '^(chr15|#)' "${GENCODE_GTF}" > "${FC_CHR15_GTF}"
    echo "  chr15 GTF: $(grep -vc '^#' "${FC_CHR15_GTF}") features -> ${FC_CHR15_GTF}"
fi

# ---- Create housekeeping panel GTF (for FC_SCOPE=PANEL) ----
# ~100 housekeeping and signaling genes spanning diverse biotypes:
# constitutive (ribosomal, glycolytic, chaperones), chromatin, cell-cycle,
# tumor suppressors, and key pathway markers for cross-cohort QC.
FC_PANEL_GTF="${FC_PANEL_GTF:-${REF_DIR}/housekeeping_panel.gtf}"
if [[ -f "${GENCODE_GTF}" && ! -f "${FC_PANEL_GTF}" ]]; then
    echo "[$(date -Is)] Creating housekeeping panel GTF..."
    PANEL_GENES=(
        # Constitutive / structural
        ACTB ACTG1 TUBB TUBA1A TUBA1B
        # Glycolysis
        GAPDH LDHA LDHB ALDOA ENO1 ENO2 PKM TPI1 PGAM1 PGK1
        # Ribosomes
        RPLP0 RPLP1 RPLP2 RPL13A RPL27A RPL32 RPS18 RPS6 RPS13
        # Translation initiation / elongation
        EEF1A1 EEF2 EIF4A1 EIF4E EIF4G1
        # Housekeeping enzymes
        HPRT1 HMBS SDHA UBC TBP PPIA B2M YWHAZ
        # RNA binding / splicing
        HNRNPA1 HNRNPC PCBP1 SFPQ SRSF4 PABPC1 PABPC4 PAIP1
        # Chaperones / UPR
        CANX CALR HSP90AA1 HSP90AB1 HSPA1A HSPA1B HSPA5 HSPA8
        # Chromatin / DNA repair
        VCP NPM1 NONO XRCC5 XRCC6 PCNA
        # HIF / hypoxia
        HIF1A EPAS1 VEGFA
        # Proliferation / cell cycle
        CCND1 CCND3 CDK4 CDK6 CDKN1A CDKN1B MKI67
        # Tumor suppressors
        TP53 MDM2 BRCA1 BRCA2
        # MYC family
        MYC MYCN
        # Oncogenic signaling
        EGFR ERBB2 MET KRAS NRAS BRAF MAP2K1
        AKT1 MTOR PIK3CA
        STAT3 STAT5A JAK1 JAK2
        CTNNB1 APC AXIN1
        NOTCH1 HES1
        # Nutrient transporters
        TFRC SLC2A1
        # Stress / UPR transcription factors
        ATF4 ATF6 XBP1 DDIT3
    )
    # Build grep pattern matching gene_name "SYMBOL" in GTF attribute field
    PATTERN=$(printf 'gene_name "%s"\\|' "${PANEL_GENES[@]}" | sed 's/\\|$//')
    { grep '^#' "${GENCODE_GTF}"; grep -E "${PATTERN}" "${GENCODE_GTF}"; } > "${FC_PANEL_GTF}"
    echo "  Panel GTF: ${#PANEL_GENES[@]} genes -> $(grep -vc '^#' "${FC_PANEL_GTF}") features -> ${FC_PANEL_GTF}"
fi

echo "[$(date -Is)] Setup complete"
