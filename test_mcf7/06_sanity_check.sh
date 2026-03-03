#!/bin/bash
#SBATCH --job-name=check_mcf7
#SBATCH --partition=common
#SBATCH --account=parrishlab
#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH --time=00:10:00
set -euo pipefail

source "${PIPELINE_CONFIG}"

VENV_PYTHON="${VENV_DIR}/bin/python3"

echo "=================================================================="
echo " L1-NEDD4 MCF7 Sanity Check (positive + negative controls)"
echo " $(date)"
echo "=================================================================="
echo " Positive controls (untreated):       ERR973734, ERR973735"
echo " Negative controls (ORF1 shRNA sh1085): ERR973728, ERR973729"
echo "   Rationale: shRNA targets L1 ORF1; RISC degrades L1 transcript"
echo "   including the read-through region, so RT CPM should be lower."
echo "=================================================================="

RT_CSV="${RESULTS_DIR}/readthrough_signal.csv"
if [[ ! -f "$RT_CSV" ]]; then
    echo "FAIL: readthrough_signal.csv not found at ${RT_CSV}" >&2
    exit 1
fi

${VENV_PYTHON} - <<'PYEOF'
import os, sys
import pandas as pd

csv = os.environ["RESULTS_DIR"] + "/readthrough_signal.csv"
df = pd.read_csv(csv)

POSITIVE_CTRL = ["ERR973734", "ERR973735"]   # MCF-7 untreated
NEGATIVE_CTRL = ["ERR973728", "ERR973729"]   # MCF-7 ORF1 shRNA sh1085

print(df[["srr", "rt_mean_depth", "total_mapped", "rt_cpm"]].to_string(index=False))
print()

fails = []
warns = []

# ── Check 1: all 4 expected samples present ─────────────────────────────────
for acc in POSITIVE_CTRL + NEGATIVE_CTRL:
    if acc not in df["srr"].values:
        fails.append(f"FAIL: {acc} missing from results CSV")

if fails:                          # can't do downstream comparisons if rows missing
    for msg in fails:
        print(msg, file=sys.stderr)
    print(f"\nRESULT: FAIL ({len(fails)} issue(s))", file=sys.stderr)
    sys.exit(1)

pos_df = df[df["srr"].isin(POSITIVE_CTRL)]
neg_df = df[df["srr"].isin(NEGATIVE_CTRL)]

# ── Check 2: adequate mapping depth for all samples ─────────────────────────
for _, row in df.iterrows():
    if row["total_mapped"] < 1_000_000:
        fails.append(
            f"FAIL: {row['srr']} has only {row['total_mapped']:,} mapped reads "
            f"(expected >1 M)"
        )

# ── Check 3: positive controls show detectable RT CPM ───────────────────────
pos_cpm_mean = pos_df["rt_cpm"].mean()
if pos_cpm_mean == 0:
    fails.append(
        "FAIL: rt_cpm is 0 for ALL positive controls — no read-through signal "
        "detected in untreated MCF-7. Check STRANDEDNESS setting "
        f"(current: {os.environ.get('STRANDEDNESS','?')}) and coordinate window."
    )
elif pos_cpm_mean < 0.001:
    warns.append(
        f"WARN: positive controls have very low mean RT CPM ({pos_cpm_mean:.6f}) "
        "— signal present but unexpectedly weak. Verify STRANDEDNESS and "
        "RT_WINDOW coordinates."
    )

# ── Check 4: negative controls are not higher than positive controls ─────────
neg_cpm_mean = neg_df["rt_cpm"].mean()
if pos_cpm_mean > 0 and neg_cpm_mean >= pos_cpm_mean:
    warns.append(
        f"WARN: negative controls (mean rt_cpm={neg_cpm_mean:.6f}) are >= "
        f"positive controls (mean rt_cpm={pos_cpm_mean:.6f}). The ORF1 shRNA "
        "knockdown may be incomplete, or consider whether ORF1 shRNA is "
        "sufficient to suppress this read-through locus. Scramble shRNA "
        "controls (ERR973731-733) would disambiguate."
    )

# ── Report ───────────────────────────────────────────────────────────────────
print("── Positive controls (untreated) ──")
print(pos_df[["srr", "rt_cpm"]].to_string(index=False))
print(f"  Mean RT CPM: {pos_cpm_mean:.6f}")
print()
print("── Negative controls (ORF1 shRNA sh1085) ──")
print(neg_df[["srr", "rt_cpm"]].to_string(index=False))
print(f"  Mean RT CPM: {neg_cpm_mean:.6f}")
print()

if pos_cpm_mean > 0 and neg_cpm_mean < pos_cpm_mean:
    print(f"  Fold-reduction (pos/neg): {pos_cpm_mean/neg_cpm_mean:.2f}x")

for msg in warns:
    print(msg)

if fails:
    print()
    for msg in fails:
        print(msg, file=sys.stderr)
    print(f"\nRESULT: FAIL ({len(fails)} issue(s))", file=sys.stderr)
    sys.exit(1)
else:
    status = "PASS" if not warns else "PASS (with warnings)"
    print(f"\nRESULT: {status}")
    print(f"  Positive controls show RT CPM > negative controls: "
          f"{'YES' if pos_cpm_mean > neg_cpm_mean else 'NO'}")
PYEOF

echo "Check complete."
