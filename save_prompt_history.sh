#!/usr/bin/env bash
# ============================================================================
# save_prompt_history.sh — Save AI prompt history to a dated markdown file
# ============================================================================
# Writes a markdown log of all user prompts (verbatim) and one-line AI action
# summaries to: ${DEST_DIR}/${DATE}_AI_prompt_hx.md
#
# Usage (called by AI agent at end of session):
#   bash save_prompt_history.sh
#
# The script reads from stdin or a heredoc. The AI agent is responsible for
# supplying the content — each entry is a fenced block of the user's verbatim
# input followed by a single-line "AI action" summary with model name.
#
# Format of each entry fed to this script:
#   USER_PROMPT_START
#   <verbatim user text>
#   USER_PROMPT_END
#   AI_ACTION: <one-line summary> | Model: <model name>
#
# Example invocation by the AI agent:
#
#   bash save_prompt_history.sh <<'HISTORY'
#   USER_PROMPT_START
#   please update me on the progress of the array jobs
#   USER_PROMPT_END
#   AI_ACTION: Checked SLURM queue, diagnosed DependencyNeverSatisfied on trim job, installed cutadapt, resubmitted pipeline chain | Model: Claude Opus 4
#   HISTORY
# ============================================================================
set -euo pipefail

DATE="${DATE:-$(date +%Y_%m_%d)}"
DEST_DIR="${DEST_DIR:-/hpc/home/nfp8/copilot/2026_03_03}"
OUTFILE="${DEST_DIR}/${DATE}_AI_prompt_hx.md"

mkdir -p "$DEST_DIR"

# If the file doesn't exist, write the header
if [[ ! -f "$OUTFILE" ]]; then
    cat > "$OUTFILE" <<HEADER
# AI Prompt History — ${DATE}

> Auto-generated log of user prompts (verbatim) and AI action summaries.
> Each entry records exactly what the user typed and a one-line description
> of the AI's response, including the model used.

---

HEADER
    echo "Created ${OUTFILE}"
fi

# Read entries from stdin and append
ENTRY_NUM=$(grep -c '^### Prompt ' "$OUTFILE" 2>/dev/null || true)
ENTRY_NUM=${ENTRY_NUM:-0}
# Strip whitespace to avoid arithmetic errors
ENTRY_NUM=$(echo "$ENTRY_NUM" | tr -d '[:space:]')
TIMESTAMP=""
IN_PROMPT=false
PROMPT_TEXT=""

while IFS= read -r line; do
    if [[ "$line" == "USER_PROMPT_START" ]]; then
        IN_PROMPT=true
        PROMPT_TEXT=""
        ENTRY_NUM=$((ENTRY_NUM + 1))
        TIMESTAMP=$(date -Is)
        continue
    fi

    if [[ "$line" == "USER_PROMPT_END" ]]; then
        IN_PROMPT=false
        # Write the prompt block
        {
            echo "### Prompt ${ENTRY_NUM}  —  ${TIMESTAMP}"
            echo ""
            echo '```'
            echo "$PROMPT_TEXT"
            echo '```'
            echo ""
        } >> "$OUTFILE"
        continue
    fi

    if $IN_PROMPT; then
        if [[ -z "$PROMPT_TEXT" ]]; then
            PROMPT_TEXT="$line"
        else
            PROMPT_TEXT="${PROMPT_TEXT}
${line}"
        fi
        continue
    fi

    # AI_ACTION line
    if [[ "$line" == AI_ACTION:* ]]; then
        {
            echo "**${line}**"
            echo ""
            echo "---"
            echo ""
        } >> "$OUTFILE"
    fi
done

echo "Appended ${ENTRY_NUM} total entries to ${OUTFILE}"
