# Copilot Agent Notes

## File Storage Policy

Only save files to the following approved locations:
- `/hpc/home/nfp8` — home directory (25 GiB total, ~24 GiB free)
- `/hpc/group/parrishlab` — parrishlab group storage (1.0 TiB total, ~974 GiB free)
- `/work/nfp8` — work scratch space (20 TiB total, ~20 TiB free)

**Do NOT save files to `/hpc/group/kolab`** — this volume is critically full (95% used, ~2 TiB remaining).

## Confirmed Permissions (as of 2026-03-02)

| Location | mkdir | read | write | Notes |
|----------|-------|------|-------|-------|
| `/hpc/home/nfp8` | ✅ | ✅ | ✅ | Accessible from login and compute nodes |
| `/hpc/group/parrishlab` | ⚠️ | ⚠️ | ⚠️ | Requires SLURM compute node; read-only on login node |
| `/work/nfp8` | ⚠️ | ⚠️ | ⚠️ | Requires SLURM compute node; read-only on login node |

> `/hpc/group/parrishlab` and `/work/nfp8` are only fully accessible from compute nodes
> (inside a SLURM job via `srun --pty bash` or `sbatch`). They appear read-only on login nodes.

## Storage Guidelines (from DCC user guide)

- `/hpc/home/nfp8` (25 GB) — personal scripts, config, environment setup. No purge.
- `/hpc/group/parrishlab` (1 TB) — persistent lab data, software, long-lived results. 7-day snapshots. No auto-purge.
- `/work/nfp8` (shared 650 TB) — large datasets under active analysis. **Files auto-purged after 75 days.**

Use `/work/nfp8` for active large-scale computation; use `/hpc/group/parrishlab` for anything that must persist longer than 75 days.

---

## Prompt History Logging (MANDATORY)

**Immediately before** summarizing the conversation (i.e., when the context window limit is approaching and a summary will be generated), the agent **MUST** write a prompt history file by running:

```bash
bash /hpc/home/nfp8/copilot/2026_03_03/save_prompt_history.sh <<'HISTORY'
USER_PROMPT_START
<verbatim text of user prompt 1>
USER_PROMPT_END
AI_ACTION: <one-line summary of what the AI did> | Model: <model name>
USER_PROMPT_START
<verbatim text of user prompt 2>
USER_PROMPT_END
AI_ACTION: <one-line summary> | Model: <model name>
...repeat for every user prompt in the session...
HISTORY
```

### Rules

1. **Every user prompt** in the session must be included **verbatim** — do not paraphrase, truncate, or omit any.
2. Each prompt is wrapped in `USER_PROMPT_START` / `USER_PROMPT_END` markers.
3. Each prompt is followed by exactly one `AI_ACTION:` line that contains:
   - A concise single-line summary of the AI's action/response
   - The model name (e.g., `Claude Opus 4`)
4. The output file is `${DATE}_AI_prompt_hx.md` in `/hpc/home/nfp8/copilot/2026_03_03/` (the `DATE` variable defaults to today's date as `YYYY_MM_DD`).
5. This step must be performed **immediately before** the conversation summary is written — not after, not at any other time.
6. If the file already exists (e.g., from a prior session segment on the same day), new entries are **appended** — the script handles this automatically.

### Why

This creates an auditable, versioned record of every instruction given to the AI agent and the action taken, independent of the ephemeral conversation context. It survives context window rollovers and can be committed to Git alongside the codebase.
