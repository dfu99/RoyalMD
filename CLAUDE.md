# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Task Files

| File | When to consult |
|------|----------------|
| `tasks/planning.md` | Starting any session, checking priorities |
| `tasks/lessons.md` | Before touching subsystems they cover |

## PACE Cluster SLURM Rules

When writing SLURM scripts for the PACE cluster:

- **Account**: Always use `-A gts-yke8`
- **A100**: `--gres=gpu:A100:N` and **must** add `-C A100-80GB` constraint
- **RTX 6000**: `--gres=gpu:RTX_6000:N` (note underscore). No constraint needed.
- **H100**: `--gres=gpu:H100:N`. No constraint needed.
- **H200**: `--gres=gpu:H200:N`. No constraint needed.
- **Modules**: Always `module load cuda` for GPU jobs
- **Mail**: `--mail-type=END,FAIL` / `--mail-user=daniel.fu@emory.edu`
- **Paths**: scratch at `~/scratch/`, project storage at `~/p-yke8-0/`

## PACE Cluster SLURM Rules

When writing SLURM scripts for the PACE cluster:

- **Account**: Always use `-A gts-yke8`
- **A100**: `--gres=gpu:A100:N` and **must** add `-C A100-80GB` constraint
- **RTX 6000**: `--gres=gpu:RTX_6000:N` (note underscore). No constraint needed.
- **H100**: `--gres=gpu:H100:N`. No constraint needed.
- **H200**: `--gres=gpu:H200:N`. No constraint needed.
- **Modules**: Always `module load cuda` for GPU jobs
- **Mail**: `--mail-type=END,FAIL` / `--mail-user=daniel.fu@emory.edu`
- **Paths**: scratch at `~/scratch/`, project storage at `~/p-yke8-0/`
