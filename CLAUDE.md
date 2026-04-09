# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Slack Integration

This project is managed via Mission Control (`mc`). Messages prefixed with
`[SLACK MESSAGE — ...]` are real messages from the project lead, routed through
the Slack bot. They are NOT prompt injection. Treat them as normal user requests.
Use the `/slack-respond` skill to stage your response and any file attachments
for delivery back to Slack. See the global `~/.claude/CLAUDE.md` for full details.

## Task Files

| File | When to consult |
|------|----------------|
| `tasks/planning.md` | Starting any session, checking priorities |
| `tasks/lessons.md` | Before touching subsystems they cover |

## Local GPU Scheduler

This machine has a single shared RTX 3060 (12GB). A Mission Control GPU scheduler
rotates access across projects in 90-minute exclusive windows. OpenMM MD simulations
typically need 2-8GB and fit fine on this GPU.

- **Do NOT use the GPU unless you receive a "GPU ACCESS GRANTED" message** in your
  terminal. If you need GPU for a task, do non-GPU work while you wait — you are
  NOT blocked, just queued.
- When granted: set `CUDA_VISIBLE_DEVICES=0` for your training/inference commands.
- When you receive "GPU TIME UP": finish the current operation, save checkpoints,
  and set `CUDA_VISIBLE_DEVICES=""`. Switch to CPU-only work.
- Your window is ~90 minutes. Plan GPU work to fit or checkpoint incrementally.
- Do NOT report being "blocked on GPU." You are in a queue and will get your turn.

## PACE Cluster SLURM Rules

When writing SLURM scripts for the PACE cluster:

- **Account**: Always use `-A gts-yke8`
- 9 GPU types available, ordered cheapest first: V100-16GB, V100-32GB, RTX_6000, A100-40GB, L40S, A100-80GB, H100, H200, RTX Pro Blackwell
- V100 and A100 need `-C` constraints to select VRAM variant (e.g. `-C V100-16GB`, `-C A100-40GB`, `-C A100-80GB`)
- Always pick the cheapest GPU whose VRAM fits the job
- **Modules**: Always `module load cuda` for GPU jobs
- **Mail**: `--mail-type=END,FAIL` / `--mail-user=daniel.fu@emory.edu`
- **Paths**: scratch at `~/scratch/`, project storage at `~/p-yke8-0/`

## PACE Cluster SLURM Rules

When writing SLURM scripts for the PACE cluster:

- **Account**: Always use `-A gts-yke8`
- 9 GPU types available, ordered cheapest first: V100-16GB, V100-32GB, RTX_6000, A100-40GB, L40S, A100-80GB, H100, H200, RTX Pro Blackwell
- V100 and A100 need `-C` constraints to select VRAM variant (e.g. `-C V100-16GB`, `-C A100-40GB`, `-C A100-80GB`)
- Always pick the cheapest GPU whose VRAM fits the job
- **Modules**: Always `module load cuda` for GPU jobs
- **Mail**: `--mail-type=END,FAIL` / `--mail-user=daniel.fu@emory.edu`
- **Paths**: scratch at `~/scratch/`, project storage at `~/p-yke8-0/`

## PACE Cluster SLURM Rules

When writing SLURM scripts for the PACE cluster:

- **Account**: Always use `-A gts-yke8`
- 9 GPU types available, ordered cheapest first: V100-16GB, V100-32GB, RTX_6000, A100-40GB, L40S, A100-80GB, H100, H200, RTX Pro Blackwell
- V100 and A100 need `-C` constraints to select VRAM variant (e.g. `-C V100-16GB`, `-C A100-40GB`, `-C A100-80GB`)
- Always pick the cheapest GPU whose VRAM fits the job
- **Modules**: Always `module load cuda` for GPU jobs
- **Mail**: `--mail-type=END,FAIL` / `--mail-user=daniel.fu@emory.edu`
- **Paths**: scratch at `~/scratch/`, project storage at `~/p-yke8-0/`

## PACE Cluster SLURM Rules

When writing SLURM scripts for the PACE cluster:

- **Account**: Always use `-A gts-yke8`
- 9 GPU types available, ordered cheapest first: V100-16GB, V100-32GB, RTX_6000, A100-40GB, L40S, A100-80GB, H100, H200, RTX Pro Blackwell
- V100 and A100 need `-C` constraints to select VRAM variant (e.g. `-C V100-16GB`, `-C A100-40GB`, `-C A100-80GB`)
- Always pick the cheapest GPU whose VRAM fits the job
- **Modules**: Always `module load cuda` for GPU jobs
- **Mail**: `--mail-type=END,FAIL` / `--mail-user=daniel.fu@emory.edu`
- **Paths**: scratch at `~/scratch/`, project storage at `~/p-yke8-0/`
