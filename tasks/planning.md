# Planning — RoyalMD

## Current Priorities

1. **Monitor ProteinTTT resubmit** — Job 5389714 queued on PACE (A100). Fixed venv install to include `openfold` from GitHub. Check logs once it starts running.

## Next Steps

- Analyze ProteinTTT results once job completes
- Compare TTT-enhanced ESMFold predictions against Protenix MSA predictions

## Recently Completed

1. **AVB3 ProteinTTT job failures diagnosed and fixed** — Root cause: broken venv missing `openfold` (required by ESMFold). Fixed submit script, nuked broken venv, resubmitted as job 5389714. (2026-03-23)
2. **MSA validation results fetched** — 10 CIF predictions (2 depths x 5 samples) pulled to `results-msa-validation/`. Confidence nearly identical between depth 0.05 and 1.00 (ranking_score ~0.733-0.738). (2026-03-23)
3. **GPU safety rule documented** — Added local GPU crash warning to `tasks/lessons.md`. All GPU work must go to PACE cluster. (2026-03-13)
