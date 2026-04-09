# Lessons — RoyalMD

_Hard-won lessons, gotchas, and things that broke before._
_This file is append-mostly. Only remove entries proven wrong._

## General

- **Never use local GPU (RTX 3060) for MD simulations or any heavy compute.** It crashes the system. Only use it for tiny tests (<100MB VRAM). For anything larger, write a SLURM script and submit to PACE cluster (account `gts-yke8`, user `dfu71`). VPN: `globalprotect connect --portal vpn.gatech.edu`.

## ProteinTTT / ESMFold

- **ESMFold requires `openfold` from GitHub.** The `fair-esm` PyPI package does NOT bundle openfold. You must install it separately: `pip install "openfold @ git+https://github.com/aqlaboratory/openfold.git"`. Without it, imports fail with `ModuleNotFoundError: No module named 'openfold'`.
- **PACE venv setup markers can hide broken installs.** If pip fails mid-install (e.g., ENOENT from corrupted site-packages), the venv may be left half-empty but the setup marker never gets written. Always verify critical imports before writing the marker. Version the marker (e.g., `.setup_complete_v2`) when changing deps.
- **AVB3 MSA depth has minimal impact on Protenix confidence.** Tested depth 0.05 (5% subsampled MSA) vs 1.00 (full MSA): ranking_score difference was <0.007. MSA depth may not be a useful lever for improving AVB3 predictions.
