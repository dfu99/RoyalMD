# 👑 Welcome to the RoyalMD

**RoyalMD** is a lightweight, educational molecular dynamics pipeline built on **OpenMM**, designed to run protein simulations on portable hardware. It automates on the fly common MD setup steps: from PDB fixing to solvation, minimization, and production runs, using a single command-line interface. The current beta version allows running molecular dynamics of any molecular system composed of proteins, nucleic acids and their complexes.

---

## ✨ Features
- Automatic **PDB fixing** (missing atoms, residues, protonation via `pdbfixer`)
- Flexible **force field selection** for protein & water models (AMBER, CHARMM)
- **Solvation and ion placement** with configurable box geometry
- Energy **minimization** and **NPT production runs**
- Flexible trajectory output and logging for benchmarks
- Designed for **Mac/Linux**

---

## 🪄 Requirements
- Python 3.12
- OpenMM
- pdbfixer
- mdtraj (for saving MD trajectories in the NetCDF format)

```bash
conda create -n MDsims python=3.12
conda activate MDsims
conda install -c conda-forge openmm pdbfixer mdtraj
conda install -c conda-forge "numpy<2.0" # <<< if you have any numpy-related errors"
```

## ⚜️ Ussage

```bash
./RoyalMD_demo.py your_structure.pdb
```

