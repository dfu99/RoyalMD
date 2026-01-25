#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import mdtraj as md


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove water molecules from a PDB before RoyalMD input.",
    )
    parser.add_argument("pdb_path", type=Path, help="Input PDB file.")
    parser.add_argument(
        "--out",
        type=Path,
        help="Output PDB file (default: <input>_nowat.pdb).",
    )
    return parser.parse_args()


def default_output_path(pdb_path: Path) -> Path:
    return pdb_path.with_name(f"{pdb_path.stem}_nowat{pdb_path.suffix}")


def main() -> int:
    args = parse_args()

    if not args.pdb_path.exists():
        print(f"PDB file not found: {args.pdb_path}", file=sys.stderr)
        return 2

    out_path = args.out or default_output_path(args.pdb_path)

    traj = md.load(str(args.pdb_path))
    water_indices = traj.topology.select("water")
    if water_indices.size == 0:
        print("No waters detected; writing original structure.")
        traj.save_pdb(str(out_path))
        print(f"Wrote {out_path}")
        return 0

    keep_indices = traj.topology.select("not water")
    stripped = traj.atom_slice(keep_indices)
    stripped.save_pdb(str(out_path))
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
