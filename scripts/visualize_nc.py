#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import mdtraj as md


DEFAULT_TOPOLOGY_CANDIDATES = (
    "equilibrated.cif",
    "solvated.pdb",
    "minimized.cif",
    "equilibration.nc",
    "heating.nc",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert an OpenMM NetCDF trajectory to VMD-friendly formats.",
    )
    parser.add_argument(
        "nc_path",
        type=Path,
        help="Path to the NetCDF trajectory (e.g., production.nc).",
    )
    parser.add_argument(
        "--top",
        type=Path,
        help="Topology file (PDB/CIF). If omitted, the script tries common filenames.",
    )
    parser.add_argument(
        "--out-prefix",
        type=Path,
        default=Path("trajectory"),
        help="Output prefix for <prefix>.pdb and <prefix>.dcd.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Stride when sampling frames for visualization.",
    )
    parser.add_argument(
        "--single-pdb",
        action="store_true",
        help="Write a multi-model PDB instead of DCD (can be large).",
    )
    return parser.parse_args()


def resolve_topology(nc_path: Path, top_arg: Path | None) -> Path:
    if top_arg:
        return top_arg

    for candidate in DEFAULT_TOPOLOGY_CANDIDATES:
        candidate_path = nc_path.parent / candidate
        if candidate_path.exists():
            return candidate_path

    raise FileNotFoundError(
        "Topology not provided and no default topology files were found."
    )


def main() -> int:
    args = parse_args()

    if args.stride < 1:
        print("--stride must be >= 1", file=sys.stderr)
        return 2

    if not args.nc_path.exists():
        print(f"NetCDF file not found: {args.nc_path}", file=sys.stderr)
        return 2

    try:
        top_path = resolve_topology(args.nc_path, args.top)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    traj = md.load(str(args.nc_path), top=str(top_path))
    if args.stride > 1:
        traj = traj[:: args.stride]

    out_prefix = args.out_prefix
    top_pdb = out_prefix.with_suffix(".pdb")
    traj[0].save_pdb(str(top_pdb))
    print(f"Wrote topology PDB to {top_pdb}")

    if args.single_pdb:
        multi_pdb = out_prefix.with_name(f"{out_prefix.name}_all.pdb")
        traj.save_pdb(str(multi_pdb))
        print(f"Wrote multi-model PDB to {multi_pdb}")
        return 0

    dcd_path = out_prefix.with_suffix(".dcd")
    traj.save_dcd(str(dcd_path))
    print(f"Wrote DCD trajectory to {dcd_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
