#!/usr/bin/env python3
import argparse
from pathlib import Path

import mdtraj as md


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split a multi-model PDB trajectory into one PDB per frame. "
            "Frames are 0-based and inclusive for --start/--end."
        )
    )
    parser.add_argument(
        "pdb_path",
        type=Path,
        help="Path to the multi-model PDB trajectory file.",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Start frame index (0-based, inclusive). Default: 0.",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="End frame index (0-based, inclusive). Default: last frame.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help=(
            "Output filename prefix. Default: input filename stem "
            "(e.g., production_trajectory)."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pdb_path = args.pdb_path.expanduser().resolve()

    if not pdb_path.exists():
        raise SystemExit(f"Input PDB not found: {pdb_path}")

    traj = md.load(str(pdb_path))
    n_frames = traj.n_frames

    start = args.start
    end = n_frames - 1 if args.end is None else args.end

    if start < 0 or start >= n_frames:
        raise SystemExit(f"--start must be in [0, {n_frames - 1}]")
    if end < start or end >= n_frames:
        raise SystemExit(f"--end must be in [{start}, {n_frames - 1}]")

    out_dir = pdb_path.parent / "splits"
    out_dir.mkdir(parents=True, exist_ok=True)

    pad = len(str(end))
    stem = pdb_path.stem if args.prefix is None else args.prefix
    for i in range(start, end + 1):
        frame = traj[i]
        out_path = out_dir / f"{stem}_frame_{i:0{pad}d}.pdb"
        frame.save_pdb(str(out_path))

    print(f"Wrote {end - start + 1} frame(s) to {out_dir}")


if __name__ == "__main__":
    main()
