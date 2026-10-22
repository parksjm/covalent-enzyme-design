"""Export actual minimization snapshots for full-system PyMOL/VMD inspection."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.minimization_trajectory import export_trajectory  # noqa: E402


def main() -> None:
    """Write a fresh compressed PDB trajectory and source-hash report."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("topology", "trajectory", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    report = export_trajectory(args.topology, args.trajectory, args.output)
    logging.info(
        "Exported %d frames, %d atoms per frame",
        report["frames"],
        report["atoms_per_frame"],
    )


if __name__ == "__main__":
    main()
