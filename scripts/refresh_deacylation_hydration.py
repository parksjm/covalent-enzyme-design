#!/usr/bin/env python3
"""Refresh accepted hydration-segment water analysis and its portable movie."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
import subprocess
import sys

LOGGER = logging.getLogger(__name__)


def refresh(base: Path, project_root: Path) -> None:
    """Run ordered analysis and export, propagating either failure.

    Args:
        base: Existing hydration campaign with accepted native segments.
        project_root: Explicit project containing the pinned analysis/export code.
    """
    base = base.resolve(strict=True)
    project_root = project_root.resolve(strict=True)
    environment = os.environ.copy()
    environment.update(
        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    commands = [
        [
            sys.executable,
            str(project_root / "scripts/analyze_deacylation_hydration.py"),
            "--base",
            str(base),
            "--output",
            str(base / "analysis"),
        ],
        [
            sys.executable,
            str(project_root / "scripts/export_hydration_vmd.py"),
            "--source",
            str(base),
            "--output",
            str(base / "vmd_progress"),
            "--water-summary",
            str(base / "analysis/per_water_summary.json"),
            "--cpptraj",
            "/data/jpc/code/ambertools26/venv/bin/cpptraj",
        ],
    ]
    for command in commands:
        LOGGER.info("Run hydration refresh: %s", command)
        subprocess.run(command, cwd=project_root, env=environment, check=True)


def main() -> None:
    """Parse the explicit campaign and project used by a pinned MD hook."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    refresh(args.base, args.project_root)


if __name__ == "__main__":
    main()
