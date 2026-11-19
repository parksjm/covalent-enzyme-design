"""Extract and plot endpoint potential energies from a coordinate-driving run."""

from __future__ import annotations

import argparse
import hashlib
import logging
import os
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.path_energies import (  # noqa: E402
    driven_path_energies,
    write_driven_path_energy_data,
)

LOGGER = logging.getLogger(__name__)


def _write_manifest(output: Path) -> None:
    files = sorted(output.glob("potential_energy.*"))
    lines = [
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}"
        for path in files
    ]
    (output / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")


def main() -> None:
    """Generate hashed JSON/CSV energy data and a PNG plot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--window-glob", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--plot-python", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    windows = sorted(path.name for path in args.run.glob(args.window_glob))
    rows = driven_path_energies(args.run, windows)
    args.output.mkdir(parents=True, exist_ok=True)
    json_path, _ = write_driven_path_energy_data(rows, args.output)
    plot_script = Path(__file__).with_name("plot_reaction_energies.py")
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    subprocess.run(
        [
            str(args.plot_python),
            str(plot_script),
            str(json_path),
            str(args.output / "potential_energy.png"),
            "--title",
            args.title,
            "--subtitle",
            "Reported Sander Etot and Etot minus reported RESTRAINT.\n"
            "First-window reference; no MM positional restraints; "
            "construction scan, not a validated barrier.",
        ],
        env=environment,
        check=True,
    )
    _write_manifest(args.output)
    LOGGER.info("Wrote %d-window energy profile to %s", len(rows), args.output)


if __name__ == "__main__":
    main()
