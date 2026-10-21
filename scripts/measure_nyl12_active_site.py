"""Measure active-site geometry from an explicit Amber topology/restart pair."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys
from typing import Sequence

from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.active_site_geometry import (  # noqa: E402
    coordinates_from_amber,
    measure_active_site,
)


TARGETS = {"pa6_l2": "pa6", "poly_xama_maxa": "poly"}
LOGGER = logging.getLogger(__name__)


def _sha256(path: Path) -> str:
    """Hash a file without changing it."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_parser() -> argparse.ArgumentParser:
    """Build the explicit-input command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--target", required=True, choices=sorted(TARGETS))
    parser.add_argument("--state", required=True, choices=("ti1", "es"))
    parser.add_argument("--output", required=True, type=Path)
    return parser


def measure_files(
    topology_path: Path,
    coordinates_path: Path,
    target: str,
    state: str,
    output_path: Path,
) -> dict[str, object]:
    """Measure an Amber model and create a provenance-rich JSON report.

    Existing output paths are rejected. Input paths must name regular files.
    """
    if target not in TARGETS:
        raise ValueError(f"Unsupported target: {target}")
    if state not in {"ti1", "es"}:
        raise ValueError("state must be explicitly set to 'ti1' or 'es'")
    topology_file = topology_path.resolve(strict=True)
    coordinates_file = coordinates_path.resolve(strict=True)
    if not topology_file.is_file() or not coordinates_file.is_file():
        raise ValueError("topology and coordinates must be regular files")
    destination = output_path.resolve()
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite existing report: {destination}")

    model = AmberParm(str(topology_file), xyz=str(coordinates_file))
    coordinates = coordinates_from_amber(model, TARGETS[target])
    report = measure_active_site(coordinates, TARGETS[target], state)
    report["target"] = target
    report["input_paths"] = {
        "topology": str(topology_file),
        "coordinates": str(coordinates_file),
    }
    report["input_hashes"] = {
        "topology_sha256": _sha256(topology_file),
        "coordinates_sha256": _sha256(coordinates_file),
    }
    report["source_hashes"] = {
        "measurement_script_sha256": _sha256(Path(__file__).resolve(strict=True)),
        "active_site_geometry_sha256": _sha256(
            ROOT / "src" / "nyl12_covalent" / "active_site_geometry.py"
        ),
        "reference_math_sha256": report["geometry_source_sha256"],
    }

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    LOGGER.info("Wrote %s active-site geometry report to %s", target, destination)
    return report


def main(argv: Sequence[str] | None = None) -> int:
    """Run the active-site measurement CLI."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    measure_files(
        args.topology,
        args.coordinates,
        args.target,
        args.state,
        args.output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
