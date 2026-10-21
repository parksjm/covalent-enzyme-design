"""Write a read-only bound on relay-water oxygen steric feasibility."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.relay_feasibility import (  # noqa: E402
    GRID_SETTINGS,
    assess_relay_feasibility,
)


def _sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_new_json(path: Path, report: dict[str, object]) -> None:
    """Create a JSON report without replacing an existing artifact."""
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")


def main() -> None:
    """Load topology and coordinates, then write the bounded diagnostic."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    topology = AmberParm(str(args.topology), xyz=str(args.coordinates))
    coordinates = np.asarray(topology.coordinates, dtype=float)
    report = assess_relay_feasibility(topology, coordinates)
    relay_source = ROOT / "src/nyl12_covalent/relay_water.py"
    feasibility_source = ROOT / "src/nyl12_covalent/relay_feasibility.py"
    report["provenance"] = {
        str(path.resolve()): _sha256(path)
        for path in (
            args.topology,
            args.coordinates,
            Path(__file__),
            relay_source,
            feasibility_source,
        )
    }
    report["grid_settings"] = GRID_SETTINGS
    _write_new_json(args.output, report)


if __name__ == "__main__":
    main()
