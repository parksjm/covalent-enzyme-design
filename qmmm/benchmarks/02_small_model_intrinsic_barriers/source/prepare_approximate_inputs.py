"""Convert Gaussian formatted checkpoints into matched SQM and xTB inputs."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import re

BOHR_TO_ANGSTROM = 0.529177210903
LOGGER = logging.getLogger(__name__)
ARRAY_HEADER = re.compile(r"^(.+?)\s+[IR]\s+N=\s*(\d+)\s*$")
SYMBOLS = {
    1: "H",
    6: "C",
    7: "N",
    8: "O",
}


def _arrays(path: Path) -> dict[str, list[float]]:
    """Read numeric arrays needed from a Gaussian formatted checkpoint."""
    lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    arrays: dict[str, list[float]] = {}
    line_index = 0
    while line_index < len(lines):
        match = ARRAY_HEADER.match(lines[line_index])
        if match is None:
            line_index += 1
            continue
        name, count_text = match.groups()
        count = int(count_text)
        values: list[float] = []
        line_index += 1
        while len(values) < count:
            values.extend(
                float(value.replace("D", "E")) for value in lines[line_index].split()
            )
            line_index += 1
        arrays[name.strip()] = values[:count]
    return arrays


def read_geometry(path: Path) -> tuple[list[int], list[tuple[float, float, float]]]:
    """Read atomic numbers and Cartesian coordinates from an FCHK file."""
    arrays = _arrays(path)
    atomic_numbers = [int(value) for value in arrays["Atomic numbers"]]
    flat = arrays["Current cartesian coordinates"]
    if len(flat) != 3 * len(atomic_numbers):
        raise ValueError(f"Coordinate count mismatch in {path}")
    coordinates = [
        tuple(value * BOHR_TO_ANGSTROM for value in flat[index : index + 3])
        for index in range(0, len(flat), 3)
    ]
    unsupported = sorted(set(atomic_numbers) - SYMBOLS.keys())
    if unsupported:
        raise ValueError(f"Unsupported atomic numbers in {path}: {unsupported}")
    return atomic_numbers, coordinates


def _write_xyz(
    path: Path,
    atomic_numbers: list[int],
    coordinates: list[tuple[float, float, float]],
) -> None:
    """Write an XYZ geometry for standalone GFN2-xTB."""
    lines = [str(len(atomic_numbers)), "DFT-optimized stationary-point geometry"]
    lines.extend(
        f"{SYMBOLS[number]:2s} {x: .10f} {y: .10f} {z: .10f}"
        for number, (x, y, z) in zip(atomic_numbers, coordinates, strict=True)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_sqm(
    path: Path,
    atomic_numbers: list[int],
    coordinates: list[tuple[float, float, float]],
    charge: int,
) -> None:
    """Write a DFTB3/3OB SQM single-point input."""
    lines = [
        "DFT-optimized intrinsic-barrier stationary point",
        " &qmmm",
        "  qm_theory='DFTB3',",
        f"  qmcharge={charge},",
        "  scfconv=1.0d-10, errconv=1.0d-6, tight_p_conv=1,",
        "  itrmax=500, maxcyc=0, verbosity=0,",
        " /",
    ]
    lines.extend(
        f" {number:2d} {SYMBOLS[number]:2s} {x: .10f} {y: .10f} {z: .10f}"
        for number, (x, y, z) in zip(atomic_numbers, coordinates, strict=True)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare(root: Path) -> None:
    """Create matched approximate-method inputs for every stationary point.

    Args:
        root: Generated barrier-benchmark input/output root.
    """
    manifest = json.loads((root / "benchmark_manifest.json").read_text())
    for record in manifest["records"]:
        calculation = root / record["directory"]
        atomic_numbers, coordinates = read_geometry(calculation / "state.fchk")
        _write_xyz(calculation / "optimized.xyz", atomic_numbers, coordinates)
        _write_sqm(
            calculation / "dftb3.in",
            atomic_numbers,
            coordinates,
            record["charge"],
        )
    LOGGER.info(
        "Prepared approximate-method inputs for %d states", len(manifest["records"])
    )


def main() -> None:
    """Parse the benchmark root and prepare inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    prepare(args.root.resolve(strict=True))


if __name__ == "__main__":
    main()
