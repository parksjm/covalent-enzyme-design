"""Prepare unconstrained intrinsic proton-transfer barrier benchmark inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

LOGGER = logging.getLogger(__name__)
TARGET_DISTANCE_A = 2.7
STATE_POINTS = {
    "reactant": 0,
    "transition_state": 4,
    "product": 8,
}


def _read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    """Read symbols and coordinates from an XYZ file."""
    lines = path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[0])
    atoms = [line.split() for line in lines[2 : 2 + atom_count]]
    if len(atoms) != atom_count:
        raise ValueError(f"Incomplete XYZ file: {path}")
    symbols = [atom[0] for atom in atoms]
    coordinates = np.asarray([[float(value) for value in atom[1:4]] for atom in atoms])
    return symbols, coordinates


def _reactive_indices(coordinates: np.ndarray, donor_h_a: float) -> dict[str, int]:
    """Recover one-based reactive indices from the matched scan construction."""
    targets = {
        "donor": np.asarray([0.0, 0.0, 0.0]),
        "acceptor": np.asarray([TARGET_DISTANCE_A, 0.0, 0.0]),
        "proton": np.asarray([donor_h_a, 0.0, 0.0]),
    }
    indices: dict[str, int] = {}
    for name, target in targets.items():
        distances = np.linalg.norm(coordinates - target, axis=1)
        index = int(np.argmin(distances))
        if distances[index] > 1.0e-5:
            raise ValueError(
                f"Unable to recover {name} index; nearest={distances[index]:.3g}"
            )
        indices[name] = index + 1
    if len(set(indices.values())) != 3:
        raise ValueError(f"Reactive atom indices are not unique: {indices}")
    return indices


def _coordinates_block(symbols: list[str], coordinates: np.ndarray) -> list[str]:
    """Render a Gaussian Cartesian coordinate block."""
    return [
        f"{symbol:2s} {x: .10f} {y: .10f} {z: .10f}"
        for symbol, (x, y, z) in zip(symbols, coordinates, strict=True)
    ]


def _render_gaussian(
    symbols: list[str],
    coordinates: np.ndarray,
    charge: int,
    state: str,
    pair: str,
) -> str:
    """Render optimization/frequency and triple-zeta single-point links."""
    if state == "transition_state":
        optimization = "Opt=(TS,CalcFC,NoEigenTest,MaxCycles=200)"
    else:
        optimization = "Opt=(CalcFC,MaxCycles=200)"
    lines = [
        "%mem=8GB",
        "%nprocshared=4",
        "%chk=state.chk",
        (
            f"#p wB97XD/def2SVP {optimization} Freq "
            "SCF=(XQC,Tight) Integral=UltraFine NoSymm"
        ),
        "",
        f"{pair} {state}: unconstrained intrinsic proton-transfer benchmark",
        "",
        f"{charge} 1",
        *_coordinates_block(symbols, coordinates),
        "",
        "--Link1--",
        "%oldchk=state.chk",
        "%chk=state_sp.chk",
        "%mem=8GB",
        "%nprocshared=4",
        (
            "#p wB97XD/aug-cc-pVTZ SP Geom=AllCheck Guess=Read "
            "SCF=(XQC,Tight) Integral=UltraFine NoSymm"
        ),
        "",
    ]
    return "\n".join(lines) + "\n"


def _render_irc(direction: str) -> str:
    """Render an IRC input that reads the accepted TS Hessian checkpoint."""
    return (
        "\n".join(
            [
                "%oldchk=state.chk",
                f"%chk=irc_{direction}.chk",
                "%mem=8GB",
                "%nprocshared=4",
                (
                    "#p wB97XD/def2SVP "
                    f"IRC=(RCFC,{direction.title()},MaxPoints=40,StepSize=5) "
                    "Geom=AllCheck Guess=Read SCF=(XQC,Tight) Integral=UltraFine NoSymm"
                ),
                "",
            ]
        )
        + "\n"
    )


def prepare(scan_root: Path, output: Path) -> None:
    """Prepare three-state inputs from the validated fixed-geometry scan.

    Args:
        scan_root: Existing 81-point matched fixed-geometry benchmark.
        output: New barrier-benchmark input directory.
    """
    scan_root = scan_root.resolve(strict=True)
    source_manifest = json.loads((scan_root / "benchmark_manifest.json").read_text())
    output.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, Any]] = []
    systems = [system["name"] for system in source_manifest["systems"]]
    for pair in systems:
        pair_records = {
            record["point"]: record
            for record in source_manifest["records"]
            if record["pair"] == pair
            and np.isclose(record["donor_acceptor_A"], TARGET_DISTANCE_A)
        }
        if set(pair_records) != set(range(9)):
            raise ValueError(f"Incomplete R={TARGET_DISTANCE_A:.1f} scan for {pair}")
        common_indices: dict[str, int] | None = None
        for state, point in STATE_POINTS.items():
            source = pair_records[point]
            source_xyz = scan_root / source["directory"] / "geometry.xyz"
            symbols, coordinates = _read_xyz(source_xyz)
            indices = _reactive_indices(coordinates, source["donor_H_A"])
            if common_indices is None:
                common_indices = indices
            elif indices != common_indices:
                raise ValueError(f"Reactive mapping changed for {pair}: {indices}")
            calculation = output / pair / state
            calculation.mkdir(parents=True)
            geometry = calculation / "initial.xyz"
            geometry.write_bytes(source_xyz.read_bytes())
            gaussian = calculation / "opt_freq_sp.gjf"
            gaussian.write_text(
                _render_gaussian(symbols, coordinates, source["charge"], state, pair),
                encoding="utf-8",
            )
            (calculation / "charge.txt").write_text(
                f"{source['charge']}\n", encoding="utf-8"
            )
            if state == "transition_state":
                for direction in ("forward", "reverse"):
                    (calculation / f"irc_{direction}.gjf").write_text(
                        _render_irc(direction), encoding="utf-8"
                    )
            records.append(
                {
                    "pair": pair,
                    "state": state,
                    "charge": source["charge"],
                    "directory": str(calculation.relative_to(output)),
                    "source_geometry": str(source_xyz),
                    "source_point": point,
                    "reactive_atom_indices_1based": indices,
                    "sha256": {
                        "initial_xyz": hashlib.sha256(
                            geometry.read_bytes()
                        ).hexdigest(),
                        "gaussian_input": hashlib.sha256(
                            gaussian.read_bytes()
                        ).hexdigest(),
                    },
                }
            )
    manifest = {
        "purpose": "Unconstrained intrinsic gas-phase proton-transfer barriers",
        "geometry_method": "omegaB97X-D/def2-SVP",
        "reference_single_point": "omegaB97X-D/aug-cc-pVTZ",
        "approximate_methods": ["DFTB3/3OB", "GFN2-xTB"],
        "donor_acceptor_seed_A": TARGET_DISTANCE_A,
        "systems": systems,
        "states": list(STATE_POINTS),
        "records": records,
    }
    (output / "benchmark_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info("Prepared %d stationary-point guesses in %s", len(records), output)


def main() -> None:
    """Parse explicit source and destination paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    prepare(args.scan_root, args.output)


if __name__ == "__main__":
    main()
