"""Analyze intrinsic proton-transfer stationary points and barrier estimates."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
import re
from typing import Any

import numpy as np

from prepare_approximate_inputs import read_geometry

HARTREE_TO_KCAL_MOL = 627.5094740631
LOGGER = logging.getLogger(__name__)
DFT_PATTERN = re.compile(r"SCF Done:\s+E\([^)]+\)\s+=\s*([-+0-9.Ee]+)")
DFTB_PATTERN = re.compile(r"Total SCF energy\s*=\s*([-+0-9.Ee]+)\s+kcal/mol")
XTB_PATTERN = re.compile(r"\| TOTAL ENERGY\s+([-+0-9.Ee]+)\s+Eh")
ZPE_PATTERN = re.compile(r"Zero-point correction=\s*([-+0-9.Ee]+)")
FREQUENCY_PATTERN = re.compile(r"Frequencies --\s+([^\n]+)")
METHODS = ("dftb3", "gfn2_xtb", "wb97xd_augccpvtz")


def _read_output_energy(path: Path, pattern: re.Pattern[str], marker: str) -> float:
    """Read the final energy and require an engine termination marker."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if marker not in text:
        raise RuntimeError(f"Missing termination marker in {path}")
    matches = pattern.findall(text)
    if not matches:
        raise ValueError(f"Energy not found in {path}")
    return float(matches[-1])


def _gaussian_evidence(path: Path, state: str) -> dict[str, Any]:
    """Extract stationary-point, frequency, energy, and ZPE evidence."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if text.count("Normal termination of Gaussian 16") != 2:
        raise RuntimeError(f"Expected two normal Gaussian links in {path}")
    if "Stationary point found" not in text:
        raise RuntimeError(f"No optimized stationary point in {path}")
    frequencies = [
        float(value) for row in FREQUENCY_PATTERN.findall(text) for value in row.split()
    ]
    imaginary = [value for value in frequencies if value < -20.0]
    expected = 1 if state == "transition_state" else 0
    if len(imaginary) != expected:
        raise RuntimeError(
            f"Expected {expected} significant imaginary frequencies in {path}; "
            f"found {imaginary}"
        )
    zpe = ZPE_PATTERN.findall(text)
    if len(zpe) != 1:
        raise RuntimeError(f"Expected one ZPE correction in {path}; found {zpe}")
    energies = DFT_PATTERN.findall(text)
    if not energies:
        raise ValueError(f"No DFT energy in {path}")
    return {
        "electronic_energy_kcal_mol": float(energies[-1]) * HARTREE_TO_KCAL_MOL,
        "zpe_kcal_mol": float(zpe[0]) * HARTREE_TO_KCAL_MOL,
        "imaginary_frequencies_cm-1": imaginary,
    }


def _distance(
    coordinates: list[tuple[float, float, float]], first: int, second: int
) -> float:
    """Return a distance between one-based atom indices in angstroms."""
    return float(
        np.linalg.norm(np.asarray(coordinates[first - 1]) - coordinates[second - 1])
    )


def _irc_evidence(root: Path, record: dict[str, Any]) -> dict[str, Any]:
    """Require opposite proton identities at the two IRC endpoints."""
    calculation = root / record["directory"]
    indices = record["reactive_atom_indices_1based"]
    endpoints: dict[str, Any] = {}
    identities: set[str] = set()
    for direction in ("forward", "reverse"):
        log = calculation / f"irc_{direction}.log"
        if "Normal termination of Gaussian 16" not in log.read_text(errors="replace"):
            raise RuntimeError(f"IRC did not terminate normally: {log}")
        _, coordinates = read_geometry(calculation / f"irc_{direction}.fchk")
        donor_h = _distance(coordinates, indices["donor"], indices["proton"])
        acceptor_h = _distance(coordinates, indices["acceptor"], indices["proton"])
        identity = "donor" if donor_h < acceptor_h else "acceptor"
        identities.add(identity)
        endpoints[direction] = {
            "donor_H_A": donor_h,
            "acceptor_H_A": acceptor_h,
            "proton_identity": identity,
        }
    if identities != {"donor", "acceptor"}:
        raise RuntimeError(
            f"IRC endpoints do not connect both proton identities: {endpoints}"
        )
    return endpoints


def _stationary_geometry(root: Path, record: dict[str, Any]) -> dict[str, Any]:
    """Measure and validate proton identity in an optimized stationary point."""
    calculation = root / record["directory"]
    _, coordinates = read_geometry(calculation / "state.fchk")
    indices = record["reactive_atom_indices_1based"]
    donor_h = _distance(coordinates, indices["donor"], indices["proton"])
    acceptor_h = _distance(coordinates, indices["acceptor"], indices["proton"])
    identity = "donor" if donor_h < acceptor_h else "acceptor"
    expected = {"reactant": "donor", "product": "acceptor"}.get(record["state"])
    if expected is not None and identity != expected:
        raise RuntimeError(
            f"{record['pair']} {record['state']} optimized to the wrong proton "
            f"identity: {identity}; d(D-H)={donor_h:.3f}, d(A-H)={acceptor_h:.3f} A"
        )
    return {
        "donor_H_A": donor_h,
        "acceptor_H_A": acceptor_h,
        "proton_identity": identity,
    }


def analyze(root: Path, output: Path) -> None:
    """Validate stationary points and write matched barrier estimates.

    Args:
        root: Completed intrinsic-barrier benchmark directory.
        output: New analysis directory.
    """
    root = root.resolve(strict=True)
    manifest = json.loads((root / "benchmark_manifest.json").read_text())
    output.mkdir(parents=True, exist_ok=False)
    states: dict[str, dict[str, dict[str, Any]]] = {}
    for record in manifest["records"]:
        calculation = root / record["directory"]
        pair = record["pair"]
        state = record["state"]
        gaussian = _gaussian_evidence(calculation / "opt_freq_sp.log", state)
        states.setdefault(pair, {})[state] = {
            "wb97xd_augccpvtz": gaussian["electronic_energy_kcal_mol"],
            "zpe": gaussian["zpe_kcal_mol"],
            "dftb3": _read_output_energy(
                calculation / "dftb3.out", DFTB_PATTERN, "Calculation Completed"
            ),
            "gfn2_xtb": _read_output_energy(
                calculation / "gfn2.out", XTB_PATTERN, "normal termination of xtb"
            )
            * HARTREE_TO_KCAL_MOL,
            "imaginary_frequencies_cm-1": gaussian["imaginary_frequencies_cm-1"],
            "geometry": _stationary_geometry(root, record),
        }

    ts_records = [
        record
        for record in manifest["records"]
        if record["state"] == "transition_state"
    ]
    irc = {record["pair"]: _irc_evidence(root, record) for record in ts_records}
    results: dict[str, Any] = {}
    table: list[dict[str, Any]] = []
    for pair, values in states.items():
        pair_result: dict[str, Any] = {
            "irc_endpoints": irc[pair],
            "stationary_geometries": {
                state: values[state]["geometry"] for state in values
            },
            "methods": {},
        }
        for method in METHODS:
            reactant = values["reactant"][method]
            transition_state = values["transition_state"][method]
            product = values["product"][method]
            method_result = {
                "forward_barrier_kcal_mol": transition_state - reactant,
                "reverse_barrier_kcal_mol": transition_state - product,
                "reaction_energy_kcal_mol": product - reactant,
            }
            if method == "wb97xd_augccpvtz":
                method_result.update(
                    {
                        "forward_barrier_zpe_kcal_mol": (
                            transition_state
                            + values["transition_state"]["zpe"]
                            - reactant
                            - values["reactant"]["zpe"]
                        ),
                        "reaction_energy_zpe_kcal_mol": (
                            product
                            + values["product"]["zpe"]
                            - reactant
                            - values["reactant"]["zpe"]
                        ),
                    }
                )
            pair_result["methods"][method] = method_result
            table.append({"pair": pair, "method": method, **method_result})
        reference = pair_result["methods"]["wb97xd_augccpvtz"]
        pair_result["barrier_errors_kcal_mol"] = {
            method: {
                "forward": pair_result["methods"][method]["forward_barrier_kcal_mol"]
                - reference["forward_barrier_kcal_mol"],
                "reverse": pair_result["methods"][method]["reverse_barrier_kcal_mol"]
                - reference["reverse_barrier_kcal_mol"],
            }
            for method in METHODS[:2]
        }
        results[pair] = pair_result

    (output / "barrier_metrics.json").write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )
    fieldnames = sorted({key for row in table for key in row})
    with (output / "barrier_metrics.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(table)
    LOGGER.info("Validated and analyzed %d intrinsic barriers", len(results))


def main() -> None:
    """Parse explicit completed-root and output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    analyze(args.root, args.output)


if __name__ == "__main__":
    main()
