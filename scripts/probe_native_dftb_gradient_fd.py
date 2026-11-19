"""Prepare and assess finite-difference checks of a native DFTB/MM gradient."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

if __package__:
    from scripts.prepare_dimer_continuation import (
        parse_option,
        read_xyz_frames,
        write_restart,
    )
else:
    from prepare_dimer_continuation import parse_option, read_xyz_frames, write_restart

LOGGER = logging.getLogger(__name__)
HARTREE_TO_KCAL_MOL = 627.5095
BOHR_TO_ANGSTROM = 0.529177210903


def _record(path: Path) -> dict[str, Any]:
    """Return source path, SHA-256 and byte size."""
    data = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
    }


def displace_active(
    coordinates: np.ndarray, active_indices: np.ndarray, mode: np.ndarray, step: float
) -> np.ndarray:
    """Displace only active atoms by a normalized mode step.

    Args:
        coordinates: Full (natom, 3) starting coordinates in Angstrom.
        active_indices: Zero-based active topology indices in XYZ file order.
        mode: Final active-coordinate vector from path_tsmode.xyz.
        step: Signed displacement norm in Angstrom.

    Returns:
        A copied full coordinate array with only the active region displaced.
    """
    if mode.shape != (len(active_indices), 3) or not np.isfinite(mode).all():
        raise ValueError("Mode array must be finite and match the active atom count")
    norm = float(np.linalg.norm(mode))
    if norm <= 1.0e-12 or not np.isfinite(norm) or not np.isfinite(step):
        raise ValueError("Mode norm and displacement must be finite and nonzero")
    displaced = np.array(coordinates, dtype=float, copy=True)
    if displaced.shape[1:] != (3,) or not np.isfinite(displaced).all():
        raise ValueError("Starting coordinates must be finite Cartesian vectors")
    displaced[active_indices] += step * mode / norm
    return displaced


def prepare_fd(
    case_base: Path,
    force_calculation: Path,
    output: Path,
    step: float = 0.1,
    center_offset: float = 0.0,
) -> dict[str, Any]:
    """Write paired restarts around an offset along the saved dimer mode."""
    if output.exists():
        raise FileExistsError(output)
    if not np.isfinite(step) or step <= 0 or not np.isfinite(center_offset):
        raise ValueError("Finite positive step and finite center offset required")
    case_base = case_base.resolve(strict=True)
    force_calculation = force_calculation.resolve(strict=True)
    topology_path = case_base / "system.prmtop"
    input_mdin = force_calculation / "dimer.mdin"
    candidate = case_base / "calculation/optimized.ncrst"
    mode_path = case_base / "calculation/path_tsmode.xyz"
    force_path = force_calculation / "path_force.xyz"
    topology = AmberParm(str(topology_path))
    active_mask = parse_option(input_mdin.read_text(), "active")
    selected = np.asarray(AmberMask(topology, active_mask).Selection(), dtype=bool)
    active_indices = np.flatnonzero(selected)
    if len(active_indices) != 4165:
        raise ValueError(
            f"Active mask selects {len(active_indices)} atoms; expected 4165"
        )
    mode = read_xyz_frames(mode_path, len(active_indices))[-1]
    if not np.isfinite(mode).all() or np.linalg.norm(mode) <= 1.0e-12:
        raise ValueError("Final path_tsmode frame is not a finite nonzero vector")
    force_rows = read_xyz_frames(force_path, len(active_indices))[-1]
    labels = [
        line.split()[0]
        for line in force_path.read_text().splitlines()[2:]
        if line.strip()
    ]
    expected_labels = [
        str(topology.atoms[index].name).strip()[0].upper() for index in active_indices
    ]
    if labels != expected_labels:
        raise ValueError(
            "path_force row labels/order differ from topology active-mask order"
        )
    reference = Rst7.open(str(candidate))
    xyz = np.asarray(reference.coordinates, dtype=float).reshape(-1, 3)
    if len(xyz) != len(topology.atoms):
        raise ValueError("Candidate restart and topology atom counts differ")
    output.mkdir(parents=True, exist_ok=False)
    restart_paths: dict[str, Path] = {}
    offsets = {"plus": center_offset + step, "minus": center_offset - step}
    for label, offset in offsets.items():
        restart_path = output / f"{label}.ncrst"
        displaced = displace_active(xyz, active_indices, mode, offset)
        write_restart(restart_path, displaced, reference)
        restart_paths[label] = restart_path
    manifest = {
        "step_angstrom": step,
        "center_offset_angstrom": center_offset,
        "offsets_from_ts_angstrom": offsets,
        "mode_norm": float(np.linalg.norm(mode)),
        "active_atom_count": int(len(active_indices)),
        "active_indices_zero_based": active_indices.tolist(),
        "path_force_rows_match_active_topology_order": True,
        "path_force_file_finite": bool(np.isfinite(force_rows).all()),
        "force_component_max_abs": float(np.max(np.abs(force_rows))),
        "sources": [
            _record(path)
            for path in (topology_path, candidate, mode_path, force_path, input_mdin)
        ],
        "prepared_restarts": {
            label: _record(path) for label, path in restart_paths.items()
        },
        "preparer": _record(Path(__file__)),
        "energy_evaluations": {"plus": None, "minus": None},
    }
    (output / "fd_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared ±%.3f A mode displacements in %s", step, output)
    return manifest


def read_final_energy(path: Path) -> float:
    """Read a final Hartree energy from native SANDER mdout."""
    pattern = (
        r"(?:Final converged energy:|Energy calculation finished, energy:)\s*"
        r"([-+\d.Ee]+)"
    )
    matches = re.findall(pattern, path.read_text())
    if not matches:
        raise ValueError(f"No final converged energy in {path}")
    energy = float(matches[-1])
    if not np.isfinite(energy):
        raise ValueError(f"Nonfinite final energy in {path}")
    return energy


def assess_fd(
    manifest_path: Path,
    plus_mdout: Path,
    minus_mdout: Path,
    center_energy_hartree: float,
    force_path: Path,
) -> dict[str, Any]:
    """Compare central finite differences with projected path-force values."""
    manifest = json.loads(manifest_path.read_text())
    step = float(manifest["step_angstrom"])
    plus, minus = read_final_energy(plus_mdout), read_final_energy(minus_mdout)
    derivative = (plus - minus) * HARTREE_TO_KCAL_MOL / (2.0 * step)
    center = center_energy_hartree
    force_rows = read_xyz_frames(force_path, manifest["active_atom_count"])[-1]
    case_base = Path(manifest["sources"][0]["path"]).parent
    topology = AmberParm(str(case_base / "system.prmtop"))
    active_mask = parse_option(
        Path(manifest["sources"][-1]["path"]).read_text(), "active"
    )
    active = np.asarray(AmberMask(topology, active_mask).Selection(), dtype=bool)
    indices = np.flatnonzero(active)
    mode_path = Path(manifest["sources"][2]["path"])
    mode = read_xyz_frames(mode_path, len(indices))[-1]
    direction = mode / np.linalg.norm(mode)
    projection = float(np.sum(force_rows * direction))
    candidates = {
        "raw_kcal_per_mol_angstrom": 1.0,
        "hartree_per_angstrom_to_kcal": HARTREE_TO_KCAL_MOL,
        "hartree_per_bohr_to_kcal_per_angstrom": HARTREE_TO_KCAL_MOL / BOHR_TO_ANGSTROM,
    }
    comparisons = {
        unit: {
            sign: projection * scale * sign,
            "absolute_error_kcal_per_mol_angstrom": abs(
                derivative - projection * scale * sign
            ),
        }
        for unit, scale in candidates.items()
        for sign in (1, -1)
    }
    result = {
        "plus_energy_hartree": plus,
        "minus_energy_hartree": minus,
        "center_energy_hartree": center,
        "plus_minus_energy_difference_hartree": plus - minus,
        "finite_difference_dE_ds_kcal_per_mol_angstrom": derivative,
        "path_force_projected_raw": projection,
        "path_force_projection_unit_hypotheses": comparisons,
        "energy_mean_offset_from_center_kcal_per_mol": ((plus + minus) / 2.0 - center)
        * HARTREE_TO_KCAL_MOL,
        "inputs": [
            _record(path)
            for path in (manifest_path, plus_mdout, minus_mdout, force_path)
        ],
    }
    result_path = manifest_path.parent / "finite_difference_result.json"
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result


def main() -> None:
    """Prepare displacements or assess completed finite-difference probes."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prepare = sub.add_parser("prepare")
    prepare.add_argument("--case-base", type=Path, required=True)
    prepare.add_argument("--force-calculation", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--step", type=float, default=0.1)
    prepare.add_argument("--center-offset", type=float, default=0.0)
    assess = sub.add_parser("assess")
    assess.add_argument("--manifest", type=Path, required=True)
    assess.add_argument("--plus-mdout", type=Path, required=True)
    assess.add_argument("--minus-mdout", type=Path, required=True)
    assess.add_argument("--center-energy-hartree", type=float, required=True)
    assess.add_argument("--force-path", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if args.action == "prepare":
        prepare_fd(
            args.case_base,
            args.force_calculation,
            args.output,
            args.step,
            args.center_offset,
        )
    else:
        print(
            json.dumps(
                assess_fd(
                    args.manifest,
                    args.plus_mdout,
                    args.minus_mdout,
                    args.center_energy_hartree,
                    args.force_path,
                ),
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
