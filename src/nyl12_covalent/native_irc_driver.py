"""Reproducible native-Sander numerical downhill paths from a dimer TS."""

from __future__ import annotations

import hashlib
import json
import logging
import re
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from nyl12_covalent import native_sander_cli
from nyl12_covalent.mass_weighted_descent import mass_weighted_descent_step
from nyl12_covalent.native_sander_cli import (
    FORCE_SOURCE_FILES,
    FORCE_SOURCE_REFERENCES,
    HARTREE_TO_KCAL_PER_MOL,
    ANGSTROM_PER_BOHR,
    PATH_FORCE_TO_GRADIENT_SCALE,
    NativeSanderEvaluator,
    load_calibration,
    prepare_native_input,
    read_distance_pairs,
)
from scripts.prepare_dimer_continuation import read_xyz_frames, write_restart

LOGGER = logging.getLogger(__name__)
FloatArray = np.ndarray


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_two_branches(
    *,
    topology_path: Path,
    restart_path: Path,
    mdin_path: Path,
    path_xyz: Path,
    mode_path: Path,
    sander: Path,
    output_dir: Path,
    calibration_path: Path,
    distance_pairs_path: Path,
    accepted_steps: int = 20,
    initial_step_size: float = 0.02,
    minimum_step_size: float = 1e-5,
    gradient_rms_stop: float | None = None,
    energy_plateau_window: int = 0,
    energy_plateau_tolerance_hartree: float = 1e-8,
) -> Path:
    """Run bounded positive/negative mass-weighted descent branches.

    This function is the execution entry point. It refuses to start until the
    native force convention is backed by an accepted finite-difference record.
    """
    if accepted_steps < 1:
        raise ValueError("accepted_steps must be positive")
    if gradient_rms_stop is not None and gradient_rms_stop <= 0.0:
        raise ValueError("gradient_rms_stop must be positive")
    if energy_plateau_window < 0 or energy_plateau_tolerance_hartree <= 0.0:
        raise ValueError("Invalid energy plateau criteria")
    scale, calibration = load_calibration(calibration_path)
    paths = [
        topology_path,
        restart_path,
        mdin_path,
        path_xyz,
        mode_path,
        sander,
        calibration_path,
        distance_pairs_path,
        Path(__file__),
        Path(__file__).with_name("mass_weighted_descent.py"),
        Path(native_sander_cli.__file__),
        Path(__file__).resolve().parents[2]
        / "scripts/run_native_mass_weighted_path.py",
        Path(__file__).resolve().parents[2] / "scripts/prepare_dimer_continuation.py",
        *FORCE_SOURCE_FILES,
    ]
    continuing = output_dir.exists()
    topology = AmberParm(str(topology_path))
    restart = Rst7.open(str(restart_path))
    xyz = np.asarray(restart.coordinates, dtype=float).reshape((-1, 3))
    if len(xyz) != len(topology.atoms):
        raise ValueError("Restart and topology atom counts differ")
    source_text = mdin_path.read_text(encoding="utf-8")
    active_mask_match = re.search(r"\bactive\s*=\s*'([^']+)'", source_text, re.I)
    if active_mask_match is None:
        raise ValueError("Source input lacks quoted active mask")
    active = np.asarray(
        AmberMask(topology, active_mask_match.group(1)).Selection(), dtype=bool
    )
    if int(active.sum()) != 4165:
        raise ValueError(f"Active mask selects {active.sum()} atoms; expected 4165")
    path_frames = read_xyz_frames(path_xyz, int(active.sum()))
    mode_frames = read_xyz_frames(mode_path, int(active.sum()))
    path_error = float(np.max(np.abs(path_frames[-1] - xyz[active])))
    if path_error > 0.0025:
        raise ValueError(
            f"Last path frame differs from TS restart by {path_error:.6g} A"
        )
    mode = mode_frames[-1]
    mode_q = np.zeros_like(xyz)
    masses = np.asarray([atom.mass for atom in topology.atoms], dtype=float)
    mode_q[active] = np.sqrt(masses[active, None]) * mode
    if not np.isfinite(mode_q).all() or np.linalg.norm(mode_q) <= 1e-12:
        raise ValueError("Final dimer mode is invalid")
    pairs = read_distance_pairs(distance_pairs_path, len(xyz))
    mdin_text = prepare_native_input(source_text)
    hashes = {str(path.resolve()): sha256(path) for path in paths}
    manifest = {
        "description": "Numerical mass-weighted downhill path pilot; not a formal IRC",
        "accepted_steps_per_branch": accepted_steps,
        "initial_step_size_mass_weighted": initial_step_size,
        "minimum_step_size_mass_weighted": minimum_step_size,
        "stop_criteria": {
            "maximum_accepted_steps_per_branch": accepted_steps,
            "gradient_rms_kcal_per_mol_angstrom_sqrt_amu": gradient_rms_stop,
            "energy_plateau_window": energy_plateau_window,
            "energy_plateau_tolerance_hartree": energy_plateau_tolerance_hartree,
        },
        "temperature_kelvin": 0.0,
        "dimer": False,
        "fixed_exterior": True,
        "active_mask_atom_count": int(active.sum()),
        "path_frame_count": len(path_frames),
        "mode_frame_count": len(mode_frames),
        "last_path_max_abs_error_angstrom": path_error,
        "force_calibration": calibration,
        "force_to_gradient_scale": scale,
        "force_conversion": {
            "formula": "gradient[kcal mol-1 A-1] = -path_force / ANGSTROM_PER_BOHR^2 * HARTREE_TO_KCAL_PER_MOL",
            "angstrom_per_bohr": ANGSTROM_PER_BOHR,
            "hartree_to_kcal_per_mol": HARTREE_TO_KCAL_PER_MOL,
            "path_force_to_gradient_scale": PATH_FORCE_TO_GRADIENT_SCALE,
            "source_references": FORCE_SOURCE_REFERENCES,
            "source_sha256": {
                str(path.resolve()): sha256(path) for path in FORCE_SOURCE_FILES
            },
        },
        "source_sha256": hashes,
        "prepared_mdin_sha256": hashlib.sha256(mdin_text.encode()).hexdigest(),
        "distance_pairs": pairs,
        "branches": ["mode_plus", "mode_minus"],
        "path_classification": "provisional numerical descent; not Hessian-based IRC",
    }
    manifest_path = output_dir / "manifest.json"
    if continuing:
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            existing.get("source_sha256") != hashes
            or existing.get("prepared_mdin_sha256") != manifest["prepared_mdin_sha256"]
        ):
            raise ValueError(
                "Existing run inputs differ; resume in a new output directory"
            )
        if (
            existing.get("initial_step_size_mass_weighted") != initial_step_size
            or existing.get("minimum_step_size_mass_weighted") != minimum_step_size
        ):
            raise ValueError("Resume step sizes must match the original run")
        old_stop = dict(existing.get("stop_criteria", {}))
        new_stop = dict(manifest["stop_criteria"])
        old_stop.pop("maximum_accepted_steps_per_branch", None)
        new_stop.pop("maximum_accepted_steps_per_branch", None)
        if old_stop != new_stop:
            raise ValueError("Resume stop criteria must match the original run")
        previous_target = int(existing["accepted_steps_per_branch"])
        if accepted_steps < previous_target:
            raise ValueError("Resume target cannot be lower than the recorded target")
        for branch in manifest["branches"]:
            branch_dir = output_dir / branch
            if not (branch_dir / "steps.jsonl").is_file():
                raise ValueError(f"Cannot resume incomplete branch record: {branch}")
            if not (branch_dir / "checkpoint.ncrst").is_file():
                raise ValueError(f"Cannot resume without branch checkpoint: {branch}")
        manifest["resumed_from_manifest_sha256"] = sha256(manifest_path)
        manifest["previous_accepted_steps_target"] = previous_target
        existing["accepted_steps_per_branch"] = accepted_steps
        existing["stop_criteria"] = manifest["stop_criteria"]
        existing["resumed_from_manifest_sha256"] = manifest[
            "resumed_from_manifest_sha256"
        ]
        existing["previous_accepted_steps_target"] = previous_target
        manifest = existing
    else:
        output_dir.mkdir(parents=True, exist_ok=False)
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    for branch, sign in (("mode_plus", 1.0), ("mode_minus", -1.0)):
        branch_dir = output_dir / branch
        if not continuing:
            branch_dir.mkdir()
        prior_records = []
        journal_path = branch_dir / "steps.jsonl"
        if continuing:
            prior_records = [
                json.loads(line)
                for line in journal_path.read_text().splitlines()
                if line.strip()
            ]
        accepted = sum(record["accepted"] is True for record in prior_records)
        checkpoint = branch_dir / "checkpoint.ncrst"
        if continuing and prior_records:
            recorded_checkpoint = prior_records[-1].get("checkpoint_sha256")
            if recorded_checkpoint != sha256(checkpoint):
                raise ValueError(
                    f"Checkpoint does not match last step record: {branch}"
                )
        if accepted >= accepted_steps:
            continue
        current = (
            np.asarray(
                Rst7.open(str(branch_dir / "checkpoint.ncrst")).coordinates, dtype=float
            ).reshape((-1, 3))
            if continuing
            else xyz.copy()
        )
        if not continuing:
            write_restart(checkpoint, current, restart)
        direction = sign * mode_q
        step_size = (
            prior_records[-1]["next_step_size_mass_weighted"]
            if prior_records
            else initial_step_size
        )
        attempts = len(prior_records)
        stop_reason: str | None = None
        journal = journal_path.open("a" if continuing else "x", encoding="utf-8")
        evaluator = NativeSanderEvaluator(
            topology_path.resolve(),
            restart,
            active,
            mdin_text,
            sander.resolve(),
            branch_dir,
            scale,
        )
        try:
            while accepted < accepted_steps:
                attempts += 1
                if attempts > accepted_steps * 12:
                    raise RuntimeError(f"{branch}: acceptance limit exhausted")
                evaluated: list[float] = []
                gradients: list[FloatArray] = []

                def evaluate(coords: FloatArray) -> tuple[float, FloatArray]:
                    energy, gradient = evaluator(coords)
                    evaluated.append(energy)
                    gradients.append(gradient.copy())
                    return energy, gradient

                result = mass_weighted_descent_step(
                    current,
                    masses,
                    evaluate,
                    step_size,
                    fixed_atoms=~active,
                    initial_direction_q=direction if accepted == 0 else None,
                    minimum_step_size=minimum_step_size,
                )
                if result.accepted:
                    current = result.coordinates.copy()
                    accepted += 1
                    direction = None
                    step_size = result.step_size
                    write_restart(checkpoint, current, restart)
                else:
                    step_size = result.step_size
                    if step_size <= minimum_step_size and not result.accepted:
                        raise RuntimeError(f"{branch}: minimum step rejected")
                record = {
                    "branch": branch,
                    "attempt": attempts,
                    "accepted_step": accepted if result.accepted else None,
                    "accepted": result.accepted,
                    "step_size_mass_weighted": result.step_size,
                    "next_step_size_mass_weighted": step_size,
                    "evaluation_energies_hartree": evaluated,
                    "evaluation_evidence": evaluator.evidence[-len(evaluated) :],
                    "energy_hartree": result.energy,
                    "gradient_rms_kcal_per_mol_angstrom_sqrt_amu": float(
                        np.sqrt(
                            np.mean(
                                (
                                    gradients[2 if result.accepted else 0][active]
                                    / np.sqrt(masses[active, None])
                                )
                                ** 2
                            )
                        )
                    ),
                    "checkpoint_sha256": sha256(checkpoint),
                    "distances_angstrom": {
                        pair["label"]: float(
                            np.linalg.norm(current[pair["i"]] - current[pair["j"]])
                        )
                        for pair in pairs
                    },
                    "coordinates_sha256": hashlib.sha256(current.tobytes()).hexdigest(),
                    "source_sha256": hashes,
                }
                journal.write(json.dumps(record, allow_nan=False) + "\n")
                journal.flush()
                accepted_records = [
                    item for item in (*prior_records, record) if item["accepted"]
                ]
                recent = (
                    accepted_records[-energy_plateau_window:]
                    if energy_plateau_window
                    else []
                )
                if (
                    gradient_rms_stop is not None
                    and result.accepted
                    and record["gradient_rms_kcal_per_mol_angstrom_sqrt_amu"]
                    <= gradient_rms_stop
                ):
                    stop_reason = "gradient_rms"
                if energy_plateau_window and len(recent) == energy_plateau_window:
                    delta = max(item["energy_hartree"] for item in recent) - min(
                        item["energy_hartree"] for item in recent
                    )
                    if delta <= energy_plateau_tolerance_hartree:
                        stop_reason = "energy_plateau"
                prior_records.append(record)
                if stop_reason:
                    break
            if stop_reason is None and accepted >= accepted_steps:
                stop_reason = "maximum_accepted_steps"
            manifest.setdefault("branch_stop_reasons", {})[branch] = (
                stop_reason or "interrupted"
            )
            manifest_path.write_text(
                json.dumps(manifest, indent=2, allow_nan=False) + "\n"
            )
            endpoint = {
                "branch": branch,
                "terminal_status": stop_reason or "interrupted",
                "accepted_steps": accepted,
                "coordinates_sha256": hashlib.sha256(current.tobytes()).hexdigest(),
                "checkpoint_sha256": sha256(branch_dir / "checkpoint.ncrst"),
                "distances_angstrom": {
                    pair["label"]: float(
                        np.linalg.norm(current[pair["i"]] - current[pair["j"]])
                    )
                    for pair in pairs
                },
                "interpretation": "endpoint candidate for later chemical acceptance; numerical path is not a formal IRC",
                "source_sha256": hashes,
            }
            (branch_dir / "endpoint_candidate.json").write_text(
                json.dumps(endpoint, indent=2, allow_nan=False) + "\n", encoding="utf-8"
            )
            LOGGER.info(
                "Stopped %s after %d accepted steps: %s", branch, accepted, stop_reason
            )
        finally:
            journal.close()
    if {str(path.resolve()): sha256(path) for path in paths} != hashes:
        raise RuntimeError("A source input changed during the path run")
    return output_dir
