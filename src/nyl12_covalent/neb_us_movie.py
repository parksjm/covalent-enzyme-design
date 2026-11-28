"""Export accepted native two-CV umbrella samples for movie rendering."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm
from scipy.io import netcdf_file

from nyl12_covalent.path_inspection import dcd_counts

LOGGER = logging.getLogger(__name__)
_BULK_MASK = ":WAT,Na+,Cl-,NA,CL"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _sample_frames(frame_count: int, frames_per_block: int) -> list[int]:
    """Return evenly spaced, unique 1-based frame indices including the end."""
    if frame_count < 1 or frames_per_block < 1 or frames_per_block > frame_count:
        raise ValueError("Frame counts must satisfy 1 <= frames_per_block <= frames")
    frames = [
        math.ceil((i + 1) * frame_count / frames_per_block)
        for i in range(frames_per_block)
    ]
    if len(set(frames)) != frames_per_block:
        raise ValueError("Frame plan contains duplicate source frames")
    return frames


def _validate_block(
    block: Path,
    topology_hash: str,
    topology_atoms: int,
    sampled_frames_1based: list[int],
) -> tuple[Path, dict[str, Any]]:
    process = _load_json(block / "process.json")
    acceptance = _load_json(block / "acceptance.json")
    if type(process.get("exit_status")) is not int or process["exit_status"] != 0:
        raise ValueError(f"Production process did not succeed: {block}")
    if "every_integration_step_SCC_converged" in acceptance:
        history = acceptance.get("electronic_history")
        scc_accepted = (
            acceptance.get("accepted") is True
            and acceptance["every_integration_step_SCC_converged"] is True
            and isinstance(history, dict)
            and history.get("passed") is True
            and history.get("verified_scc_history") is True
            and history.get("unresolved_trials") == []
            and history.get("fatal_markers") == []
        )
    elif "SCC_error_detected" in acceptance:
        scc_accepted = acceptance["SCC_error_detected"] is False
    else:
        raise ValueError(f"Unknown production acceptance schema: {block}")
    temperature = acceptance.get("final_temperature_K")
    temperature_gate = acceptance.get("temperature_gate_K")
    if (
        acceptance.get("completed_steps") != 10000
        or acceptance.get("finite_energy_temperature") is not True
        or not scc_accepted
        or not isinstance(temperature, (int, float))
        or not isinstance(temperature_gate, list)
        or len(temperature_gate) != 2
        or not temperature_gate[0] <= temperature <= temperature_gate[1]
        or acceptance.get("frozen_outer_max_displacement_A") != 0.0
    ):
        raise ValueError(f"Production acceptance failed: {block}")
    trajectory = block / "trajectory.nc"
    observables = _load_json(block / "observables" / "observables.json")
    hashes = observables.get("source_sha256", {})
    if hashes.get(str(trajectory)) != _sha256(trajectory):
        raise ValueError(f"Trajectory hash mismatch: {trajectory}")
    matched_topologies = [
        digest
        for path, digest in hashes.items()
        if Path(path).name in {"system.prmtop", "topology.prmtop"}
    ]
    if topology_hash not in matched_topologies:
        raise ValueError(f"Observable topology hash mismatch: {block}")
    if observables.get("frames") != 100:
        raise ValueError(f"Expected 100 observed frames: {block}")
    with netcdf_file(trajectory, "r", mmap=True) as nc:
        coordinates = nc.variables.get("coordinates")
        time = nc.variables.get("time")
        if (
            coordinates is None
            or coordinates.shape != (100, topology_atoms, 3)
            or time is None
            or time.shape != (100,)
        ):
            raise ValueError(f"Expected 100 NetCDF frames: {trajectory}")
        for frame in sampled_frames_1based:
            xyz = np.asarray(coordinates[frame - 1, :, :]).copy()
            if (
                xyz.shape != (topology_atoms, 3)
                or not np.isfinite(xyz).all()
                or float(np.max(np.abs(xyz))) < 1.0
                or float(np.ptp(xyz, axis=0).max()) < 1.0
            ):
                raise ValueError(
                    f"Sampled trajectory coordinates are invalid in frame {frame}: "
                    f"{trajectory}"
                )
        time_values = time[:].copy().tolist()
        time_units = getattr(time, "units", b"")
        del coordinates, time
        if time_units not in {
            b"picosecond",
            b"picoseconds",
            "picosecond",
            "picoseconds",
        }:
            raise ValueError(
                f"Unsupported NetCDF time units {time_units!r}: {trajectory}"
            )
    controls = {
        str(path): _sha256(path)
        for path in (
            block / "process.json",
            block / "acceptance.json",
            block / "observables" / "observables.json",
        )
    }
    return trajectory, {
        "observables": observables,
        "roles_1based": observables.get("roles_1based"),
        "times": time_values,
        "time_units": "ps",
        "trajectory_sha256": hashes[str(trajectory)],
        "control_sha256": controls,
    }


def _cpptraj(script: Path, binary: Path, output_log: Path) -> None:
    environment = os.environ.copy()
    library_dir = str(binary.resolve().parent.parent / "lib")
    existing_library_path = environment.get("LD_LIBRARY_PATH")
    environment["LD_LIBRARY_PATH"] = (
        f"{library_dir}:{existing_library_path}"
        if existing_library_path
        else library_dir
    )
    result = subprocess.run(
        [str(binary.resolve()), "-i", str(script.resolve())],
        capture_output=True,
        text=True,
        check=False,
        cwd=script.parent,
        env=environment,
    )
    output_log.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"cpptraj failed ({result.returncode}); see {output_log}")


def export_part(
    campaign: Path,
    windows: list[str],
    cpptraj: Path,
    output: Path,
    frames_per_block: int = 5,
    production_blocks: int = 6,
    prepared_manifest: Path | None = None,
) -> dict[str, Any]:
    """Export selected production frames from a native two-CV NEB-US campaign.

    Args:
        campaign: Campaign root containing ``prepared.json``.
        windows: Requested window labels. Output follows prepared order.
        cpptraj: AmberTools cpptraj executable.
        output: New output directory; existing paths are rejected.
        frames_per_block: Number of evenly spaced frames sampled from each block.
        production_blocks: Number of completed production blocks selected from
            the start of each window's declared production stages.
        prepared_manifest: Explicit continuation preparation; defaults to the
            original campaign manifest. Native inputs are never rewritten.

    Returns:
        Manifest mapping source frames to the exported DCD frames.
    """
    campaign = campaign.resolve()
    prepared_path = (prepared_manifest or campaign / "prepared.json").resolve(
        strict=True
    )
    config = _load_json(prepared_path)
    if type(production_blocks) is not int or production_blocks < 1:
        raise ValueError("Production block count must be a positive integer")
    topology = (campaign / config["topology"]).resolve(strict=True)
    topology_hash = _sha256(topology)
    parm = AmberParm(str(topology))
    natom = len(parm.atoms)
    declared = config["windows"]
    by_label = {item["label"]: item for item in declared}
    unknown = set(windows) - set(by_label)
    if unknown:
        raise ValueError(f"Unknown campaign windows: {sorted(unknown)}")
    ordered = [item for item in declared if item["label"] in set(windows)]
    if not ordered or len(ordered) != len(set(windows)):
        raise ValueError("At least one unique requested window is required")
    if not cpptraj.is_file():
        raise FileNotFoundError(cpptraj)
    cpptraj = cpptraj.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    qm_ids: list[int] = []
    expected_roles: dict[str, int] | None = None
    block_records: list[dict[str, Any]] = []
    source_hashes: dict[str, str] = {
        str(topology): topology_hash,
        str(prepared_path): _sha256(prepared_path),
    }
    plan: list[dict[str, Any]] = []
    for window in ordered:
        available_stages = [
            stage for stage in window["stages"] if stage.get("kind") == "production"
        ]
        if len(available_stages) < production_blocks:
            raise ValueError(
                f"Not enough declared production blocks: {window['label']}"
            )
        stages = available_stages[:production_blocks]
        for stage in stages:
            block = campaign / window["label"] / stage["label"]
            frames = _sample_frames(100, frames_per_block)
            trajectory, metadata = _validate_block(block, topology_hash, natom, frames)
            observables = metadata["observables"]
            ids = observables.get("qm_atom_ids")
            if (
                not isinstance(ids, list)
                or not ids
                or any(not isinstance(i, int) for i in ids)
            ):
                raise ValueError(f"Invalid QM atom IDs: {block}")
            if qm_ids and qm_ids != ids:
                raise ValueError("QM atom IDs differ between production blocks")
            qm_ids = ids
            roles = metadata["roles_1based"]
            required_roles = {"carbon", "leaving_n", "og1", "hg1", "nalpha"}
            if (
                not isinstance(roles, dict)
                or required_roles - roles.keys()
                or any(
                    type(value) is not int or not 1 <= value <= natom
                    for value in roles.values()
                )
                or not set(roles.values()).issubset(qm_ids)
            ):
                raise ValueError(f"Invalid observable role mapping: {block}")
            if expected_roles is not None and expected_roles != roles:
                raise ValueError("Observable role mappings differ between blocks")
            expected_roles = roles
            trajectory_hash = metadata["trajectory_sha256"]
            source_hashes[str(trajectory)] = trajectory_hash
            source_hashes.update(metadata["control_sha256"])
            block_records.append({"trajectory": trajectory, "frames": frames})
            feature_names = window["bias"].get("feature_names")
            centers = window["bias"].get("centers_A")
            if (
                not isinstance(feature_names, list)
                or len(feature_names) != 2
                or len(set(feature_names)) != 2
                or not isinstance(centers, list)
                or len(centers) != 2
                or not set(feature_names).issubset({"q_acyl", "q_act", "q_del"})
                or "q_acyl" not in feature_names
                or not any(name in feature_names for name in ("q_act", "q_del"))
            ):
                raise ValueError(
                    f"Invalid two-CV relay bias metadata: {window['label']}"
                )
            bias_centers = dict(zip(feature_names, map(float, centers), strict=True))
            for frame in frames:
                time_values = metadata["times"]
                time_ps = float(time_values[frame - 1])
                plan.append(
                    {
                        "output_frame0based": len(plan),
                        "window": window["label"],
                        "chart": window["chart"],
                        "block": stage["label"],
                        "source": str(trajectory),
                        "source_frame1based": frame,
                        "window_time_ps": time_ps,
                        "time_basis": "NetCDF time (ps)",
                        "bias_centers_A": bias_centers,
                    }
                )

    if any(atom_id < 1 or atom_id > natom for atom_id in qm_ids):
        raise ValueError("QM atom ID outside topology")
    qm_mask = ",".join(map(str, sorted(set(qm_ids))))
    strip_mask = f"({_BULK_MASK})&!@{qm_mask}"
    kept = [
        atom.idx + 1
        for atom in parm.atoms
        if not (
            atom.residue.name in {"WAT", "Na+", "Cl-", "NA", "CL"}
            and atom.idx + 1 not in qm_ids
        )
    ]
    kept_set = set(kept)
    kept_names = [
        (atom.name, atom.residue.name)
        for atom in parm.atoms
        if atom.idx + 1 in kept_set
    ]
    if not set(qm_ids).issubset(kept):
        raise ValueError("QM atoms were removed by the bulk-solvent mask")

    parm_script = output / "write_topology.in"
    parm_script.write_text(
        f"parm {topology}\nparmstrip {strip_mask}\n"
        "parmwrite out part.parm7\nrun\nquit\n"
    )
    _cpptraj(parm_script, cpptraj, output / "write_topology.log")
    dcd_script = output / "write_trajectory.in"
    commands = [f"parm {topology}"]
    for record in block_records:
        frames = record["frames"]
        stride = frames[1] - frames[0] if len(frames) > 1 else 1
        if frames == list(range(frames[0], frames[-1] + 1, stride)):
            commands.append(
                f"trajin {record['trajectory']} {frames[0]} {frames[-1]} {stride}"
            )
        else:
            commands.extend(
                f"trajin {record['trajectory']} {frame} {frame} 1" for frame in frames
            )
    commands.extend([f"strip {strip_mask}", "trajout part.dcd dcd", "run", "quit"])
    dcd_script.write_text("\n".join(commands) + "\n")
    _cpptraj(dcd_script, cpptraj, output / "write_trajectory.log")

    out_topology = AmberParm(str(output / "part.parm7"))
    out_names = [(atom.name, atom.residue.name) for atom in out_topology.atoms]
    if out_names != kept_names:
        raise ValueError(
            "Exported topology atom order/names differ from selected source atoms"
        )
    dcd_frame_count, dcd_atom_count = dcd_counts(output / "part.dcd")
    if (dcd_frame_count, dcd_atom_count) != (len(plan), len(kept)):
        raise ValueError("Exported DCD dimensions differ from frame plan/topology")

    manifest = {
        "windows": [item["label"] for item in ordered],
        "window_charts": {item["label"]: item["chart"] for item in ordered},
        "roles_1based": expected_roles,
        "frames": len(plan),
        "atoms": len(kept),
        "original_atom_ids_1based": kept,
        "qm_atom_ids_1based": sorted(set(qm_ids)),
        "topology_sha256": topology_hash,
        "source_sha256": source_hashes,
        "frame_plan": plan,
        "part.dcd": "part.dcd",
        "part.parm7": "part.parm7",
        "part_sha256": {
            "part.dcd": _sha256(output / "part.dcd"),
            "part.parm7": _sha256(output / "part.parm7"),
        },
        "coordinates_aligned_or_rewrapped": False,
        "strip_mask": strip_mask,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info(
        "Exported %d frames from %d windows to %s", len(plan), len(ordered), output
    )
    return manifest
