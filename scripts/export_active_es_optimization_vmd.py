#!/usr/bin/env python3
"""Export the active unrestrained E-S optimization for inspection in VMD."""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm, NetCDFTraj

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.vmd_export_io import sha, validate_dcd  # noqa: E402

LOGGER = logging.getLogger(__name__)
EXPECTED_TOPOLOGY_SHA256 = (
    "ef49ad15041019b439fc20f51562192724a341b37bf3089aa18dd5a31852eb55"
)


def _snapshot(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    shutil.copy2(source, destination)


def _read_restart(topology: Path, coordinates: Path) -> tuple[np.ndarray, np.ndarray]:
    parm = AmberParm(str(topology), xyz=str(coordinates))
    xyz = np.asarray(parm.coordinates)
    box = np.asarray(parm.box)
    if xyz.shape != (103251, 3) or box.shape != (6,):
        raise ValueError(f"Unexpected restart dimensions: {coordinates}")
    return xyz, box


def _read_trajectory(
    path: Path, fixed_box: np.ndarray
) -> list[tuple[np.ndarray, np.ndarray]]:
    trajectory = NetCDFTraj.open_old(str(path))
    coordinates = np.asarray(trajectory.coordinates)
    boxes = np.asarray(trajectory.box)
    trajectory.close()
    if coordinates.ndim != 3 or coordinates.shape[1:] != (103251, 3):
        raise ValueError(f"Unexpected trajectory dimensions: {coordinates.shape}")
    if boxes.shape != (coordinates.shape[0], 6):
        raise ValueError(f"Unexpected trajectory box dimensions: {boxes.shape}")
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for xyz, box in zip(coordinates, boxes, strict=True):
        if not np.isfinite(xyz).all() or not np.isfinite(box).all():
            raise ValueError("Nonfinite optimizer coordinates or box")
        if np.all(box == 0):
            box = fixed_box.copy()
        elif not np.allclose(box, fixed_box, atol=1e-10, rtol=0):
            raise ValueError("Optimizer box differs from the fixed reference box")
        states.append((xyz, box))
    return states


def _write_viewer(output: Path, frame_count: int) -> None:
    common_source = ROOT / "output_vmd_inspection_20261024/nyl12_es_vmd_common.tcl"
    common = common_source.read_text(encoding="utf-8")
    expected = (
        "proc optimization {} {nyl12_load stripped.parm7 optimization.dcd "
        "48 $::nyl12_qm_indices $::nyl12_stripped_atom_count}"
    )
    replacement = (
        "proc optimization {} {nyl12_load stripped.parm7 optimization.dcd "
        f"{frame_count} $::nyl12_qm_indices $::nyl12_stripped_atom_count}}"
    )
    if common.count(expected) != 1:
        raise ValueError("Canonical VMD optimization procedure changed")
    (output / "nyl12_es_vmd_common.tcl").write_text(
        common.replace(expected, replacement), encoding="utf-8"
    )
    (output / "view_optimization.vmd").write_text(
        "# Stripped active unrestrained E-S optimization; QM atoms use Licorice.\n"
        "set ::nyl12_dir [file dirname [file normalize [info script]]]\n"
        "foreach name {qm_indices.tcl nyl12_es_vmd_common.tcl} {\n"
        '    if {![file isfile [file join $::nyl12_dir $name]]} {error "Missing viewer input: $name"}\n'
        "    source [file join $::nyl12_dir $name]\n"
        "}\n"
        "optimization\n"
        f'puts "{frame_count} optimizer states: initial, saved cycles, latest restart."\n',
        encoding="utf-8",
    )


def _gzip_copy(source: Path, destination: Path) -> None:
    with source.open("rb") as incoming, destination.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as outgoing:
            shutil.copyfileobj(incoming, outgoing)


def _write_manifest(
    output: Path,
    source_paths: dict[str, Path],
    states: list[tuple[np.ndarray, np.ndarray]],
    kept: np.ndarray,
    validation: dict[str, Any],
) -> None:
    manifest = {
        "display_frames": len(states),
        "frame_order": [
            "starting.ncrst",
            "all complete frames in the optimization.nc snapshot",
            "latest optimized.ncrst snapshot",
        ],
        "optimizer_states_not_md_times": True,
        "sources": {
            name: {
                "path": path.as_posix(),
                "sha256_at_export": sha(path),
                "bytes_at_export": path.stat().st_size,
            }
            for name, path in source_paths.items()
        },
        "stripped_atom_count": len(kept),
        "stripped_to_original_atom_indices_0based": kept.tolist(),
        "retained_water_topology_residues_1based": [3328, 10325, 20216],
        "selection": (
            "All protein and LIG residues plus the fixed three-water union selected "
            "within 4 A of a QM heavy atom in the validated E-S trajectory export"
        ),
        "validation": validation,
        "virtual_link_caps_inserted": 0,
    }
    (output / "frame_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def export(run_dir: Path, output: Path, cpptraj: Path) -> None:
    """Create and independently validate a stripped VMD trajectory.

    Args:
        run_dir: Active unrestrained minimization directory.
        output: New VMD output directory.
        cpptraj: AmberTools cpptraj executable.

    Raises:
        FileExistsError: If the destination already exists.
        ValueError: If source identities or exported coordinates fail validation.
    """
    if output.exists():
        raise FileExistsError(output)
    topology = run_dir.parent / "system/md_solvated.prmtop"
    starting = run_dir / "starting.ncrst"
    active_trajectory = run_dir / "optimization.nc"
    active_endpoint = run_dir / "optimized.ncrst"
    if sha(topology) != EXPECTED_TOPOLOGY_SHA256:
        raise ValueError("Active topology does not match the validated VMD topology")
    canonical = ROOT / "output_vmd_inspection_20261024"
    mapping = json.loads((canonical / "frame_manifest.json").read_text())
    kept = np.asarray(mapping["stripped_to_original_atom_indices_0based"], dtype=int)
    if len(kept) != 19976:
        raise ValueError("Unexpected canonical stripped atom count")

    output.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix="nyl12_es_vmd_", dir="/tmp") as temp_name:
        temporary = Path(temp_name)
        trajectory_snapshot = temporary / "optimization_snapshot.nc"
        endpoint_snapshot = temporary / "optimized_snapshot.ncrst"
        _snapshot(active_trajectory, trajectory_snapshot)
        _snapshot(active_endpoint, endpoint_snapshot)
        initial_xyz, fixed_box = _read_restart(topology, starting)
        saved_states = _read_trajectory(trajectory_snapshot, fixed_box)
        endpoint_xyz, endpoint_box = _read_restart(topology, endpoint_snapshot)
        if not np.allclose(endpoint_box, fixed_box, atol=1e-10, rtol=0):
            raise ValueError("Latest restart box differs from the fixed reference box")
        states = (
            [(initial_xyz, fixed_box)] + saved_states + [(endpoint_xyz, endpoint_box)]
        )
        cpptraj_input = output / "export.cpptraj"
        box_labels = ("x", "y", "z", "alpha", "beta", "gamma")
        box_command = "box " + " ".join(
            f"{label} {value:.16g}"
            for label, value in zip(box_labels, fixed_box, strict=True)
        )
        cpptraj_input.write_text(
            "\n".join(
                [
                    f"parm {topology}",
                    f"trajin {starting}",
                    f"trajin {trajectory_snapshot}",
                    f"trajin {endpoint_snapshot}",
                    "strip !:1-1313,3328,10325,20216",
                    box_command,
                    f"trajout {output / 'optimization.dcd'} dcd namdcell",
                    "run",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        environment = os.environ.copy()
        amber_lib = cpptraj.parent.parent / "lib"
        environment["LD_LIBRARY_PATH"] = str(amber_lib)
        with (output / "cpptraj.log").open("w", encoding="utf-8") as log:
            subprocess.run(
                [str(cpptraj), "-i", str(cpptraj_input)],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                env=environment,
            )
        selected_states = [(xyz[kept], box) for xyz, box in states]
        validation = validate_dcd(output / "optimization.dcd", selected_states)
        source_paths = {
            "topology": topology,
            "starting_coordinates": starting,
            "optimization_snapshot": trajectory_snapshot,
            "latest_restart_snapshot": endpoint_snapshot,
        }
        _write_manifest(output, source_paths, states, kept, validation)

    _gzip_copy(canonical / "stripped.parm7", output / "stripped.parm7.gz")
    shutil.copy2(canonical / "qm_indices.tcl", output / "qm_indices.tcl")
    _write_viewer(output, len(states))
    (output / "README.md").write_text(
        "# Active unrestrained E-S optimization\n\n"
        "In VMD, run: `source /absolute/path/view_optimization.vmd`.\n\n"
        "The protein NewCartoon representation is loaded but hidden. The 98 physical "
        "QM atoms and full substrate use Licorice; non-QM substrate carbons are orange. "
        "Three complete nearby TIP3P waters are retained; bulk water and ions are removed. "
        "The final frame is the latest restart snapshot and may not be converged.\n",
        encoding="utf-8",
    )
    payload = [
        "README.md",
        "cpptraj.log",
        "export.cpptraj",
        "frame_manifest.json",
        "nyl12_es_vmd_common.tcl",
        "optimization.dcd",
        "qm_indices.tcl",
        "stripped.parm7.gz",
        "view_optimization.vmd",
    ]
    (output / "MANIFEST.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in payload),
        encoding="utf-8",
    )
    LOGGER.info("Exported and verified %d frames to %s", len(states), output)


def main() -> int:
    """Parse arguments and export the active optimization.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpptraj", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    export(args.run_dir.resolve(), args.output.resolve(), args.cpptraj.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
