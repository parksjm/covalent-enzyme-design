#!/usr/bin/env python3
"""Extract and audit late Tier A/B WT snapshots for modular-PA6 QM/MM."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
import parmed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.catalytic_alignment import align_catalytic_site  # noqa: E402
from nyl12_covalent.modular_es_snapshot_selection import (  # noqa: E402
    classify_geometry,
    read_geometry_timeseries,
    select_modular_es_frames,
)

LOGGER = logging.getLogger(__name__)


def source_record(path: Path) -> dict[str, Any]:
    """Return a stable size, timestamp, and SHA-256 source record."""
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source changed while hashing: {path}")
    return {
        "path": str(path.resolve()),
        "bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256": digest.hexdigest(),
    }


def cpptraj_text(topology: Path, trajectory: Path, frame: int, output: Path) -> str:
    """Return exact cpptraj input for extracting one imaged solvated frame."""
    if frame < 1:
        raise ValueError("Trajectory frame must be positive")
    return (
        f"parm {topology.resolve()}\n"
        f"trajin {trajectory.resolve()} {frame} {frame} 1\n"
        "autoimage anchor :1313-1316\n"
        f"trajout {output.resolve()} netcdf\n"
        "run\n"
    )


def explicit_native_torsion_selections(
    geometry_csv: Path,
    torsion_csv: Path,
    frames: list[int],
) -> list[dict[str, Any]]:
    """Validate explicit frames against Tier A/B and native-torsion records."""
    if (
        not frames
        or len(frames) != len(set(frames))
        or any(frame < 1 for frame in frames)
    ):
        raise ValueError("Explicit frames must be unique positive integers")
    geometry = {
        int(row["frame"]): row for row in read_geometry_timeseries(geometry_csv)
    }
    torsions: dict[int, dict[str, str]] = {}
    with torsion_csv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            frame = int(float(row["frame"]))
            if frame in torsions:
                raise ValueError(f"Duplicate torsion frame {frame}")
            torsions[frame] = row
    selected = []
    for frame in frames:
        if frame not in geometry or frame not in torsions:
            raise ValueError(f"Explicit frame absent from source tables: {frame}")
        metrics = geometry[frame]
        tier = classify_geometry(metrics)
        torsion = torsions[frame]
        if tier not in {"A", "B"} or torsion["native_torsions"] != "True":
            raise ValueError(
                f"Explicit frame fails joint geometry/torsion gate: {frame}"
            )
        if torsion["tier"] != tier:
            raise ValueError(f"Tier disagreement at explicit frame {frame}")
        selected.append(
            {
                **metrics,
                "frame": frame,
                "tier": tier,
                "psi_delta": float(torsion["psi_delta"]),
                "chi1_delta": float(torsion["chi1_delta"]),
            }
        )
    return selected


def _extract_candidate(
    topology: Path,
    trajectory: Path,
    selection: dict[str, Any],
    target: Path,
    cpptraj: Path,
    reference_cif: Path,
) -> dict[str, Any]:
    target.mkdir()
    netcdf = target / "solvated.nc"
    cpptraj_input = target / "extract.cpptraj"
    cpptraj_input.write_text(
        cpptraj_text(topology, trajectory, selection["frame"], netcdf),
        encoding="utf-8",
    )
    environment = dict(os.environ)
    amberhome = cpptraj.parent.parent
    environment.update(AMBERHOME=str(amberhome), OMP_NUM_THREADS="1")
    library = str(amberhome / "lib")
    environment["LD_LIBRARY_PATH"] = (
        f"{library}:{environment['LD_LIBRARY_PATH']}"
        if environment.get("LD_LIBRARY_PATH")
        else library
    )
    with (target / "cpptraj.log").open("w", encoding="utf-8") as log:
        subprocess.run(
            [str(cpptraj), "-i", str(cpptraj_input)],
            check=True,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    structure = parmed.load_file(str(topology), xyz=str(netcdf))
    coordinates = np.asarray(structure.coordinates, dtype=float)
    if (
        coordinates.shape != (len(structure.atoms), 3)
        or not np.isfinite(coordinates).all()
    ):
        raise ValueError(f"Invalid extracted coordinates: {target}")
    restart = target / "solvated.rst7"
    structure.save(str(restart), overwrite=False)
    structure[":1-1316"].save(str(target / "protein_substrate.pdb"), overwrite=False)
    alignment = align_catalytic_site(
        reference_cif,
        topology,
        restart,
        restart,
        target / "native_alignment",
    )
    return {
        "selection": selection,
        "directory": str(target.resolve()),
        "coordinates": source_record(restart),
        "native_alignment": alignment,
    }


def prepare(config_path: Path, output: Path) -> dict[str, Any]:
    """Extract configured candidates and write a source-pinned audit manifest."""
    config = json.loads(config_path.read_text(encoding="utf-8"))
    topology = Path(config["topology"])
    cpptraj = Path(config["cpptraj"])
    reference_cif = Path(config["reference_cif"])
    source_paths = [config_path, Path(__file__), topology, cpptraj, reference_cif]
    replicas = config["replicas"]
    selection_mode = config.get("selection_mode")
    if selection_mode not in {"geometry_midpoint", "explicit_native_torsion"}:
        raise ValueError("Explicit supported selection_mode is required")
    for replica in replicas:
        source_paths.extend(
            [
                Path(replica["trajectory"]),
                Path(replica["geometry_csv"]),
                Path(replica["summary_json"]),
            ]
        )
        if selection_mode == "explicit_native_torsion":
            source_paths.append(Path(replica["torsion_csv"]))
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing configured inputs: {missing}")
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)

    topology_record = source_record(topology)
    windows = [tuple(window) for window in config["windows_ns"]]
    candidates = []
    for replica in replicas:
        summary = json.loads(Path(replica["summary_json"]).read_text(encoding="utf-8"))
        if summary["topology_sha256"] != topology_record["sha256"]:
            raise ValueError(f"Topology hash mismatch for {replica['name']}")
        if selection_mode == "geometry_midpoint":
            selections = select_modular_es_frames(
                Path(replica["geometry_csv"]), windows
            )
        else:
            selections = explicit_native_torsion_selections(
                Path(replica["geometry_csv"]),
                Path(replica["torsion_csv"]),
                replica["frames"],
            )
        for window_number, selection in enumerate(selections, start=1):
            label = (
                f"{replica['name']}_window_{window_number:02d}_"
                f"frame_{selection['frame']:04d}_{selection['time_ns']:.2f}ns"
            )
            LOGGER.info("Extracting %s", label)
            record = _extract_candidate(
                topology,
                Path(replica["trajectory"]),
                selection,
                output / label,
                cpptraj,
                reference_cif,
            )
            record["replica"] = replica["name"]
            candidates.append(record)
    manifest = {
        "schema_version": 1,
        "purpose": "late WT modular-PA6 Tier A/B candidates for QM/MM screening",
        "selection_mode": selection_mode,
        "selection_rule": config["selection_rule"],
        "sources": [source_record(path) for path in source_paths],
        "topology": topology_record,
        "candidates": candidates,
        "scientific_acceptance": {
            "snapshot_extraction_complete": True,
            "qmmm_hamiltonian_validated": False,
            "catalytic_es_after_qmmm_validated": False,
        },
    }
    (output / "snapshot_selection_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    """Parse command-line arguments and prepare the snapshot ensemble."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(message)s",
    )
    prepare(args.config.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
