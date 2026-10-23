#!/usr/bin/env python3
"""Export an immutable Gaussian optimization log as a reviewable trajectory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from datetime import datetime, timezone
import sys
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import GetDihedralDeg

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
from scripts.gaussian_log_to_xyz import Frame, parse_log, trajectory_text  # noqa: E402
from nyl12_covalent.gaussian_checks import validate_geometry  # noqa: E402

LOGGER = logging.getLogger(__name__)
EXPECTED_ATOMS = 39
EXPECTED_HYDROGENS = 23
ORIGIN_NAMES = {"THR:HG1": "H3"}


def _distance(frame: Frame, first: int, second: int) -> float:
    return float(
        np.linalg.norm(np.asarray(frame.coordinates[first]) - frame.coordinates[second])
    )


def _load_inputs(
    graph: Path, atom_map: Path, policy: Path
) -> tuple[Chem.Mol, dict[str, int], dict[str, Any]]:
    molecules = list(Chem.SDMolSupplier(str(graph), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Graph SDF must contain exactly one valid molecule")
    molecule = molecules[0]
    mapping_data = json.loads(atom_map.read_text(encoding="utf-8"))
    origins = mapping_data["atom_index_by_origin"]
    mapping_atoms = mapping_data["atoms"]
    policy_data = json.loads(policy.read_text(encoding="utf-8"))
    if molecule.GetNumAtoms() != EXPECTED_ATOMS or len(mapping_atoms) != EXPECTED_ATOMS:
        raise ValueError("Expected exactly 39 graph and atom-map entries")
    if (
        sum(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms())
        != EXPECTED_HYDROGENS
    ):
        raise ValueError("Expected exactly 23 explicit hydrogens")
    if any(origin.startswith("WAT:") for origin in origins):
        raise ValueError("Water-free atom map contains a WAT origin")
    if len(origins) != EXPECTED_ATOMS:
        raise ValueError("Atom map origin inventory is incomplete")
    for idx, atom in enumerate(mapping_atoms):
        if atom["index"] != idx or origins.get(atom["origin_id"]) != idx:
            raise ValueError(
                "Atom-map indexes are not contiguous or mutually consistent"
            )
        if atom["element"] != molecule.GetAtomWithIdx(idx).GetSymbol():
            raise ValueError(f"Atom-map element mismatch at atom {idx + 1}")
    return molecule, origins, policy_data


def _pdb_atom_name(origin: str) -> str:
    return ORIGIN_NAMES.get(origin, origin.rsplit(":", 1)[-1])


def _write_pdb(
    path: Path, frames: tuple[Frame, ...], molecule: Chem.Mol, origins: dict[str, int]
) -> None:
    """Write complete models with explicit input connectivity and formal charges."""
    origin_by_index = {index: origin for origin, index in origins.items()}
    connectivity = []
    for bond in molecule.GetBonds():
        first, second = bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1
        order = int(bond.GetBondTypeAsDouble())
        if order not in (1, 2, 3):
            raise ValueError("Unsupported nonintegral PDB bond order")
        connectivity.append(f"CONECT{first:5d}" + f"{second:5d}" * order + "\n")
        connectivity.append(f"CONECT{second:5d}" + f"{first:5d}" * order + "\n")
    lines = [
        "REMARK Fixed input graph in CONECT; inspect metrics for actual bond distances.\n"
    ]
    for model, frame in enumerate(frames, 1):
        lines.append(f"MODEL     {model:4d}\n")
        for idx, (atom, point) in enumerate(
            zip(molecule.GetAtoms(), frame.coordinates, strict=True), 1
        ):
            origin = origin_by_index[idx - 1]
            residue = origin.split(":", 1)[0]
            resname = {"THR": "THR", "NME": "NME", "PA6": "TIL"}[residue]
            chain, resid = {"THR": ("E", 1), "NME": ("E", 2), "PA6": ("I", 1)}[residue]
            name = _pdb_atom_name(origin)
            x, y, z = point
            charge = atom.GetFormalCharge()
            charge_code = (
                f"{abs(charge)}{'+' if charge > 0 else '-'}" if charge else "  "
            )
            lines.append(
                f"HETATM{idx:5d} {name:>4s} {resname:>3s} {chain}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {atom.GetSymbol():>2s}{charge_code}\n"
            )
        lines.extend(connectivity)
        lines.append("ENDMDL\n")
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def _write_viewer(output: Path, count: int) -> None:
    commands = [
        "delete waterfree_ti1_opt",
        "load optimization_trajectory.pdb, waterfree_ti1_opt, discrete=0",
        "hide everything, waterfree_ti1_opt",
        "show sticks, waterfree_ti1_opt",
        "set valence, 0",
        "util.cnc waterfree_ti1_opt",
        'label (waterfree_ti1_opt and resn THR and name OG1), "OG1"',
        'label (waterfree_ti1_opt and resn TIL and name C54), "C54"',
        'label (waterfree_ti1_opt and resn TIL and name O49), "O49"',
        'label (waterfree_ti1_opt and resn TIL and name N59), "N59"',
        "set movie_fps, 3",
        f"mset 1 -{count}",
        "zoom waterfree_ti1_opt",
    ]
    (output / "view_optimization.pml").write_text(
        "\n".join(commands) + "\n", encoding="utf-8"
    )


def export(
    log_path: Path,
    graph_path: Path,
    atom_map_path: Path,
    policy_path: Path,
    output: Path,
) -> dict[str, Any]:
    """Snapshot and export all complete input-orientation frames without overwriting.

    Args:
        log_path: Live Gaussian log, read once as bytes.
        graph_path: Explicit-hydrogen SDF defining the fixed input graph.
        atom_map_path: JSON mapping atom origins to Gaussian order.
        policy_path: Pilot policy with frozen atoms and prospective thresholds.
        output: New directory for the immutable snapshot and viewer files.

    Returns:
        Frame inventory and provenance report; this does not claim acceptance.
    """
    if output.exists():
        raise FileExistsError(f"Output directory already exists: {output}")
    source_bytes = log_path.read_bytes()
    captured_utc = datetime.now(timezone.utc).isoformat()
    source_hash = hashlib.sha256(source_bytes).hexdigest()
    text = source_bytes.decode("utf-8", errors="strict")
    molecule, origins, policy = _load_inputs(graph_path, atom_map_path, policy_path)
    parsed = parse_log(text, orientation="input", allow_incomplete=True)
    expected_elements = [atom.GetAtomicNum() for atom in molecule.GetAtoms()]
    for frame in parsed.frames:
        if list(frame.atomic_numbers) != expected_elements:
            raise ValueError(
                f"39-atom element order differs from input graph at log line {frame.line}"
            )
    output.mkdir(parents=True)
    (output / "optimization.log.txt").write_bytes(source_bytes)
    xyz_text = trajectory_text(parsed, log_path.name, source_hash)
    (output / "optimization_trajectory.xyz").write_text(xyz_text, encoding="utf-8")
    _write_pdb(output / "optimization_trajectory.pdb", parsed.frames, molecule, origins)

    xyz_seed_path = graph_path.with_name("native_seed.xyz")
    seed_lines = xyz_seed_path.read_text(encoding="utf-8").splitlines()
    seed = np.asarray(
        [[float(value) for value in line.split()[1:4]] for line in seed_lines[2:]]
    )
    if seed.shape != (EXPECTED_ATOMS, 3):
        raise ValueError("Seed XYZ does not contain 39 coordinates")
    idx = {
        name: origins[name]
        for name in ("THR:OG1", "PA6:C54", "PA6:O49", "PA6:N59", "PA6:C62")
    }
    frozen = policy["frozen_atom_indices_0based"]
    chi_atoms = [int(i) - 1 for i in policy["chi1_constraint_indices_1based"]]
    start_conf = Chem.Conformer(EXPECTED_ATOMS)
    for i, point in enumerate(seed):
        start_conf.SetAtomPosition(i, tuple(float(v) for v in point))
    start_chi = GetDihedralDeg(start_conf, *chi_atoms)
    hydrogen = [
        atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 1
    ]
    metrics_path = output / "trajectory_metrics.csv"
    columns = [
        "frame",
        "original_log_line",
        "OG1_C54_A",
        "C54_O49_A",
        "C54_N59_A",
        "C54_C62_A",
        "minimum_O49_H_A",
        "max_frozen_displacement_A",
        "chi1_error_deg",
        "graph_geometry_screen_pass",
    ]
    with metrics_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for frame_number, frame in enumerate(parsed.frames, 1):
            coordinates = np.asarray(frame.coordinates)
            conformer = Chem.Conformer(EXPECTED_ATOMS)
            for i, point in enumerate(coordinates):
                conformer.SetAtomPosition(i, tuple(float(v) for v in point))
            delta = GetDihedralDeg(conformer, *chi_atoms) - start_chi
            chi_error = abs((delta + 180.0) % 360.0 - 180.0)
            geometry = validate_geometry(
                molecule, list(frame.atomic_numbers), coordinates
            )
            writer.writerow(
                {
                    "frame": frame_number,
                    "original_log_line": frame.line,
                    "OG1_C54_A": _distance(frame, idx["THR:OG1"], idx["PA6:C54"]),
                    "C54_O49_A": _distance(frame, idx["PA6:C54"], idx["PA6:O49"]),
                    "C54_N59_A": _distance(frame, idx["PA6:C54"], idx["PA6:N59"]),
                    "C54_C62_A": _distance(frame, idx["PA6:C54"], idx["PA6:C62"]),
                    "minimum_O49_H_A": float(
                        np.min(
                            np.linalg.norm(
                                coordinates[hydrogen] - coordinates[idx["PA6:O49"]],
                                axis=1,
                            )
                        )
                    ),
                    "max_frozen_displacement_A": float(
                        np.max(
                            np.linalg.norm(coordinates[frozen] - seed[frozen], axis=1)
                        )
                    ),
                    "chi1_error_deg": chi_error,
                    "graph_geometry_screen_pass": geometry["accepted"],
                }
            )
    termination_lines = text.rstrip().splitlines()
    source_dir = output / "source"
    source_dir.mkdir()
    source_paths = [
        graph_path,
        atom_map_path,
        policy_path,
        xyz_seed_path,
        Path(__file__),
        ROOT / "scripts/gaussian_log_to_xyz.py",
        ROOT / "src/nyl12_covalent/gaussian_checks.py",
    ]
    source_hashes = {}
    for number, source in enumerate(source_paths):
        data = source.read_bytes()
        target = source_dir / f"{number}_{source.name}"
        target.write_bytes(data)
        source_hashes[str(source.resolve())] = hashlib.sha256(data).hexdigest()
    report: dict[str, Any] = {
        "captured_utc": captured_utc,
        "source_path": str(log_path.resolve()),
        "input_source_sha256": source_hashes,
        "coordinate_precision_note": "XYZ prints ten decimals; source Gaussian tables typically contain six decimals in angstrom. PDB is rounded to0.001A.",
        "frame_selection_note": "Every complete Input orientation table is retained, including trial/repeated snapshots; frame numbers are not asserted to be accepted optimizer steps or MD time.",
        "source_basename": log_path.name,
        "source_sha256": source_hash,
        "orientation": "Input",
        "orientation_counts": parsed.counts,
        "ignored_trailing_incomplete_counts": parsed.skipped_trailing,
        "frames": len(parsed.frames),
        "atoms_per_frame": EXPECTED_ATOMS,
        "explicit_H_per_frame": EXPECTED_HYDROGENS,
        "frame_source_lines": [frame.line for frame in parsed.frames],
        "optimization_completed": "Optimization completed." in text,
        "normal_termination": bool(
            termination_lines
            and termination_lines[-1]
            .strip()
            .startswith("Normal termination of Gaussian 16")
        ),
        "geometry_acceptance_claimed": False,
        "snapshot_note": "Frames are complete Input orientation tables from the snapshotted log. This export does not assert an optimized minimum, MD trajectory, or preserved bond connectivity.",
        "connectivity_note": "PDB CONECT records describe only the fixed input SDF graph and do not prove bonds remain intact in each frame.",
    }
    (output / "trajectory_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    _write_viewer(output, len(parsed.frames))
    (output / "README.md").write_text(
        "# Water-free TI1 optimization trajectory\n\n"
        f"Snapshot: {len(parsed.frames)} geometry frames, 39 atoms and 23 explicit H per frame.\n\n"
        "On your laptop, change PyMOL's working directory to the downloaded folder "
        "containing this README (normally nyl12-inspection-20261021/qm/waterfree_ti1_optimization). "
        "Then run:\n\n```pml\n@view_optimization.pml\nmplay\n```\n\n"
        "Use mstop or the frame slider to inspect individual states. VMD can load "
        "optimization_trajectory.xyz as a multi-frame XYZ file.\n\n"
        "The movie contains complete Input orientation tables from the immutable log snapshot, "
        "including any trial/repeated geometries. These are not MD times or guaranteed accepted "
        "optimizer steps. It is not an optimized minimum. All hydrogens are explicit; water is "
        "intentionally absent from this39-atom model.\n\n"
        "PDB CONECT records retain the initial chemical graph for display and do not prove "
        "bonds survive. Check trajectory_metrics.csv for the four reactive bond distances, "
        "minimum oxyanion-H distance, frozen-atom shifts, chi1 and graph screens.\n\n"
        "The exporter can be rerun with a new --output directory while Gaussian runs. "
        "Each export reads the log once, keeps that snapshot and never overwrites previous exports.\n",
        encoding="utf-8",
    )
    return report


def main() -> None:
    """Parse the required paths and export the trajectory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--graph", type=Path, required=True)
    parser.add_argument("--atom-map", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = export(args.log, args.graph, args.atom_map, args.policy, args.output)
    LOGGER.info(
        "Exported %d input-orientation frames (%s); no geometry acceptance claimed",
        report["frames"],
        report["source_sha256"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
