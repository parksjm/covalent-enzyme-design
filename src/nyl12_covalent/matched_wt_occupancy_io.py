"""Bounded-memory read-only streaming of matched WT production trajectories."""

from __future__ import annotations

from collections import Counter
import csv
import hashlib
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from scipy.io import netcdf_file

from nyl12_covalent.catalytic_alignment import SITE_RESIDUES, _cif_atoms
from nyl12_covalent.leaving_group_relay_mapping import resolve_atom
from nyl12_covalent.leaving_group_relays import catalytic_geometry
from nyl12_covalent.matched_wt_occupancy import block_index, fractions, frame_flags
from nyl12_covalent.modular_es_snapshot_selection import classify_geometry
from nyl12_covalent.nalpha_water_relay import cell_matrix, native_geometry


def digest(path: Path) -> str:
    """Return a file SHA-256 using bounded reads.

    Args:
        path: Existing file.

    Returns:
        Hexadecimal SHA-256.
    """
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def resolve_geometry(topology: Any, descriptors: dict[str, dict]) -> dict[str, int]:
    """Resolve case-owned mapping and verify the retained scissile graph.

    Args:
        topology: Independent ParmEd topology.
        descriptors: Explicit one-based residue/name atom descriptors.

    Returns:
        Zero-based catalytic indices.
    """
    mapped = {key: resolve_atom(topology, desc) for key, desc in descriptors.items()}
    if len(set(mapped.values())) != len(mapped):
        raise ValueError("Geometry atoms must be distinct")
    carbon = topology.atoms[mapped["carbon"]]
    for key in ("oxygen", "nleave", "adjacent_c"):
        if topology.atoms[mapped[key]] not in carbon.bond_partners:
            raise ValueError(f"Scissile carbon missing declared bond to {key}")
    return mapped


def verify_microstate(topology: Any, state: str) -> None:
    """Validate all four declared Thr termini and Asp-pair proton inventories.

    Args:
        topology: Independently loaded matched-state topology.
        state: Canonical matched WT state name.
    """
    expected_nh = 3 if state.startswith("nh3_") else 2
    for thr, asp272, asp274 in (
        (929, 968, 970),
        (1025, 1064, 1066),
        (1121, 1160, 1162),
        (1217, 1256, 1258),
    ):
        nalpha = topology.atoms[
            resolve_atom(
                topology, {"residue_index1": thr, "resname": "THR", "atom": "N"}
            )
        ]
        if sum(atom.atomic_number == 1 for atom in nalpha.bond_partners) != expected_nh:
            raise ValueError("Thr terminal proton inventory differs from state")
        for residue_index, protonated in (
            (asp272, False),
            (asp274, state.endswith("ash274")),
        ):
            residue = topology.residues[residue_index - 1]
            if residue.name != ("ASH" if protonated else "ASP"):
                raise ValueError("Asp residue identity differs from state")
            protons = {
                a.name
                for a in residue.atoms
                if a.atomic_number == 1 and a.name.startswith("HD")
            }
            if protons != ({"HD2"} if protonated else set()):
                raise ValueError("Asp carboxyl proton inventory differs from state")


def analyze_case(case: dict, config: dict, output: Path) -> dict[str, Any]:
    """Stream every saved production frame into geometry/native/Tier evidence.

    Args:
        case: State-owned topology, trajectories and geometry mappings.
        config: Explicit cadence, time offset, reference CIF and native gates.
        output: Existing new analysis directory.

    Returns:
        Replica counts, 10 ns block occupancies and source hashes/stat guards.
    """
    paths = [
        Path(case["topology"]).resolve(strict=True),
        *[Path(p).resolve(strict=True) for p in case["trajectories"]],
        Path(config["reference_cif"]).resolve(strict=True),
    ]
    stats = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in paths}
    hashes = {str(p): digest(p) for p in paths}
    topology = parmed.load_file(str(paths[0]))
    verify_microstate(topology, case["state"])
    mapped = resolve_geometry(topology, case["geometry_atoms"])
    reference = _cif_atoms(paths[-1], "A", {auth for _, auth, _ in SITE_RESIDUES})
    counts: Counter[str] = Counter()
    blocks: dict[int, Counter[str]] = {}
    frame = 0
    previous_time = None
    csv_path = output / f"{case['name']}.csv"
    with csv_path.open("x", newline="") as handle:
        writer = None
        for trajectory in paths[1:-1]:
            with netcdf_file(str(trajectory), mmap=True) as source:
                shape = source.variables["coordinates"].shape
                if shape[1:] != (len(topology.atoms), 3):
                    raise ValueError("Trajectory atom inventory differs from topology")
                for local in range(shape[0]):
                    frame += 1
                    if frame > case["expected_frames"]:
                        raise ValueError("Too many production frames")
                    time = float(source.variables["time"][local])
                    production_time = frame * config["dt_ps"]
                    if not np.isclose(
                        time,
                        production_time + config["time_offset_ps"],
                        atol=0.05,
                        rtol=0,
                    ):
                        raise ValueError(
                            "Production time/cadence differs from manifest"
                        )
                    if previous_time is not None and time <= previous_time:
                        raise ValueError("Non-increasing trajectory time")
                    previous_time = time
                    xyz = np.array(source.variables["coordinates"][local], dtype=float)
                    if not np.isfinite(xyz).all():
                        raise ValueError("Nonfinite coordinates")
                    cell = cell_matrix(
                        np.array(source.variables["cell_lengths"][local]),
                        np.array(source.variables["cell_angles"][local]),
                    )
                    metrics = catalytic_geometry(xyz, mapped, cell)
                    native = native_geometry(
                        xyz, topology, reference, cell, mapped["nalpha"]
                    )
                    if not all(np.isfinite(v) for v in native.values()):
                        raise ValueError("Nonfinite native geometry")
                    flags = frame_flags(metrics, bool(native["native_pass"]))
                    row = {
                        "case": case["name"],
                        "frame": frame,
                        "time_ps": time,
                        "production_time_ps": production_time,
                        **metrics,
                        "tier": classify_geometry(metrics) or "none",
                        **native,
                    }
                    if writer is None:
                        writer = csv.DictWriter(handle, fieldnames=list(row))
                        writer.writeheader()
                    writer.writerow(row)
                    counts.update({k: int(v) for k, v in flags.items()})
                    block = blocks.setdefault(
                        block_index(frame, config["dt_ps"]), Counter()
                    )
                    block.update({k: int(v) for k, v in flags.items()})
                    block["frames"] += 1
    if frame != case["expected_frames"]:
        raise ValueError("Missing production frames")
    if any((p.stat().st_size, p.stat().st_mtime_ns) != stats[str(p)] for p in paths):
        raise ValueError("Source changed during analysis")
    return {
        "name": case["name"],
        "state": case["state"],
        "frames": frame,
        "counts": dict(counts),
        "fractions": fractions(counts, frame),
        "blocks_10ns": {
            str(k): {"frames": v["frames"], "fractions": fractions(v, v["frames"])}
            for k, v in sorted(blocks.items())
        },
        "geometry_mapping": mapped,
        "input_sha256": hashes,
        "input_size_mtime_ns": stats,
        "frame_csv": str(csv_path.resolve()),
    }
