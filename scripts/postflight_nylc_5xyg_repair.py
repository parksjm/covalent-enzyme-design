#!/usr/bin/env python3
"""Read-only source mapping, geometry, and new-atom clash audit for 5XYG."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
SOURCE = ROOT / "analysis/negoro_2023_protonation/table1_job/panel/5XYG_deposited.pdb"
RUN = (
    ROOT
    / "analysis/negoro_2023_protonation/pdbfixer_repair_20261120/repair_run_20261121"
)
EXPECTED_RANGES = {"A": (1, 266), "B": (1, 266), "C": (267, 355), "D": (267, 355)}
LOGGER = logging.getLogger(__name__)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peptide_audit(
    atoms: dict[tuple[str, int, str, str], tuple[float, float, float]],
) -> dict[str, Any]:
    """Measure every expected adjacent peptide C-N distance."""
    distances = []
    missing = []
    for chain, (start, end) in EXPECTED_RANGES.items():
        for residue in range(start, end):
            c = (chain, residue, " ", "C")
            n = (chain, residue + 1, " ", "N")
            if c not in atoms or n not in atoms:
                missing.append([chain, residue, residue + 1])
                continue
            distances.append(
                {
                    "chain": chain,
                    "residue_pair": [residue, residue + 1],
                    "distance_angstrom": math.dist(atoms[c], atoms[n]),
                }
            )
    bad = [row for row in distances if not 1.1 <= row["distance_angstrom"] <= 1.8]
    return {
        "peptide_bond_count": len(distances),
        "missing_pairs": missing,
        "out_of_range_count": len(bad),
        "out_of_range": bad,
        "minimum_angstrom": min(
            (row["distance_angstrom"] for row in distances), default=None
        ),
        "maximum_angstrom": max(
            (row["distance_angstrom"] for row in distances), default=None
        ),
    }


def clash_audit(
    candidate: Path, modeled: set[tuple[str, int, str, str]]
) -> dict[str, Any]:
    """Report nonbonded heavy-atom distances below 1.5 A involving modeled atoms."""
    from openmm import unit
    from openmm.app import PDBFile

    pdb = PDBFile(str(candidate))
    atoms = list(pdb.topology.atoms())
    positions = pdb.positions.value_in_unit(unit.angstrom)
    keys = {
        atom.index: (
            atom.residue.chain.id,
            int(atom.residue.id),
            atom.residue.insertionCode or " ",
            atom.name,
        )
        for atom in atoms
    }
    coords = {keys[i]: tuple(map(float, positions[i])) for i in range(len(atoms))}
    heavy = {
        a.index
        for a in atoms
        if a.element is not None and a.element.symbol not in {"H", "D"}
    }
    bonds: dict[int, set[int]] = defaultdict(set)
    for a, b in pdb.topology.bonds():
        bonds[a.index].add(b.index)
        bonds[b.index].add(a.index)
    excluded = {}
    for i in heavy:
        near = set(bonds[i])
        for j in bonds[i]:
            near.update(bonds[j])
        excluded[i] = near
    cell_size = 1.5
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for i in heavy:
        cells[tuple(math.floor(float(x) / cell_size) for x in positions[i])].append(i)
    clashes = []
    for cell, indices in cells.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    other_cell = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    for i in indices:
                        for j in cells.get(other_cell, ()):
                            if j <= i or j in excluded[i]:
                                continue
                            k1, k2 = keys[i], keys[j]
                            if k1 not in modeled and k2 not in modeled:
                                continue
                            distance = math.dist(coords[k1], coords[k2])
                            if distance < 1.5:
                                clashes.append(
                                    {
                                        "atom1": k1,
                                        "atom2": k2,
                                        "distance_angstrom": distance,
                                    }
                                )
    clashes.sort(key=lambda row: row["distance_angstrom"])
    return {
        "heavy_atom_count": len(heavy),
        "newly_modeled_heavy_atom_count": len(modeled),
        "cutoff_angstrom": 1.5,
        "excluded_graph_neighbors": "1-2 and 1-3",
        "clash_count": len(clashes),
        "clashes": clashes,
    }


def run_postflight(
    candidate: Path, status_path: Path, report_path: Path, acceptance_path: Path | None
) -> dict[str, Any]:
    """Map original atoms and report candidate counts, peptide geometry, and clashes."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from repair_nylc_5xyg_pdbfixer import read_atoms, select_protein_atoms

    command = json.loads((RUN / "command.json").read_text())
    native_status = json.loads(status_path.read_text())
    source_digest = sha256(SOURCE)
    if source_digest != command["input_sha256"]:
        raise ValueError("Original structure hash differs from launch receipt")
    selected, choices = select_protein_atoms(SOURCE.read_text().splitlines())
    source_map = {}
    for chain, lines in selected.items():
        for line in lines:
            source_map[(chain, int(line[22:26]), line[26:27], line[12:16].strip())] = (
                line
            )
    staged = read_atoms(RUN / "result/5XYG_split_protein_input.pdb")
    output = read_atoms(candidate)
    if set(source_map) != set(staged):
        raise ValueError("Split input differs from selected original heavy atoms")
    missing = sorted(set(source_map) - set(output))
    mapping = []
    maximum_drift = 0.0
    for key, line in source_map.items():
        if key in output:
            xyz = tuple(float(line[start : start + 8]) for start in (30, 38, 46))
            drift = math.dist(xyz, output[key])
            maximum_drift = max(maximum_drift, drift)
            mapping.append(
                {
                    "source_auth_chain": "A" if key[0] in {"A", "C"} else "B",
                    "source_auth_residue": key[1],
                    "source_atom": key[3],
                    "output_chain": key[0],
                    "output_residue": key[1],
                    "output_atom": key[3],
                    "coordinate_displacement_angstrom": drift,
                }
            )
    mapping_path = report_path.with_name(report_path.stem + "_atom_mapping.json")
    mapping_path.write_text(
        json.dumps(
            {
                "source": str(SOURCE),
                "source_sha256": source_digest,
                "selected_altlocs": choices,
                "retained_original_heavy_atom_count": len(mapping),
                "missing_original_atoms": missing,
                "mapping": mapping,
            },
            indent=2,
        )
        + "\n"
    )
    residues: dict[str, set[int]] = defaultdict(set)
    for chain, residue, _, _ in output:
        residues[chain].add(residue)
    modeled = set(output) - set(source_map)
    acceptance = (
        json.loads(acceptance_path.read_text())
        if acceptance_path and acceptance_path.exists()
        else None
    )
    report = {
        "state": "accepted_postflight"
        if acceptance and native_status.get("exit_status") == 0
        else "candidate_diagnostic",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "candidate": str(candidate),
        "native_status": native_status,
        "source_sha256": source_digest,
        "mapped_original_heavy_atom_count": len(mapping),
        "missing_original_atoms": missing,
        "modeled_heavy_atom_count": len(modeled),
        "repaired_atom_count": len(output),
        "residues_per_chain": {key: len(values) for key, values in residues.items()},
        "maximum_original_heavy_displacement_angstrom": maximum_drift,
        "atom_mapping": str(mapping_path),
        "atom_mapping_sha256": sha256(mapping_path),
        "peptide_geometry": peptide_audit(output),
        "clash_audit": clash_audit(candidate, modeled),
        "acceptance_receipt": acceptance,
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info(
        "Postflight state %s; %d modeled-atom clashes reported",
        report["state"],
        report["clash_audit"]["clash_count"],
    )
    return report


def main() -> int:
    """Wait boundedly for status, then diagnose the specified candidate."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate-path",
        type=Path,
        default=RUN / "result/5XYG_repaired_protein_softmin.pdb",
    )
    parser.add_argument(
        "--status-path", type=Path, default=RUN / "softmin_continuation_status.json"
    )
    parser.add_argument(
        "--acceptance-path",
        type=Path,
        default=RUN / "result/softmin_scientific_acceptance.json",
    )
    parser.add_argument(
        "--report-path", type=Path, default=RUN / "softmin_postflight.json"
    )
    parser.add_argument("--wait-seconds", type=int, default=7200)
    parser.add_argument("--poll-seconds", type=int, default=30)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    deadline = time.monotonic() + args.wait_seconds
    while time.monotonic() <= deadline:
        if args.status_path.exists():
            if not args.candidate_path.exists():
                args.report_path.write_text(
                    json.dumps(
                        {
                            "state": "no_candidate_after_worker_exit",
                            "worker_status": json.loads(args.status_path.read_text()),
                        },
                        indent=2,
                    )
                    + "\n"
                )
                return 2
            run_postflight(
                args.candidate_path,
                args.status_path,
                args.report_path,
                args.acceptance_path,
            )
            return 0
        time.sleep(args.poll_seconds)
    args.report_path.write_text(
        json.dumps(
            {"state": "timeout_waiting_for_worker", "wait_seconds": args.wait_seconds},
            indent=2,
        )
        + "\n"
    )
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
