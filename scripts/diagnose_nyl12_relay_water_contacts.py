"""Report threshold-violating contacts for saved relay-water candidates."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.relay_water import _contact_limits  # noqa: E402


def _atom_label(atom: object) -> dict[str, object]:
    residue = atom.residue
    topology_residue = residue.idx + 1
    if topology_residue <= 928:
        chain_offset = (topology_residue - 1) // 232
        chain = "ABCD"[chain_offset]
        canonical_residue = (topology_residue - 1) % 232 + 1
        chain_residue = canonical_residue
    elif topology_residue <= 1312:
        chain_offset = (topology_residue - 929) // 96
        chain = "EFGH"[chain_offset]
        chain_residue = (topology_residue - 929) % 96 + 1
        canonical_residue = chain_residue + 232
    elif topology_residue == 1313:
        chain, chain_residue, canonical_residue = "I", 1, 1
    else:
        chain, chain_residue, canonical_residue = "J", 1, 1
    parents = [
        {"name": parent.name, "element": parent.atomic_number}
        for parent in atom.bond_partners
    ] if atom.atomic_number == 1 else []
    return {
        "atom_index_zero_based": atom.idx,
        "atom_name": atom.name,
        "element": atom.atomic_number,
        "residue_name": residue.name,
        "residue_index_one_based": residue.idx + 1,
        "parmed_residue_number_raw": residue.number,
        "parmed_chain_raw": getattr(residue, "chain", ""),
        "topology_residue_one_based": topology_residue,
        "canonical_chain": chain,
        "chain_residue_number": chain_residue,
        "canonical_residue_number": canonical_residue,
        "hydrogen_parents": parents,
        "is_n_h": atom.atomic_number == 1 and any(p["element"] == 7 for p in parents),
    }


def diagnose(
    topology: AmberParm,
    candidate_data: dict[str, object],
    provenance: dict[str, str],
) -> dict[str, object]:
    """Calculate saved candidate contacts using the unchanged geometric gates.

    Args:
        topology: Fully protonated Amber topology with coordinates loaded.
        candidate_data: Parsed water search JSON containing candidate coordinates.
        provenance: Mapping of absolute input and code paths to SHA-256 hashes.

    Returns:
        JSON-serializable contact report with canonical residue labels.
    """
    coords = np.asarray(topology.coordinates, dtype=float)
    thr = {atom.name: atom for atom in topology.residues[928].atoms}
    midpoint = (coords[thr["N"].idx] + coords[thr["OG1"].idx]) / 2
    selected = [
        atom.idx for atom in topology.atoms
        if atom.residue.idx != 1313 and np.linalg.norm(coords[atom.idx] - midpoint) < 9.0
    ]
    oxygen_limits, hydrogen_limits = _contact_limits(topology, selected)
    limits = {index: (oxygen_limits[pos], hydrogen_limits[pos]) for pos, index in enumerate(selected)}
    water_indices = candidate_data["water_atom_indices_zero_based"]
    reports = []
    for rank, candidate in enumerate(candidate_data["candidates"], 1):
        water_xyz = np.asarray(candidate["water_xyz_A"], dtype=float)
        contacts = []
        for water_name, point in zip(("O", "H1", "H2"), water_xyz, strict=True):
            for index in selected:
                atom = topology.atoms[index]
                distance = float(np.linalg.norm(coords[index] - point))
                limit = float(limits[index][0 if water_name == "O" else 1])
                if distance < limit:
                    label = _atom_label(atom)
                    contacts.append({
                        "water_atom": water_name,
                        **label,
                        "distance_A": distance,
                        "threshold_A": limit,
                        "shortfall_A": limit - distance,
                    })
        contacts.sort(key=lambda row: row["shortfall_A"], reverse=True)
        reports.append({
            "rank": rank,
            "search_accepted": candidate["accepted"],
            "search_max_shortfall_A": candidate["maximum_contact_shortfall_A"],
            "water_xyz_A": candidate["water_xyz_A"],
            "offending_contact_count": len(contacts),
            "max_shortfall_A": max((row["shortfall_A"] for row in contacts), default=0.0),
            "n_h_contacts": sum(bool(row["is_n_h"]) for row in contacts),
            "contacts": contacts,
        })
    return {
        "source_method": candidate_data["method"],
        "thresholds_source": "nyl12_covalent.relay_water._contact_limits; unchanged",
        "water_atom_indices_zero_based": water_indices,
        "input_sha256": provenance,
        "local_environment_atom_count": len(selected),
        "contact_rule": "same 9.0 A midpoint environment and per-element O/H limits used by search_relay_water; report strict distance < threshold",
        "candidates": reports,
    }


def _verify_search_inputs(
    search: dict[str, object], topology_path: Path, coordinates_path: Path
) -> None:
    """Ensure candidate coordinates were generated from these exact inputs."""
    hashes = search.get("input_sha256")
    if not isinstance(hashes, dict):
        raise ValueError("search JSON has no input_sha256 provenance")
    for path in (topology_path, coordinates_path):
        key = str(path.resolve())
        expected = hashes.get(key)
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if expected != actual:
            raise ValueError(f"search input hash mismatch for {key}")


def _write_new_report(path: Path, report: dict[str, object]) -> None:
    """Write a report without replacing a prior diagnostic artifact."""
    payload = json.dumps(report, indent=2) + "\n"
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        stream.write(payload)


def main() -> None:
    """Load an Amber topology and write a JSON contact report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--search-json", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    _verify_search_inputs(
        json.loads(args.search_json.read_text(encoding="utf-8")),
        args.topology,
        args.coordinates,
    )
    search = json.loads(args.search_json.read_text(encoding="utf-8"))
    topology = AmberParm(str(args.topology), xyz=str(args.coordinates))
    relay_module = ROOT / "src/nyl12_covalent/relay_water.py"
    report = diagnose(
        topology,
        search,
        {
            str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (args.topology, args.coordinates, args.search_json, Path(__file__), relay_module)
        },
    )
    _write_new_report(args.output, report)


if __name__ == "__main__":
    main()
