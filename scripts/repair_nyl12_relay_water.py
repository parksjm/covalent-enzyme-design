"""Search and independently screen reproducible Nyl12 TI1 relay-water seeds."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import sys
from typing import Any

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.active_site_geometry import (  # noqa: E402
    coordinates_from_amber,
    measure_active_site,
)
from nyl12_covalent.inspection import write_inspection_pdb  # noqa: E402
from nyl12_covalent.relay_repair import search_repair  # noqa: E402
from nyl12_covalent.relay_repair_geometry import component_after_cut  # noqa: E402
from nyl12_covalent.ti1_checks import check_geometry  # noqa: E402

LOGGER = logging.getLogger(__name__)
TARGETS = {"pa6_l4": ("pa6_l2", "pa6"), "maxa_maxa": ("poly_xama_maxa", "poly")}


def write_json(path: Path, value: Any) -> None:
    """Create a JSON artifact without overwriting an existing file."""
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def check_invariance(model: Any, before: np.ndarray, after: np.ndarray) -> dict:
    """Require fixed backbone/substrate and preserved non-water bond lengths."""
    asn = {atom.name: atom for atom in model.residues[185].atoms}
    thr = {atom.name: atom for atom in model.residues[928].atoms}
    allowed = set(component_after_cut(asn["CB"], asn["CA"]))
    allowed.update(a.idx for a in thr["N"].bond_partners if a.atomic_number == 1)
    water = {a.idx for a in model.residues[1313].atoms}
    allowed.update(water)
    fixed = sorted(set(range(len(model.atoms))) - allowed)
    fixed_change = float(np.max(np.abs(after[fixed] - before[fixed])))
    bond_changes = []
    for bond in model.bonds:
        first, second = bond.atom1.idx, bond.atom2.idx
        if first in water or second in water:
            continue
        old = np.linalg.norm(before[first] - before[second])
        new = np.linalg.norm(after[first] - after[second])
        bond_changes.append(abs(float(new - old)))
    maximum_bond_change = max(bond_changes, default=0.0)
    changed = np.where(np.linalg.norm(after - before, axis=1) > 1e-8)[0]
    return {
        "accepted": fixed_change <= 1e-8 and maximum_bond_change <= 1e-8,
        "fixed_atom_maximum_coordinate_change_A": fixed_change,
        "nonwater_maximum_bond_length_change_A": maximum_bond_change,
        "allowed_atom_indices_zero_based": sorted(allowed),
        "changed_atom_indices_zero_based": changed.tolist(),
    }


def repair_files(
    topology: Path,
    coordinates: Path,
    target: str,
    output: Path,
    step_degrees: int,
) -> dict:
    """Create reports and export only a geometrically accepted starting seed.

    Args:
        topology: Product-gated TI1 Amber topology.
        coordinates: Matching original Amber coordinates.
        target: Canonical substrate identifier.
        output: New result directory, never an existing artifact directory.
        step_degrees: Explicit Asn and ammonium rotation grid increment.

    Returns:
        Scientific acceptance, independently of process completion.
    """
    legacy, kind = TARGETS[target]
    if output.exists():
        raise FileExistsError(output)
    topology, coordinates = (
        topology.resolve(strict=True),
        coordinates.resolve(strict=True),
    )
    model = AmberParm(str(topology), xyz=str(coordinates))
    original = np.asarray(model.coordinates, dtype=float).copy()
    initial = check_geometry(topology, coordinates, legacy)
    initial_metrics = measure_active_site(
        coordinates_from_amber(model, kind), kind, "ti1"
    )
    output.mkdir(parents=True)
    provenance = {
        "target": target,
        "archival_target_alias": legacy,
        "inputs": {
            label: {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for label, path in (("topology", topology), ("coordinates", coordinates))
        },
        "step_degrees": step_degrees,
        "purpose": "Geometric TI1 starting seed; not an optimized minimum or validated force field",
        "source_snapshot": str(ROOT),
    }
    write_json(output / f"{target}_ti1_provenance.json", provenance)
    write_json(output / "original_preflight.json", initial)
    write_json(output / "original_active_site.json", initial_metrics)
    report, candidate = search_repair(
        model, original, legacy, step_degrees=step_degrees
    )
    write_json(output / "search.json", report)
    acceptance = {"accepted": False, "candidate_found": candidate is not None}
    if candidate is not None:
        invariant = check_invariance(model, original, candidate)
        write_json(output / "invariance.json", invariant)
        if not invariant["accepted"]:
            write_json(
                output / "acceptance.json", acceptance | {"invariance": invariant}
            )
            raise ValueError("Repair altered fixed atoms or non-water bond lengths")
        model.coordinates = candidate
        restart = output / "candidate.rst7"
        model.save(str(restart), format="rst7", overwrite=False)
        preflight = check_geometry(topology, restart, legacy)
        write_json(output / f"{target}_ti1_audit.json", preflight)
        final_metrics = measure_active_site(
            coordinates_from_amber(model, kind), kind, "ti1"
        )
        write_json(output / "active_site.json", final_metrics)
        acceptance.update(accepted=bool(preflight["accepted"]), preflight=preflight)
        if preflight["accepted"]:
            mapping = write_inspection_pdb(
                model, candidate, output / f"{target}_ti1_starting.pdb"
            )
            write_json(output / f"{target}_residue_map.json", mapping)
    write_json(output / "acceptance.json", acceptance)
    LOGGER.info("%s geometric seed accepted: %s", target, acceptance["accepted"])
    return acceptance


def main() -> None:
    """Run the explicit-input geometric repair CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--target", required=True, choices=sorted(TARGETS))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--step-degrees", required=True, type=int)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    repair_files(
        args.topology, args.coordinates, args.target, args.output, args.step_degrees
    )


if __name__ == "__main__":
    main()
