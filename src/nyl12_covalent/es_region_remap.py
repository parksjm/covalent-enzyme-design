"""Remap an audited E-S QM region onto a coordinate-compatible topology."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from parmed.amber import AmberParm

from nyl12_covalent.es_qm_region import audit_es_qm_region
from nyl12_covalent.es_system_preparation import sha256


def _identity(atom: Any) -> dict[str, Any]:
    return {
        "index": atom.idx + 1,
        "residue": atom.residue.idx + 1,
        "resname": atom.residue.name,
        "name": atom.name,
        "element": int(atom.atomic_number),
        "mass": float(atom.mass),
        "type": atom.type,
        "charge": float(atom.charge),
    }


def remap_region(
    template: Path, topology: Path, coordinates: Path, output: Path
) -> dict[str, Any]:
    """Remap a reviewed atom-index template and independently audit the result.

    The operation permits coordinates and solvent count to change, but every
    selected atom identity and every crossing QM/MM bond must remain identical.
    """
    for path in (template, topology, coordinates):
        if not path.is_file():
            raise FileNotFoundError(path)
    payload = json.loads(template.read_text())
    parm = AmberParm(str(topology), xyz=str(coordinates))
    payload["topology"] = str(topology.resolve())
    payload["coordinates"] = str(coordinates.resolve())
    payload["remapped_from"] = {
        "path": str(template.resolve()),
        "sha256": sha256(template),
    }
    audits: list[dict[str, Any]] = []
    for variant in payload["variants"]:
        variant["topology_sha256"] = sha256(topology)
        variant["coordinates_sha256"] = sha256(coordinates)
        selected: set[int] = set()
        for group in variant["groups"]:
            for row in group["atoms"]:
                atom = parm.atoms[row["index"] - 1]
                current = _identity(atom)
                expected = (
                    row["residue"],
                    row["resname"],
                    row["name"],
                    row["element"],
                )
                observed = (
                    current["residue"],
                    current["resname"],
                    current["name"],
                    current["element"],
                )
                if observed != expected:
                    raise ValueError(f"QM identity changed at atom {row['index']}")
                row.clear()
                row.update(current)
                selected.add(current["index"])
        actual_boundaries = sorted(
            (number, partner.idx + 1)
            for number in selected
            for partner in parm.atoms[number - 1].bond_partners
            if partner.idx + 1 not in selected
        )
        expected_boundaries = sorted(
            (row["qm_index"], row["mm_index"]) for row in variant["boundaries"]
        )
        if actual_boundaries != expected_boundaries:
            raise ValueError("QM/MM boundary graph changed during remapping")
        audits.append(audit_es_qm_region(topology, coordinates, variant))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    audit_path = output.with_name(output.stem + "_audit.json")
    audit_path.write_text(json.dumps(audits, indent=2, allow_nan=False) + "\n")
    return {"config": str(output.resolve()), "audit": str(audit_path.resolve())}
