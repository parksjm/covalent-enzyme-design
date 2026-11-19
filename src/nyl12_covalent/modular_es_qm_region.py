"""Build the reviewed E-S QM chemistry for modular PA6 topologies."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

from parmed.amber import AmberParm

from nyl12_covalent.es_qm_region import audit_es_qm_region
from nyl12_covalent.es_system_preparation import sha256

LIGAND_KEYS = (
    (1314, "C4"),
    (1314, "H41"),
    (1314, "H42"),
    (1314, "C5"),
    (1314, "H51"),
    (1314, "H52"),
    (1314, "C"),
    (1314, "O"),
    (1315, "N"),
    (1315, "HN"),
    (1315, "C1"),
    (1315, "H11"),
    (1315, "H12"),
    (1315, "C2"),
    (1315, "H21"),
    (1315, "H22"),
)


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


def _atom(parm: Any, residue: int, name: str) -> Any:
    matches = [atom for atom in parm.residues[residue - 1].atoms if atom.name == name]
    if len(matches) != 1:
        raise ValueError(f"Expected one atom at residue {residue} named {name}")
    return matches[0]


def build_modular_es_region(
    template: Path,
    topology: Path,
    coordinates: Path,
    output: Path,
) -> dict[str, Any]:
    """Build and audit the 98-real-atom modular PA6 E-S QM region.

    Protein groups are preserved atom-for-atom from the reviewed whole-ligand
    template. The reacting modular amide is the carbonyl of P6I residue 1314,
    the amide nitrogen of P6I residue 1315, and two methylenes on either side.
    """
    for path in (template, topology, coordinates):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output.exists():
        raise FileExistsError(output)
    source = json.loads(template.read_text(encoding="utf-8"))
    variants = [
        row for row in source["variants"] if row["name"] == "primary_106_nuclei"
    ]
    if len(variants) != 1:
        raise ValueError("Exactly one reviewed primary QM template is required")
    variant = copy.deepcopy(variants[0])
    parm = AmberParm(str(topology), xyz=str(coordinates))

    ligand_groups = [
        group for group in variant["groups"] if group["name"] == "PA6_reacting_amide"
    ]
    if len(ligand_groups) != 1:
        raise ValueError("Reviewed template must contain one PA6 reacting-amide group")
    for group in variant["groups"]:
        if group is ligand_groups[0]:
            group["atoms"] = [_identity(_atom(parm, *key)) for key in LIGAND_KEYS]
            continue
        for row in group["atoms"]:
            atom = parm.atoms[row["index"] - 1]
            observed = (
                atom.residue.idx + 1,
                atom.residue.name,
                atom.name,
                atom.atomic_number,
            )
            expected = (row["residue"], row["resname"], row["name"], row["element"])
            if observed != expected:
                raise ValueError(f"Reviewed protein QM identity changed: {expected}")
            row.clear()
            row.update(_identity(atom))

    ligand_roles = {
        "scissile_carbon": _atom(parm, 1314, "C").idx + 1,
        "scissile_oxygen": _atom(parm, 1314, "O").idx + 1,
        "leaving_nitrogen": _atom(parm, 1315, "N").idx + 1,
    }
    variant["roles"].update(ligand_roles)
    role_names = {
        "thr_n_terminus": ("N", 7, 2),
        "nucleophile": ("OG1", 8, 1),
        "oxyanion_donor_backbone": ("N", 7, 1),
        "oxyanion_donor_sidechain": ("ND2", 7, 2),
        "asp272_od1": ("OD1", 8, 0),
        "asp272_od2": ("OD2", 8, 0),
        "asp274_od1": ("OD1", 8, 0),
        "asp274_od2": ("OD2", 8, 0),
        "scissile_carbon": ("C", 6, 0),
        "scissile_oxygen": ("O", 8, 0),
        "leaving_nitrogen": ("N", 7, 1),
        "donor_carbonyl_c": ("C", 6, 0),
        "donor_carbonyl_o": ("O", 8, 0),
    }
    variant["role_identities"] = {
        role: {"name": name, "element": element, "bonded_hydrogens": hydrogens}
        for role, (name, element, hydrogens) in role_names.items()
    }
    selected = {row["index"] for group in variant["groups"] for row in group["atoms"]}
    boundaries = sorted(
        (number, partner.idx + 1)
        for number in selected
        for partner in parm.atoms[number - 1].bond_partners
        if partner.idx + 1 not in selected
    )
    variant["boundaries"] = [
        {"qm_index": qm_index, "mm_index": mm_index}
        for qm_index, mm_index in boundaries
    ]
    variant["topology_sha256"] = sha256(topology)
    variant["coordinates_sha256"] = sha256(coordinates)
    payload = {
        "case": "Nyl12 WT modular RESP/GAFF2 PA6 L4 E-S",
        "topology": str(topology.resolve()),
        "coordinates": str(coordinates.resolve()),
        "template_source": {
            "path": str(template.resolve()),
            "sha256": sha256(template),
        },
        "modular_reacting_amide": {
            "carbonyl_residue": 1314,
            "leaving_nitrogen_residue": 1315,
            "full_substrate_residues": [1313, 1314, 1315, 1316],
        },
        "variants": [variant],
    }
    audit = audit_es_qm_region(topology, coordinates, variant)
    if (audit["real_atom_count"], audit["link_H_count"], audit["total_QM_nuclei"]) != (
        98,
        8,
        106,
    ):
        raise ValueError("Modular E-S region must retain 98 real atoms and 8 links")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    audit_path = output.with_name(output.stem + "_audit.json")
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {"config": str(output.resolve()), "audit": str(audit_path.resolve())}
