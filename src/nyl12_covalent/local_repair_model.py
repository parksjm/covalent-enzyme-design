"""Prepare the explicit alpha/beta PA6 TI1 diagnostic without reparameterizing."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm, AmberMask

from .af3_postprocess import read_atom_sites
from .local_hydrogens import transfer_attached_coordinates
from .local_water_seed import seed_water

QM_MASK = ":233,329,330|(:234&!@C,O)"


def identity(residue_index: int) -> tuple[str, int]:
    """Return canonical chain and residue identities for the 330-residue test."""
    if 0 <= residue_index < 232:
        return "A", residue_index + 1
    if 232 <= residue_index < 328:
        return "E", residue_index - 231
    if residue_index in (328, 329):
        return ("I" if residue_index == 328 else "J"), 1
    raise ValueError("Unexpected dimer residue index")


def export_pdb(model: AmberParm, path: Path) -> None:
    """Write inspection coordinates with canonical A/E/I/J residue numbering."""
    lines = []
    previous = None
    elements = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
    for atom, xyz in zip(model.atoms, model.coordinates, strict=True):
        chain, number = identity(atom.residue.idx)
        if previous is not None and previous != chain:
            lines.append("TER")
        previous = chain
        x, y, z = xyz
        lines.append(
            f"ATOM  {atom.idx + 1:5d} {atom.name:>4s} {atom.residue.name:>3s} "
            f"{chain}{number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00          {elements[atom.atomic_number]:>2s}"
        )
    path.write_text("\n".join([*lines, "TER", "END", ""]))


def prepare_model(topology: Path, coordinates: Path, cif: Path, output: Path) -> dict:
    """Extract the covalent dimer, map AF3 heavy atoms and construct hydrogens.

    Args:
        topology: Original product-gated 1314-residue TI1 Amber topology.
        coordinates: Matching original coordinates used for bonded H frames.
        cif: Two-restraint AF3 alpha/beta prediction.
        output: Existing fresh attempt directory.

    Returns:
        Explicit mapping, water parameters and hydrogen construction provenance.
    """
    source = AmberParm(str(topology), xyz=str(coordinates))
    if len(source.residues) != 1314:
        raise ValueError("Expected original 1314-residue topology")
    source.strip("!(:1-232,929-1024,1313-1314)")
    if len(source.residues) != 330:
        raise ValueError("Dimer extraction failed")
    old = np.array(source.coordinates, copy=True)
    sites = read_atom_sites(cif)
    mapped = {(a.chain, a.sequence, a.name): a for a in sites}
    if len(mapped) != len(sites):
        raise ValueError("Duplicate AF3 atom identity")
    heavy = {}
    missing = []
    used = set()
    for atom in source.atoms:
        if atom.atomic_number == 1:
            continue
        chain, residue = identity(atom.residue.idx)
        key = chain, residue, atom.name
        if key not in mapped:
            missing.append(atom)
            continue
        site = mapped[key]
        expected = {
            "HIE": "HIS",
            "HID": "HIS",
            "HIP": "HIS",
            "LIG": "PA6TI1",
            "WAT": "HOH",
        }.get(atom.residue.name, atom.residue.name)
        if site.residue != expected:
            raise ValueError(f"Residue mismatch at {key}: {site.residue} != {expected}")
        heavy[atom.idx] = np.array(site.xyz)
        used.add(key)
    if set(mapped) != used:
        raise ValueError(f"AF3 atoms absent from topology: {set(mapped) - used}")
    rebuilt = []
    for atom in missing:
        if atom.name != "OXT" or atom.residue.idx not in (231, 327):
            raise ValueError(f"Unmapped heavy atom {atom.residue.idx + 1}:{atom.name}")
        group = {a.name: a.idx for a in atom.residue.atoms}
        center = group["C"]
        anchors = [group["CA"], group["O"]]
        a = old[anchors] - old[center]
        b = np.array([heavy[i] - heavy[center] for i in anchors])
        u, _, vt = np.linalg.svd(a.T @ b)
        correction = np.diag([1.0, 1.0, np.linalg.det(u @ vt)])
        rotation = u @ correction @ vt
        heavy[atom.idx] = heavy[center] + (old[atom.idx] - old[center]) @ rotation
        rebuilt.append(
            {
                "residue": atom.residue.idx + 1,
                "atom": "OXT",
                "method": "original C-OXT vector fitted to AF3 CA/C/O frame",
            }
        )
    water = source.residues[329]
    water_start = water.atoms[0].idx
    if [a.name for a in water.atoms] != ["O", "H1", "H2"]:
        raise ValueError("Expected TIP3P WAT O/H1/H2 order")
    xyz, hydrogen_report = transfer_attached_coordinates(
        old[:water_start],
        {i: p for i, p in heavy.items() if i < water_start},
        [a.atomic_number for a in source.atoms[:water_start]],
        [
            (b.atom1.idx, b.atom2.idx)
            for b in source.bonds
            if b.atom1.idx < water_start and b.atom2.idx < water_start
        ],
    )
    thr = {a.name: a for a in source.residues[232].atoms}
    wat_xyz, seed_report = seed_water(
        heavy[water_start],
        xyz[thr["N"].idx],
        xyz[thr["OG1"].idx],
        xyz,
        np.array([a.atomic_number for a in source.atoms[:water_start]]),
        preserve_oxygen=True,
    )
    source.coordinates = np.vstack([xyz, wat_xyz])
    water_parameters = [
        {
            "name": a.name,
            "type": a.type,
            "charge": a.charge,
            "rmin_half_A": a.rmin,
            "epsilon_kcal_mol": a.epsilon,
        }
        for a in water.atoms
    ]
    if (
        water.name != "WAT"
        or water.atoms[0].type != "OW"
        or water.atoms[0].epsilon <= 0
    ):
        raise ValueError("Missing expected TIP3P oxygen Lennard-Jones parameters")
    if abs(sum(a.charge for a in water.atoms)) > 1e-6:
        raise ValueError("Water is not neutral")
    source.write_parm(str(output / "model.prmtop"))
    source.save(str(output / "prepared.rst7"), format="rst7", overwrite=False)
    export_pdb(source, output / "prepared.pdb")
    selected = list(AmberMask(source, QM_MASK).Selected())
    selected_set = set(selected)
    cuts = [
        (b.atom1, b.atom2)
        for b in source.bonds
        if (b.atom1.idx in selected_set) != (b.atom2.idx in selected_set)
    ]
    if (
        len(selected) != 109
        or len(cuts) != 1
        or {(a.residue.idx, a.name) for a in cuts[0]} != {(233, "CA"), (233, "C")}
    ):
        raise ValueError("QM boundary differs from the audited chemical selection")
    report = {
        "system": "PA6 L4 TI1 alpha A / beta E / ligand I / water J",
        "original_residue_selection": "1-232,929-1024,1313-1314",
        "source_parameters_changed": False,
        "heavy_atoms_mapped_from_AF3": len(used),
        "rebuilt_terminal_atoms": rebuilt,
        "hydrogen_transfer": hydrogen_report,
        "water_seed": seed_report,
        "water_parameters": water_parameters,
        "qmmask": QM_MASK,
        "qmcharge": 0,
        "qm_atoms": len(selected),
        "qm_electrons_with_link_H": sum(source.atoms[i].atomic_number for i in selected)
        + 1,
        "accepted_minimum": False,
    }
    (output / "preparation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
