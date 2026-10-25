"""File-based adversarial checks for explicit ES QM core accounting."""

from __future__ import annotations

import copy
import hashlib
from pathlib import Path

import numpy as np
import parmed
import pytest
from parmed.topologyobjects import Atom, AtomType, Bond, BondType

from nyl12_covalent import es_qm_region as module


@pytest.fixture
def qm_fixture(tmp_path: Path) -> tuple[Path, Path, dict]:
    """Build six chemically capped components with real Amber files and links."""
    structure = parmed.Structure()
    groups, roles, boundaries = [], {}, []
    group_atoms: list[int] = []
    bonds: list[tuple[int, int]] = []

    def add(name: str, z: int, residue: str, number: int, selected: bool = True) -> int:
        mass = {1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999}[z]
        typ = {"CA": "XC", "CB": "2C"}.get(name, {1: "H", 6: "CT", 7: "N", 8: "O"}[z])
        atype = AtomType(typ, z, mass, z)
        atype.set_lj_params(0.1, 1.5)
        atom = Atom(name=name, type=typ, atomic_number=z, mass=mass)
        atom.atom_type = atype
        structure.add_atom(atom, residue, number)
        if selected:
            group_atoms.append(atom.idx + 1)
        return atom.idx

    def bond(left: int, right: int) -> None:
        bonds.append((left, right))

    def hydrogen(parent: int, name: str) -> None:
        atom = structure.atoms[parent]
        bond(parent, add(name, 1, atom.residue.name, atom.residue.number))

    def finish(name: str, charge: int) -> None:
        groups.append(
            {"name": name, "formal_charge": charge, "indices": group_atoms.copy()}
        )
        group_atoms.clear()

    # Neutral native Thr plus an N-methyl cap, one saturated boundary.
    n = add("N", 7, "THR", 1)
    roles["thr_n_terminus"] = n + 1
    hydrogen(n, "H1")
    hydrogen(n, "H2")
    ca = add("CA", 6, "THR", 1)
    bond(n, ca)
    hydrogen(ca, "HA")
    cb = add("CB", 6, "THR", 1)
    bond(ca, cb)
    hydrogen(cb, "HB")
    cg = add("CG2", 6, "THR", 1)
    bond(cb, cg)
    for k in range(3):
        hydrogen(cg, f"HG2{k}")
    og = add("OG1", 8, "THR", 1)
    roles["nucleophile"] = og + 1
    bond(cb, og)
    hydrogen(og, "HG1")
    c = add("C", 6, "THR", 1)
    o = add("O", 8, "THR", 1)
    bond(ca, c)
    bond(c, o)
    capn = add("N", 7, "NME", 2)
    bond(c, capn)
    hydrogen(capn, "H")
    capc = add("CA", 6, "NME", 2)
    bond(capn, capc)
    hydrogen(capc, "H1")
    hydrogen(capc, "H2")
    boundaries.append(capc)
    finish("Thr", 0)
    # N-methylacetamide representing the complete Phe backbone donor.
    methyl = add("CA", 6, "ACE", 3)
    for k in range(3):
        hydrogen(methyl, f"H{k}")
    c = add("C", 6, "ACE", 3)
    o = add("O", 8, "ACE", 3)
    bond(methyl, c)
    bond(c, o)
    roles["donor_carbonyl_c"], roles["donor_carbonyl_o"] = c + 1, o + 1
    n = add("N", 7, "PHE", 4)
    bond(c, n)
    hydrogen(n, "H")
    roles["oxyanion_donor_backbone"] = n + 1
    ca = add("CA", 6, "PHE", 4)
    bond(n, ca)
    hydrogen(ca, "H1")
    hydrogen(ca, "H2")
    boundaries.append(ca)
    finish("Phe", 0)
    for number, label, asn in [
        (5, "Asn", True),
        (6, "Asp272", False),
        (7, "Asp274", False),
    ]:
        resname = "ASN" if asn else "ASP"
        cb = add("CB", 6, resname, number)
        hydrogen(cb, "HB2")
        hydrogen(cb, "HB3")
        c = add("CG", 6, resname, number)
        bond(cb, c)
        o = add("OD1", 8, resname, number)
        bond(c, o)
        if asn:
            n = add("ND2", 7, resname, number)
            bond(c, n)
            hydrogen(n, "HD21")
            hydrogen(n, "HD22")
            roles["oxyanion_donor_sidechain"] = n + 1
        else:
            other = add("OD2", 8, resname, number)
            bond(c, other)
            roles[label.lower() + "_od1"], roles[label.lower() + "_od2"] = (
                o + 1,
                other + 1,
            )
            structure.atoms[o].charge = structure.atoms[other].charge = -0.5
        boundaries.append(cb)
        finish(label, 0 if asn else -1)
    carbon = add("C54", 6, "LIG", 8)
    oxygen = add("O49", 8, "LIG", 8)
    bond(carbon, oxygen)
    n = add("N59", 7, "LIG", 8)
    bond(carbon, n)
    hydrogen(n, "H59")
    roles.update(
        scissile_carbon=carbon + 1, scissile_oxygen=oxygen + 1, leaving_nitrogen=n + 1
    )
    methyl = add("C62", 6, "LIG", 8)
    bond(carbon, methyl)
    for k in range(3):
        hydrogen(methyl, f"H62{k}")
    nc = add("C65", 6, "LIG", 8)
    bond(n, nc)
    hydrogen(nc, "H651")
    hydrogen(nc, "H652")
    boundaries.append(nc)
    finish("PA6", 0)
    boundary_pairs = []
    for parent in boundaries:
        external = add(f"CX{parent}", 6, "ENV", 9, False)
        bond(parent, external)
        for k in range(3):
            h = add(f"H{parent}{k}", 1, "ENV", 9, False)
            bond(external, h)
        boundary_pairs.append({"qm_index": parent + 1, "mm_index": external + 1})
    bond_type = BondType(100, 1.5)
    structure.bond_types.append(bond_type)
    for left, right in bonds:
        structure.bonds.append(
            Bond(structure.atoms[left], structure.atoms[right], type=bond_type)
        )
    structure.coordinates = np.random.default_rng(1).normal(
        size=(len(structure.atoms), 3)
    )
    parm = parmed.amber.AmberParm.from_structure(structure)
    topology, coordinates = tmp_path / "es.prmtop", tmp_path / "es.rst7"
    parm.save(str(topology))
    parm.save(str(coordinates), format="rst7")
    reloaded = parmed.amber.AmberParm(str(topology), xyz=str(coordinates))
    for group in groups:
        group["atoms"] = [
            module._identity(reloaded.atoms[i - 1]) for i in group.pop("indices")
        ]
    spec = {
        "name": "synthetic",
        "qmcharge": -2,
        "multiplicity": 1,
        "allow_multiple_link_parents": False,
        "groups": groups,
        "roles": roles,
        "boundaries": boundary_pairs,
        "topology_sha256": hashlib.sha256(topology.read_bytes()).hexdigest(),
        "coordinates_sha256": hashlib.sha256(coordinates.read_bytes()).hexdigest(),
    }
    return topology, coordinates, spec


def test_real_amber_accounting(qm_fixture: tuple) -> None:
    topology, coordinates, spec = qm_fixture
    before = (topology.read_bytes(), coordinates.read_bytes())
    result = module.audit_es_qm_region(topology, coordinates, spec)
    assert before == (topology.read_bytes(), coordinates.read_bytes())
    assert result["qmcharge"] == -2 and result["electron_count_with_links"] % 2 == 0
    assert result["link_H_count"] == 6
    assert result["total_QM_nuclei"] == result["real_atom_count"] + 6
    assert len(result["proton_inventory"]["thr_n_terminus"]["hydrogen_indices"]) == 2
    assert result["source_charge_integer_residual"] == pytest.approx(0, abs=1e-6)
    assert result["embedding_and_link_distance_validated"] is False


@pytest.mark.parametrize(
    "mutation",
    [
        "hash",
        "duplicate",
        "identity",
        "mass",
        "boolindex",
        "charge",
        "multiplicity",
        "formal",
        "boundary",
        "role",
        "missing_h",
    ],
)
def test_rejects_invalid_explicit_specs(qm_fixture: tuple, mutation: str) -> None:
    topology, coordinates, original = qm_fixture
    spec = copy.deepcopy(original)
    if mutation == "hash":
        spec["topology_sha256"] = "bad"
    elif mutation == "duplicate":
        spec["groups"][0]["atoms"].append(spec["groups"][0]["atoms"][0])
    elif mutation == "identity":
        spec["groups"][0]["atoms"][0]["name"] = "WRONG"
    elif mutation == "mass":
        spec["groups"][0]["atoms"][0]["mass"] = float("nan")
    elif mutation == "boolindex":
        spec["groups"][0]["atoms"][0]["index"] = True
    elif mutation == "charge":
        spec["qmcharge"] = -1
    elif mutation == "multiplicity":
        spec["multiplicity"] = True
    elif mutation == "formal":
        spec["groups"][0]["formal_charge"] = -1
    elif mutation == "boundary":
        spec["boundaries"].pop()
    elif mutation == "role":
        spec["roles"]["nucleophile"] = spec["roles"]["thr_n_terminus"]
    else:
        spec["groups"][0]["atoms"] = [
            row for row in spec["groups"][0]["atoms"] if row["name"] != "HG1"
        ]
    with pytest.raises(ValueError):
        module.audit_es_qm_region(topology, coordinates, spec)


def test_missing_input(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        module.audit_es_qm_region(tmp_path / "missing", tmp_path / "other", {})
