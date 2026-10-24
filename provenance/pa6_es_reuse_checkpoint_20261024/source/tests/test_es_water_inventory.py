"""Synthetic file-based coverage for the descriptive water inventory."""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import parmed
import pytest
from parmed.topologyobjects import Atom, AtomType, Bond, BondType

from nyl12_covalent import es_water_inventory as module


@pytest.fixture
def water_fixture(tmp_path: Path) -> tuple[Path, Path, dict]:
    """Write a small actual Amber topology and restart with two waters."""
    structure = parmed.Structure()
    specification = [
        ("N", 7, "THR", 1, (0, 0, 0)),
        ("H1", 1, "THR", 1, (0, 1, 0)),
        ("OG1", 8, "THR", 1, (2, 0, 0)),
        ("HG1", 1, "THR", 1, (2, 1, 0)),
        ("OD1", 8, "ASP", 2, (0, 4, 0)),
        ("OD2", 8, "ASP", 2, (1, 4, 0)),
        ("OD1", 8, "ASP", 3, (2, 4, 0)),
        ("OD2", 8, "ASP", 3, (3, 4, 0)),
        ("O", 8, "WAT", 4, (1, 2, 0)),
        ("H1", 1, "WAT", 4, (1, 1, 0)),
        ("H2", 1, "WAT", 4, (2, 2, 0)),
        ("O", 8, "WAT", 5, (20, 20, 20)),
        ("H1", 1, "WAT", 5, (20, 21, 20)),
        ("H2", 1, "WAT", 5, (21, 20, 20)),
    ]
    masses = {1: 1.008, 7: 14.007, 8: 15.999}
    for name, element, residue, number, _ in specification:
        atom_type = AtomType(f"T{element}", element, masses[element], element)
        atom_type.set_lj_params(0.1, 1.5)
        atom = Atom(
            name=name,
            type=f"T{element}",
            atomic_number=element,
            mass=masses[element],
            charge=0,
        )
        atom.atom_type = atom_type
        structure.add_atom(atom, residue, number)
    bond_type = BondType(100, 1)
    structure.bond_types.append(bond_type)
    for left, right in [(0, 1), (2, 3), (8, 9), (8, 10), (11, 12), (11, 13)]:
        structure.bonds.append(
            Bond(structure.atoms[left], structure.atoms[right], type=bond_type)
        )
    structure.coordinates = np.asarray([row[4] for row in specification], dtype=float)
    parm = parmed.amber.AmberParm.from_structure(structure)
    topology, coordinates = tmp_path / "es.prmtop", tmp_path / "es.rst7"
    parm.save(str(topology))
    parm.save(str(coordinates), format="rst7")
    roles = {
        "thr_n_terminus": {"residue": 1, "resname": "THR", "atom": "N"},
        "nucleophile": {"residue": 1, "resname": "THR", "atom": "OG1"},
        "asp272_od1": {"residue": 2, "resname": "ASP", "atom": "OD1"},
        "asp272_od2": {"residue": 2, "resname": "ASP", "atom": "OD2"},
        "asp274_od1": {"residue": 3, "resname": "ASP", "atom": "OD1"},
        "asp274_od2": {"residue": 3, "resname": "ASP", "atom": "OD2"},
    }
    return topology, coordinates, roles


def test_file_inventory(water_fixture: tuple) -> None:
    topology, coordinates, roles = water_fixture
    before = [path.read_bytes() for path in (topology, coordinates)]
    result = module.inventory_es_waters(topology, coordinates, roles)
    assert before == [path.read_bytes() for path in (topology, coordinates)]
    assert result["local_water_count"] == result["geometric_only_candidate_count"] == 1
    water = result["waters"][0]
    assert water["oxygen"]["index"] == 9
    assert water["residue_index"] == 4
    assert water["hydrogens"][0]["index"] == 10
    assert water["oxygen_distances_A"]["nucleophile"] == pytest.approx(5**0.5)
    assert water["nearest_asp_oxygen_distance_A"] == pytest.approx(2)
    geometry = water["donor_h_acceptor_geometry"]
    assert geometry["thr_n_h_water_o"]["best_angle_degrees"] == pytest.approx(135)
    assert geometry["thr_n_h_water_o"]["angle_at_least_120_degrees_count"] == 1
    assert geometry["water_o_h_thr_n"]["explicit_hydrogen_count"] == 2


@pytest.mark.parametrize("mutation", ["missing", "bool", "name", "atom", "row", "pair"])
def test_bad_mapping(water_fixture: tuple, mutation: str) -> None:
    topology, coordinates, original = water_fixture
    roles = copy.deepcopy(original)
    if mutation == "missing":
        del roles["nucleophile"]
    elif mutation == "row":
        roles["nucleophile"] = "invalid"
    elif mutation == "pair":
        roles["asp274_od2"]["residue"] = 2
    else:
        field, value = {
            "bool": ("residue", True),
            "name": ("resname", "SER"),
            "atom": ("atom", "O"),
        }[mutation]
        roles["nucleophile"][field] = value
    with pytest.raises(ValueError):
        module.inventory_es_waters(topology, coordinates, roles)


@pytest.mark.parametrize(
    "mutation",
    ["missing_h", "ambiguous", "nonfinite", "shape", "coincident", "no_thr_h", "far"],
)
def test_invalid_or_changed_geometry(
    water_fixture: tuple, monkeypatch: pytest.MonkeyPatch, mutation: str
) -> None:
    topology, coordinates, roles = water_fixture
    parm = parmed.amber.AmberParm(str(topology), xyz=str(coordinates))
    if mutation in {"missing_h", "no_thr_h"}:
        pair = {8, 9} if mutation == "missing_h" else {2, 3}
        bond = next(b for b in parm.bonds if {b.atom1.idx, b.atom2.idx} == pair)
        bond.delete()
        parm.bonds.remove(bond)
    elif mutation == "ambiguous":
        parm.atoms[9].atomic_number = 8
    else:
        xyz = parm.coordinates.copy()
        if mutation == "nonfinite":
            xyz[0, 0] = np.nan
        elif mutation == "shape":
            monkeypatch.setattr(
                type(parm), "coordinates", property(lambda self: np.zeros((1, 3)))
            )
        elif mutation == "coincident":
            xyz[1] = xyz[0]
        elif mutation == "far":
            xyz[8:11] += 20
        if mutation != "shape":
            parm.coordinates = xyz
    monkeypatch.setattr(module, "AmberParm", lambda *args, **kwargs: parm)
    if mutation in {"no_thr_h", "far"}:
        result = module.inventory_es_waters(topology, coordinates, roles)
        if mutation == "far":
            assert result["local_water_count"] == 0
        else:
            assert (
                result["waters"][0]["donor_h_acceptor_geometry"]["thr_og1_h_water_o"][
                    "best_angle_degrees"
                ]
                is None
            )
    else:
        with pytest.raises(ValueError):
            module.inventory_es_waters(topology, coordinates, roles)


def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        module.inventory_es_waters(tmp_path / "absent", tmp_path / "other", {})
