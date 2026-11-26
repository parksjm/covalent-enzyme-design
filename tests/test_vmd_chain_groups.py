"""Validate native and reconstructed display chains without touching topology."""

from types import SimpleNamespace

import pytest

from nyl12_covalent.vmd_chain_groups import display_chain_groups


class _Topology:
    def __init__(self, residues: list[SimpleNamespace]) -> None:
        self.residues = residues

    @property
    def atoms(self) -> list[SimpleNamespace]:
        return [atom for residue in self.residues for atom in residue.atoms]


def topology_fixture() -> _Topology:
    """Build minimal mature Nyl12 residues with distinct physical C/N atoms."""
    residues = []
    for index in range(1312):
        beta_start = index in {928, 1024, 1120, 1216}
        residue = SimpleNamespace(
            idx=index, name="THR" if beta_start else "ALA", chain="", atoms=[]
        )
        for offset, name in enumerate(("N", "C")):
            residue.atoms.append(
                SimpleNamespace(
                    idx=index * 2 + offset, name=name, bond_partners=[], residue=residue
                )
            )
        residues.append(residue)
    return _Topology(residues)


def test_reconstructs_mature_groups_without_mutating_native_labels() -> None:
    top = topology_fixture()
    mapping = {index + 1: index for index in range(2624)}
    groups = display_chain_groups(top, mapping)
    assert set(groups) == set("ABCDEFGH")
    assert len(groups["A"]) == 464 and len(groups["E"]) == 192
    assert groups["E"][0] == 1856 and groups["H"][-1] == 2623
    assert all(residue.chain == "" for residue in top.residues)


def test_preserves_explicit_chains_and_stripped_mapping() -> None:
    top = topology_fixture()
    top.residues = top.residues[:2]
    top.residues[0].chain = "A"
    top.residues[1].chain = "B"
    assert display_chain_groups(top, {1: 9, 4: 2}) == {"A": [9], "B": [2]}


@pytest.mark.parametrize(
    "issue", ["partial", "count", "sequence", "thr", "bond", "atom"]
)
def test_rejects_unjustified_chain_assignments(issue: str) -> None:
    top = topology_fixture()
    if issue == "partial":
        top.residues[0].chain = "A"
    elif issue == "count":
        top.residues.pop()
    elif issue == "sequence":
        top.residues[232].name = "GLY"
    elif issue == "thr":
        for index in (928, 1024, 1120, 1216):
            top.residues[index].name = "ALA"
    elif issue == "bond":
        top.residues[231].atoms[1].bond_partners.append(top.residues[232].atoms[0])
    else:
        top.residues[231].atoms.pop()
    with pytest.raises(ValueError):
        display_chain_groups(top, {1: 0})


def test_rejects_empty_protein_or_empty_display() -> None:
    with pytest.raises(ValueError, match="No protein"):
        display_chain_groups(_Topology([]), {})
    with pytest.raises(ValueError, match="No mapped"):
        display_chain_groups(topology_fixture(), {})
