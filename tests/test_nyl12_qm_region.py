from __future__ import annotations

from pathlib import Path

import parmed as pmd
import pytest

from nyl12_covalent.qm_region import audit_qm_region


def _topology() -> pmd.Structure:
    structure = pmd.Structure()

    def atom(
        residue: int, residue_name: str, name: str, z: int, charge: float = 0.0
    ) -> pmd.Atom:
        result = pmd.Atom(name=name, type=name, atomic_number=z, charge=charge)
        structure.add_atom(result, residue_name, residue)
        return result

    for number in range(1, 929):
        atom(number, "DUM", "X", 6)

    n929 = atom(929, "THR", "N", 7, -0.4)
    h929_1 = atom(929, "THR", "H1", 1, 0.2)
    h929_2 = atom(929, "THR", "H2", 1, 0.2)
    structure.bonds.extend((pmd.Bond(n929, h929_1), pmd.Bond(n929, h929_2)))

    # Deliberately put MM C before QM CA to test direction-independent boundary logic.
    c930 = atom(930, "THR", "C", 6)
    ca930 = atom(930, "THR", "CA", 6)
    n930 = atom(930, "THR", "N", 7)
    atom(930, "THR", "O", 8)
    structure.bonds.extend((pmd.Bond(c930, ca930), pmd.Bond(ca930, n930)))

    for number in range(931, 1313):
        atom(number, "DUM", "X", 6)

    lig_c = atom(1313, "LIG", "C1", 6, -0.1)
    lig_h = atom(1313, "LIG", "H1", 1, 0.1)
    structure.bonds.append(pmd.Bond(lig_c, lig_h))

    water_o = atom(1314, "WAT", "O", 8)
    water_h1 = atom(1314, "WAT", "H1", 1)
    water_h2 = atom(1314, "WAT", "H2", 1)
    structure.bonds.extend((pmd.Bond(water_o, water_h1), pmd.Bond(water_o, water_h2)))
    return structure


def _audit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    charge: int = 0,
) -> dict[str, object]:
    topology_path = tmp_path / "ti1.prmtop"
    topology_path.write_text("topology fixture\n", encoding="utf-8")
    source_path = tmp_path / "input.pdb"
    source_path.write_text("source fixture\n", encoding="utf-8")
    topology = _topology()
    monkeypatch.setattr("nyl12_covalent.qm_region.AmberParm", lambda _: topology)
    return audit_qm_region(topology_path, [source_path], charge)


def test_audit_reports_selection_charge_boundary_and_hashes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    report = _audit(tmp_path, monkeypatch)

    assert report["qmmask"] == ":929,1313,1314|(:930&!@C,O)"
    assert report["qmcharge"] == 0
    assert report["selected_partial_charge"] == pytest.approx(0.0)
    assert report["electron_count_parity"] == "even"
    assert report["link_boundary"]["qm_atom"]["name"] == "CA"
    assert report["link_boundary"]["mm_atom"]["name"] == "C"
    assert report["link_boundary"]["link_hydrogen"] == "H"
    assert len(report["selected_atoms"]) == 10
    assert report["topology_sha256"]
    assert report["sources"][0]["sha256"]
    assert report["qm_formal_charge_is_independent_of_partial_charge"] is True


def test_residue_indices_ignore_external_residue_numbers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    topology = _topology()
    for residue in topology.residues:
        residue.number += 100
    from nyl12_covalent import qm_region

    monkeypatch.setattr(qm_region, "AmberParm", lambda _: topology)
    topology_path = tmp_path / "ti1.prmtop"
    topology_path.write_text("fixture\n", encoding="utf-8")
    assert (
        qm_region.audit_qm_region(topology_path)["qmmask"]
        == ":929,1313,1314|(:930&!@C,O)"
    )


def test_rejects_hydrogen_bonded_to_unselected_atom(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    topology = _topology()
    hydrogen = topology.atoms[929]
    excluded_oxygen = next(
        atom for atom in topology.atoms if atom.name == "O" and atom.residue.idx == 929
    )
    topology.bonds.append(pmd.Bond(hydrogen, excluded_oxygen))
    from nyl12_covalent import qm_region

    monkeypatch.setattr(qm_region, "AmberParm", lambda _: topology)
    topology_path = tmp_path / "ti1.prmtop"
    topology_path.write_text("fixture\n", encoding="utf-8")
    with pytest.raises(ValueError, match="bonded to an unselected atom"):
        qm_region.audit_qm_region(topology_path)


def test_rejects_odd_electron_count_for_singlet(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    topology = _topology()
    topology.atoms[929].atomic_number = 2
    with pytest.raises(ValueError, match="odd electron count"):
        from nyl12_covalent import qm_region

        monkeypatch.setattr(qm_region, "AmberParm", lambda _: topology)
        topology_path = tmp_path / "ti1.prmtop"
        topology_path.write_text("fixture\n", encoding="utf-8")
        qm_region.audit_qm_region(topology_path)


@pytest.mark.parametrize("charge", [True, 0.0, "0", -1, 1])
def test_rejects_non_target_formal_charge_or_type(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    charge: object,
) -> None:
    with pytest.raises(ValueError, match="requires integer qmcharge=0"):
        _audit(tmp_path, monkeypatch, charge=charge)  # type: ignore[arg-type]


def test_does_not_force_partial_charge_to_equal_formal_charge(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    topology = _topology()
    topology.atoms[928].charge += 0.37
    from nyl12_covalent import qm_region

    monkeypatch.setattr(qm_region, "AmberParm", lambda _: topology)
    topology_path = tmp_path / "ti1.prmtop"
    topology_path.write_text("fixture\n", encoding="utf-8")
    report = qm_region.audit_qm_region(topology_path)
    assert report["selected_partial_charge"] == pytest.approx(0.37)
    assert report["qmcharge"] == 0
