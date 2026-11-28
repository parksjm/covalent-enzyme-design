"""Focused scientific input and original-coordinate retention checks."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parents[1] / "scripts/refine_rebuilt_protein_atoms.py"
SPEC = importlib.util.spec_from_file_location("refinement", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_core_preserves_atom_identity_despite_reordering() -> None:
    """An added atom and arbitrary dictionary order cannot change the core map."""
    original = {("A", 1, " ", "CA"): (1.0, 2.0, 3.0)}
    observed = {("B", 2, " ", "CA"): (8.0, 8.0, 8.0), **original}
    assert MODULE.assert_core_retained(original, observed, 0.002) == 0.0


def test_core_rejects_missing_experimental_atom() -> None:
    """Missing original atoms must fail rather than be remapped by list order."""
    with pytest.raises(ValueError, match="missing"):
        MODULE.assert_core_retained({"CA": (0, 0, 0)}, {"CB": (0, 0, 0)}, 0.002)


def test_core_rejects_motion_beyond_gate() -> None:
    """The experimental core retention gate cannot be silently relaxed."""
    with pytest.raises(ValueError, match="moved"):
        MODULE.assert_core_retained({"CA": (0, 0, 0)}, {"CA": (0.003, 0, 0)}, 0.002)


def test_sha256_matches_standard_reference(tmp_path: Path) -> None:
    """Provenance hashes identify exact file contents."""
    source = tmp_path / "input"
    source.write_bytes(b"abc")
    assert MODULE.sha256(source) == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_geometry_audit_rejects_bad_peptide_and_reports_clash(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Invalid modeled peptides and clashes must both reject a complete sequence."""
    import sys
    from types import ModuleType, SimpleNamespace

    chain = SimpleNamespace(id="A")
    left = SimpleNamespace(id="1", name="ALA", insertionCode="", chain=chain)
    right = SimpleNamespace(id="2", name="GLY", insertionCode="", chain=chain)
    left.atoms = lambda: [SimpleNamespace(name="N")]
    right.atoms = lambda: [SimpleNamespace(name="OXT")]
    chain.residues = lambda: iter([left, right])
    topology = SimpleNamespace(chains=lambda: iter([chain]), bonds=lambda: iter([]))
    fake_app = ModuleType("openmm.app")
    fake_app.PDBFile = lambda _: SimpleNamespace(topology=topology)
    monkeypatch.setitem(sys.modules, "openmm", ModuleType("openmm"))
    monkeypatch.setitem(sys.modules, "openmm.app", fake_app)
    fake_postflight = ModuleType("postflight_nylc_5xyg_repair")
    fake_postflight.clash_audit = lambda *_: {"clash_count": 1}
    monkeypatch.setitem(sys.modules, fake_postflight.__name__, fake_postflight)
    fake_repair = ModuleType("repair_nylc_5xyg_pdbfixer")
    original = {("A", 1, " ", "C"): (0, 0, 0), ("A", 2, " ", "N"): (2, 0, 0)}
    fake_repair.read_atoms = lambda _: original
    monkeypatch.setitem(sys.modules, fake_repair.__name__, fake_repair)
    result = MODULE.geometry_audit(
        Path("candidate.pdb"), original, {"A": ["ALA", "GLY"]}, {"A": {"n": 1, "c": 2}}
    )
    assert result["failures"] == ["peptide_geometry_or_bond", "modeled_atom_clashes"]
    assert result["bad_peptides"][0]["distance_A"] == 2


def test_atom_key_retains_chain_number_and_insertion() -> None:
    """Author numbering and insertion codes define frozen atom identity."""
    from types import SimpleNamespace

    residue = SimpleNamespace(
        chain=SimpleNamespace(id="D"), id="267", insertionCode="A"
    )
    assert MODULE.atom_key(SimpleNamespace(residue=residue, name="N")) == (
        "D",
        267,
        "A",
        "N",
    )
