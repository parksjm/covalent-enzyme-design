"""Strict native serialization checks, including identical-copy permutations."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import ipsae2_native_mapping as mapping


def test_atom_order_and_confidence_binding() -> None:
    expected = [("A", 1, "ALA", "CA", 0), ("B", 1, "ALA", "CA", 1)]
    observed = [("A", 1, "ALA", "CA", 80.0), ("B", 1, "ALA", "CA", 90.0)]
    mapping.compare_native_order(expected, observed, np.array([0.8, 0.9]))
    with pytest.raises(ValueError, match="atom order"):
        mapping.compare_native_order(expected, observed[::-1], np.array([0.8, 0.9]))
    # Even if identical copies are renamed back to A/B, token confidence exposes
    # the identity permutation. The same check applies to identical ligands.
    renamed_swap = [(*observed[0][:4], 90.0), (*observed[1][:4], 80.0)]
    with pytest.raises(ValueError, match="token binding"):
        mapping.compare_native_order(expected, renamed_swap, np.array([0.8, 0.9]))
    with pytest.raises(ValueError, match="atom count"):
        mapping.compare_native_order(expected, [], np.array([0.8, 0.9]))
    with pytest.raises(ValueError, match="Nonfinite"):
        mapping.compare_native_order(expected, observed, np.array([np.nan, 0.9]))


@pytest.mark.parametrize("bfactor", [np.nan, np.inf, -np.inf])
def test_nonfinite_cif_confidence_rejected(bfactor: float) -> None:
    expected = [("A", 1, "ALA", "CA", 0)]
    observed = [("A", 1, "ALA", "CA", bfactor)]
    with pytest.raises(ValueError, match="Nonfinite native CIF B factor"):
        mapping.compare_native_order(expected, observed, np.array([0.8]))


@pytest.fixture
def native_fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict:
    chain_dtype = [
        (key, kind)
        for key, kind in [
            ("name", "U1"),
            ("asym_id", "i4"),
            ("mol_type", "i4"),
            ("entity_id", "i4"),
            ("sym_id", "i4"),
            ("atom_num", "i4"),
            ("res_num", "i4"),
            ("res_idx", "i4"),
        ]
    ]
    chains = np.array(
        [("A", 0, 0, 0, 0, 1, 1, 0), ("I", 1, 3, 1, 0, 1, 1, 1)], dtype=chain_dtype
    )
    residues = np.array(
        [("ALA", 0, 0, 1), ("LIG1", 0, 1, 1)],
        dtype=[
            ("name", "U4"),
            ("res_idx", "i4"),
            ("atom_idx", "i4"),
            ("atom_num", "i4"),
        ],
    )
    atoms = np.array([("CA",), ("C1",)], dtype=[("name", "U2")])
    tokens = np.array(
        [(0, 0, 1, 0, "ALA", 0), (1, 1, 1, 0, "LIG1", 1)],
        dtype=[
            ("token_idx", "i4"),
            ("atom_idx", "i4"),
            ("atom_num", "i4"),
            ("res_idx", "i4"),
            ("res_name", "U4"),
            ("asym_id", "i4"),
        ],
    )
    structure = SimpleNamespace(chains=chains, residues=residues, atoms=atoms)
    monkeypatch.setattr(mapping, "_native_arrays", lambda *_: (structure, tokens))
    observed = [("A", 1, "ALA", "CA", 80.0), ("I", 1, "LIG1", "C1", 90.0)]
    monkeypatch.setattr(mapping, "_cif_atoms", lambda _: observed)
    for file in ("raw.cif", "structure.npz"):
        (tmp_path / file).write_text("fixture")
    np.savez(tmp_path / "pae.npz", pae=np.zeros((2, 2)))
    np.savez(tmp_path / "plddt.npz", plddt=np.array([0.8, 0.9]))
    (tmp_path / "confidence.json").write_text(
        '{"pair_chains_iptm":{"0":{"0":1,"1":1},"1":{"0":1,"1":1}},"chains_ptm":{"0":1,"1":1}}'
    )
    for source in (
        "data/tokenize/boltz2.py",
        "data/types.py",
        "data/write/writer.py",
        "data/write/mmcif.py",
        "model/loss/confidence.py",
        "model/modules/confidencev2.py",
    ):
        path = tmp_path / source
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("source")
    return dict(
        cif_path=tmp_path / "raw.cif",
        pae_path=tmp_path / "pae.npz",
        confidence_path=tmp_path / "confidence.json",
        plddt_path=tmp_path / "plddt.npz",
        structure_path=tmp_path / "structure.npz",
        boltz_source_root=tmp_path,
    )


def test_receipt_contains_exact_artifact_and_token_mapping(
    native_fixture: dict,
) -> None:
    receipt = mapping.qualify_native_mapping(**native_fixture)
    assert receipt["accepted"]
    assert receipt["protein_chains"] == ["A"]
    assert receipt["ligand_chains"] == ["I"]
    assert receipt["tokens"][1]["atom_names"] == ["C1"]
    assert receipt["native_chain_map"][1]["token_start"] == 1
    assert len(receipt["artifacts"]["cif"]["sha256"]) == 64
    assert receipt["physical_groups_qualified"] is False


@pytest.mark.parametrize(
    "failure", ["shape", "nonfinite", "negative", "confidence", "chain_ptm", "order"]
)
def test_receipt_rejects_invalid_inputs(
    native_fixture: dict, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    if failure == "shape":
        np.savez(native_fixture["pae_path"], pae=np.zeros((3, 3)))
    elif failure in {"nonfinite", "negative"}:
        pae = np.zeros((2, 2))
        pae[0, 1] = np.nan if failure == "nonfinite" else -1
        np.savez(native_fixture["pae_path"], pae=pae)
    elif failure in {"confidence", "chain_ptm"}:
        import json

        path = native_fixture["confidence_path"]
        confidence = json.loads(path.read_text())
        confidence["pair_chains_iptm" if failure == "confidence" else "chains_ptm"] = {}
        path.write_text(json.dumps(confidence))
    else:
        monkeypatch.setattr(
            mapping,
            "_cif_atoms",
            lambda _: [("I", 1, "LIG1", "C1", 90.0), ("A", 1, "ALA", "CA", 80.0)],
        )
    receipt = mapping.qualify_native_mapping(**native_fixture)
    assert not receipt["accepted"]
    assert receipt["errors"]


@pytest.mark.parametrize(
    "element,residue,model,error",
    [
        ("C", "ALA", "1", None),
        ("H", "ALA", "1", "hydrogen"),
        ("C", "WAT", "1", "water"),
        ("C", "ALA", "2", "single native"),
    ],
)
def test_cif_reader(
    monkeypatch: pytest.MonkeyPatch,
    element: str,
    residue: str,
    model: str,
    error: str | None,
) -> None:
    import sys

    rows = [["A", "1", residue, "CA", "80.0", element, model]]
    block = SimpleNamespace(find=lambda *args: rows)
    document = SimpleNamespace(sole_block=lambda: block)
    monkeypatch.setitem(
        sys.modules,
        "gemmi",
        SimpleNamespace(cif=SimpleNamespace(read_file=lambda _: document)),
    )
    if error:
        with pytest.raises(ValueError, match=error):
            mapping._cif_atoms(Path("fake.cif"))
    else:
        assert mapping._cif_atoms(Path("fake.cif")) == [("A", 1, "ALA", "CA", 80.0)]
        rows.clear()
        with pytest.raises(ValueError, match="No native CIF"):
            mapping._cif_atoms(Path("fake.cif"))
