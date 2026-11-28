"""Verify genuine four-copy crystal templates without substrate coordinates."""

from pathlib import Path
import json

import pytest

from nyl12_covalent import crystal_protein_templates as templates

SOURCE = Path("cif/n12c.cif")
CANONICAL = Path(
    "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/canonical_protection.json"
)


def test_actual_crystal_mature_templates_preserve_coordinates_and_gaps(
    tmp_path: Path,
) -> None:
    """Keep all four actual chains, every observed protein atom and ordinary mismatch."""
    query = json.loads(CANONICAL.read_text())["source_sequences"]
    result = templates.prepare_crystal_templates(SOURCE, query, tmp_path / "templates")
    assert set(result["chains"]) == set("ABCDEFGH")
    assert result["ligand_coordinates"] == result["waters"] == result["additives"] == 0
    original = templates.crystal_records(SOURCE)
    combined = templates.crystal_records(
        tmp_path / "templates/assembly_eight_mature_protein_chains.cif"
    )
    assert len(combined) == len(original)
    for chain, info in result["chains"].items():
        rows = templates.crystal_records(tmp_path / "templates" / info["template_file"])
        assert {r["auth_asym_id"] for r in rows} == {chain}
        assert all(r["group_PDB"] == "ATOM" for r in rows)
        expected = {
            (r["auth_seq_id"], r["label_atom_id"]): r
            for r in original
            if r["auth_asym_id"] == info["source_auth_chain"]
            and (
                int(r["auth_seq_id"]) <= 232
                if info["role"] == "alpha"
                else int(r["auth_seq_id"]) >= 233
            )
        }
        assert len(rows) == len(expected)
        for r in rows:
            precursor = int(r["auth_seq_id"]) + (232 if chain in "EFGH" else 0)
            source = expected[str(precursor), r["label_atom_id"]]
            assert all(
                r[k] == source[k]
                for k in (
                    "Cartn_x",
                    "Cartn_y",
                    "Cartn_z",
                    "label_comp_id",
                    "type_symbol",
                )
            )
        assert len(info["queryIndices"]) == len(info["templateIndices"])
        assert set(info["queryIndices"]).isdisjoint(info["missing_query_indices0"])
        assert sorted(info["queryIndices"] + info["missing_query_indices0"]) == list(
            range(len(info["query_sequence"]))
        )
        assert len(info["template_sequence"]) == len(info["templateIndices"])
    mutation = result["chains"]["C"]["substitutions"]
    assert mutation == [
        {
            "query_index0": 227,
            "query_residue": "P",
            "template_residue": "A",
            "source_precursor_residue": 228,
        }
    ]
    for chain in "EFGH":
        info = result["chains"][chain]
        assert info["template_sequence"][0] == "T"
        assert info["missing_query_indices0"] == [93, 94, 95]
    assert not any(
        "struct_conn" in f.read_text() for f in (tmp_path / "templates").glob("*.cif")
    )
    with pytest.raises(FileExistsError):
        templates.write_template(original, tmp_path / "templates/chain_A.cif")
    with pytest.raises(ValueError, match="Empty"):
        templates.write_template([], tmp_path / "empty.cif")


def test_invalid_crystal_and_duplicate_atoms_fail(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Invalid or duplicated source identities fail instead of fabricating geometry."""
    f = tmp_path / "invalid.cif"
    f.write_text("data_empty\n#\n")
    with pytest.raises(ValueError, match="No original"):
        templates.crystal_records(f)
    rows = templates.crystal_records(SOURCE)
    query = json.loads(CANONICAL.read_text())["source_sequences"]
    monkeypatch.setattr(templates, "crystal_records", lambda _: rows + [rows[0]])
    with pytest.raises(ValueError, match="Duplicate"):
        templates.prepare_crystal_templates(SOURCE, query, tmp_path / "duplicate")
