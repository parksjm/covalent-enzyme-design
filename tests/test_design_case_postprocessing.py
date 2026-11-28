"""Check actual-case export orchestration without launching native engines."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from nyl12_covalent import design_case_postprocessing as target
from nyl12_covalent.relay_any_amide import load_graph

ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/wt_protein_template_controls_20261127/inspection_export_20261127"
)
LEGACY = SOURCE.parent / "amber_relaxation_20261127/postprocess_template_relaxation.py"


@pytest.fixture
def inputs(tmp_path: Path) -> dict[str, Any]:
    """Copy actual validated map/graph and create explicit native placeholders."""
    kwargs: dict[str, Any] = {
        "case": "relay_variant_alpha",
        "target_phi_sign": -1,
        "output_dir": tmp_path / "publication",
    }
    for name in [
        "topology",
        "before_restart",
        "after_restart",
        "conversion",
        "legacy_processor",
    ]:
        path = tmp_path / name
        path.write_text(name)
        kwargs[name] = path
    for name, filename in [
        ("whole_context_mapping", "whole_l4_context_mapping.json"),
        ("source_mol2", "whole_l4_source.mol2"),
    ]:
        path = tmp_path / filename
        path.write_bytes((SOURCE / filename).read_bytes())
        kwargs[name] = path
    return kwargs


def geometry_rows(source_mol2: Path) -> list[dict[str, Any]]:
    """Create all48 site/ligand rows with actual graph-derived amide identities."""
    amides = load_graph(source_mol2).amides
    return [
        {
            "site": site,
            "ligand": ligand,
            "carbonyl_c": amide.carbon,
            "carbonyl_o": amide.oxygen,
            "leaving_n": amide.nitrogen,
            "d1": 3.0,
            "herschlag_phi": 85.0,
            "tier": "A",
        }
        for site in ["AE", "BF", "CG", "DH"]
        for ligand in "IJKL"
        for amide in amides
    ]


def install_exporter(
    monkeypatch: pytest.MonkeyPatch, rows: list[dict[str, Any]], calls: list[Any]
) -> None:
    """Replace dynamic module execution with a captured, complete native export."""

    def process(args: Any) -> None:
        calls.append(args)
        args.output_dir.mkdir()
        (args.output_dir / "README.md").write_text("legacy WT evidence unchanged")
        for state in ["before", "after"]:
            folder = args.output_dir / f"{args.case}_{state}"
            folder.mkdir()
            (folder / "geometry_any_amide.json").write_text(
                json.dumps({"all_amide_site_rows": rows})
            )

    module = SimpleNamespace(process=process)
    spec = SimpleNamespace(loader=SimpleNamespace(exec_module=lambda _: None))
    monkeypatch.setattr(
        target.importlib.util, "spec_from_file_location", lambda *_: spec
    )
    monkeypatch.setattr(target.importlib.util, "module_from_spec", lambda _: module)


def test_actual_case_all48_negative_selection_and_immutable_evidence(
    inputs: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    rows = geometry_rows(inputs["source_mol2"])
    # The positive A on this pair is closer; a different negative B must survive.
    rows[0]["d1"] = 2.5
    rows[1].update(tier="B", herschlag_phi=-85.0, d1=3.1)
    original = deepcopy(rows)
    legacy_hash = hashlib.sha256(LEGACY.read_bytes()).hexdigest()
    calls: list[Any] = []
    install_exporter(monkeypatch, rows, calls)
    report = target.process_design_case(**inputs)
    assert calls[0].case == "relay_variant_alpha"
    assert calls[0].before_restart == inputs["before_restart"]
    assert calls[0].after_restart == inputs["after_restart"]
    assert calls[0].output_dir == inputs["output_dir"] / "canonical_before_after"
    assert report["case"] == calls[0].case
    for state in ["before", "after"]:
        gate = report["snapshots"][state]
        assert gate["input_row_count"] == 48
        assert gate["assigned_ab_count"] == 1
        assert gate["assigned"][0]["tier"] == "B"
        assert gate["all_amide_site_rows"][0]["tier"] == "A"
        saved = json.loads(
            (inputs["output_dir"] / f"pose_gate_{state}.json").read_text()
        )
        assert saved == gate
        assert {tuple(r["monomer_boundary"]) for r in gate["all_amide_site_rows"]} == {
            (1, 2),
            (2, 3),
            (3, 4),
        }
    assert rows == original
    assert hashlib.sha256(LEGACY.read_bytes()).hexdigest() == legacy_hash
    assert "relay_variant_alpha" in (inputs["output_dir"] / "README.md").read_text()
    assert (
        calls[0].output_dir / "README.md"
    ).read_text() == "legacy WT evidence unchanged"
    assert len(report["input_sha256"]) == 7
    assert (inputs["output_dir"] / "case_postprocessing.json").exists()


def test_exact_real_source_registers(inputs: dict[str, Any]) -> None:
    rows = geometry_rows(inputs["source_mol2"])
    original = deepcopy(rows)
    annotated = target.annotate_amide_registers(
        rows,
        whole_context_mapping=inputs["whole_context_mapping"],
        source_mol2=inputs["source_mol2"],
    )
    expected = {
        "C55": ([3, 4], [3, 1]),
        "C54": ([2, 3], [2, 2]),
        "C52": ([1, 2], [1, 3]),
    }
    for row in annotated:
        assert (
            row["monomer_boundary"],
            row["hypothetical_fragment_lengths"],
        ) == expected[row["carbonyl_c"]]
        assert row["register_is_product_hypothesis"] is True
    assert rows == original


@pytest.mark.parametrize(
    "mutation,match",
    [
        ("duplicate", "Duplicate"),
        ("missing", "cover every"),
        ("unit", "four integer"),
        ("element", "element mismatch"),
        ("oxygen", "share C unit"),
        ("nitrogen", "share C unit"),
    ],
)
def test_bad_mapping_fails_before_export(
    inputs: dict[str, Any], mutation: str, match: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = json.loads(inputs["whole_context_mapping"].read_text())
    assignments = data["blocks"][0]["assignments"]
    if mutation == "duplicate":
        assignments.append(dict(assignments[0]))
    elif mutation == "missing":
        assignments.pop()
    elif mutation == "unit":
        assignments[0]["residue_index_1based"] = 0
    elif mutation == "element":
        assignments[0]["element"] = "Xe"
    else:
        from rdkit import Chem

        mol = Chem.MolFromMol2File(str(inputs["source_mol2"]), removeHs=False)
        name = "O50" if mutation == "oxygen" else "N57"
        index = next(
            a.GetIdx() + 1
            for a in mol.GetAtoms()
            if a.GetProp("_TriposAtomName") == name
        )
        next(a for a in assignments if a["atom_index_1based"] == index)[
            "residue_index_1based"
        ] = 1
    inputs["whole_context_mapping"].write_text(json.dumps(data))
    calls: list[Any] = []
    install_exporter(monkeypatch, geometry_rows(inputs["source_mol2"]), calls)
    with pytest.raises(ValueError, match=match):
        target.process_design_case(**inputs)
    assert calls == []
    assert not inputs["output_dir"].exists()


def test_unknown_amide_row_fails(inputs: dict[str, Any]) -> None:
    rows = geometry_rows(inputs["source_mol2"])
    rows[0]["carbonyl_c"] = "C_fake"
    with pytest.raises(ValueError, match="absent from audited graph"):
        target.annotate_amide_registers(
            rows,
            whole_context_mapping=inputs["whole_context_mapping"],
            source_mol2=inputs["source_mol2"],
        )


@pytest.mark.parametrize("case", ["", "../case", ".", ".."])
def test_invalid_case_fails(inputs: dict[str, Any], case: str) -> None:
    inputs["case"] = case
    with pytest.raises(ValueError, match="actual case identity"):
        target.process_design_case(**inputs)


def test_existing_output_never_overwritten(inputs: dict[str, Any]) -> None:
    inputs["output_dir"].mkdir()
    with pytest.raises(FileExistsError):
        target.process_design_case(**inputs)


def test_missing_input_fails_before_output(inputs: dict[str, Any]) -> None:
    inputs["before_restart"].unlink()
    with pytest.raises(FileNotFoundError):
        target.process_design_case(**inputs)
    assert not inputs["output_dir"].exists()


def test_missing_dynamic_import_spec_fails(
    inputs: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        target.importlib.util, "spec_from_file_location", lambda *_: None
    )
    with pytest.raises(ImportError):
        target.process_design_case(**inputs)


def test_exporter_input_mutation_is_rejected(
    inputs: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    rows = geometry_rows(inputs["source_mol2"])
    calls: list[Any] = []
    install_exporter(monkeypatch, rows, calls)
    original_write = target._write

    def write_and_mutate(path: Path, value: Any) -> None:
        original_write(path, value)
        inputs["topology"].write_text("unexpected mutation")

    monkeypatch.setattr(target, "_write", write_and_mutate)
    with pytest.raises(ValueError, match="input changed"):
        target.process_design_case(**inputs)
    assert not (inputs["output_dir"] / "case_postprocessing.json").exists()


def test_explicit_positive_target_publishes_generic_case_gate(
    inputs: dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    inputs["target_phi_sign"] = +1
    rows = geometry_rows(inputs["source_mol2"])
    original = deepcopy(rows)
    calls: list[Any] = []
    install_exporter(monkeypatch, rows, calls)
    result = target.process_design_case(**inputs)
    assert result["target_phi_sign"] == +1
    for state in ["before", "after"]:
        gate = result["snapshots"][state]
        assert gate["target_phi_sign"] == +1
        assert gate["eligible_ab_count"] == 48
        assert gate["assigned_ab_count"] == 4
        assert gate["passes_pose_gate"] is True
        assert gate["signed_ab_counts"] == {"positive": 48, "negative": 0}
        assert (inputs["output_dir"] / f"pose_gate_{state}.json").exists()
        assert not (inputs["output_dir"] / f"negative_pose_gate_{state}.json").exists()
    assert rows == original
    assert "target phi sign +1" in (inputs["output_dir"] / "README.md").read_text()


@pytest.mark.parametrize("sign", [0, True, False, 1.0])
def test_invalid_postprocessor_target_fails_before_export(
    inputs: dict[str, Any], sign: Any
) -> None:
    inputs["target_phi_sign"] = sign
    with pytest.raises(ValueError, match="target_phi_sign"):
        target.process_design_case(**inputs)
    assert not inputs["output_dir"].exists()
