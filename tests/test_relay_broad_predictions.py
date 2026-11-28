"""Scientific-input and bounded-native-execution checks for eight-model comparison."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from rdkit import Chem

from nyl12_covalent import relay_broad_graph as graph
from nyl12_covalent import relay_broad_predictions as runner
from nyl12_covalent import relay_broad_predictor_runtime as runtime
from nyl12_covalent.relay_broad_prediction_io import audit_prediction


@pytest.fixture
def design(tmp_path: Path) -> Path:
    """Prepare a visibly synthetic one-mutation technical sequence fixture."""
    native = json.loads(runner.DEFAULT_PROTECTION.read_text())["source_sequences"]
    alpha = native["alpha"]
    index = 20
    alpha = alpha[:index] + ("A" if alpha[index] != "A" else "V") + alpha[index + 1 :]
    path = tmp_path / "technical_design.json"
    path.write_text(
        json.dumps(
            dict(
                candidate_id="TECHNICAL_FIXTURE_NOT_SELECTED",
                alpha_sequence=alpha,
                beta_sequence=native["beta"],
            )
        )
    )
    return path


def model_cif(path: Path, job: dict) -> None:
    """Write a complete heavyatom synthetic software fixture, without geometricutility."""
    columns = [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_comp_id",
        "auth_asym_id",
        "auth_seq_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "B_iso_or_equiv",
        "pdbx_PDB_model_num",
        "label_alt_id",
    ]
    lines = ["data_technical", "#", "loop_", *["_atom_site." + x for x in columns]]
    serial = 0
    for chain in "ABCDEFGH":
        seq = job["alpha_sequence"] if chain in "ABCD" else job["beta_sequence"]
        mol = Chem.MolFromSequence(seq)
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            name = info.GetName().strip()
            if name == "OXT":
                continue
            serial += 1
            lines.append(
                f"ATOM {serial} {atom.GetSymbol()} {name} {info.GetResidueName()} {chain} {info.GetResidueNumber()} 0 0 0 80 1 ."
            )
    for chain in "IJKL":
        for name, source in job["ligand_atom_maps"][chain].items():
            serial += 1
            lines.append(
                f"HETATM {serial} {source[0]} {name} LIG {chain} 1 0 0 0 80 1 ."
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([*lines, "#"]) + "\n")


def test_authentic_es_graph_and_native_name_mapping(tmp_path: Path) -> None:
    """FreshCCDretains79physicalatoms,threeamidecarbonylsandzwitterionicgraph."""
    mol = graph.source_molecule(graph.CONFIG / "whole_l4_source.mol2")
    table = json.loads(runner.DEFAULT_ATOMS.read_text())
    maps = graph.ligand_maps(mol, table)
    assert len(maps["af3"]) == 33 and len(maps["boltz"]) == 33
    assert set(maps["af3"]) == set(maps["boltz"].values())
    ccd = tmp_path / "fresh.cif"
    graph.write_es_ccd(mol, ccd)
    text = ccd.read_text()
    assert graph.COMPONENT in text and "TI1" not in text
    assert "_chem_comp_bond.value_order" in text
    assert sum(a.GetFormalCharge() == 1 for a in mol.GetAtoms()) == 1
    assert sum(a.GetFormalCharge() == -1 for a in mol.GetAtoms()) == 1
    with pytest.raises(FileExistsError):
        graph.write_es_ccd(mol, ccd)
    bad = dict(table, smiles="CC")
    with pytest.raises(ValueError):
        graph.ligand_maps(mol, bad)
    with pytest.raises(FileNotFoundError):
        graph.source_molecule(tmp_path / "missing.mol2")


def test_eight_unguided_inputs_and_idempotence(tmp_path: Path, design: Path) -> None:
    """Botharchitecturesusequery-onlyMSAandemptytemplates;onlyeightone-samplejobs."""
    out = tmp_path / "comparison"
    manifest = runner.prepare_comparison(out, design)
    assert len(manifest["jobs"]) == 8
    assert {(j["role"], j["predictor"], j["seed"]) for j in manifest["jobs"]} == {
        (r, p, s)
        for r in ["WT", "DESIGN"]
        for p in ["af3", "boltz"]
        for s in runner.SEEDS
    }
    for job in manifest["jobs"]:
        data = json.loads(Path(job["input"]).read_text())
        command = job["command"]
        assert job["samples"] == 1 and job["environment"]["CUDA_VISIBLE_DEVICES"] == "1"
        assert "--use_potentials" not in command
        assert "bondedAtomPairs" not in data
        if job["predictor"] == "af3":
            assert "--contact_guidance_weight=0" in command
            assert (
                "--num_diffusion_samples=1" in command and "--buckets=1536" in command
            )
            for entry in data["sequences"][:2]:
                p = entry["protein"]
                assert p["templates"] == [] and p["pairedMsa"] == ""
                assert p["unpairedMsa"].count(">") == 1
        else:
            assert data["templates"] == [] and data["constraints"] == []
            assert command[command.index("--num_workers") + 1] == "0"
            assert command[command.index("--diffusion_samples") + 1] == "1"
            for entry in data["sequences"][:2]:
                assert Path(entry["protein"]["msa"]).read_text().count(">") == 1
    assert runner.prepare_comparison(out, design) == manifest
    Path(manifest["jobs"][0]["input"]).write_text("{}")
    with pytest.raises(ValueError, match="inputchanged"):
        runner.prepare_comparison(out, design)


def test_sequential_execution_and_graph_only_acceptance(
    tmp_path: Path, design: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exactly8nativecommandsrunsequentially; geometryisnotusedasasoftwaregate."""
    out = tmp_path / "comparison"
    manifest = runner.prepare_comparison(out, design)
    calls = []

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        if command in [j["command"] for j in manifest["jobs"]]:
            job = next(j for j in manifest["jobs"] if j["command"] == command)
            calls.append(job["job_id"])
            model_cif(Path(job["expected_cif"]), job)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(runtime.subprocess, "run", fake_run)
    monkeypatch.setattr(
        runtime, "export_geometry", lambda *args: dict(geometry_json="technicalfixture")
    )
    result = runner.run_comparison(manifest, out)
    assert result["model_count"] == 8 and result["software_acceptance"]
    assert len(calls) == 8
    assert all(m["finite_coordinates"] for m in result["models"])
    runner.run_comparison(manifest, out)
    assert len(calls) == 8
    first = manifest["jobs"][0]
    path = Path(first["expected_cif"])
    path.write_text(path.read_text().replace(" 0 0 0 80 ", " nan 0 0 80 ", 1))
    with pytest.raises(ValueError, match="Nonfinite"):
        audit_prediction(
            path,
            first["alpha_sequence"],
            first["beta_sequence"],
            first["ligand_atom_maps"]["I"],
        )
    with pytest.raises(ValueError, match="changed"):
        runner.run_comparison(manifest, out)


def test_invalid_inputs_fail_explicitly(
    tmp_path: Path, design: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Rejectunsupportedstatesandstopafterthenativefirstfailurewithoutfallback."""
    with pytest.raises(ValueError):
        runner._sequences(dict(alpha_sequence="A", beta_sequence="T"))
    with pytest.raises(ValueError):
        runner._sequences(dict(alpha_sequence="A" * 232, beta_sequence="A" * 96))
    with pytest.raises(ValueError):
        runner.commands(dict(predictor="unknown"))
    out = tmp_path / "comparison"
    manifest = runner.prepare_comparison(out, design)

    def fail(command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        return subprocess.CompletedProcess(
            command, 7 if command == manifest["jobs"][0]["command"] else 0
        )

    monkeypatch.setattr(runtime.subprocess, "run", fail)
    with pytest.raises(RuntimeError, match="Nativepredictionfailed"):
        runner.run_comparison(manifest, out)
    with pytest.raises(ValueError):
        runner.run_comparison(dict(jobs=[]), out)


def test_geometry_export_uses_explicit_atom_names(
    tmp_path: Path, design: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Export48rowswithoutguessingbondsorforcingWTcontactgeometry."""
    import csv
    from nyl12_covalent import relay_any_amide
    from nyl12_covalent.relay_broad_prediction_io import export_geometry

    out = tmp_path / "comparison"
    manifest = runner.prepare_comparison(out, design)
    job = manifest["jobs"][0]
    cif = Path(job["expected_cif"])
    model_cif(cif, job)
    calls = []

    def evaluate(atoms: dict, maps: dict, graph_value: object) -> dict:
        calls.append((atoms, maps))
        rows = [
            dict(ligand=chain, site=site, carbonyl_c=f"C{amide}", tier="reject")
            for chain in "IJKL"
            for site in "ABCD"
            for amide in range(3)
        ]
        return dict(
            all_amide_site_rows=rows,
            assigned=[rows[i] for i in [0, 12, 24, 36]],
            assigned_tier_ab_count=0,
        )

    monkeypatch.setattr(relay_any_amide, "evaluate_model", evaluate)
    result = export_geometry(cif, job["ligand_atom_maps"], cif.parent)
    with Path(result["geometry_csv"]["all_amide_site_rows"]).open() as handle:
        assert len(list(csv.DictReader(handle))) == 48
    assert len(calls) == 1 and ("E", 1, "OG1") in calls[0][0]
    assert calls[0][1] == job["ligand_atom_maps"]
    assert json.loads(Path(result["geometry_json"]).read_text())[
        "software_acceptance_separate"
    ]
