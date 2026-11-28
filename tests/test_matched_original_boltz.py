"""Test bounded matched search/staging without native GPU dependencies."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from nyl12_covalent import matched_original_boltz as m


def plan(tmp_path: Path) -> dict:
    """Build exact 328-residue synthetic queries and explicit database inputs."""
    db = tmp_path / "db"
    db.mkdir()
    for name in m.DATABASE_NAMES:
        (db / name).write_text(name)
    executable = tmp_path / "colabfold_search"
    executable.write_text("native")
    return {
        "target_phi_sign": 1,
        "search": {"database": str(db), "executable": str(executable)},
        "candidates": [
            {
                "candidate": "design06",
                "sequences": {"alpha": "A" * 232, "beta": "G" * 96},
                "gpu_uuid": "GPU-test",
            }
        ],
        "seeds": [20261126, 20261127],
    }


def test_search_native_header_and_exact_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One query uses its native header output and never query-swaps a cache."""
    p = plan(tmp_path)
    monkeypatch.setattr(m, "journal", lambda *args: None)
    monkeypatch.setattr(m, "require_idle", lambda gpu: None)
    calls = []

    def run(argv: list, **kwargs: object) -> SimpleNamespace:
        calls.append(argv)
        output = Path(argv[-1])
        output.mkdir()
        (output / "design06.a3m").write_text(
            ">design06\n"
            + "A" * 232
            + "G" * 96
            + "\n>homolog\n"
            + "C" * 232
            + "T" * 96
            + "\n"
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(m.subprocess, "run", run)
    source = m.search(p, tmp_path / "out", "design06")
    assert source.name == "design06.a3m"
    assert m.search(p, tmp_path / "out", "design06") == source
    assert len(calls) == 1
    p["candidates"][0]["sequences"]["alpha"] = "C" * 232
    with pytest.raises(ValueError, match="cache identity"):
        m.search(p, tmp_path / "out", "design06")


def test_search_failure_retains_process(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Nonzero native status cannot publish scientific acceptance."""
    p = plan(tmp_path)
    monkeypatch.setattr(m, "journal", lambda *args: None)
    monkeypatch.setattr(m, "require_idle", lambda gpu: None)
    monkeypatch.setattr(
        m.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=9)
    )
    with pytest.raises(RuntimeError, match="Search failed"):
        m.search(p, tmp_path / "out", "design06")
    work = tmp_path / "out/searches/design06"
    assert json.loads((work / "process.json").read_text())["exit_status"] == 9
    assert not (work / "source_acceptance.json").exists()


def test_prepare_full_mature_tetramer_without_guidance(tmp_path: Path) -> None:
    """Staging keeps eight mature chains and four unguided ligand copies."""
    p = plan(tmp_path)
    source = tmp_path / "source.a3m"
    source.write_text(
        ">q\n" + "A" * 232 + "G" * 96 + "\n>h\n" + "C" * 232 + "T" * 96 + "\n"
    )
    accepted = tmp_path / "accepted.json"
    accepted.write_text(
        json.dumps({"status": "accepted", "precursor_a3m_sha256": m.sha(source)})
    )
    p["candidates"][0].update(
        precursor_a3m=str(source), source_acceptance=str(accepted)
    )
    graph = tmp_path / "graph.json"
    graph.write_text(
        json.dumps(
            {
                "heavy_atoms": 33,
                "formal_charge": 0,
                "covalent_ES_bonds": [],
                "smiles": "PA6",
            }
        )
    )
    reference = tmp_path / "reference.json"
    reference.write_text("{}")
    p.update(
        ligand_graph_manifest=str(graph), approved_reference_manifest=str(reference)
    )
    cases = m.prepare(p, tmp_path / "out", "design06")
    assert len(cases) == 2
    case = json.loads((cases[0] / "case.json").read_text())
    data = json.loads(Path(case["input"]).read_text())
    assert set(data) == {"version", "sequences"}
    assert data["sequences"][0]["protein"]["id"] == list("ABCD")
    assert data["sequences"][1]["protein"]["id"] == list("EFGH")
    assert data["sequences"][2]["ligand"]["id"] == list("IJKL")


@pytest.mark.parametrize("sign", [-1, 0])
def test_reject_wrong_target(tmp_path: Path, sign: int) -> None:
    """Historical negative-only configuration cannot advance."""
    p = plan(tmp_path)
    p["target_phi_sign"] = sign
    with pytest.raises(ValueError, match="target_phi_sign"):
        m.candidate(p, "design06")


def test_native_prediction_completion_receipt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A native exit0 publishes a hash-bound CIF receipt for independent Amber work."""
    from nyl12_covalent import matched_original_native as n

    p = plan(tmp_path)
    executable, checkpoint = tmp_path / "boltz", tmp_path / "checkpoint"
    executable.write_text("native")
    checkpoint.write_text("weights")
    p["engine"] = {
        "executable": str(executable),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": m.sha(checkpoint),
        "cache": str(tmp_path),
        "env_overrides": {},
    }
    p.update(
        ligand_graph_manifest="graph.json", approved_reference_manifest="approved.json"
    )
    work = tmp_path / "out/cases/design06_seed_20261126"
    work.mkdir(parents=True)
    source, featurizer = work / "input.yaml", work / "featurizer.py"
    source.write_text("{}")
    featurizer.write_text("native")
    msas = work / "msas"
    msas.mkdir()
    for name in [
        "alpha_paired.csv",
        "beta_paired.csv",
        "precursor_pairing_receipt.json",
    ]:
        (msas / name).write_text("binding")
    case = {
        "case": work.name,
        "input": str(source),
        "msas": str(msas),
        "remote_output": str(work / "runs"),
    }
    (work / "case.json").write_text(json.dumps(case))
    record = work / f"runs/boltz_results_{work.name}/processed/records/{work.name}.json"
    record.parent.mkdir(parents=True)
    record.write_text("native record")
    case["source_sha256"] = {str(source): m.sha(source)}
    (work / "case.json").write_text(json.dumps(case))
    acceptance = {
        "record_sha256": m.sha(record),
        "paired_source_receipt_sha256": m.sha(msas / "precursor_pairing_receipt.json"),
        "input_sha256": m.sha(source),
        "native_featurizer_source_path": str(featurizer),
        "native_featurizer_source_sha256": m.sha(featurizer),
    }
    (work / "cpu_acceptance.json").write_text(json.dumps(acceptance))
    monkeypatch.setattr(n, "audit", lambda *args: acceptance)
    monkeypatch.setattr(n, "journal", lambda *args: None)
    monkeypatch.setattr(n, "require_idle", lambda *args: None)

    def run(argv: list, **kwargs: object) -> SimpleNamespace:
        cif = (
            work
            / f"runs/boltz_results_{work.name}/predictions/{work.name}/{work.name}_model_0.cif"
        )
        cif.parent.mkdir(parents=True, exist_ok=True)
        cif.write_text("data_native")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(n.subprocess, "run", run)
    n.predict_candidate(p, tmp_path / "out", "design06", seed=20261126)
    receipt = json.loads((work / "native_prediction.json").read_text())
    assert receipt["status"] == "native_prediction_completed"
    assert receipt["prediction_cif_sha256"] == m.sha(Path(receipt["prediction_cif"]))
    assert receipt["scientific_acceptance"] == "Not evaluated; process completion only"
    assert "--no_kernels" in receipt["command"]


@pytest.mark.parametrize("exit_status", [0, 4])
def test_native_audit_process_boundary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, exit_status: int
) -> None:
    """Native audit requires both process success and an explicit accepted receipt."""
    from nyl12_covalent import matched_original_native as n

    p = {
        "engine": {
            "python": "/explicit/native/python",
            "cache": "/cache",
            "env_overrides": {},
        }
    }
    (tmp_path / "case.json").write_text("{}")
    (tmp_path / "cpu_acceptance.json").write_text(
        json.dumps({"status": "accepted", "stock_native_entry_only": True})
    )
    commands = []

    def run(argv: list, **kwargs: object) -> SimpleNamespace:
        commands.append((argv, kwargs))
        return SimpleNamespace(returncode=exit_status)

    monkeypatch.setattr(n.subprocess, "run", run)
    if exit_status:
        with pytest.raises(RuntimeError, match="CPU audit failed"):
            n.audit(p, tmp_path)
    else:
        assert n.audit(p, tmp_path)["status"] == "accepted"
    assert commands[0][0][0] == "/explicit/native/python"
    assert commands[0][1]["env"]["CUDA_VISIBLE_DEVICES"] == ""


def test_native_rejects_seed_outside_matched_pair(tmp_path: Path) -> None:
    """Single-case queue scheduling cannot create an extra prediction seed."""
    from nyl12_covalent import matched_original_native as n

    p = plan(tmp_path)
    executable = tmp_path / "boltz"
    executable.write_text("stock")
    checkpoint = tmp_path / "weights"
    checkpoint.write_text("weights")
    p["engine"] = {
        "executable": str(executable),
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": m.sha(checkpoint),
    }
    with pytest.raises(ValueError, match="Selected seed"):
        n.predict_candidate(p, tmp_path / "out", "design06", seed=42)
