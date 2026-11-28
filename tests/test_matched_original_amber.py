"""Verify matched Amber control flow without invoking scientific engines."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from nyl12_covalent import matched_original_amber as target


def save(path: Path, data: Any) -> Path:
    path.write_text(json.dumps(data))
    return path


@pytest.fixture
def plan(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    cfg = {}
    for key in ("engine", "tleap", "mapping", "mol2", "amberhome"):
        path = tmp_path / key
        path.write_text(key)
        cfg[key] = str(path)
    release = tmp_path / "parameters"
    release.mkdir()
    (release / "release.json").write_text("release")
    cfg.update(release=str(release), gpu="GPU-reserved")
    adapter = tmp_path / "adapter.py"
    adapter.write_text(
        "def reference_atoms(): return {}\ndef pairing(a,b): return ({}, {})\ndef canonical_cif(a,b,c):\n b.write_text('cif')\n return {}\n"
    )
    processor = tmp_path / "processor.py"
    processor.write_text("processor")
    reference = save(tmp_path / "approved.json", {"target_phi_sign": 1})
    graph = save(tmp_path / "graph.json", {"ligand_atom_maps": {"boltz": {"I": {}}}})
    cases = []
    for candidate in ("relay_prot101_06", "relay_prot101_03", "WT"):
        for seed in (20261126, 20261127):
            case = f"{candidate}_seed_{seed}"
            inp = {
                "case": case,
                "candidate": candidate,
                "seed": seed,
                "sequences": {"alpha": "A", "beta": "T"},
            }
            case_input = save(tmp_path / f"{case}.json", inp)
            raw = tmp_path / f"{case}.cif"
            raw.write_text("original")
            receipt = {
                "case": case,
                "candidate": candidate,
                "seed": seed,
                "status": "native_prediction_completed",
                "exit_status": 0,
                "prediction_cif": str(raw),
                "prediction_cif_sha256": target.sha256(raw),
                "approved_reference_manifest": str(reference),
                "target_phi_sign": 1,
                "ligand_graph_manifest": str(graph),
                "input_sha256": {str(case_input): target.sha256(case_input)},
            }
            native = save(tmp_path / f"{case}_native.json", receipt)
            cases.append(
                {
                    **{k: inp[k] for k in ("case", "candidate", "seed")},
                    "case_input": str(case_input),
                    "native_receipt": str(native),
                }
            )
    pins = {str(p): target.sha256(p) for p in tmp_path.rglob("*") if p.is_file()}
    monkeypatch.setattr(target, "ENGINE", Path(cfg["engine"]))
    monkeypatch.setattr(target, "ENGINE_SHA256", target.sha256(Path(cfg["engine"])))
    monkeypatch.setattr(target.socket, "gethostname", lambda: "violet")
    monkeypatch.setattr(target, "_journal", lambda *_: None)
    return {
        "schema_version": 1,
        "target_phi_sign": 1,
        "cases": cases,
        "relaxation_config": cfg,
        "canonical_adapter": str(adapter),
        "legacy_processor": str(processor),
        "approved_reference_manifest": str(reference),
        "input_sha256": pins,
    }


def test_exact_six_case_protocol(plan: dict[str, Any]) -> None:
    target.validate_plan(plan)
    assert {c["candidate"] for c in plan["cases"]} == {
        "relay_prot101_06",
        "relay_prot101_03",
        "WT",
    }


@pytest.mark.parametrize(
    "fault",
    [
        "sign",
        "bool",
        "schema",
        "cases",
        "name",
        "engine",
        "pin",
        "changed",
        "reference",
    ],
)
def test_invalid_plan_fails(plan: dict[str, Any], fault: str) -> None:
    if fault == "sign":
        plan["target_phi_sign"] = -1
    elif fault == "bool":
        plan["target_phi_sign"] = True
    elif fault == "schema":
        plan["schema_version"] = 2
    elif fault == "cases":
        plan["cases"].pop()
    elif fault == "name":
        plan["cases"][0]["case"] = "../bad"
    elif fault == "engine":
        plan["relaxation_config"]["engine"] = "/tmp/wrong"
    elif fault == "pin":
        plan["input_sha256"].pop(plan["canonical_adapter"])
    elif fault == "changed":
        Path(plan["legacy_processor"]).write_text("changed")
    else:
        p = Path(plan["approved_reference_manifest"])
        save(p, {"target_phi_sign": -1})
        plan["input_sha256"][str(p)] = target.sha256(p)
    with pytest.raises(ValueError):
        target.validate_plan(plan)


def test_interchain_clashes_are_exact_and_bounded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    atoms = {("A", 1, str(i)): ("C", (0, 0, 0)) for i in range(22)}
    atoms.update(
        {
            ("B", 1, "C"): ("C", (0.5, 0, 0)),
            ("C", 1, "H"): ("H", (0, 0, 0)),
            ("D", 1, "C"): ("C", (1.5, 0, 0)),
        }
    )
    monkeypatch.setattr(target, "pdb_atoms", lambda _: atoms)
    report = target.severe_clashes(tmp_path / "source.pdb")
    assert report["interchain_heavy_pairs_below_threshold"] == 22
    assert len(report["examples"]) == 20
    assert report["accepted"] is False
    atoms.clear()
    with pytest.raises(ValueError, match="Empty"):
        target.severe_clashes(tmp_path)
    atoms[("A", 1, "C")] = ("C", (float("nan"), 0, 0))
    with pytest.raises(ValueError, match="Non-finite"):
        target.severe_clashes(tmp_path)


def native_stubs(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(target, "audit_prediction", lambda *_: {"accepted": True})
    monkeypatch.setattr(target, "read_atom_sites", lambda _: [])
    monkeypatch.setattr(target, "export_control", lambda *args: {"all48": []})
    monkeypatch.setattr(target, "annotate_amide_registers", lambda rows, **_: rows)
    monkeypatch.setattr(
        target, "evaluate_signed_amide_gate", lambda rows, **_: {"input_row_count": 48}
    )


def test_native_preservation_and_roundtrip(
    plan: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    native_stubs(monkeypatch)
    out = tmp_path / "native"
    out.mkdir()
    source = target._canonicalize(plan["cases"][0], plan, out)
    assert source == out / "canonical/canonical_mature_local.pdb"
    assert (out / "original_native_prediction.cif").read_text() == "original"
    assert (
        json.loads((out / "native_all48_positive_gate.json").read_text())[
            "input_row_count"
        ]
        == 48
    )


@pytest.mark.parametrize(
    "fault", ["identity", "status", "hash", "reference", "input", "graph", "case_input"]
)
def test_native_receipt_rejection(
    plan: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fault: str
) -> None:
    native_stubs(monkeypatch)
    case = plan["cases"][0]
    path = Path(case["native_receipt"])
    receipt = json.loads(path.read_text())
    if fault == "identity":
        receipt["seed"] = 7
    elif fault == "status":
        receipt["exit_status"] = 1
    elif fault == "hash":
        receipt["prediction_cif_sha256"] = "wrong"
    elif fault == "reference":
        receipt["approved_reference_manifest"] = "/tmp/wrong"
    elif fault == "input":
        receipt["input_sha256"][case["case_input"]] = "wrong"
    elif fault == "graph":
        plan["input_sha256"].pop(receipt["ligand_graph_manifest"])
    else:
        case["case_input"] = str(save(tmp_path / "wrong.json", {"case": "wrong"}))
    save(path, receipt)
    out = tmp_path / "native"
    out.mkdir()
    with pytest.raises(ValueError):
        target._canonicalize(case, plan, out)


def test_run_case_rejects_assembly_without_amber(
    plan: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(target, "_canonicalize", lambda *_: tmp_path / "preserved.pdb")
    monkeypatch.setattr(target, "severe_clashes", lambda _: {"accepted": False})
    monkeypatch.setattr(
        target, "relax_model", lambda *_: pytest.fail("Must not launch Amber")
    )
    report = target.run_case(plan["cases"][0], plan, tmp_path / "case")
    assert report["status"] == "failed_assembly_precheck"
    assert not report["amber_launched"]


def test_run_case_canonicalizes_real_native_paths(
    plan: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(target, "_canonicalize", lambda *_: tmp_path / "source.pdb")
    monkeypatch.setattr(target, "severe_clashes", lambda _: {"accepted": True})
    monkeypatch.setattr(
        target,
        "relax_model",
        lambda *_: {"final_rms_gradient": 0.047399, "force_converged": True},
    )
    calls = []

    def process(**kwargs: Any) -> dict[str, Any]:
        calls.append(kwargs)
        return {"snapshots": {s: {"assigned_ab_count": 4} for s in ("before", "after")}}

    monkeypatch.setattr(target, "process_design_case", process)
    out = tmp_path / "accepted"
    result = target.run_case(plan["cases"][0], plan, out)
    assert result["positive_ab_after"] == 4
    assert (
        calls[0]["topology"] == out / "relaxation/minimization/endpoint_complex.prmtop"
    )
    assert calls[0]["before_restart"] == out / "relaxation/system/dry.rst7"
    assert calls[0]["target_phi_sign"] == 1


def test_campaign_failure_census_continues_other_cases(
    plan: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    first = plan["cases"][0]["case"]

    def run(case: dict[str, Any], *_: Any) -> dict[str, Any]:
        if case["case"] == first:
            raise RuntimeError("native assembly evidence retained")
        return {"case": case["case"], "status": "completed"}

    monkeypatch.setattr(target, "run_case", run)
    out = tmp_path / "campaign"
    result = target.run_campaign(plan, out)
    assert result["failure_count"] == 1
    assert len(result["cases"]) == 6
    assert json.loads((out / "failure_census.json").read_text()) == result
    with pytest.raises(FileExistsError):
        target.run_campaign(plan, out)


def test_single_case_worker_and_host_guard(
    plan: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        target,
        "run_case",
        lambda case, *_: {"case": case["case"], "status": "completed"},
    )
    name = plan["cases"][0]["case"]
    result = target.run_campaign(plan, tmp_path / "worker", case_name=name)
    assert list(result["cases"]) == [name]
    with pytest.raises(ValueError, match="Unknown"):
        target.run_campaign(plan, tmp_path, case_name="unknown")
    monkeypatch.setattr(target.socket, "gethostname", lambda: "azalea")
    with pytest.raises(ValueError, match="Violet"):
        target.run_campaign(plan, tmp_path)
