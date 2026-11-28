"""Focused strict-contract tests for the bounded broad continuation."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import relay_broad_runtime as runtime


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, Path("scripts") / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_immutable_write_hash_export(tmp_path):
    path = tmp_path / "data.json"
    runtime.write_new(path, {"accepted": True})
    assert len(runtime.digest(path)) == 64
    with pytest.raises(FileExistsError):
        runtime.write_new(path, {})
    source = tmp_path / "source.pdb"
    source.write_text("ATOM      1  CA  THR E 233    unchanged\n")
    runtime.export_local(source, tmp_path / "local.pdb")
    assert "E   1" in (tmp_path / "local.pdb").read_text()
    runtime.export_local(tmp_path / "local.pdb", tmp_path / "again.pdb")
    assert (tmp_path / "again.pdb").read_bytes() == (
        tmp_path / "local.pdb"
    ).read_bytes()
    source.write_text("ATOM      1  CA  THR E 999    unchanged\n")
    with pytest.raises(ValueError, match="numbering"):
        runtime.export_local(source, tmp_path / "invalid.pdb")


def make_parent(tmp_path):
    parent = {
        "sequence_scope": "native_WT",
        "parent_id": "parent01_t2_seed01",
        "parent_source_id": "medoid01",
        "noise_A": 2,
    }
    for field in (
        "pdb_path",
        "canonical_mapping",
        "pa6_graph_receipt",
        "scientific_acceptance",
    ):
        path = tmp_path / field
        path.write_text(json.dumps({"accepted": True}))
        parent[field] = str(path)
        parent[field + "_sha256"] = runtime.digest(path)
    return parent


def test_parent_receipts_selection_and_hash_rejection(tmp_path):
    parent = make_parent(tmp_path)
    runtime.verify_parent(parent)
    pilot = load_script("run_relay_broad_pilot")
    assert len(pilot.selected_remodels([{"parents": [parent, parent]}])) == 1
    parent["sequence_scope"] = "designed"
    with pytest.raises(ValueError, match="native-sequence"):
        runtime.verify_parent(parent)
    parent["sequence_scope"] = "native_WT"
    Path(parent["scientific_acceptance"]).write_text('{"accepted":false}')
    parent["scientific_acceptance_sha256"] = runtime.digest(
        Path(parent["scientific_acceptance"])
    )
    with pytest.raises(ValueError, match="rejected"):
        runtime.verify_parent(parent)
    parent["pdb_path_sha256"] = "bad"
    with pytest.raises(ValueError, match="hash mismatch"):
        runtime.verify_parent(parent)


def test_native_process_receipt_journal_and_nonzero(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime, "journal", lambda *args: None)
    code = runtime.native_run(
        tmp_path,
        [sys.executable, "-c", "raise SystemExit(4)"],
        tmp_path / "attempt",
        {"OMP_NUM_THREADS": "2"},
    )
    assert code == 4
    assert json.loads((tmp_path / "attempt.process.json").read_text())["exit_code"] == 4
    monkeypatch.setattr(
        runtime.subprocess, "run", lambda *args, **kwargs: SimpleNamespace()
    )
    runtime.journal(tmp_path, "intent", "test", "reason")


def test_preselection_local_scores_nonwt_and_tiebreaks():
    def candidate(name, mutations, score):
        return {
            "candidate_id": name,
            "mutation_count": mutations,
            "scores": [{"state": "TS", "local_AE_design_nll": score}],
        }

    values = [candidate("WT", 0, -100), candidate("a", 10, 1), candidate("b", 2, 1)]
    assert [x["candidate_id"] for x in runtime.rank_broad(values)] == ["b", "a"]
    with pytest.raises(ValueError, match="non-WT"):
        runtime.rank_broad([candidate("WT", 0, 1)])
    with pytest.raises(ValueError, match="finite"):
        runtime.rank_broad([candidate("bad", 1, float("nan"))])
    assert runtime.rank_broad([candidate("ok", 1, 3)])[0]["local_ts_worst_nll"] == 3


def test_parser_audit_actual_contract(monkeypatch, tmp_path):
    names = [(c, p) for c in "ABCDEFGH" for p in range(1, 233 if c in "ABCD" else 97)]
    protein = {
        "chain_letters": [x[0] for x in names],
        "R_idx": np.array([x[1] for x in names]),
        "mask": np.ones(1312),
        "Y_m": np.ones(33),
        "xyz_37_m": np.ones((1312, 37)),
    }
    monkeypatch.setitem(
        sys.modules, "data_utils", SimpleNamespace(parse_PDB=lambda *a, **k: [protein])
    )
    pdb = tmp_path / "input.pdb"
    pdb.write_text("test")
    assert runtime.parser_audit(pdb, tmp_path / "run.py")["accepted"]
    protein["Y_m"] = np.ones(32)
    with pytest.raises(ValueError, match="33 ligand"):
        runtime.parser_audit(pdb, tmp_path / "run.py")
    protein["Y_m"] = np.ones(33)
    protein["xyz_37_m"][112, 0] = 0
    with pytest.raises(ValueError, match="fixed sidechain"):
        runtime.parser_audit(pdb, tmp_path / "run.py")
    protein["mask"][0] = 0
    with pytest.raises(ValueError, match="canonical"):
        runtime.parser_audit(pdb, tmp_path / "run.py")


def test_gate_wait_accept_failure_timeout_and_hash(tmp_path, monkeypatch):
    controller = load_script("run_relay_broad_controller")
    parent = make_parent(tmp_path)
    gate, terminal = tmp_path / "gate.json", tmp_path / "terminal.json"
    gate.write_text(json.dumps({"accepted": True, "parents": [parent, parent]}))
    assert controller.wait_for_gate(gate, terminal, float("inf"))["accepted"]
    gate.write_text('{"accepted":true,"parents":[]}')
    with pytest.raises(ValueError, match="exactly two"):
        controller.wait_for_gate(gate, terminal, float("inf"))
    gate.write_text('{"accepted":false}')
    terminal.write_text('{"terminal":true,"accepted":false}')
    with pytest.raises(RuntimeError, match="Final qualification"):
        controller.wait_for_gate(gate, terminal, float("inf"))
    with pytest.raises(TimeoutError):
        controller.wait_for_gate(gate, terminal, -1)
    controller.check_hashes({str(gate): runtime.digest(gate)})
    with pytest.raises(ValueError, match="changed"):
        controller.check_hashes({str(gate): "bad"})


def test_distinct_global_local_summaries(monkeypatch):
    monkeypatch.setattr(
        runtime,
        "summarize_score",
        lambda stats, seq, mask: {"AE_design_nll": mask["n"]},
    )
    assert runtime.score_summary({}, "", {"n": 2}, {"n": 3}) == {
        "AE_design_nll": 2,
        "local_AE_design_nll": 3,
    }


def test_full_bounded_driver_mocked_native_engines(tmp_path, monkeypatch):
    pilot = load_script("run_relay_broad_pilot")
    native_sequences = {c: "A" for c in "ABCDEFGH"}
    templates = [
        {"id": "p1", "pdb": "pdb", "state": "TS", "raw_samples": 6},
        {"id": "p2", "pdb": "pdb", "state": "TS", "raw_samples": 6},
        {"id": "ES", "pdb": "pdb", "state": "ES"},
    ]
    monkeypatch.setattr(pilot, "prepare_templates", lambda *a: templates)
    monkeypatch.setattr(pilot, "journal", lambda *a: None)
    monkeypatch.setattr(pilot, "native_run", lambda *a: 0)
    monkeypatch.setattr(pilot, "_pdb_sequences", lambda *a: native_sequences)
    monkeypatch.setattr(
        pilot,
        "_fasta",
        lambda *a: [
            ("native", ":".join(native_sequences.values())),
            *[("sample", "mutant")] * 6,
        ],
    )
    monkeypatch.setattr(
        pilot,
        "audit_sequence",
        lambda seq, *a: {
            "sequence": seq,
            "mutation_count": 1,
            "alpha_sequence": "A",
            "beta_sequence": "A",
        },
    )
    monkeypatch.setattr(pilot, "score_summary", lambda *a: {"local_AE_design_nll": 2})
    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(load=lambda *a, **k: {"decoding_order": np.zeros((2, 1312))}),
    )
    config = {
        "source_sha256": {},
        "request": {
            "generation": {
                "temperature": 0.1,
                "native_identity_bias": 0,
                "max_raw_samples": 24,
                "mask": {},
                "redesigned_residues": "A2",
                "symmetry_residues": "A2,B2,C2,D2",
                "symmetry_weights": "1,0,0,0",
            }
        },
        "native": {
            "environment": {},
            "python": "python",
            "generation_engine": "run.py",
            "score_engine": "score.py",
            "score_wrapper": "wrapper.py",
            "checkpoint": "weights",
        },
        "local_pocket_mask": {},
    }
    output = tmp_path / "run"
    pilot.run(config, output)
    receipt = json.loads((output / "scientific_acceptance.json").read_text())
    assert receipt["raw_count"] == 12
    assert receipt["unique_count"] == 1
    assert (output / "preselected_candidate.json").exists()
    config["request"]["generation"]["temperature"] = 0.2
    with pytest.raises(ValueError, match="Unauthorized"):
        pilot.run(config, tmp_path / "bad")


def test_controller_expansion_to_predictions_mocked(tmp_path, monkeypatch):
    controller = load_script("run_relay_broad_controller")
    output = tmp_path / "controller"
    output.mkdir()
    runner = tmp_path / "rfd3.py"
    runner.write_text("source")
    expansion = tmp_path / "expansion_manifest.json"
    expansion.write_text(json.dumps({"jobs": [{}] * 6, "runner": str(runner)}))
    release = tmp_path / "release.json"
    release.write_text(
        json.dumps(
            {
                "execution_release_checked": True,
                "qualification_gate_verified": True,
                "output_scientific_acceptance": False,
                "manifest": str(expansion),
                "manifest_sha256": runtime.digest(expansion),
                "runner_sha256": runtime.digest(runner),
            }
        )
    )
    expanded = tmp_path / "expanded.json"
    expanded.write_text('{"status":"scientific_failed","parents":[]}')
    native = tmp_path / "native.json"
    native.write_text('{"python":"python","environment":{}}')
    request = tmp_path / "request.json"
    request.write_text("{}")
    monkeypatch.setattr(
        controller, "wait_for_gate", lambda *a: {"accepted": True, "parents": []}
    )
    monkeypatch.setattr(controller, "journal", lambda *a: None)
    commands = []

    def fake_run(root, argv, out, env):
        commands.append(argv)
        if out.name == "broad_pilot":
            generation = output / "generation_run01"
            generation.mkdir()
            (generation / "scientific_acceptance.json").write_text('{"accepted":true}')
        return 1 if out.name == "expansion" else 0

    monkeypatch.setattr(controller, "native_run", fake_run)
    (tmp_path / "gate.json").write_text("{}")
    config = {
        "output": str(output),
        "source_sha256": {},
        "gate": str(tmp_path / "gate.json"),
        "terminal": "terminal",
        "max_wait_seconds": 1,
        "expansion_release": str(release),
        "rfd3_python": "python",
        "rfd3_environment": {},
        "release_only_argv": ["release.py", "--release-only"],
        "expanded_parents": str(expanded),
        "native_manifest": str(native),
        "request": str(request),
        "local_pocket_mask": {},
        "broad_runner": "broad.py",
        "prediction_argv": ["predict.py", "{design_json}"],
        "prediction_environment": {},
    }
    controller.execute(config)
    assert len(commands) == 4
    assert commands[-1][-1].endswith("preselected_candidate.json")
    assert (output / "continuation_scientific_acceptance.json").exists()


def test_prepare_templates_exports_originals_and_es(tmp_path, monkeypatch):
    pilot = load_script("run_relay_broad_pilot")
    source = tmp_path / "parent.pdb"
    source.write_text("END\n")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(pilot, "parser_audit", lambda *a: {"accepted": True})
    manifest = {
        "request": {
            "original_parent_ids": ["p1", "p2"],
            "original_parent_pdbs": {"p1": str(source), "p2": str(source)},
        },
        "accepted_manifests": [{"parents": []}],
        "native": {
            "generation_engine": "run.py",
            "ES_template": {"id": "ES", "pdb": str(source)},
        },
    }
    templates = pilot.prepare_templates(manifest, out)
    assert len(templates) == 3
    assert templates[0]["raw_samples"] == 6
    assert templates[-1]["state"] == "ES"


@pytest.mark.parametrize("failed", [False, True])
def test_controller_main_terminal_receipts(tmp_path, monkeypatch, failed):
    controller = load_script("run_relay_broad_controller")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"output": str(tmp_path)}))
    monkeypatch.setattr(sys, "argv", ["controller", str(config)])
    monkeypatch.setattr(controller, "journal", lambda *a: None)

    def execute(configuration):
        if failed:
            raise RuntimeError("scientific gate rejected")

    monkeypatch.setattr(controller, "execute", execute)
    if failed:
        with pytest.raises(RuntimeError, match="scientific gate"):
            controller.main()
    else:
        controller.main()
    terminal = json.loads((tmp_path / "continuation_terminal.json").read_text())
    assert terminal["terminal"] is True
    assert terminal["accepted"] is not failed
