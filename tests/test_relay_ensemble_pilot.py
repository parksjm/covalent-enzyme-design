"""Exercise complete-sequence auditing and bounded orchestration without GPUs."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent.relay_ensemble_sequences import (
    audit_sequence,
    pool_samples,
    rank_candidates,
    summarize_score,
)

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str) -> object:
    """Import an existing entry point for isolated orchestration tests."""
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def inputs() -> tuple[dict, dict, str]:
    """Return a canonical native assembly and single-site mutation mask."""
    native = {**dict.fromkeys("ABCD", "A" * 232), **dict.fromkeys("EFGH", "A" * 96)}
    chains = {c: "C" + seq[1:] if c in "ABCD" else seq for c, seq in native.items()}
    return native, {"alpha": [1], "beta": []}, ":".join(chains.values())


def score_stats(sequence: str) -> dict:
    """Create known finite per-residue log probabilities for independent gathering."""
    alphabet = list("ACDEFGHIKLMNPQRSTVWYX")
    seq = sequence.replace(":", "")
    return {
        "residue_names": dict(
            enumerate(
                f"{c}{p}"
                for c in "ABCDEFGH"
                for p in range(1, (232 if c in "ABCD" else 96) + 1)
            )
        ),
        "alphabet": alphabet,
        "native_sequence": np.array([alphabet.index(a) for a in seq]),
        "log_probs": np.full((2, 1312, 21), -2.0),
        "decoding_order": np.array([np.arange(1312), np.arange(1312)[::-1]]),
        "mask": np.ones(1312),
    }


def test_audit_pool_and_rank(inputs: tuple) -> None:
    """Preserve canonical mutations/lineages and independently assess Pareto order."""
    native, mask, sequence = inputs
    sample = audit_sequence(sequence, native, mask)
    assert sample["mutation_count"] == 1
    assert sample["mutations"][0]["merged_position"] == 1
    pooled = pool_samples([{**sample, "sample_id": "a"}, {**sample, "sample_id": "b"}])
    assert pooled[0]["lineage"] == ["a", "b"]
    assert pool_samples([]) == []
    with pytest.raises(ValueError, match="budget"):
        pool_samples([{"sequence": str(i), "sample_id": str(i)} for i in range(25)])
    broken = sequence.split(":")
    broken[1] = "A" * 232
    with pytest.raises(ValueError, match="ties"):
        audit_sequence(":".join(broken), native, mask)
    with pytest.raises(ValueError, match="untied/fixed"):
        audit_sequence(sequence, native, {"alpha": [], "beta": []})
    protected = {
        c: s[:112] + "C" + s[113:] if c in "ABCD" else s for c, s in native.items()
    }
    with pytest.raises(ValueError, match="protected"):
        audit_sequence(
            ":".join(protected.values()), native, {"alpha": [113], "beta": []}
        )
    rows = [
        {
            "candidate_id": str(i),
            "mutation_count": 1,
            "scores": [{"state": "TS", "AE_design_nll": v} for v in values],
        }
        for i, values in enumerate(([1, 2, 3], [2, 3, 4], [0, 0, 5]))
    ]
    ranked = rank_candidates(rows)
    assert [r["candidate_id"] for r in ranked] == ["2", "0", "1"]
    assert [r["pareto"] for r in ranked] == [True, True, False]
    with pytest.raises(ValueError, match="three"):
        rank_candidates([{"scores": []}])


def test_complete_candidate_score_gather(inputs: tuple) -> None:
    """Reject native-sequence scoring, incomplete maps and nonfinite native output."""
    _, mask, sequence = inputs
    stats = score_stats(sequence)
    assert summarize_score(stats, sequence, mask)["AE_design_nll"] == 2.0
    for key, value, message in (
        ("native_sequence", np.zeros(1312), "complete candidate"),
        ("residue_names", {0: "Z1"}, "order"),
        ("log_probs", np.zeros((1, 12, 21)), "dimensions"),
        ("log_probs", np.full((1, 1312, 21), np.nan), "Nonfinite"),
        ("mask", np.zeros(1312), "missing"),
    ):
        with pytest.raises(ValueError, match=message):
            summarize_score({**stats, key: value}, sequence, mask)
    with pytest.raises(ValueError, match="Empty"):
        summarize_score(stats, sequence, {"alpha": [], "beta": []})


def test_runner_success_and_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, inputs: tuple
) -> None:
    """Run all orchestration paths using native-shaped fake engine artifacts."""
    runner = load_script("run_relay_ensemble_pilot")
    native, mask, sequence = inputs
    pdbs = [tmp_path / f"parent_{i}.pdb" for i in range(3)]
    for p in pdbs:
        p.write_text("fixture")
    files = {
        str(p): __import__("hashlib").sha256(p.read_bytes()).hexdigest() for p in pdbs
    }
    required = ["weights", "run", "score", "wrapper", runner.__file__]
    for filename in required:
        p = Path(filename) if filename == runner.__file__ else tmp_path / filename
        if not p.exists():
            p.write_text("fixture")
        files[str(p)] = __import__("hashlib").sha256(p.read_bytes()).hexdigest()
    manifest = {
        "approved_for_launch": True,
        "samples_per_template": 8,
        "templates": [{"id": f"ts{i}", "pdb": str(p)} for i, p in enumerate(pdbs)],
        "mask": mask,
        "python": sys.executable,
        "checkpoint": str(tmp_path / "weights"),
        "generation_engine": str(tmp_path / "run"),
        "score_engine": str(tmp_path / "score"),
        "score_wrapper": str(tmp_path / "wrapper"),
        "support": "support",
        "environment": {**runner.ENVIRONMENT, "PYTHONPATH": "support"},
        "source_sha256": files,
        "symmetry_residues": "A1,B1,C1,D1",
        "symmetry_weights": "1,0,0,0",
        "redesigned_residues": "A1 B1 C1 D1",
        "parser_audit": {"accepted": True},
        "ES_control": {"status": "pending"},
    }
    monkeypatch.setattr(runner, "_pdb_sequences", lambda p: native)
    runner.validate(manifest)
    for key, bad in (
        ("approved_for_launch", False),
        ("samples_per_template", 9),
        ("environment", {}),
        ("symmetry_weights", "0.25,0.25,0.25,0.25"),
        ("symmetry_residues", "A2,B2,C2,D2"),
        ("redesigned_residues", "A2"),
        ("source_sha256", {}),
        ("parser_audit", {"accepted": False}),
    ):
        with pytest.raises(ValueError):
            runner.validate({**manifest, key: bad})
    loaded_sequence = {"value": sequence}
    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(load=lambda *a, **k: score_stats(loaded_sequence["value"])),
    )

    def fake_native(argv: list, output: Path, env: dict) -> None:
        if argv[1] == manifest["score_wrapper"]:
            loaded_sequence["value"] = json.loads(Path(argv[2]).read_text())["sequence"]
        if "--number_of_batches" in argv:
            out = Path(argv[argv.index("--out_folder") + 1])
            (out / "seqs").mkdir(parents=True)
            pdb = Path(argv[argv.index("--pdb_path") + 1])
            records = [("native", ":".join(native.values()))] + [
                ("design", sequence)
            ] * 8
            (out / "seqs" / f"{pdb.stem}.fa").write_text(
                "".join(f">{h}\n{s}\n" for h, s in records)
            )

    monkeypatch.setattr(runner, "native_run", fake_native)
    rows = runner.run(manifest, tmp_path / "success")
    assert len(rows) == 1 and rows[0]["mutation_count"] == 1
    assert len(rows[0]["lineage"]) == 24 and len(rows[0]["scores"]) == 3
    assert json.loads((tmp_path / "success/scientific_acceptance.json").read_text())[
        "accepted"
    ]
    generated = runner.run(manifest, tmp_path / "generation_only", generation_only=True)
    assert len(generated) == 1 and "scores" not in generated[0]
    es_manifest = {
        **manifest,
        "ES_template": {"id": "ES", "pdb": str(pdbs[0])},
        "ES_control": {"status": "accepted_exact_mapping"},
    }
    rescored = runner.run(
        es_manifest,
        tmp_path / "rescored",
        generation_source=tmp_path / "generation_only",
    )
    assert len(rescored[0]["scores"]) == 4
    assert rescored[0]["scores"][-1]["state"] == "ES"
    with pytest.raises(FileExistsError):
        runner.run(manifest, tmp_path / "success")
    monkeypatch.setattr(runner, "pool_samples", lambda samples: [])
    assert runner.run(manifest, tmp_path / "empty") == []
    assert (
        json.loads((tmp_path / "empty/scientific_acceptance.json").read_text())[
            "outcome"
        ]
        == "empty_pool"
    )


def test_native_process_receipts(tmp_path: Path) -> None:
    """A failed engine retains its exit record and raises instead of falling back."""
    runner = load_script("run_relay_ensemble_pilot")
    runner.native_run([sys.executable, "-c", "pass"], tmp_path / "ok", {})
    with pytest.raises(__import__("subprocess").CalledProcessError):
        runner.native_run(
            [sys.executable, "-c", "raise SystemExit(7)"], tmp_path / "bad", {}
        )
    assert json.loads((tmp_path / "bad.process.json").read_text())["exit_code"] == 7


def test_score_wrapper_injects_sequence_without_context_mutation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, inputs: tuple
) -> None:
    """Native parsing retains coordinate objects while every candidate AA is injected."""
    wrapper = load_script("score_relay_ensemble_sequence")
    native, mask, sequence = inputs

    class Tensor(np.ndarray):
        def clone(self) -> np.ndarray:
            return self.copy()

    alphabet = list("ACDEFGHIKLMNPQRSTVWYX")
    coordinates = object()
    protein = {
        "S": np.zeros(1312, dtype=int).view(Tensor),
        "chain_letters": [c for c in "ABCDEFGH" for _ in range(len(native[c]))],
        "R_idx": np.array(
            [p for c in "ABCDEFGH" for p in range(1, len(native[c]) + 1)]
        ),
        "xyz_37": coordinates,
    }
    data = SimpleNamespace(
        alphabet=alphabet, parse_PDB=lambda *a, **k: (protein, 1, 2, 3, 4)
    )
    torch = SimpleNamespace(
        cuda=SimpleNamespace(is_available=lambda: True), set_num_threads=lambda n: None
    )
    monkeypatch.setitem(sys.modules, "data_utils", data)
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setattr(wrapper, "_pdb_sequences", lambda p: native)
    request = tmp_path / "request.json"
    request.write_text(
        json.dumps(
            {
                "sequence": sequence,
                "mask": mask,
                "score_engine": str(tmp_path / "score.py"),
                "pdb": "parent.pdb",
                "argv": ["--use_sequence", "1"],
            }
        )
    )
    monkeypatch.setattr(sys, "argv", ["wrapper", str(request)])

    def native_main(path: str, run_name: str) -> None:
        parsed = data.parse_PDB("parent.pdb", parse_all_atoms=True)
        assert parsed[0]["xyz_37"] is coordinates
        assert parsed[1:] == (1, 2, 3, 4)
        observed = "".join(alphabet[i] for i in parsed[0]["S"])
        assert observed == sequence.replace(":", "")
        assert sys.argv[1:] == ["--use_sequence", "1"]

    monkeypatch.setattr(wrapper.runpy, "run_path", native_main)
    wrapper.main()
    torch.cuda.is_available = lambda: False
    monkeypatch.setattr(sys, "argv", ["wrapper", str(request)])
    with pytest.raises(RuntimeError, match="no CPU fallback"):
        wrapper.main()
