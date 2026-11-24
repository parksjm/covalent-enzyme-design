"""Verify portable source closure and metadata-only conditional preparation."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import prepare_window033_relay_path as prep


def fixture_sources(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "repo"
    scripts = root / "scripts"
    package = root / "src/nyl12_covalent"
    scripts.mkdir(parents=True)
    package.mkdir(parents=True)
    for name in (
        "watch_window033_relay.py",
        "run_window033_relay_path.py",
        "prepare_window033_relay_path.py",
        "record_agent_action.py",
    ):
        (scripts / name).write_text("from nyl12_covalent.helper import value\n")
    (package / "__init__.py").write_text("")
    (package / "helper.py").write_text("from .sibling import value\n")
    (package / "sibling.py").write_text("value = 1\n")
    monkeypatch.setattr(prep, "ROOT", root)
    monkeypatch.setattr(
        prep, "__file__", str(scripts / "prepare_window033_relay_path.py")
    )
    return root


def test_static_import_closure_and_missing_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = fixture_sources(tmp_path, monkeypatch)
    target = tmp_path / "snapshot"
    files = prep.freeze_sources(target)
    assert len(files) == 7
    assert (target / "src/nyl12_covalent/sibling.py").read_text() == "value = 1\n"
    assert all(prep.digest(Path(p)) == h for p, h in files.items())
    (root / "src/nyl12_covalent/sibling.py").unlink()
    with pytest.raises(FileNotFoundError, match="Unresolved"):
        prep.freeze_sources(tmp_path / "another")


def test_prepare_preserves_exact_inputs_without_quantum_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture_sources(tmp_path, monkeypatch)
    seed = tmp_path / "seed"
    frozen = seed / "minimum_verification/01_frozen"
    (frozen / "inputs").mkdir(parents=True)
    (frozen / "inputs/system.prmtop").write_text("exact topology")
    (frozen / "inputs/optimized.ncrst").write_text("exact coordinates")
    (frozen / "manifest.json").write_text("{}")
    (frozen / "release_acceptance.json").write_text(
        json.dumps(
            {
                "accepted": True,
                "force_converged": True,
                "final_energy_hartree": -2.0,
            }
        )
    )
    viewer = seed / "vmd_progress_verified"
    viewer.mkdir()
    for name in ("stripped.parm7", "view_relay.vmd"):
        (viewer / name).write_text("original viewer")
    top_hash = prep.digest(frozen / "inputs/system.prmtop")
    (viewer / "manifest.json").write_text(json.dumps({"topology_sha256": top_hash}))
    (tmp_path / "prepared.json").write_text(json.dumps({"parameter_sha256": {}}))
    monkeypatch.setattr(prep, "verify_frozen", lambda *a: None)
    monkeypatch.setattr(
        prep,
        "NativeState",
        lambda root: SimpleNamespace(
            topology=root / "inputs/system.prmtop",
            sander=Path("/not-executed/sander"),
        ),
    )
    real_digest = prep.digest
    monkeypatch.setattr(
        prep,
        "digest",
        lambda path: (
            "binary-hash"
            if str(path).startswith(("/not-executed", "/data/jpc/code/ambertools26"))
            else real_digest(path)
        ),
    )
    monkeypatch.setattr(prep.subprocess, "check_output", lambda *a, **k: "revision\n")
    output = tmp_path / "campaign"
    config = prep.prepare(seed, output)
    assert (output / "model/inputs/optimized.ncrst").read_text() == "exact coordinates"
    assert config["source_topology_sha256"] == top_hash
    assert config["reference_energy_hartree"] == -2
    assert (
        json.loads((output / "preparation.json").read_text())[
            "native_calculations_launched"
        ]
        is False
    )
    with pytest.raises(FileExistsError):
        prep.prepare(seed, output)
    doc = json.loads((frozen / "release_acceptance.json").read_text())
    doc["accepted"] = False
    (frozen / "release_acceptance.json").write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="not force-accepted"):
        prep.prepare(seed, tmp_path / "bad")
