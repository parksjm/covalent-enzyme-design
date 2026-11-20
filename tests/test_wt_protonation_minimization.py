"""Acceptance and isolated preparation/queue ordering without running an engine."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nyl12_covalent import wt_protonation_minimization as support
from scripts import run_wt_protonation_campaign as runner


def native(rms: float = 0.05) -> str:
    """Return a native minimization completion summary."""
    return (
        f" NSTEP ENERGY RMS GMAX NAME NUMBER\n10000 -100.0 {rms} 0.9 N 1\n5. TIMINGS\n"
    )


def test_native_gates_distinguish_preparation_and_force_convergence() -> None:
    assert support.native_acceptance(native(), 0)["MD_preparation_accepted"]
    assert not support.native_acceptance(native(), 0)["force_converged"]
    assert support.native_acceptance(native(0.01), 0)["force_converged"]
    for text, code in [
        (native(0.101), 0),
        (native(), 1),
        (native() + "NaN", 0),
        ("TIMINGS", 0),
        ("", 0),
    ]:
        with pytest.raises(ValueError):
            support.native_acceptance(text, code)


def fixture_plan(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Prepare tiny pinned inputs and a harmless placeholder executable."""
    root = tmp_path / "campaign"
    root.mkdir()
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / "sample.py").write_text("# frozen\n")
    engine = tmp_path / "pmemd.cuda"
    engine.write_text("placeholder")
    engine.chmod(0o700)
    monkeypatch.setattr(runner, "SOURCES", ("sample.py",))
    monkeypatch.setattr(runner, "ENGINE", engine)
    support.write_json(
        root / "manifest.json", {"states": dict.fromkeys(support.STATES, {})}
    )
    for state in support.STATES:
        folder = root / state
        folder.mkdir()
        (folder / "topology.prmtop").write_text("topology")
        (folder / "start.rst7").write_text("start")
        support.write_json(
            folder / "preparation_acceptance.json",
            dict(
                preparation_accepted=True,
                process_complete=True,
                topology_sha256=support.digest(folder / "topology.prmtop"),
                coordinates_sha256=support.digest(folder / "start.rst7"),
            ),
        )
    return root, repo


def test_freeze_minimum_inputs_and_all_minima_precede_queue(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root, repo = fixture_plan(tmp_path, monkeypatch)
    runner.prepare(root, repo)
    assert (root / "frozen_protocol/sample.py").read_text() == "# frozen\n"
    plan = json.loads((root / "orchestration_plan.json").read_text())
    calls = []
    assert [r["gpu"] for r in plan["states"]] == ["0", "1", "0"]
    monkeypatch.setattr(runner, "require_gpu", lambda gpu: dict(index=gpu))
    monkeypatch.setattr(
        runner,
        "coordinate_acceptance",
        lambda top, rst, state: dict(
            restart_readable=True,
            topology_sha256=support.digest(top),
            coordinates_sha256=support.digest(rst),
        ),
    )

    def fake_run(argv: list[str], **kwargs: object) -> SimpleNamespace:
        calls.append(argv)
        if "-o" in argv:
            Path(argv[argv.index("-o") + 1]).write_text(native())
            Path(argv[argv.index("-r") + 1]).write_text("minimized")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    runner.execute(root)
    assert len(calls) == 5 and all("-o" in c for c in calls[:3])
    assert [c[2] for c in calls[3:]] == ["plan", "launch"]
    dataset = json.loads((root / "datasets.json").read_text())["datasets"]
    assert len(dataset) == 3 and all(
        len(d["protected_residues"]) == 24 for d in dataset.values()
    )
    rows = [
        r
        for p in (root / "md_configs").glob("*.json")
        for r in json.loads(p.read_text())["replicas"]
    ]
    assert (
        len(rows) == 12
        and len({r["seed"] + 10000 * i for r in rows for i in range(6)}) == 72
    )
    assert [r["gpu"] for r in rows].count("0") == 6
    with pytest.raises(FileExistsError):
        runner.execute(root)


def test_invalid_state_and_changed_frozen_source_fail(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    with pytest.raises(ValueError):
        support.protected_rows("other")
    root, repo = fixture_plan(tmp_path, monkeypatch)
    runner.prepare(root, repo)
    (root / "frozen_protocol/sample.py").write_text("changed")
    with pytest.raises(ValueError, match="Frozen source"):
        runner.execute(root)


def test_coordinate_and_proton_geometry_gates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import parmed.amber
    import prepared_md_support

    names = ("N", "H1", "H2", "OG1", "HG1")
    atoms = [SimpleNamespace(name=n, idx=i) for i, n in enumerate(names)]
    top = SimpleNamespace(
        atoms=atoms,
        residues=[SimpleNamespace(name="THR", atoms=atoms)],
        coordinates=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 2, 0], [0, 3, 0]]),
        box=np.array([30, 30, 30, 90, 90, 90]),
    )
    monkeypatch.setattr(parmed.amber, "AmberParm", lambda *a, **k: top)
    monkeypatch.setattr(prepared_md_support, "validate_protected", lambda *a: None)
    monkeypatch.setattr(
        support,
        "protected_rows",
        lambda state: [dict(amber_residue=1, sequence_residue=233, site="AE")],
    )
    path = tmp_path / "input"
    path.write_text("test")
    assert support.coordinate_acceptance(path, path, support.STATES[0])[
        "proton_inventory_accepted"
    ]
    top.coordinates[1] = [4, 0, 0]
    with pytest.raises(ValueError, match="Proton bond"):
        support.coordinate_acceptance(path, path, support.STATES[0])
    top.coordinates = np.array([[float("nan")] * 3] * 5)
    with pytest.raises(ValueError, match="Invalid topology"):
        support.coordinate_acceptance(path, path, support.STATES[0])
