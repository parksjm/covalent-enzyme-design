"""Preparation, native numerical gates and uncensored proton diagnostic regression."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.io import netcdf_file

from nyl12_covalent import neb_umbrella_campaign as campaign
from nyl12_covalent import neb_umbrella_runtime as runtime
from nyl12_covalent import neb_umbrella_diagnostics as diagnostics


def rows() -> list[dict]:
    return [
        dict(
            image=i + 1,
            q_acyl=float(i) / 10,
            q_h=float(i) / 20,
            full_restart=f"image{i}",
            sha256="test",
        )
        for i in range(24)
    ]


def test_centers_and_priority() -> None:
    original = rows()
    result = campaign.dense_centers(original, 0.1)
    points = np.array([[r["q_acyl_center"], r["q_h_center"]] for r in result])
    assert np.linalg.norm(np.diff(points, axis=0), axis=1).max() <= 0.1 + 1e-12
    assert [r["original_neb_image"] for r in result if r["original_neb_image"]] == list(
        range(1, 25)
    )
    assert {r["source_restart"] for r in result} <= {
        r["full_restart"] for r in original
    }
    priority = campaign.priority_order(result, 8)
    assert len(set(priority)) == len(result)
    assert priority[:2] == [result[0]["label"], result[-1]["label"]]
    assert (
        next(r["label"] for r in result if r["original_neb_image"] == 8) == priority[2]
    )


@pytest.mark.parametrize("step", [0, -1, float("nan")])
def test_bad_centers(step: float) -> None:
    with pytest.raises(ValueError):
        campaign.dense_centers(rows(), step)


def fixture_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    source = tmp_path / "source.py"
    source.write_text("source")
    monkeypatch.setattr(campaign, "ROOT", tmp_path)
    monkeypatch.setattr(campaign, "source_paths", lambda: [source])
    monkeypatch.setattr(campaign, "build_mdin", lambda *args: "native input")
    monkeypatch.setattr(campaign, "amber_restraints", lambda *args: "native restraint")
    image = tmp_path / "restart"
    image.write_text("full restart")
    top = tmp_path / "top"
    top.write_text("topology")
    old = tmp_path / "old"
    old.write_text("old source")
    audit = dict(
        passed=True,
        windows=[
            dict(r, full_restart=str(image), sha256=campaign.digest(image))
            for r in rows()
        ],
        active_atoms=4165,
        qm_atoms=98,
        qm_charge=-2,
        roles_1based={"og1": 1, "hg1": 2, "nalpha": 3},
        topology=str(top),
        source_input=str(old),
        active_mask=":1",
        source_sha256={str(top): campaign.digest(top)},
    )
    auditpath = tmp_path / "audit.json"
    campaign.write_json(auditpath, audit)
    config = tmp_path / "request.json"
    campaign.write_json(
        config,
        dict(
            seed_audit=str(auditpath),
            output=str(tmp_path / "out"),
            sander=str(source),
            python=str(source),
            workers=8,
            seed_base=202800001,
            temperature_K=300,
            temperature_gate_K=[250, 350],
            heating_temperature_gate_K=[25, 350],
            maximum_cv_step_A=0.1,
            k_acyl=100,
            k_h=100,
            source_peak_image=8,
        ),
    )
    return config


def test_prepare_freeze(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    request = fixture_config(tmp_path, monkeypatch)
    prepared = campaign.prepare(request)
    cfg = campaign.frozen(prepared)
    assert len(cfg["windows"]) == 47
    assert len({s["seed"] for r in cfg["windows"] for s in r["stages"]}) == 47 * 8
    assert len(cfg["priority_schedule"]) == 47
    with pytest.raises(FileExistsError):
        campaign.prepare(request)
    Path(cfg["windows"][0]["seed"]).write_text("changed")
    with pytest.raises(ValueError, match="changed"):
        campaign.frozen(prepared)
    prepared.write_text("changed")
    with pytest.raises(ValueError, match="configchanged"):
        campaign.frozen(prepared)


@pytest.mark.parametrize(
    "key,value", [("workers", 9), ("seed_base", 0), ("temperature_K", 290)]
)
def test_invalid_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, key: str, value: object
) -> None:
    path = fixture_config(tmp_path, monkeypatch)
    cfg = json.loads(path.read_text())
    cfg[key] = value
    campaign.write_json(path, cfg)
    with pytest.raises(ValueError):
        campaign.prepare(path)


def native_output(
    stage: Path, tail: str = "", temp: float = 300, steps: int = 10000
) -> None:
    (stage / "run.mdout").write_text(
        f"NSTEP = {steps} TIME(PS) = 5.0 TEMP(K) = {temp}\n EPtot = -123.0\n" + tail
    )


def test_stage_acceptance(tmp_path: Path) -> None:
    native_output(tmp_path)
    assert (
        runtime.stage_acceptance(tmp_path, "production", [250, 350])["completed_steps"]
        == 10000
    )
    for tail in ["NaN", "SCC convergence failed", "ERROR: failure", "EPtot = Infinity"]:
        native_output(tmp_path, tail)
        with pytest.raises(ValueError):
            runtime.stage_acceptance(tmp_path, "production", [250, 350])
    native_output(tmp_path, temp=400)
    with pytest.raises(ValueError, match="temperature"):
        runtime.stage_acceptance(tmp_path, "production", [250, 350])
    native_output(tmp_path, steps=50)
    with pytest.raises(ValueError, match="timestep"):
        runtime.stage_acceptance(tmp_path, "production", [250, 350])


def test_proton_switch_retained() -> None:
    roles = dict(
        carbon=1, leaving_n=2, thr_og=3, proton=4, thr_n=5, nalpha_h1=6, nalpha_h2=7
    )
    xyz = np.array(
        [
            [0.0, 0, 0],
            [2, 0, 0],
            [-2, 0, 0],
            [-1.1, 0, 0],
            [0, 5, 0],
            [0, 6, 0],
            [0, 4, 0],
        ]
    )
    top = SimpleNamespace(
        atoms=[
            SimpleNamespace(idx=i, atomic_number=1 if i in (3, 5, 6) else 7)
            for i in range(7)
        ]
    )
    normal = runtime.observables(xyz, top, roles)
    assert not normal["unexpected_proton_switch_flag"]
    xyz[3] = [0, 5.9, 0]
    switched = runtime.observables(xyz, top, roles)
    assert switched["unexpected_proton_switch_flag"]
    assert switched["direct_route_interpretation_blocked"]
    assert "retained" in switched["sample_retention"]


def test_trajectory_cv(tmp_path: Path) -> None:
    path = tmp_path / "trajectory.nc"
    with netcdf_file(path, "w") as stream:
        stream.createDimension("frame", 2)
        stream.createDimension("atom", 4)
        stream.createDimension("spatial", 3)
        var = stream.createVariable("coordinates", "d", ("frame", "atom", "spatial"))
        var[:] = np.array(
            [
                [[0, 0, 0], [2, 0, 0], [-1, 0, 0], [0, 0, 0]],
                [[0, 0, 0], [3, 0, 0], [-1, 0, 0], [-1, 0, 0]],
            ]
        )
    assert np.allclose(
        diagnostics.trajectory_cvs(
            path, dict(carbon=1, leaving_n=2, thr_og=3, proton=4)
        ),
        [[1, -1], [2, -4]],
    )


def test_run_window_success_and_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:

    roles = dict(
        carbon=1, leaving_n=2, thr_og=3, proton=4, thr_n=5, nalpha_h1=6, nalpha_h2=7
    )
    xyz = np.zeros((4166, 3))
    source = SimpleNamespace(
        coordinates=xyz,
        box=np.array([30, 30, 30, 90, 90, 90]),
        velocities=xyz,
        hasvels=True,
    )
    top = SimpleNamespace(atoms=[])
    monkeypatch.setattr(runtime.parmed, "load_file", lambda p: top)
    monkeypatch.setattr(runtime.Rst7, "open", lambda p: source)
    monkeypatch.setattr(
        runtime, "AmberMask", lambda *a: SimpleNamespace(Selected=lambda: range(4165))
    )
    monkeypatch.setattr(runtime, "observables", lambda *a: {"all_retained": True})
    monkeypatch.setattr(diagnostics, "block_diagnostics", lambda *a: None)
    monkeypatch.setattr(diagnostics, "export_view", lambda *a: None)
    prepared = tmp_path / "prepared.json"
    prepared.write_text("{}")
    topfile = tmp_path / "top"
    topfile.write_text("top")
    seed = tmp_path / "seed"
    seed.write_text("seed")
    row = dict(
        label="window_001",
        seed=str(seed),
        stages=[dict(label="stage1", kind="production")],
    )
    cfg = dict(
        windows=[row],
        topology=str(topfile),
        active_mask=":1",
        sander="/amber/bin/sander",
        roles_1based=roles,
        temperature_gate_K=[250, 350],
        heating_temperature_gate_K=[25, 350],
    )
    monkeypatch.setattr(runtime, "frozen", lambda p: cfg)
    stage = tmp_path / "window_001/stage1"
    stage.mkdir(parents=True)
    (stage / "run.mdin").write_text("input")
    (stage / "umbrella.RST").write_text("restraint")

    def run(stage: Path, *args: object) -> Path:
        native_output(stage)
        path = stage / "restart.ncrst"
        path.write_text("checkpoint")
        return path

    monkeypatch.setattr(runtime, "_run_stage", run)
    runtime.run_window(prepared, "window_001")
    assert (
        json.loads((stage.parent / "window_status.json").read_text())["status"]
        == "completed"
    )
    with pytest.raises(FileExistsError):
        runtime.run_window(prepared, "window_001")
    # Explicit second attempt directory; numerical errors create a retained failed receipt.
    row["label"] = "window_002"
    stage2 = tmp_path / "window_002/stage1"
    stage2.mkdir(parents=True)
    (stage2 / "run.mdin").write_text("input")
    (stage2 / "umbrella.RST").write_text("restraint")

    def fail(*args: object) -> Path:
        raise RuntimeError("SCC diagnostic")

    monkeypatch.setattr(runtime, "_run_stage", fail)
    with pytest.raises(RuntimeError):
        runtime.run_window(prepared, "window_002")
    assert (
        json.loads((stage2.parent / "window_status.json").read_text())["status"]
        == "failed"
    )


def test_block_diagnostics_exact_cross_bias(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from nyl12_covalent import neb_umbrella_observables as helper

    monkeypatch.setattr(
        "parmed.amber.mask.AmberMask",
        lambda *a: SimpleNamespace(Selected=lambda: [0, 1]),
    )

    def extract(*args: object) -> dict:
        folder = args[-1]
        folder.mkdir()
        (folder / "observables.csv").write_text("q_acyl,q_h\n0.01,0.02\n0.03,0.04\n")
        return {}

    monkeypatch.setattr(helper, "extract_observables", extract)
    windows = [
        dict(
            label=f"window_{i:03d}",
            q_acyl_center=float(i),
            q_h_center=float(i),
            bias=dict(
                q_acyl_center=float(i),
                q_h_center=float(i),
                k_acyl=100,
                k_h=100,
                outer_extent=10,
            ),
        )
        for i in (1, 2)
    ]
    cfg = dict(
        roles_1based={},
        topology="top",
        windows=windows,
        temperature_K=300,
        bias_convention="Amber rk",
    )
    for row in windows:
        stage = tmp_path / row["label"] / "03_production_block_01_5ps"
        stage.mkdir(parents=True)
        (stage / "run.mdin").write_text("qmmask='@1,2'")
        (stage / "trajectory.nc").write_text("trajectory")
        diagnostics.block_diagnostics(cfg, row, stage, object())
        with np.load(stage / "samples_and_cross_bias.npz") as data:
            assert data["bias_kcal_mol"].shape == (2, 2)
            assert data["bias_kcal_mol"][0, 0] == pytest.approx(
                100 * (0.99**2 + 0.98**2)
            )
        assert (stage / "cv_coverage.png").stat().st_size > 100
    receipt = json.loads((stage / "overlap_checkpoint.json").read_text())
    assert receipt["available_windows"] == 2
    assert receipt["overlaps"][0]["histogram_intersection"] == pytest.approx(1)


def test_export_portable_endpoints(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import parmed
    from nyl12_covalent import relay_progress_vmd as renderer

    top = parmed.Structure()
    elements = [6, 7, 8, 1, 7, 1, 1]
    for i, element in enumerate(elements):
        top.add_atom(
            parmed.Atom(name=f"A{i}", atomic_number=element),
            "ALA" if i == 6 else "P6N",
            2 if i == 6 else 1,
            chain="A",
        )
    for a, b in [(0, 1), (0, 2), (2, 3), (4, 5), (4, 6)]:
        top.bonds.append(parmed.Bond(top.atoms[a], top.atoms[b]))
    coords = np.arange(21, dtype=float).reshape(7, 3)
    monkeypatch.setattr(
        diagnostics.Rst7, "open", lambda *a: SimpleNamespace(coordinates=coords)
    )
    monkeypatch.setattr(
        "parmed.amber.mask.AmberMask",
        lambda *a: SimpleNamespace(Selected=lambda: range(7)),
    )
    monkeypatch.setattr(renderer, "substrate_indices", lambda *a: set(range(7)))
    folder = tmp_path / "window_001"
    folder.mkdir()
    cfg = dict(roles_1based=dict(carbon=1, leaving_n=2, thr_og=3))
    row = dict(label="window_001", seed="seed", bias={}, stages=[dict(label="heat")])
    (folder / "heat").mkdir()
    (folder / "heat/run.mdin").write_text("qmmask='@1-7'")
    for i in (1, 2):
        stage = folder / f"03_production_block_{i:02d}_5ps"
        stage.mkdir()
        (stage / "acceptance.json").write_text("{}")
    diagnostics.export_view(cfg, row, folder, top)
    output = folder / "progress/latest"
    assert (output / "endpoints.xyz").is_file()
    assert "DynamicBonds 1.4 0.1 12" in (output / "view_umbrella.vmd").read_text()
    assert (
        "::nyl12_movie_policy::apply $molid"
        in (output / "view_umbrella.vmd").read_text()
    )
    assert json.loads((output / "manifest.json").read_text())["frames"] == 2
    for line in (output / "MANIFEST.sha256").read_text().splitlines():
        checksum, name = line.split()
        assert campaign.digest(output / name) == checksum


def test_queue_launch_and_priority_worker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import scripts.run_neb_umbrella as driver

    prepared = tmp_path / "prepared.json"
    prepared.write_text("{}")
    windows = [dict(label=f"window_{i:03d}") for i in range(4)]
    cfg = dict(
        windows=windows,
        priority_schedule=[r["label"] for r in reversed(windows)],
        workers=2,
        sander="/amber/bin/sander",
        python="/python",
        runtime_script="script.py",
        session_prefix="test",
        prepared_input_sha256={},
    )
    monkeypatch.setattr(driver, "frozen", lambda p: cfg)
    calls = []

    def launch(*args: object, **kwargs: object) -> None:
        calls.append(args)
        job = args[1]
        job.mkdir(parents=True)
        (job / "status.json").write_text('{"ended_utc":"now","exit_status":0}\n')
        if "window_" in job.parent.name:
            (job.parent / "window_status.json").write_text('{"status":"completed"}')

    monkeypatch.setattr(driver, "launch_job", launch)
    driver.launch(prepared)
    assert len(calls) == 2
    driver.worker(prepared, 0, poll=0)
    assert calls[2][3][-1] == "window_003"
    assert calls[3][3][-1] == "window_001"
    assert (
        json.loads((tmp_path / "workers/worker_01/status.json").read_text())["status"]
        == "completed"
    )
    with pytest.raises(ValueError):
        driver.worker(prepared, 3)


def test_heat_230_explicit_gate_accepted(tmp_path: Path) -> None:
    native_output(tmp_path, temp=230, steps=2000)
    assert (
        runtime.stage_acceptance(tmp_path, "heat", [25, 350])["final_temperature_K"]
        == 230
    )


def test_production_230_rejected(tmp_path: Path) -> None:
    native_output(tmp_path, temp=230, steps=10000)
    with pytest.raises(ValueError, match="temperature"):
        runtime.stage_acceptance(tmp_path, "production", [250, 350])


def test_native_short_filenames(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stage = tmp_path / ("long_stage_directory_" * 8)
    stage.mkdir()
    topology = tmp_path / "source.prmtop"
    topology.write_text("topology")
    start = tmp_path / "source.ncrst"
    start.write_text("restart")

    def execute(command: list[str], **kwargs: object) -> SimpleNamespace:
        assert kwargs["cwd"] == stage
        assert all(len(arg) < 80 and "/" not in arg for arg in command[1:])
        assert (stage / "input.ncrst").read_text() == "restart"
        assert (stage / "topology.prmtop").read_text() == "topology"
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runtime.subprocess, "run", execute)
    assert (
        runtime._run_stage(stage, Path("/amber/bin/sander"), topology, start, {})
        == stage / "restart.ncrst"
    )
    assert json.loads((stage / "process.json").read_text())["exit_status"] == 0
