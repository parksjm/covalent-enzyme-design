"""Mocked force-consistency checks; these tests perform no QM calculation."""

import copy
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from nyl12_covalent import es_dftb_engine as engine


@pytest.fixture
def config(tmp_path: Path) -> dict[str, Any]:
    paths = {}
    for name in ("python_executable", "topology", "coordinates", "parameters"):
        path = tmp_path / name
        path.write_text(name)
        paths[name] = str(path)
    return {
        **{key: value for key, value in paths.items() if key != "parameters"},
        "parameter_files": [paths["parameters"]],
        "environment": {
            "AMBERHOME": "/amber",
            "LD_LIBRARY_PATH": "/lib",
            "OMP_NUM_THREADS": "1",
        },
        "mm_options": dict(ntb=1, igb=0, ntc=1, ntf=1, ntr=0, ifqnt=1, cut=8.0),
        "qm_options": dict(
            qm_theory="DFTB3",
            qmcharge=0,
            qmmask=":1",
            adjust_q=0,
            qm_ewald=1,
            qm_pme=1,
            scfconv=1e-10,
            qmmm_int=1,
            qmcut=8.0,
            lnk_dis=1.09,
            lnk_atomic_no=1,
            lnk_method=1,
        ),
        "force_absolute_tolerance": 1e-7,
        "force_relative_tolerance": 1e-7,
        "energy_reproducibility_tolerance": 1e-9,
        "step_sizes_angstrom": [0.001, 0.002],
        "components": [{"atom": 1, "axis": 0}],
        "expected_inventory": dict(
            real_atoms=1, link_atoms=1, electron_count=7, boundary_pairs=[[1, 2]]
        ),
    }


def _save(path: Path, value: Any) -> Path:
    path.write_text(json.dumps(value))
    return path


def test_valid_config_with_explicit_extra_environment(config: dict[str, Any]) -> None:
    config["environment"]["CUDA_VISIBLE_DEVICES"] = ""
    engine._validate(config)


@pytest.mark.parametrize(
    "key", ["python_executable", "topology", "coordinates", "parameter_files"]
)
def test_missing_prerequisites(config: dict[str, Any], key: str) -> None:
    del config[key]
    with pytest.raises(ValueError, match="prerequisite"):
        engine._validate(config)


@pytest.mark.parametrize(
    "section,key,value",
    [
        ("environment", "OMP_NUM_THREADS", "2"),
        ("mm_options", "ntc", 2),
        ("mm_options", "cut", float("nan")),
        ("mm_options", "cut", 0),
        ("qm_options", "qm_theory", "PM6"),
        ("qm_options", "qm_ewald", 0),
        ("qm_options", "qmcut", 9.0),
        ("qm_options", "lnk_dis", 1.1),
        ("qm_options", "scfconv", 1e-5),
        ("qm_options", "scfconv", float("nan")),
    ],
)
def test_invalid_options(
    config: dict[str, Any], section: str, key: str, value: Any
) -> None:
    config[section][key] = value
    with pytest.raises(ValueError):
        engine._validate(config)


def test_missing_required_environment_and_inventory(config: dict[str, Any]) -> None:
    bad = copy.deepcopy(config)
    del bad["environment"]["AMBERHOME"]
    with pytest.raises(ValueError, match="environment"):
        engine._validate(bad)
    del config["expected_inventory"]["electron_count"]
    with pytest.raises(ValueError, match="inventory"):
        engine._validate(config)


@pytest.mark.parametrize(
    "steps", [[], [0.01], [0.01, 0.01], [-1, 1], [1, float("inf")]]
)
def test_invalid_steps(config: dict[str, Any], steps: list[float]) -> None:
    config["step_sizes_angstrom"] = steps
    with pytest.raises(ValueError, match="steps"):
        engine._validate(config)


@pytest.mark.parametrize(
    "key,value",
    [
        ("force_absolute_tolerance", -1),
        ("force_relative_tolerance", float("nan")),
        ("energy_reproducibility_tolerance", float("inf")),
    ],
)
def test_invalid_tolerances(config: dict[str, Any], key: str, value: float) -> None:
    config[key] = value
    with pytest.raises(ValueError, match="tolerance"):
        engine._validate(config)


@pytest.mark.parametrize(
    "components",
    [
        [],
        [{"atom": 0, "axis": 0}],
        [{"atom": True, "axis": 0}],
        [{"atom": 1, "axis": 3}],
    ],
)
def test_invalid_components(
    config: dict[str, Any], components: list[dict[str, Any]]
) -> None:
    config["components"] = components
    with pytest.raises(ValueError):
        engine._validate(config)


@pytest.fixture
def fake_api(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, config: dict[str, Any]
) -> Any:
    atoms = [
        SimpleNamespace(idx=0, atomic_number=6, charge=0.2),
        SimpleNamespace(idx=1, atomic_number=1, charge=-0.2),
    ]
    structure = SimpleNamespace(
        coordinates=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        box=np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0]),
        atoms=atoms,
        bonds=[SimpleNamespace(atom1=atoms[0], atom2=atoms[1])],
    )
    structure.get_coordinates = lambda: np.array([structure.coordinates])
    monkeypatch.setattr(engine, "AmberParm", lambda *args, **kwargs: structure)
    monkeypatch.setattr(
        engine, "AmberMask", lambda *args: SimpleNamespace(Selection=lambda: [1, 0])
    )
    api_file = tmp_path / "mock_sander.py"
    api_file.write_text("Mock only; no scientific engine")
    api = SimpleNamespace(
        __file__=str(api_file),
        _pys=SimpleNamespace(__file__=str(api_file)),
        bad_energy=False,
        bad_force=False,
        change_box=False,
    )
    api.pme_input = lambda: SimpleNamespace(**config["mm_options"])
    api.QmInputOptions = lambda: SimpleNamespace(**config["qm_options"])

    class Context:
        def __init__(self, xyz: Any, box: Any) -> None:
            self.xyz = xyz
            self.box = box.copy()

        def __enter__(self) -> Any:
            return self

        def __exit__(self, *args: Any) -> None:
            return None

        def energy_forces(self) -> Any:
            if api.change_box:
                self.box[0] += 1
            energy = float("nan") if api.bad_energy else float(np.sum(self.xyz**2))
            forces = np.full(self.xyz.shape, np.nan) if api.bad_force else -2 * self.xyz
            return SimpleNamespace(tot=energy), forces

    api.setup = lambda topology, xyz, box, mm, qm: Context(xyz, box)
    original_import = engine.importlib.import_module
    monkeypatch.setattr(
        engine.importlib,
        "import_module",
        lambda name: api if name == "sander" else original_import(name),
    )
    return SimpleNamespace(api=api, structure=structure)


def test_worker_analytic_force_and_perturbation(
    tmp_path: Path, config: dict[str, Any], fake_api: Any
) -> None:
    request = _save(
        tmp_path / "request.json",
        {"config": config, "perturbation": {"atom": 1, "axis": 0, "delta": 0.25}},
    )
    result = engine.evaluate_request(request)
    assert result["forces"] == [-2.5]
    assert result["inventory"]["boundary_pairs"] == [[1, 2]]
    assert np.load(tmp_path / "analytic_forces.npy").shape == (2, 3)
    assert json.loads((tmp_path / "result.json").read_text()) == result


@pytest.mark.parametrize(
    "problem",
    [
        "bad_energy",
        "bad_force",
        "change_box",
        "inventory",
        "coordinates",
        "box",
        "frame_count",
    ],
)
def test_worker_rejects_bad_results(
    tmp_path: Path, config: dict[str, Any], fake_api: Any, problem: str
) -> None:
    if problem in ("bad_energy", "bad_force", "change_box"):
        setattr(fake_api.api, problem, True)
    elif problem == "inventory":
        config["expected_inventory"]["real_atoms"] = 2
    elif problem == "coordinates":
        fake_api.structure.coordinates[0, 0] = np.nan
    elif problem == "box":
        fake_api.structure.box[0] = np.inf
    else:
        fake_api.structure.get_coordinates = lambda: np.zeros((2, 2, 3))
    request = _save(tmp_path / "request.json", {"config": config, "perturbation": None})
    with pytest.raises(ValueError):
        engine.evaluate_request(request)


def test_unknown_api_option(config: dict[str, Any], fake_api: Any) -> None:
    fake_api.api.pme_input = lambda: SimpleNamespace()
    with pytest.raises(ValueError, match="Unsupported"):
        engine._options(fake_api.api, config)


def _mock_subprocess(
    monkeypatch: pytest.MonkeyPatch, mode: str, config: dict[str, Any]
) -> None:
    def run(command: list[str], **kwargs: Any) -> Any:
        directory = Path(kwargs["cwd"])
        request = json.loads(Path(command[-1]).read_text())
        delta = (
            0 if request["perturbation"] is None else request["perturbation"]["delta"]
        )
        energy = (1 + delta) ** 2
        if mode == "repeat" and directory.name == "baseline_repeat":
            energy += 1
        if mode == "source" and directory.name == "baseline_repeat":
            Path(config["parameter_files"][0]).write_text("changed")
        if mode == "failure":
            return SimpleNamespace(returncode=3)
        if mode == "broken_environment":
            raise TypeError("mock subprocess rejects nonstring environment value")
        _save(
            directory / "result.json",
            {
                "energy_kcal_per_mol": energy,
                "forces": [2.0 if mode == "wrong_force" else -2.0],
            },
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(subprocess, "run", run)


@pytest.mark.parametrize(
    "mode,accepted",
    [("correct", True), ("wrong_force", False), ("repeat", False), ("source", False)],
)
def test_central_difference_and_independent_gates(
    tmp_path: Path,
    config: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
    accepted: bool,
) -> None:
    _mock_subprocess(monkeypatch, mode, config)
    result = engine.check_es_dftb_forces(
        _save(tmp_path / "config.json", config), tmp_path / "evidence"
    )
    assert result["accepted"] is accepted
    assert len(result["components"]) == 2
    assert result["components"][0]["finite_difference_force"] == pytest.approx(-2.0)
    assert result["sources_unchanged"] is (mode != "source")
    assert (tmp_path / "evidence/source/es_dftb_engine.py").is_file()
    assert (
        json.loads((tmp_path / "evidence/baseline/process_status.json").read_text())[
            "exit_status"
        ]
        == 0
    )


@pytest.mark.parametrize(
    "mode,error",
    [("failure", subprocess.CalledProcessError), ("broken_environment", TypeError)],
)
def test_controller_records_failures(
    tmp_path: Path,
    config: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
    error: type[Exception],
) -> None:
    _mock_subprocess(monkeypatch, mode, config)
    if mode == "broken_environment":
        config["environment"]["EXTRA"] = None
    with pytest.raises(error):
        engine.check_es_dftb_forces(
            _save(tmp_path / "config.json", config), tmp_path / "evidence"
        )
    assert (
        json.loads((tmp_path / "evidence/scientific_acceptance.json").read_text())[
            "accepted"
        ]
        is False
    )


def test_missing_file_and_existing_output(
    tmp_path: Path, config: dict[str, Any]
) -> None:
    config_path = _save(tmp_path / "config.json", config)
    with pytest.raises(FileExistsError):
        engine.check_es_dftb_forces(config_path, tmp_path)
    Path(config["topology"]).unlink()
    with pytest.raises(FileNotFoundError):
        engine.check_es_dftb_forces(config_path, tmp_path / "new")
