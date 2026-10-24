"""Small-file tests for source provenance, production gates and seed exports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from scripts import audit_nyl12_pa6_md_reuse as audit


def _production(directory: Path, **overrides: Any) -> None:
    values = {
        "imin": 0,
        "irest": 1,
        "nstlim": 100,
        "dt": 0.002,
        "ntwx": 10,
        "ig": 123,
        "ntr": 0,
        "ntb": 2,
    }
    values.update(overrides)
    directory.mkdir(exist_ok=True)
    fields = ", ".join(f"{key}={value}" for key, value in values.items())
    (directory / "07_prod_100ns.in").write_text(f"Test\n&cntrl\n{fields},\n/\n")
    (directory / "07_prod_100ns.out").write_text(
        "NSTEP = 0\nNSTEP = 100\nFinal Performance Info\n"
    )


def test_production_uses_actual_parmed_mdin(tmp_path: Path) -> None:
    _production(tmp_path)
    result = audit._production_record(tmp_path)
    assert result["normal_completion"] is True
    assert result["duration_ns"] == pytest.approx(0.0002)
    assert result["frame_dt_ns"] == pytest.approx(0.00002)
    assert result["declared"]["ig"] == 123


@pytest.mark.parametrize(
    "change",
    [
        "marker",
        "steps",
        "incomplete",
        "restrained",
        "minimize",
        "restart",
        "periodicity",
        "zero_timestep",
        "negative_timestep",
        "zero_output_stride",
        "zero_duration",
    ],
)
def test_production_rejects_corruption(tmp_path: Path, change: str) -> None:
    overrides = {
        "restrained": {"ntr": 1},
        "minimize": {"imin": 1},
        "restart": {"irest": 0},
        "periodicity": {"ntb": 1},
        "zero_timestep": {"dt": 0.0},
        "negative_timestep": {"dt": -0.002},
        "zero_output_stride": {"ntwx": 0},
        "zero_duration": {"nstlim": 0},
    }
    _production(tmp_path, **overrides.get(change, {}))
    output = tmp_path / "07_prod_100ns.out"
    if change in {"marker", "steps", "incomplete"}:
        output.write_text(
            {
                "marker": "NSTEP = 100\n",
                "steps": "Final Performance Info\n",
                "incomplete": "NSTEP = 99\nFinal Performance Info\n",
            }[change]
        )
    with pytest.raises(ValueError):
        audit._production_record(tmp_path)


def test_source_hash(tmp_path: Path) -> None:
    path = tmp_path / "source.dat"
    path.write_bytes(b"original source\n")
    record = audit._source_record(path)
    assert record["bytes"] == 16
    assert record["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record["path"] == str(path.resolve())


def test_changed_during_hash(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "source.dat"
    path.write_bytes(b"data")
    original = Path.stat
    calls = 0

    def changing_stat(self: Path, *args: Any, **kwargs: Any) -> Any:
        nonlocal calls
        if self == path:
            calls += 1
            return SimpleNamespace(st_size=4, st_mtime_ns=calls)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Path, "stat", changing_stat)
    with pytest.raises(ValueError, match="Source changed"):
        audit._source_record(path)


@pytest.fixture
def native() -> tuple[Any, list[dict[str, Any]]]:
    residues, rows = [], []
    for number, (label, name, atoms) in enumerate(
        [
            ("thr233", "THR", ["N", "OG1"]),
            ("asp272", "ASP", ["OD1", "OD2"]),
            ("asp274", "ASP", ["OD1", "OD2"]),
        ],
        1,
    ):
        residues.append(
            SimpleNamespace(
                chain="A",
                number=number,
                name=name,
                atoms=[
                    SimpleNamespace(
                        name=atom, xx=float(number), xy=0.0, xz=0.0, other_locations={}
                    )
                    for atom in atoms
                ],
            )
        )
        rows.append({"label": label, "chain": "A", "number": number, "resname": name})
    crystal = SimpleNamespace(
        residues=residues, get_coordinates=lambda: np.zeros((1, 6, 3))
    )
    return crystal, rows


@pytest.mark.parametrize(
    "corruption", [None, "models", "identity", "duplicate", "alternate", "nonfinite"]
)
def test_native_mapping(
    native: tuple[Any, list[dict[str, Any]]],
    monkeypatch: pytest.MonkeyPatch,
    corruption: str | None,
) -> None:
    crystal, rows = native
    if corruption == "models":
        crystal.get_coordinates = lambda: np.zeros((2, 6, 3))
    elif corruption == "identity":
        crystal.residues[0].name = "SER"
    elif corruption == "duplicate":
        crystal.residues[0].atoms.append(crystal.residues[0].atoms[0])
    elif corruption == "alternate":
        crystal.residues[0].atoms[0].other_locations = {"B": object()}
    elif corruption == "nonfinite":
        crystal.residues[0].atoms[0].xx = float("nan")
    monkeypatch.setattr(audit.CIFFile, "parse", lambda *a, **k: crystal)
    if corruption:
        with pytest.raises(ValueError):
            audit._native_mapping(Path("native.cif"), rows)
    else:
        result = audit._native_mapping(Path("native.cif"), rows)
        assert len(result["distances_angstrom"]) == 8
        assert all(row["native_identity_verified"] for row in result["rows"])


@pytest.mark.parametrize("corruption", [None, "missing", "models", "nonfinite"])
def test_extract_frame(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    corruption: str | None,
) -> None:
    class Structure:
        atoms = [object(), object()]
        residues = [SimpleNamespace(name="WAT"), SimpleNamespace(name="THR")]
        coordinates = np.full((2, 3), np.nan if corruption == "nonfinite" else 0.0)

        def get_coordinates(self) -> np.ndarray:
            return np.zeros((2 if corruption == "models" else 1, 2, 3))

        def save(self, path: str, **kwargs: Any) -> None:
            Path(path).write_text("small export")

        def __getitem__(self, mask: str) -> Any:
            assert mask == ":1-1313"
            return self

    def run(*args: Any, **kwargs: Any) -> None:
        assert kwargs["check"] is True
        if corruption != "missing":
            (tmp_path / "replica_01_frame_0007" / "solvated.nc").write_bytes(b"frame")

    monkeypatch.setattr(audit.subprocess, "run", run)
    monkeypatch.setattr(audit, "AmberParm", lambda *a, **k: Structure())
    if corruption:
        with pytest.raises(ValueError):
            audit._extract_frame(
                tmp_path, "replica_01", {"full_frame": 7}, tmp_path, tmp_path
            )
    else:
        result = audit._extract_frame(
            tmp_path, "replica_01", {"full_frame": 7}, tmp_path, tmp_path
        )
        assert result["water_count"] == 1
        assert result["atom_count"] == 2
        assert Path(result["directory"], "solvated.rst7").is_file()


@pytest.mark.parametrize("corruption", [None, "seed", "identity"])
def test_main_small_sources(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    corruption: str | None,
) -> None:
    config = {
        "system": str(tmp_path),
        "amberhome": str(tmp_path),
        "replicas": ["replica_01", "replica_02"],
        "atom_map": {},
        "native_mapping": [{"topology_residue": 1, "resname": "THR"}],
        "native_reference": str(tmp_path / "native.cif"),
        "geometry_csv": str(tmp_path / "geometry.csv"),
        "model": "wt",
        "source_files": ["topology"],
        "replica_files": ["trajectory"],
        "windows_ns": [[20, 45]],
        "geometry_stride": 10,
    }
    for name in ["native.cif", "geometry.csv"]:
        (tmp_path / name).write_text("source")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(config))
    output = tmp_path / "output"
    monkeypatch.setattr(
        audit.sys,
        "argv",
        ["audit", "--manifest", str(manifest), "--output", str(output)],
    )

    def production(directory: Path) -> dict[str, Any]:
        seed = 1 if corruption == "seed" else int(directory.name[-1])
        return {
            "directory": str(directory),
            "declared": {"ig": seed},
            "frame_dt_ns": 0.02,
        }

    monkeypatch.setattr(audit, "_production_record", production)
    monkeypatch.setattr(audit, "audit_es_topology", lambda *a: {"checked": True})
    monkeypatch.setattr(audit, "_native_mapping", lambda *a: {"checked": True})
    monkeypatch.setattr(
        audit,
        "AmberParm",
        lambda *a: SimpleNamespace(
            residues=[
                SimpleNamespace(name="SER" if corruption == "identity" else "THR")
            ]
        ),
    )
    monkeypatch.setattr(audit, "_source_record", lambda path: {"path": str(path)})
    monkeypatch.setattr(audit, "select_es_frames", lambda *a: [{"full_frame": 1}])
    monkeypatch.setattr(audit, "_extract_frame", lambda *a: {"export": True})
    if corruption:
        with pytest.raises(ValueError):
            audit.main()
    else:
        audit.main()
        result = json.loads((output / "reuse_audit.json").read_text())
        assert len(result["selected_frames"]) == 2
        assert result["acceptance"]["reuse_audit_complete"] is True
        assert result["acceptance"]["qmmm_hamiltonian_validated"] is False
        assert (output / "source" / "manifest.json").is_file()
