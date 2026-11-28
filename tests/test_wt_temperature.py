"""Check the bounded WT E-S thermal protocol and actual accepted start preparation."""

from __future__ import annotations

from pathlib import Path
import shutil
from unittest.mock import Mock

import numpy as np
import pytest

from nyl12_covalent import wt_temperature_protocol as protocol
from nyl12_covalent.hydration_followon_gates import read
from scripts import prepare_wt_temperature as preparation


def test_native_protocol() -> None:
    inputs = protocol.mdin_files(34801000)
    assert len(inputs) == 3
    for index, (name, steps) in enumerate(protocol.STAGES):
        text = inputs[name]
        protocol.validate_input(text, name, 34801000 + index)
        assert f"nstlim={steps}" in text and "irest=1" in text
        assert "ntr=0, ibelly=0, ifqnt=0" in text
    assert "value1=300.0, value2=348.15" in inputs["01_ramp_nvt"]
    for seed in (-1, 0, 2147483646, True):
        with pytest.raises(ValueError):
            protocol.mdin_files(seed)


@pytest.mark.parametrize(
    "stage,old,new",
    [
        ("01_ramp_nvt", "ntr=0", "ntr=1"),
        ("01_ramp_nvt", "value1=300.0", "value1=0.0"),
        ("01_ramp_nvt", "value2=348.15", "value2=400"),
        ("03_check_npt", "taup=2.0", "taup=4.0"),
        ("03_check_npt", "/", "/\nDISANG='hidden.RST'"),
    ],
)
def test_unintended_changes_rejected(stage: str, old: str, new: str) -> None:
    inputs = protocol.mdin_files(1)
    seed = 1 + [name for name, _ in protocol.STAGES].index(stage)
    with pytest.raises(ValueError):
        protocol.validate_input(inputs[stage].replace(old, new), stage, seed)


@pytest.mark.parametrize(
    "text,passes",
    [
        ("NSTEP = 100000 TIME(PS) = 100700.0 TEMP(K) = 348.0\n EPtot = -123.4", True),
        ("NSTEP = 100000 TIME(PS) = 100700.0 TEMP(K) = 399.0\n EPtot = -123.4", False),
        ("NSTEP = 99999 TIME(PS) = 100700.0 TEMP(K) = 348.0\n EPtot = -123.4", False),
        ("NSTEP = 100000 TIME(PS) = 100700.0 TEMP(K) = 348.0\n EPtot = 1E999", False),
        ("FATAL ERROR: NaN", False),
    ],
)
def test_thermal_native_audit(text: str, passes: bool) -> None:
    if passes:
        record = protocol.audit_output(text, 100000)
        assert record["thermal_equilibration_accepted"] is False
        assert record["sampling_converged"] is False
    else:
        with pytest.raises(ValueError):
            protocol.audit_output(text, 100000)


def test_collect_exclusive(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = tmp_path / "source"
    source.write_bytes(b"original")
    target = tmp_path / "target"
    preparation.collect("violet", source, target)
    assert target.read_bytes() == b"original"
    with pytest.raises(FileExistsError):
        preparation.collect("violet", source, target)
    mocked = Mock()
    monkeypatch.setattr(preparation.subprocess, "run", mocked)
    preparation.collect("azalea", source, tmp_path / "remote")
    assert mocked.call_args.args[0][0] == "scp"


def test_actual_preparation_without_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    staged = preparation.BASE
    original_collect = preparation.collect
    remote_files = {
        "06_production_90ns.acceptance.json": "acceptance.json",
        "extension_100ns.json": "extension_100ns.json",
        "status.json": "job_status.json",
        "06_production_90ns.in": "native_300K.in",
        "06_production_90ns.out": "native_300K.out",
        "100ns_process_complete.json": "100ns_process_complete.json",
    }

    def collect(host: str, source: Path, target: Path) -> None:
        if host == "violet":
            original_collect(host, source, target)
        elif source.name.endswith(".rst7"):
            shutil.copyfile(
                staged / "inputs/replica_02_accepted300K_100ns.rst7", target
            )
        else:
            shutil.copyfile(
                staged / "provenance/replica_02" / remote_files[source.name], target
            )

    monkeypatch.setattr(preparation, "BASE", tmp_path)
    monkeypatch.setattr(preparation, "collect", collect)
    preparation.prepare()
    plan = read(tmp_path / "plan.json")
    assert plan["duration_per_replica_ps"] == 3200
    assert plan["native_launch_authorized"] is False
    assert len(plan["replicas"]) == 2
    assert plan["replicas"][0]["start_sha256"] != plan["replicas"][1]["start_sha256"]
    assert all(rep["endpoint_identity"]["atoms"] == 130815 for rep in plan["replicas"])
    assert "blackwell" in plan["replicas"][1]["host_runtime"]["engine"]
    assert (
        plan["replicas"][0]["host_runtime"]["sha256"]
        != plan["replicas"][1]["host_runtime"]["sha256"]
    )
    with pytest.raises(FileExistsError):
        preparation.prepare()


def test_endpoint_gates(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    top = tmp_path / "topology"
    top.write_bytes(b"wrong")
    with pytest.raises(ValueError, match="topology"):
        preparation.endpoint_identity(top, tmp_path / "restart")
    monkeypatch.setattr(preparation, "digest", lambda _: preparation.TOPOLOGY_SHA)
    parm = Mock(atoms=[Mock()], residues=[], bonds=[])
    frame = Mock(
        coordinates=np.zeros((1, 3)),
        velocities=np.ones((1, 3)),
        box=np.array([1, 1, 1, 90, 90, 90]),
        time=100499.9999,
    )
    monkeypatch.setattr(preparation, "AmberParm", lambda _: parm)
    monkeypatch.setattr(preparation.Rst7, "open", lambda _: frame)
    with pytest.raises(ValueError, match="substrates missing"):
        preparation.endpoint_identity(top, top)
    residue = Mock()
    residue.name = "P6N"
    parm.residues = [residue]
    assert preparation.endpoint_identity(top, top)["atoms"] == 1
    frame.velocities = np.zeros((1, 3))
    with pytest.raises(ValueError, match="Invalid full-system"):
        preparation.endpoint_identity(top, top)
    frame.velocities = np.ones((1, 3))
    frame.time = 500
    with pytest.raises(ValueError, match="final accepted"):
        preparation.endpoint_identity(top, top)
