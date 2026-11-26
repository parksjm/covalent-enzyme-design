"""Prepare and execute the sequential hydration MD extension."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from . import deacylation_hydration as hydration

SEGMENT_NAMES = tuple(f"extend_{index:02d}" for index in range(1, 5))
STEPS = 10_000
SEEDS = (2_026_112_402, 2_026_112_403, 2_026_112_404, 2_026_112_405)


def _write(path: Path, value: dict[str, Any]) -> None:
    """Write finite JSON atomically."""
    hydration.write_json(path, value)


def prepare_extension(campaign: Path, pilot: Path, repo: Path) -> None:
    """Freeze four unbiased 5 ps segments from the accepted 6 ps endpoint.

    Args:
        campaign: Existing logical hydration reconnaissance directory.
        pilot: Completed hydration reconnaissance directory.
        repo: Repository containing source, journal, and logger.
    """
    if (campaign / "extension_20ps_config.json").exists():
        raise FileExistsError(campaign / "extension_20ps_config.json")
    pilot = pilot.resolve(strict=True)
    cfg_path = pilot / "config.json"
    cfg = json.loads(cfg_path.read_text())
    endpoint = json.loads(
        (pilot / "hydration_endpoint_status_20261124.json").read_text()
    )
    accepted = endpoint["final_acceptance"]
    if endpoint["pilot_status"]["status"] != "completed":
        raise ValueError("Pilot is not completed")
    if endpoint["source_restart"] != "equil/restart.ncrst":
        raise ValueError("Endpoint does not name the final accepted equil restart")
    start = pilot / endpoint["source_restart"]
    if hydration.digest(start) != endpoint["restart_sha256"]:
        raise ValueError("Accepted endpoint restart identity changed")
    if (
        accepted["numeric"]["completed_steps"] != STEPS
        or not accepted["numeric"]["electronic_history"]["verified_scc_history"]
        or accepted["boundary"]["frozen_max_displacement_A"] > 1e-6
        or accepted["boundary"]["box_unchanged"] is not True
        or accepted["chemistry"]["acyl_bond_retained"] is not True
        or accepted["chemistry"]["C_OG_A"] >= 1.8
    ):
        raise ValueError("Final pilot endpoint lacks required acceptance gates")
    if json.loads((pilot / "status.json").read_text())["status"] != "completed":
        raise ValueError("Pilot status is not completed")
    pilot_process = json.loads((pilot / "equil/process.json").read_text())
    pilot_acceptance = json.loads(
        (pilot / "equil/scientific_acceptance.json").read_text()
    )
    if (
        pilot_process.get("status") != "completed"
        or pilot_process.get("exit_code") != 0
    ):
        raise ValueError("Pilot equil process did not complete successfully")
    if pilot_acceptance != accepted:
        raise ValueError(
            "Pilot scientific acceptance differs from frozen endpoint status"
        )
    for name, expected in cfg["input_sha256"].items():
        if hydration.digest(pilot / name) != expected:
            raise ValueError(f"Pilot input identity changed: {name}")
    if hydration.digest(Path(cfg["audit_path"])) != cfg["audit_sha256"]:
        raise ValueError("Original endpoint audit identity changed")
    if hydration.digest(Path(cfg["sander"]["path"])) != cfg["sander"]["sha256"]:
        raise ValueError("Sander executable identity changed")

    campaign.mkdir(parents=True, exist_ok=True)
    source = (pilot / "inputs/source.mdin").read_text()
    base_record = json.loads((pilot / "config.json").read_text())
    active = base_record["active_mask"]
    for index, name in enumerate(SEGMENT_NAMES, start=1):
        segment = campaign / name
        inputs = segment / "inputs"
        inputs.mkdir(parents=True)
        (inputs / "segment.mdin").write_text(
            hydration.build_input(source, active, "equil", SEEDS[index - 1]).replace(
                "Acyl enzyme hydration only: equil",
                f"Acyl enzyme hydration extension: {name}",
                1,
            )
        )
        for src_name, dst_name in (
            ("system.prmtop", "system.prmtop"),
            ("source.mdin", "source.mdin"),
        ):
            (inputs / dst_name).write_bytes((pilot / "inputs" / src_name).read_bytes())
        if index == 1:
            (inputs / "restart.ncrst").write_bytes(start.read_bytes())
        _write(
            segment / "status.json",
            {
                "status": "prepared",
                "segment": index,
                "duration_ps": 5.0,
                "predecessor": "accepted_pilot"
                if index == 1
                else SEGMENT_NAMES[index - 2],
            },
        )

    source_files = {
        "nyl12_covalent/deacylation_hydration.py": repo
        / "src/nyl12_covalent/deacylation_hydration.py",
        "nyl12_covalent/__init__.py": repo / "src/nyl12_covalent/__init__.py",
        "nyl12_covalent/deacylation_hydration_prepare.py": repo
        / "src/nyl12_covalent/deacylation_hydration_prepare.py",
        "nyl12_covalent/deacylation_hydration_extension.py": repo
        / "src/nyl12_covalent/deacylation_hydration_extension.py",
        "nyl12_covalent/dftb_scc_validation.py": repo
        / "src/nyl12_covalent/dftb_scc_validation.py",
        "nyl12_covalent/neb_umbrella_inputs.py": repo
        / "src/nyl12_covalent/neb_umbrella_inputs.py",
        "scripts/run_deacylation_hydration_extension.py": repo
        / "scripts/run_deacylation_hydration_extension.py",
    }
    snapshot = campaign / "extension_source_snapshot"
    for rel, src in source_files.items():
        dst = snapshot / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(src.read_bytes())

    segment_inputs = {
        name: {
            path.name: hydration.digest(path)
            for path in sorted((campaign / name / "inputs").iterdir())
        }
        for name in SEGMENT_NAMES
    }
    config = {
        "scope": "20 ps unbiased NVT hydration continuation in four sequential 5 ps segments; no reaction, no position restraints, no water promotion; sampling not presumed converged",
        "pilot_dir": str(pilot),
        "journal_repo": cfg["journal_repo"],
        "pilot_config_sha256": hydration.digest(cfg_path),
        "endpoint_status_sha256": hydration.digest(
            pilot / "hydration_endpoint_status_20261124.json"
        ),
        "endpoint_restart_sha256": endpoint["restart_sha256"],
        "endpoint_restart_path": str(start),
        "topology_sha256": cfg["input_sha256"]["inputs/system.prmtop"],
        "native_input_sha256": cfg["input_sha256"]["inputs/source.mdin"],
        "original_qmmm": {"atoms": 98, "charge": -2, "links": 8},
        "active_atoms": 4165,
        "active_shell_A": 12,
        "active_mask": cfg["active_mask"],
        "qmmask": cfg["qmmask"],
        "sander": cfg["sander"],
        "dftb_parameters_sha256": cfg["dftb_parameters_sha256"],
        "segment_names": list(SEGMENT_NAMES),
        "segment_steps": STEPS,
        "dt_ps": 0.0005,
        "temperature_K": 300.0,
        "langevin_gamma_per_ps": 2.0,
        "seeds": dict(zip(SEGMENT_NAMES, SEEDS, strict=True)),
        "seed_provenance": "Use distinct explicit Langevin seeds 2026112402-2026112405 for the four segments; velocities continue through accepted restart files. Each sander process starts a new stochastic stream; this is not claimed as bitwise RNG continuation.",
        "segment_input_sha256": segment_inputs,
        "source_snapshot_sha256": {
            str(path.relative_to(snapshot)): hydration.digest(path)
            for path in sorted(path for path in snapshot.rglob("*") if path.is_file())
        },
        "postprocess_argv": None,
        "postprocess_source_sha256": {},
    }
    _write(campaign / "extension_20ps_config.json", config)
    _write(
        campaign / "extension_20ps_status.json",
        {"status": "prepared", "accepted_segments": 0},
    )


def _load_identity(
    campaign: Path,
) -> tuple[dict[str, Any], AmberParm, np.ndarray, np.ndarray, np.ndarray]:
    """Verify static identity and load original frozen-boundary reference."""
    cfg = json.loads((campaign / "extension_20ps_config.json").read_text())
    pilot = Path(cfg["pilot_dir"])
    if hydration.digest(pilot / "config.json") != cfg["pilot_config_sha256"]:
        raise ValueError("Pilot configuration identity changed")
    if (
        hydration.digest(pilot / "hydration_endpoint_status_20261124.json")
        != cfg["endpoint_status_sha256"]
    ):
        raise ValueError("Pilot endpoint acceptance identity changed")
    if (
        hydration.digest(Path(cfg["endpoint_restart_path"]))
        != cfg["endpoint_restart_sha256"]
    ):
        raise ValueError("Pilot accepted restart identity changed")
    if hydration.digest(Path(cfg["sander"]["path"])) != cfg["sander"]["sha256"]:
        raise ValueError("Sander executable identity changed")
    for path, expected in cfg["dftb_parameters_sha256"].items():
        if hydration.digest(Path(path)) != expected:
            raise ValueError(f"DFTB parameter identity changed: {path}")
    for rel, expected in cfg["source_snapshot_sha256"].items():
        if hydration.digest(campaign / "extension_source_snapshot" / rel) != expected:
            raise ValueError(f"Source snapshot identity changed: {rel}")
    if not cfg.get("postprocess_argv") or not cfg.get("postprocess_source_sha256"):
        raise ValueError(
            "A reviewed source-pinned postprocess hook is required before launch"
        )
    _verify_postprocess_sources(campaign, cfg)
    top = campaign / SEGMENT_NAMES[0] / "inputs/system.prmtop"
    if hydration.digest(top) != cfg["topology_sha256"]:
        raise ValueError("Topology identity changed")
    parm = AmberParm(str(top))
    active = np.asarray(AmberMask(parm, cfg["active_mask"]).Selection(), dtype=bool)
    qm = np.asarray(AmberMask(parm, cfg["qmmask"]).Selection(), dtype=bool)
    links = sum(bool(qm[b.atom1.idx]) != bool(qm[b.atom2.idx]) for b in parm.bonds)
    if active.sum() != 4165 or qm.sum() != 98 or links != 8 or not np.all(active[qm]):
        raise ValueError("QM or active-region identity changed")
    base = Rst7.open(cfg["endpoint_restart_path"])
    coords = np.asarray(base.coordinates).reshape(-1, 3)
    box = np.asarray(base.box)
    if coords.shape != (104922, 3) or box.shape != (6,):
        raise ValueError("Accepted pilot restart shape changed")
    return cfg, parm, active, coords, box


def _verify_predecessor(campaign: Path, index: int) -> Path:
    """Require the exact accepted endpoint from the previous segment."""
    if index == 0:
        return Path(
            json.loads((campaign / "extension_20ps_config.json").read_text())[
                "endpoint_restart_path"
            ]
        )
    previous = campaign / SEGMENT_NAMES[index - 1]
    status = json.loads((previous / "status.json").read_text())
    acceptance_path = previous / "scientific_acceptance.json"
    restart = previous / "restart.ncrst"
    if (
        status.get("status") != "accepted"
        or not acceptance_path.is_file()
        or not restart.is_file()
    ):
        raise ValueError("Predecessor segment is not scientifically accepted")
    acceptance = json.loads(acceptance_path.read_text())
    process = json.loads((previous / "process.json").read_text())
    if process.get("status") != "process_completed" or process.get("exit_code") != 0:
        raise ValueError("Predecessor native process did not complete successfully")
    if status.get("restart_sha256") != hydration.digest(restart):
        raise ValueError("Predecessor restart identity changed")
    if (
        acceptance.get("acyl_retained") is not True
        or acceptance.get("scc_verified") is not True
        or acceptance.get("finite_temperature_passed") is not True
        or acceptance.get("native_process_passed") is not True
        or acceptance.get("boundary", {}).get("box_unchanged") is not True
        or acceptance.get("boundary", {}).get("frozen_max_displacement_A", 1.0) > 1e-6
    ):
        raise ValueError("Predecessor scientific acceptance is incomplete")
    return restart


def _verify_postprocess_sources(
    campaign: Path, config: dict[str, Any]
) -> dict[str, str]:
    """Verify source hashes and retained copies for the reviewed hook."""
    copies = campaign / "extension_source_snapshot/postprocess"
    hashes = config["postprocess_source_sha256"]
    if not isinstance(hashes, dict) or not hashes:
        raise ValueError("Postprocess source hash map must be nonempty")
    observed = {}
    for source, expected in hashes.items():
        source_path = Path(source)
        if hydration.digest(source_path) != expected:
            raise ValueError(f"Postprocess source identity changed: {source}")
        retained = copies / source_path.name
        if hydration.digest(retained) != expected:
            raise ValueError(f"Retained postprocess source identity changed: {source}")
        observed[source] = expected
    return observed


def _journal(config: dict[str, Any], phase: str, action: str) -> None:
    """Append one hydration-extension event to the active action journal."""
    repo = Path(config["journal_repo"])
    subprocess.run(
        [
            sys.executable,
            str(repo / "scripts/record_agent_action.py"),
            "--log",
            str(repo / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "hydration_extension",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Root-reviewed authorized hydration extension segment",
            "--rollback",
            "Retain inputs, outputs, restart, status and acceptance; never retry in place",
        ],
        check=True,
    )


def _run_postprocess(campaign: Path, segment: Path, config: dict[str, Any]) -> None:
    """Run the reviewed hook and write a separate success/failure receipt."""
    before = _verify_postprocess_sources(campaign, config)
    replacements = {
        "campaign": str(campaign.resolve()),
        "segment_dir": str(segment.resolve()),
        "trajectory": str((segment / "trajectory.nc").resolve()),
        "restart": str((segment / "restart.ncrst").resolve()),
    }
    command = [part.format(**replacements) for part in config["postprocess_argv"]]
    receipt: dict[str, Any] = {
        "status": "running",
        "argv": command,
        "source_sha256_before": before,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    _journal(
        config,
        "intent",
        f"Run postprocess for {segment}: argv={command!r}; source_sha256={before!r}",
    )
    try:
        result = subprocess.run(
            command, cwd=segment, check=False, capture_output=True, text=True
        )
        stdout_path = segment / "postprocess.stdout.log"
        stderr_path = segment / "postprocess.stderr.log"
        stdout_path.write_text(result.stdout)
        stderr_path.write_text(result.stderr)
        after = _verify_postprocess_sources(campaign, config)
        receipt.update(
            status="completed" if result.returncode == 0 else "failed",
            exit_code=result.returncode,
            stdout_sha256=hydration.digest(stdout_path),
            stderr_sha256=hydration.digest(stderr_path),
            source_sha256_after=after,
        )
    except Exception as exc:
        receipt.update(status="failed", error=f"{type(exc).__name__}: {exc}")
    receipt["ended_utc"] = datetime.now(timezone.utc).isoformat()
    _write(segment / "postprocess.json", receipt)
    _journal(config, "outcome", f"Postprocess for {segment} receipt={receipt!r}")


def run_segment(campaign: Path, segment_name: str) -> dict[str, Any]:
    """Run one segment after verifying all identities and its predecessor.

    Args:
        campaign: Frozen four-segment campaign directory.
        segment_name: One of ``extend_01`` through ``extend_04``.

    Returns:
        Process and scientific acceptance status for this segment.
    """
    campaign = campaign.resolve(strict=True)
    if segment_name not in SEGMENT_NAMES:
        raise ValueError(f"Unknown segment: {segment_name}")
    index = SEGMENT_NAMES.index(segment_name)
    cfg, _parm, active, baseline, box = _load_identity(campaign)
    segment = campaign / segment_name
    status_path = segment / "status.json"
    status = json.loads(status_path.read_text())
    if status.get("status") != "prepared":
        raise ValueError("Segment is not in prepared state; no automatic retries")
    predecessor = _verify_predecessor(campaign, index)
    inputs = segment / "inputs"
    for name, expected in cfg["segment_input_sha256"][segment_name].items():
        if hydration.digest(inputs / name) != expected:
            raise ValueError(f"Segment input identity changed: {name}")
    start_target = inputs / "restart.ncrst"
    if index == 0:
        if (
            not start_target.is_file()
            or hydration.digest(start_target) != cfg["endpoint_restart_sha256"]
        ):
            raise ValueError("Prepared pilot restart identity changed")
    else:
        if start_target.exists():
            raise FileExistsError(
                "Dependent start restart already exists before predecessor acceptance"
            )
        start_target.write_bytes(predecessor.read_bytes())
    start_hash = hydration.digest(start_target)
    status.update(status="running", start_restart_sha256=start_hash)
    _write(status_path, status)
    command = [
        cfg["sander"]["path"],
        "-O",
        "-i",
        str(inputs / "segment.mdin"),
        "-p",
        str(inputs / "system.prmtop"),
        "-c",
        str(start_target),
        "-o",
        str(segment / "run.mdout"),
        "-r",
        str(segment / "restart.ncrst"),
        "-x",
        str(segment / "trajectory.nc"),
        "-inf",
        str(segment / "mdinfo"),
    ]
    env = dict(os.environ)
    env.update(
        AMBERHOME="/data/jpc/code/ambertools26/venv",
        LD_LIBRARY_PATH="/data/jpc/code/ambertools26/venv/lib",
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    seed = cfg["seeds"][segment_name]
    _write(
        segment / "command.json",
        {
            "argv": command,
            "cwd": str(segment),
            "start_restart_sha256": start_hash,
            "seed": seed,
        },
    )
    (segment / "engine.log").open("x").close()
    _journal(
        cfg,
        "intent",
        f"Run {segment_name}: argv={command!r}; start_restart_sha256={start_hash}",
    )
    status.update(started_utc=datetime.now(timezone.utc).isoformat())
    _write(status_path, status)
    with (segment / "engine.log").open("a") as log:
        process = subprocess.run(
            command,
            cwd=segment,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    result: dict[str, Any] = {
        "status": "process_failed" if process.returncode else "process_completed",
        "exit_code": process.returncode,
        "command": command,
        "start_restart_sha256": start_hash,
        "ended_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write(segment / "process.json", result)
    if process.returncode:
        failed_process = {**status, **result}
        _write(status_path, failed_process)
        _journal(
            cfg, "outcome", f"{segment_name} native process result={failed_process!r}"
        )
        return failed_process
    try:
        numeric = hydration.audit_output((segment / "run.mdout").read_text(), STEPS)
        xyz, boundary = hydration.audit_restart(
            segment / "restart.ncrst", baseline, active, box
        )
        distance = float(np.linalg.norm(xyz[19915] - xyz[13915]))
        acceptance = {
            "numeric": numeric,
            "boundary": boundary,
            "C_OG_A": distance,
            "finite_temperature_passed": 250.0
            <= numeric["final_temperature_K"]
            <= 350.0,
            "acyl_retained": distance < 1.8,
            "scc_verified": numeric["electronic_history"]["verified_scc_history"],
            "native_process_passed": True,
            "sampling_converged": False,
        }
        _write(segment / "scientific_acceptance.json", acceptance)
        if not acceptance["acyl_retained"]:
            stopped = {**status, **result, "status": "chemistry_stopped"}
            _write(status_path, stopped)
            _journal(
                cfg,
                "outcome",
                f"{segment_name} stopped at acyl-retention gate; acceptance={acceptance!r}; status={stopped!r}",
            )
            return stopped
        result.update(
            status="accepted",
            restart_sha256=hydration.digest(segment / "restart.ncrst"),
        )
        _write(status_path, {**status, **result})
        _journal(
            cfg,
            "outcome",
            f"{segment_name} native/scientific acceptance passed; acceptance={acceptance!r}; status={result!r}",
        )
    except Exception as exc:
        failed = {
            **status,
            **result,
            "status": "scientific_failed",
            "error": f"{type(exc).__name__}: {exc}",
        }
        _write(status_path, failed)
        _journal(
            cfg, "outcome", f"{segment_name} scientific gates failed; status={failed!r}"
        )
        return failed
    _run_postprocess(campaign, segment, cfg)
    return {
        **status,
        **result,
        "postprocess": json.loads((segment / "postprocess.json").read_text()),
    }


def run_extension(campaign: Path) -> dict[str, Any]:
    """Run four segments in order and stop on process or scientific failure."""
    outcomes = []
    _write(
        campaign / "extension_20ps_status.json",
        {"status": "running", "accepted_segments": 0, "segments": []},
    )
    for index, name in enumerate(SEGMENT_NAMES, start=1):
        result = run_segment(campaign, name)
        outcomes.append(result)
        if result["status"] != "accepted":
            _write(
                campaign / "extension_20ps_status.json",
                {"status": "stopped", "segment": name, "segments": outcomes},
            )
            return {"status": "stopped", "segment": name, "segments": outcomes}
        _write(
            campaign / "extension_20ps_status.json",
            {"status": "running", "accepted_segments": index, "segments": outcomes},
        )
    _write(
        campaign / "extension_20ps_status.json",
        {
            "status": "completed",
            "accepted_segments": 4,
            "sampling_converged": False,
            "segments": outcomes,
        },
    )
    return {
        "status": "completed",
        "accepted_segments": 4,
        "sampling_converged": False,
        "segments": outcomes,
    }
