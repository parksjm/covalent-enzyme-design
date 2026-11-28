"""Execute exactly one frozen targeted relay label with native acceptance."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from .neb_umbrella_campaign import digest, write_json
from .neb_umbrella_runtime import _run_stage
from .neb_umbrella_segments import _checkpoint
from .relay_targeted import verify_target
from .relay_umbrella_runtime import ensure_diagnostics, native_acceptance


def run_support_target(prepared: Path, label: str, expected_sha256: str) -> None:
    """Run one fresh label through 5 ps equilibration and six 5 ps blocks.

    Args:
        prepared: Immutable child preparation JSON.
        label: Exactly one declared target label.
        expected_sha256: Caller-supplied preparation hash, frozen before launch.
    """
    run_reviewed_target(prepared, label, expected_sha256, 33)


def run_reviewed_target(
    prepared: Path, label: str, expected_sha256: str, expected_count: int
) -> None:
    """Execute one immutable reviewed target with the declared campaign size."""
    prepared = prepared.resolve(strict=True)
    cfg = verify_target(prepared, expected_sha256)
    if cfg["production_horizon_ps"] != 30.0 or len(cfg["windows"]) != expected_count:
        raise ValueError("Expectedapproved33state30pssupportprotocol")
    matches = [row for row in cfg["windows"] if row["label"] == label]
    if len(matches) != 1:
        raise ValueError("Select exactly one declared target label")
    row = matches[0]
    base = prepared.parent
    folder = base / label
    status = folder / "window_status.json"
    if status.exists():
        raise FileExistsError("Target already attempted; reviewed recovery required")
    topology_path = Path(cfg["topology"])
    topology = parmed.load_file(str(topology_path))
    initial = Rst7.open(str(base / row["seed"]))
    baseline = np.asarray(initial.coordinates).reshape(-1, 3).copy()
    active = np.asarray(
        list(AmberMask(topology, cfg["system"]["active_mask"]).Selected())
    )
    qm = list(AmberMask(topology, cfg["system"]["qm_atom_ids_mask"]).Selected())
    if len(baseline) != 104922 or len(active) != 4165 or len(qm) != 98:
        raise ValueError("Audited atom inventory changed")
    movable = np.zeros(len(baseline), dtype=bool)
    movable[active] = True
    restart = base / row["seed"]
    results: list[dict[str, Any]] = []
    try:
        for spec in row["stages"]:
            verify_target(prepared, expected_sha256)
            stage = (base / spec["input"]).parent
            if any(
                (stage / name).exists()
                for name in (
                    "process.json",
                    "input_identity.json",
                    "run.mdout",
                    "restart.ncrst",
                )
            ):
                raise FileExistsError(f"Stage already attempted: {stage}")
            write_json(status, dict(status="running", stage=spec["label"]))
            identity = {
                str(path): digest(path)
                for path in (
                    prepared,
                    topology_path,
                    restart,
                    stage / "run.mdin",
                    stage / "umbrella.RST",
                )
            }
            write_json(
                stage / "input_identity.json",
                dict(
                    sha256=identity,
                    environment=cfg["targets_plan"]["environment"],
                    seed=spec["seed"],
                ),
            )
            restart = _run_stage(
                stage,
                Path(cfg["targets_plan"]["sander"]),
                topology_path,
                restart,
                cfg["targets_plan"]["environment"],
            )
            receipt = native_acceptance(stage, spec["steps"], cfg["temperature_gate_K"])
            checkpoint, _ = _checkpoint(restart, baseline, movable, initial.box)
            receipt.update(checkpoint)
            write_json(stage / "acceptance.json", receipt)
            write_json(
                stage / "sampling_role.json",
                dict(
                    kind=spec["kind"],
                    pmf_included=spec["pmf_included"],
                    sampling_role=spec["sampling_role"],
                ),
            )
            analysis = ensure_diagnostics(cfg, row, stage, topology)
            results.append(
                dict(
                    stage=spec["label"],
                    native_accepted=True,
                    pmf_included=spec["pmf_included"],
                    diagnostics_completed=analysis,
                )
            )
            if (
                sum(item["pmf_included"] for item in results) == 2
                and spec["kind"] == "production"
            ):
                write_json(
                    folder / "diagnostic_checkpoint_10ps.json",
                    dict(
                        production_ps=10.0,
                        sampling_converged=False,
                        automatic_continuation_to_ps=30.0,
                        checkpoint=str(restart),
                        checkpoint_sha256=digest(restart),
                    ),
                )
        write_json(
            status,
            dict(
                status="completed",
                production_ps=30.0,
                sampling_converged=False,
                stages=results,
            ),
        )
    except Exception as error:
        write_json(
            status,
            dict(
                status="failed",
                stages=results,
                error=f"{type(error).__name__}: {error}",
            ),
        )
        raise
