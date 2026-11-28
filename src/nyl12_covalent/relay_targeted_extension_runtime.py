"""Run one reviewed target's four-stage continuation, leaving its old status intact."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from parmed.amber import Rst7
from parmed.amber.mask import AmberMask

from .neb_umbrella_campaign import digest, write_json
from .neb_umbrella_runtime import _run_stage
from .neb_umbrella_segments import _checkpoint
from .relay_targeted_extension import verify_extension
from .relay_umbrella_runtime import ensure_diagnostics, native_acceptance


def run_extension(manifest: Path, label: str, expected_sha256: str) -> None:
    """Append20ps from the accepted04restart, preserving velocities and replica identity.

    Args:
        manifest: Frozen isolated continuation manifest.
        label: Exactly one declared existing target label.
        expected_sha256: Explicit launch-time manifest hash.
    """
    plan = verify_extension(manifest, expected_sha256)
    cfg, base = plan["cfg"], Path(plan["base"])
    rows = [r for r in cfg["windows"] if r["label"] == label]
    if len(rows) != 1:
        raise ValueError("Select exactly one reviewed continuation label")
    row = rows[0]
    status = manifest.parent / "status" / f"{label}.json"
    if status.exists():
        raise FileExistsError("Extension already attempted; explicit recovery required")
    status.parent.mkdir(exist_ok=True)
    topology_path = Path(cfg["topology"])
    topology = parmed.load_file(str(topology_path))
    initial = Rst7.open(str(base / row["seed"]))
    baseline = np.asarray(initial.coordinates).reshape(-1, 3).copy()
    active = list(AmberMask(topology, cfg["system"]["active_mask"]).Selected())
    qm = list(AmberMask(topology, cfg["system"]["qm_atom_ids_mask"]).Selected())
    if baseline.shape != (104922, 3) or len(active) != 4165 or len(qm) != 98:
        raise ValueError("Native continuation inventory differs from104922/4165/98")
    movable = np.zeros(len(baseline), dtype=bool)
    movable[active] = True
    restart = base / label / "04_production_5ps/restart.ncrst"
    previous, _ = _checkpoint(restart, baseline, movable, initial.box)
    if (
        previous["checkpoint_sha256"]
        != json.loads((restart.parent / "acceptance.json").read_text())[
            "checkpoint_sha256"
        ]
    ):
        raise ValueError("Accepted04restart changed")
    records: list[dict[str, Any]] = []
    try:
        for spec in row["extension_stages"]:
            verify_extension(manifest, expected_sha256)
            stage = (base / spec["input"]).parent
            if any(
                (stage / name).exists()
                for name in (
                    "process.json",
                    "input_identity.json",
                    "run.mdout",
                    "restart.ncrst",
                    "input.ncrst",
                )
            ):
                raise FileExistsError(f"Extension stage already attempted:{stage}")
            write_json(
                status, {"status": "running", "stage": spec["label"], "stages": records}
            )
            identity = {
                str(p): digest(p)
                for p in (
                    manifest,
                    topology_path,
                    restart,
                    stage / "run.mdin",
                    stage / "umbrella.RST",
                )
            }
            write_json(
                stage / "input_identity.json",
                {
                    "sha256": identity,
                    "environment": cfg["targets_plan"]["environment"],
                    "seed": spec["seed"],
                    "replica_id": label,
                    "velocity_policy": "nativeirest1ntx5continuation",
                },
            )
            restart = _run_stage(
                stage,
                Path(cfg["targets_plan"]["sander"]),
                topology_path,
                restart,
                cfg["targets_plan"]["environment"],
            )
            receipt = native_acceptance(stage, 10000, cfg["temperature_gate_K"])
            checkpoint, _ = _checkpoint(restart, baseline, movable, initial.box)
            receipt.update(checkpoint)
            write_json(stage / "acceptance.json", receipt)
            write_json(
                stage / "sampling_role.json",
                {
                    "kind": "production",
                    "pmf_included": True,
                    "sampling_role": "production",
                    "replica_id": label,
                    "continuation_of": "03_production_5ps+04_production_5ps",
                },
            )
            diagnostics = ensure_diagnostics(cfg, row, stage, topology)
            records.append(
                {
                    "stage": spec["label"],
                    "native_accepted": True,
                    "diagnostics_completed": diagnostics,
                    "checkpoint_sha256": digest(restart),
                }
            )
            if diagnostics is not True:
                raise RuntimeError(
                    "Full native bias/proton diagnostics failed; stop continuation"
                )
        write_json(
            status,
            {
                "status": "completed",
                "production_ps_total": 30.0,
                "extension_ps": 20.0,
                "replica_id": label,
                "stages": records,
                "native_and_diagnostics_accepted": True,
                "sampling_converged": False,
            },
        )
    except Exception as error:
        write_json(
            status,
            {
                "status": "failed",
                "stages": records,
                "replica_id": label,
                "error": f"{type(error).__name__}: {error}",
            },
        )
        raise
