"""Prepare an isolated 10-to-30 ps continuation of four existing relay targets."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
from typing import Any

from .neb_umbrella_campaign import digest, write_json
from .relay_targeted import verify_hashes, verify_target
from .relay_umbrella_prepare import native_stage_input

LABELS = (
    "activation_021_from_022",
    "activation_022_from_021",
    "delivery_027_qdel_1p50",
    "delivery_028_qdel_1p75",
)


def verify_prefix(
    cfg: dict[str, Any], base: Path, row: dict[str, Any]
) -> dict[str, str]:
    """Require the complete accepted 10 ps target prefix and unchanged last restart."""
    pins = {}
    if [s["label"] for s in row["stages"]] != [
        "01_equilibration_1ps",
        "02_equilibration_4ps",
        "03_production_5ps",
        "04_production_5ps",
    ]:
        raise ValueError(
            "Target prefix is not the reviewed 5 ps equilibration plus10ps production"
        )
    for spec in row["stages"]:
        stage = (base / spec["input"]).parent
        acceptance = json.loads((stage / "acceptance.json").read_text())
        process = json.loads((stage / "process.json").read_text())
        if acceptance["accepted"] is not True or process["exit_status"] != 0:
            raise ValueError(f"Unaccepted native prefix:{stage}")
        if (
            acceptance["completed_steps"] != spec["steps"]
            or acceptance["every_integration_step_SCC_converged"] is not True
        ):
            raise ValueError("Native prefix integration/SCC acceptance changed")
        if digest(stage / "restart.ncrst") != acceptance["checkpoint_sha256"]:
            raise ValueError("Accepted prefix checkpoint changed")
        if (
            json.loads((stage / "analysis_outcome.json").read_text())["status"]
            != "completed"
        ):
            raise ValueError("Native prefix lacks completed full diagnostics")
        for name in (
            "acceptance.json",
            "process.json",
            "restart.ncrst",
            "analysis_outcome.json",
            "input_identity.json",
            "run.mdin",
            "umbrella.RST",
        ):
            path = stage / name
            pins[str(path.resolve())] = digest(path)
    return pins


def prepare_extension(prepared: Path, output: Path) -> Path:
    """Materialize four appended 5 ps stages per target without launching jobs.

    Args:
        prepared: Original immutable additional-sampling preparation.
        output: Exactly its child extension_20261127, which must not exist.

    Returns:
        New extension manifest; original preparation and status remain unchanged.
    """
    prepared = prepared.resolve(strict=True)
    cfg = verify_target(prepared, digest(prepared))
    base = prepared.parent
    output = output.absolute()
    if output != base / "extension_20261127" or output.exists():
        raise ValueError("Use the new logical child extension_20261127 only")
    if {r["label"] for r in cfg["windows"]} != set(LABELS):
        raise ValueError("Exactly the four reviewed targets are required")
    expected = {
        "active_atoms": 4165,
        "full_atoms": 104922,
        "qm_charge": -2,
        "qm_real_atoms": 98,
        "qm_method": "DFTB3/3ob-3-1",
    }
    if (
        any(cfg["system"][k] != v for k, v in expected.items())
        or cfg["temperature_K"] != 300
    ):
        raise ValueError("Target physical inventory/method differs from authorization")
    pins = {str(prepared): digest(prepared)}
    rows, new_files = [], {}
    old_seeds = {s["seed"] for r in cfg["windows"] for s in r["stages"]}
    for index, row in enumerate(cfg["windows"]):
        pins.update(verify_prefix(cfg, base, row))
        source = Path(cfg["targets_plan"]["source_mdin"]).read_text()
        stages = []
        for number in range(5, 9):
            stage = base / row["label"] / f"{number:02d}_production_5ps"
            if stage.exists():
                raise FileExistsError(stage)
            seed = 2026122000 + index * 10 + number
            if seed in old_seeds or not 0 < seed < 2147483647:
                raise ValueError(
                    "Thermostat continuation seeds must be unique positive int32"
                )
            mdin = native_stage_input(
                source,
                "production",
                seed,
                Path("umbrella.RST"),
                cfg["system"]["active_mask"],
                steps=10000,
            )
            if "irest=1, ntx=5" not in mdin or "dt=0.0005000" not in mdin:
                raise ValueError(
                    "Continuation must retain native velocities and0.5fs timestep"
                )
            old_rst = base / row["label"] / "04_production_5ps/umbrella.RST"
            new_files[stage / "run.mdin"] = mdin
            new_files[stage / "umbrella.RST"] = old_rst.read_text()
            stages.append(
                {
                    "label": stage.name,
                    "input": str((stage / "run.mdin").relative_to(base)),
                    "steps": 10000,
                    "seed": seed,
                    "kind": "production",
                    "pmf_included": True,
                    "sampling_role": "production",
                }
            )
        rows.append(
            {
                **row,
                "stages": [*row["stages"], *stages],
                "extension_stages": stages,
                "replica_id": row["label"],
                "continuation_of_production_ps": 10.0,
            }
        )
    output.mkdir()
    for path, text in new_files.items():
        path.parent.mkdir(exist_ok=True)
        with path.open("x") as handle:
            handle.write(text)
        pins[str(path.resolve())] = digest(path)
    overlay = {**cfg, "windows": rows, "production_horizon_ps": 30.0}
    overlay["frozen_files_sha256"] = {
        **cfg["frozen_files_sha256"],
        **{str(p.relative_to(base)): digest(p) for p in new_files},
    }
    snapshot = output / "source_snapshot"
    package = Path(__file__).resolve().parent
    shutil.copytree(
        package,
        snapshot / "src/nyl12_covalent",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    script = package.parents[1] / "scripts/extend_relay_target.py"
    (snapshot / "scripts").mkdir()
    shutil.copy2(script, snapshot / "scripts" / script.name)
    sources = {str(p): digest(p) for p in snapshot.rglob("*.py")}
    manifest = output / "manifest.json"
    write_json(
        manifest,
        {
            "base": str(base),
            "source_prepared": str(prepared),
            "source_prepared_sha256": digest(prepared),
            "cfg": overlay,
            "pinned_files_sha256": pins,
            "source_sha256": sources,
            "labels": list(LABELS),
            "total_production_ps": 30.0,
            "extension_production_ps": 20.0,
            "independent_replicas_are_not_continuous": True,
            "execution_entrypoint": str(snapshot / "scripts" / script.name),
        },
    )
    manifest.with_suffix(".sha256").write_text(digest(manifest) + "\n")
    return manifest


def verify_extension(manifest: Path, expected_sha256: str) -> dict[str, Any]:
    """Require immutable scientific, restart, source and extension inputs."""
    if digest(manifest) != expected_sha256:
        raise ValueError("Extension manifest hash changed")
    plan = json.loads(manifest.read_text())
    if digest(Path(plan["source_prepared"])) != plan["source_prepared_sha256"]:
        raise ValueError("Original preparation changed")
    verify_hashes(plan["source_sha256"])
    verify_hashes(plan["pinned_files_sha256"])
    verify_hashes(plan["cfg"]["external_files_sha256"])
    return plan
