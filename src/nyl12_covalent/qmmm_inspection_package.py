"""Organize QM/MM inspection files by elementary stage and window number."""

from __future__ import annotations

import hashlib
import json
import shutil
import csv
import subprocess
from pathlib import Path
from typing import Any


def _copy(source: Path, target: Path, sources: dict[str, Any], root: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    pending = target.with_name("." + target.name + ".pending")
    shutil.copy2(source, pending)
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if hashlib.sha256(pending.read_bytes()).hexdigest() != digest:
        raise ValueError("Package copy digest mismatch")
    pending.replace(target)
    sources[str(target.relative_to(root))] = {"source": str(source), "sha256": digest}


WINDOW_FILES = {
    "acceptance.json": "geometry.json",
    "process.json": "run.json",
    "minimize.mdin": "optimize.mdin",
    "restraints.RST": "restraints.RST",
}


def _planned_windows(config: dict[str, Any], data: dict[str, Any]) -> dict[str, Any]:
    """Describe queued/running inputs without inventing optimized endpoints."""
    if "window_plan" not in config:
        return {}
    specification = config["window_plan"]
    plan = Path(specification["path"])
    stage = config["stage_names"][str(specification["stage"])]
    completed = {row["window"] for row in data["rows"] if row["accepted"]}
    windows = []
    for targets in json.loads(plan.read_text())["schedule"]:
        label = targets["label"]
        if Path(label).name != label:
            raise ValueError("Invalid planned window name")
        number = int(label.rsplit("_", 1)[1])
        folder = Path(config["run"]) / label
        files = {
            name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
            for name in WINDOW_FILES
            if (folder / name).is_file()
        }
        process = (
            json.loads((folder / "process.json").read_text())
            if "process.json" in files
            else None
        )
        if label in completed:
            state = "completed"
        elif process is not None:
            state = "running" if process.get("exit_status") is None else "not_accepted"
        else:
            state = "preparing" if folder.exists() else "queued"
        windows.append(
            {
                "original_window": label,
                "state": state,
                "directory": str(Path(stage) / f"window_{number:02d}"),
                "targets": targets,
                "available_source_sha256": files,
                "optimized_endpoint_available": label in completed,
            }
        )
    return {
        "plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
        "windows": windows,
    }


def build_qmmm_package(config_path: Path) -> Path:
    """Publish the latest verified movie, plots, scripts and numbered windows.

    Args:
        config_path: Explicit run, layout, stage naming and source file mapping.

    Returns:
        Human-readable package directory with a complete checksum manifest.
    """
    config = json.loads(config_path.read_text())
    independent = set(config.get("independent_packages", []))
    if any(Path(name).parts != (name,) or name in {".", ".."} for name in independent):
        raise ValueError("Independent packages must be named child directories")
    run = Path(config["run"])
    latest = (run / "live_inspection/latest").resolve(strict=True)
    data = json.loads((latest / "reaction_distances.json").read_text())
    planned = _planned_windows(config, data)
    destination = Path(config["destination"])
    destination.mkdir(parents=True, exist_ok=True)
    previous = destination / "package_sources.json"
    if (
        previous.exists()
        and json.loads(previous.read_text())["identity"] == data["identity"]
        and json.loads(previous.read_text())["layout_config_sha256"]
        == hashlib.sha256(config_path.read_bytes()).hexdigest()
        and json.loads(previous.read_text()).get("planned_windows", {}) == planned
    ):
        return destination
    sources: dict[str, Any] = {}
    for name, target in {
        "reaction_path.dcd": "vmd/reaction_path.dcd",
        "reaction_path.vmd": "vmd/view_path.vmd",
        "stripped.parm7.gz": "vmd/stripped.parm7.gz",
        "qm_indices.tcl": "vmd/qm_indices.tcl",
        "nyl12_es_vmd_common.tcl": "vmd/nyl12_es_vmd_common.tcl",
        "reaction_distances.png": "plots/distances.png",
        "reaction_distances.csv": "data/distances.csv",
        "reaction_distances.json": "data/distances.json",
    }.items():
        _copy(latest / name, destination / target, sources, destination)
    for target, source in config["source_files"].items():
        _copy(Path(source), destination / target, sources, destination)
    if "stage_measurements" in config:
        from nyl12_covalent.vmd_reactive_distances import append_reactive_labels

        append_reactive_labels(destination / "vmd/view_path.vmd", data["rows"], config)
        sources["vmd/view_path.vmd"]["source_sha256"] = sources["vmd/view_path.vmd"][
            "sha256"
        ]
        sources["vmd/view_path.vmd"]["transform"] = (
            "Append declared reactive distance labels"
        )
        sources["vmd/view_path.vmd"]["sha256"] = hashlib.sha256(
            (destination / "vmd/view_path.vmd").read_bytes()
        ).hexdigest()
    if "energy_plot" in config:
        from nyl12_covalent.path_energies import relay_endpoint_energies

        energies = relay_endpoint_energies(
            run, data["rows"], config["energy_distance_keys"]
        )
        energy_data = destination / "data/energies.json"
        energy_data.write_text(json.dumps(energies, indent=2) + "\n")
        columns = [key for key in energies[0] if key != "source_sha256"]
        with (destination / "data/energies.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(energies)
        subprocess.run(
            [
                config["energy_plot"]["python"],
                config["energy_plot"]["script"],
                str(energy_data),
                str(destination / "plots/energies.png"),
            ],
            check=True,
        )
    mapping = []
    for index, row in enumerate(data["rows"]):
        label = row["window"]
        stage = config["stage_names"][str(row["stage"])]
        window_number = int(label.rsplit("_", 1)[1])
        relative = Path(stage) / f"window_{window_number:02d}"
        original = run / label
        for source_name, target in {
            "acceptance.json": "geometry.json",
            "process.json": "run.json",
            "minimize.mdin": "optimize.mdin",
            "restraints.RST": "restraints.RST",
        }.items():
            _copy(
                original / source_name,
                destination / relative / target,
                sources,
                destination,
            )
        mapping.append(
            {
                "vmd_frame": index,
                "directory": str(relative),
                "original_window": label,
                "accepted_current_policy": row["accepted"],
            }
        )
    (destination / "window_order.json").write_text(json.dumps(mapping, indent=2) + "\n")
    for item in planned.get("windows", []):
        folder = destination / item["directory"]
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "status.json").write_text(json.dumps(item, indent=2) + "\n")
        for name in item["available_source_sha256"]:
            _copy(
                run / item["original_window"] / name,
                folder / WINDOW_FILES[name],
                sources,
                destination,
            )
    text = f"""# PA6: Asp272-assisted QM/MM path

Start here:
- `plots/distances.png`: all declared reaction distances versus window.
- `plots/energies.png`: relative energies in kcal/mol, including driving bias
  and a separately labeled bias-removed curve; data in `data/energies.*`.
- `vmd/view_path.vmd`: open in VMD, with matching topology and endpoint DCD.
  Forming/breaking bond distances are labeled in Å and switch by elementary stage.
- `data/distances.csv` and `.json`: measurements, jump flags and source hashes.
- `scripts/run_path.py`: exact current continuation controller used on violet.
- `scripts/initial_driver.py`: original controller, retained as history.
- `01_og1_to_asp272/window_00`, `window_01`, ...: forming OG1–C while
  transferring HG1 to Asp272.
- `02_asp272_to_nterm/window_00`, ...: subsequent transfer to terminal N;
  appears only if the first-stage chemical gate passes.

VMD frame order is in `window_order.json`; one frame is one optimized window,
not an elementary chemical event or physical time. A whole numbered stage is
an elementary-step hypothesis. Re-run the laptop pull and reload VMD for updates.
The current package contains {len(mapping)} finalized endpoint frames.
Planned/running windows have a status.json and available inputs before completion;
their coordinates enter the movie only after endpoint acceptance is recorded.

Original scientific run (unchanged): `{run}`.
Stage/window directories contain small inputs and original acceptance records;
full solvated restarts, scratch and engine trajectories stay on violet.
Original acceptance records remain intact, including earlier tracking stops.
The active policy is in `inputs/continuation_policy.json`; previous tracking
policy records remain in `inputs/previous_tracking_policy.json`.

Scripts are provenance copies with workstation-specific paths, not a portable
calculation environment. Do not rerun them on the laptop or overwrite an existing
run. The original two-step pathway is still a construction pilot, not a proven
mechanism or a measured barrier. No raw QM/MM calculation was relocated.
"""
    (destination / "README.md").write_text(text)
    record = {
        "identity": data["identity"],
        "run": str(run),
        "sources": sources,
        "layout_config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "planned_windows": planned,
    }
    # Commit the package identity only after all payloads and checksums are ready.
    pending_record = destination / ".package_sources.json.pending"
    pending_record.write_text(json.dumps(record, indent=2) + "\n")
    pending_record.replace(previous)
    lines = []
    for path in sorted(destination.rglob("*")):
        if (
            path.is_file()
            and path.name != "MANIFEST.sha256"
            and not path.name.endswith(".pending")
            and path.relative_to(destination).parts[0] not in independent
        ):
            lines.append(
                hashlib.sha256(path.read_bytes()).hexdigest()
                + "  "
                + str(path.relative_to(destination))
                + "\n"
            )
    pending_manifest = destination / ".MANIFEST.pending"
    pending_manifest.write_text("".join(lines))
    pending_manifest.replace(destination / "MANIFEST.sha256")
    return destination
