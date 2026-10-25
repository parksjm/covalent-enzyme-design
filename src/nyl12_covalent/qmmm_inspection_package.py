"""Organize QM/MM inspection files by elementary stage and window number."""

from __future__ import annotations

import hashlib
import json
import shutil
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


def build_qmmm_package(config_path: Path) -> Path:
    """Publish the latest verified movie, plots, scripts and numbered windows.

    Args:
        config_path: Explicit run, layout, stage naming and source file mapping.

    Returns:
        Human-readable package directory with a complete checksum manifest.
    """
    config = json.loads(config_path.read_text())
    run = Path(config["run"])
    latest = (run / "live_inspection/latest").resolve(strict=True)
    data = json.loads((latest / "reaction_distances.json").read_text())
    destination = Path(config["destination"])
    destination.mkdir(parents=True, exist_ok=True)
    previous = destination / "package_sources.json"
    if (
        previous.exists()
        and json.loads(previous.read_text())["identity"] == data["identity"]
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
    text = f"""# PA6: Asp272-assisted QM/MM path

Start here:
- `plots/distances.png`: all declared reaction distances versus window.
- `vmd/view_path.vmd`: open in VMD, with matching topology and endpoint DCD.
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

Original scientific run (unchanged): `{run}`.
Stage/window directories contain small inputs and original acceptance records;
full solvated restarts, scratch and engine trajectories stay on violet.
Window10 was subsequently accepted under the documented 0.20 Å tracking policy;
its original 0.15 Å failure record remains intact. See `inputs/continuation_policy.json`.

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
