"""Publish traceable N-terminal relay PES and molecular-viewer progress."""

from __future__ import annotations

import csv
import json
import logging
import math
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import parmed

from nyl12_covalent.stationary_native import write_json
from nyl12_covalent.vmd_export_io import sha, validate_dcd

LOGGER = logging.getLogger(__name__)
STRIP_NAMES = {"WAT", "Na+", "Cl-", "NA", "CL"}


def _input_paths(campaign: Path, config: dict[str, Any]) -> tuple[Path, Path, Path]:
    topology = (campaign / config["model"] / "inputs/system.prmtop").resolve(
        strict=True
    )
    template = (campaign / config["viewer_template"]).resolve(strict=True)
    cpptraj = Path(config["cpptraj"])
    if not cpptraj.is_absolute() or sha(topology) != config["source_topology_sha256"]:
        raise ValueError("Source topology hash does not match configuration")
    cpptraj = cpptraj.resolve(strict=True)
    manifest = json.loads((template / "manifest.json").read_text())
    if manifest["topology_sha256"] != config["source_topology_sha256"]:
        raise ValueError("Viewer template source topology does not match")
    return topology, template, cpptraj


def _states(topology: Path, template: Path, rows: list[dict[str, Any]]) -> list:
    parm = parmed.load_file(str(topology))
    manifest = json.loads((template / "manifest.json").read_text())
    mapping = manifest["stripped_to_original_atom_indices_0based"]
    expected = [a.idx for a in parm.atoms if a.residue.name not in STRIP_NAMES]
    if mapping != expected:
        raise ValueError("Viewer mapping differs from declared strip selection")
    stripped = parmed.load_file(str(template / "stripped.parm7"))
    if len(stripped.atoms) != len(mapping):
        raise ValueError("Viewer topology atom count differs from mapping")
    result = []
    for row in rows:
        path = Path(row["restart"])
        if not path.is_absolute():
            raise ValueError("Restart path must be absolute")
        restart = parmed.amber.Rst7.open(str(path.resolve(strict=True)))
        xyz = np.asarray(restart.coordinates).reshape(-1, 3)
        box = np.asarray(restart.box)
        if xyz.shape != (len(parm.atoms), 3) or box.shape != (6,):
            raise ValueError(
                "Restart coordinates or periodic box do not match topology"
            )
        if not np.isfinite(xyz).all() or not np.isfinite(box).all() or (box <= 0).any():
            raise ValueError("Restart has invalid coordinates or box")
        result.append((xyz[mapping], box))
    return result


def _export(
    branch: Path,
    topology: Path,
    template: Path,
    cpptraj: Path,
    rows: list[dict[str, Any]],
) -> dict[str, Any]:
    states = _states(topology, template, rows)
    viewer = branch / "vmd_progress"
    viewer.mkdir(parents=True, exist_ok=True)
    shutil.copy2(template / "stripped.parm7", viewer / "stripped.parm7")
    script = (template / "view_relay.vmd").read_text()
    if "optimization.dcd" not in script:
        raise ValueError("Viewer template has no expected DCD reference")
    (viewer / "view_path.vmd").write_text(
        script.replace("optimization.dcd", "reaction_path.dcd")
    )
    pending = viewer / "reaction_path.pending.dcd"
    for path in [topology, pending, *[Path(row["restart"]) for row in rows]]:
        if any(c in str(path) for c in ("\n", "\r", '"')):
            raise ValueError("Unsafe cpptraj path")
    commands = [f'parm "{topology}"']
    commands.extend(f'trajin "{row["restart"]}"' for row in rows)
    commands += ["strip :WAT,Na+,Cl-,NA,CL", f'trajout "{pending}" dcd namdcell', "run"]
    script_path = viewer / "export.cpptraj"
    script_path.write_text("\n".join(commands) + "\n")
    result = subprocess.run(
        [str(cpptraj), "-i", str(script_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    (viewer / "cpptraj.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"cpptraj failed with exit {result.returncode}")
    checked = validate_dcd(pending, states)
    os.replace(pending, viewer / "reaction_path.dcd")
    return checked


def _plot(branch: Path, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    directory = branch / "plots"
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    keys = [f"{h}_{a}" for h in ("H1", "H2", "HG", "HN") for a in ("N", "OG", "NL")]
    keys += ["C_OG", "C_NL", "C_O"]
    for index, row in enumerate(rows):
        d = row["geometry"]["distances_A"]
        delta = (
            {
                key: float(d[key] - rows[index - 1]["geometry"]["distances_A"][key])
                for key in keys
            }
            if index
            else {}
        )
        records.append(
            {
                "label": row["label"],
                "distance_changes_A": delta,
                "large_signed_changes_A": {
                    k: v for k, v in delta.items() if abs(v) > 0.35
                },
            }
        )
    x = np.arange(len(rows))
    baseline = float(rows[0]["physical_kcal_mol"])
    fig, ax = plt.subplots()
    for key, label in (
        ("physical_kcal_mol", "physical"),
        ("total_kcal_mol", "biased total"),
    ):
        ax.plot(x, [row[key] - baseline for row in rows], marker="o", label=label)
    ax.set(
        xlabel="Path step",
        ylabel="Energy relative to start (kcal/mol)",
        title="PES: restrained minima, not a PMF",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(directory / "energies.png")
    plt.close(fig)
    fig, axes = plt.subplots(5, 1, figsize=(10, 15), sharex=True)
    for ax, h in zip(axes[:4], ("H1", "H2", "HG", "HN"), strict=True):
        for atom in ("N", "OG", "NL"):
            key = f"{h}_{atom}"
            ax.plot(x, [r["geometry"]["distances_A"][key] for r in rows], label=key)
        ax.legend(ncol=3)
        ax.set_ylabel("Å")
    for key in ("C_OG", "C_NL", "C_O"):
        axes[4].plot(x, [r["geometry"]["distances_A"][key] for r in rows], label=key)
    axes[4].legend(ncol=3)
    for index, record in enumerate(records):
        if record["large_signed_changes_A"]:
            for ax in axes:
                ax.axvline(index, color="red", alpha=0.35, linestyle="--")
    axes[4].set(xlabel="Path step", ylabel="Å")
    fig.tight_layout()
    fig.savefig(directory / "distances.png")
    plt.close(fig)
    return records


def publish(campaign: Path, branch_root: Path, rows: list[dict[str, Any]]) -> None:
    """Publish validated restart trajectory, energies and proton-distance progress.

    Args:
        campaign: Campaign directory containing config.json and source model.
        branch_root: Destination for this branch's derived artifacts.
        rows: Ordered start and accepted or attempted native-window records.

    Raises:
        ValueError: If input identity, mapping or source coordinates fail validation.
        RuntimeError: If cpptraj cannot export the requested trajectory.
    """
    if not rows or rows[0]["label"] != "window033_start":
        raise ValueError("First row must be window033_start")
    config = json.loads((campaign / "config.json").read_text())
    if not math.isfinite(float(config["reference_energy_hartree"])):
        raise ValueError("Reference energy must be finite")
    branch_root.mkdir(parents=True, exist_ok=True)
    topology, template, cpptraj = _input_paths(campaign, config)
    validation = _export(branch_root, topology, template, cpptraj, rows)
    records = _plot(branch_root, rows)
    write_json(
        branch_root / "progress.json",
        {
            "rows": rows,
            "distance_changes": records,
            "dcd_validation": validation,
            "reference_energy_hartree": config["reference_energy_hartree"],
        },
    )
    fields = [
        "label",
        "restart",
        "accepted",
        "native_converged",
        "target_A",
        "observed_coordinate_A",
        "total_kcal_mol",
        "bias_kcal_mol",
        "physical_kcal_mol",
    ]
    csv_path = branch_root / "progress.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("Published %d relay path frames to %s", len(rows), branch_root)
