"""Publish mapped diagnostic path distances, potential energies and VMD."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

from nyl12_covalent.path_energies import (
    driven_path_energies,
    write_driven_path_energy_data,
)
from scripts.export_direct_pt_progress import export, write_checksum_manifest


def publish_progress(
    args: argparse.Namespace,
    root: Path,
    coordinate: str = "asynchronous",
    drive_directory: str = "01_coupled_coordinate_drive",
) -> Path:
    """Publish every finalized window, preserving rejected endpoints as evidence."""
    status = json.loads((args.output / "status.json").read_text())
    bundle = export(
        args.output,
        args.topology,
        args.cpptraj,
        args.output / "vmd_progress",
        coordinate,
    )
    folders = [
        args.output / drive_directory / label for label in status["completed_windows"]
    ]
    release_name = (
        "02_unrestrained_ti1_endpoint"
        if coordinate == "asynchronous"
        else "02_unrestrained_endpoint"
    )
    release = args.output / release_name
    if status["state"] == "completed" and (release / "acceptance.json").is_file():
        folders.append(release)
    rows = []
    for folder in folders:
        acceptance = json.loads((folder / "acceptance.json").read_text())
        retained = (
            acceptance["chemistry_retained"]
            if "chemistry_retained" in acceptance
            else acceptance["proton_identity_retained"]
        )
        rows.append(
            {
                "window": folder.name,
                "stage": coordinate,
                "accepted": acceptance["native_converged"] and retained,
                "distances_A": acceptance["geometry_A"],
            }
        )
    keys = [key for key in rows[0]["distances_A"] if key != "proton_difference_A"]
    data = {
        "rows": rows,
        "pairs": {key: {"group": "Mapped reaction distances"} for key in keys},
        "jumps": [
            {
                "from": a["window"],
                "to": b["window"],
                "distance": key,
                "change_A": b["distances_A"][key] - a["distances_A"][key],
            }
            for a, b in zip(rows, rows[1:])
            for key in keys
            if abs(b["distances_A"][key] - a["distances_A"][key]) > 0.35
        ],
        "jump_threshold_A": 0.35,
    }
    distance_path = bundle / "reaction_distances.json"
    distance_path.write_text(json.dumps(data, indent=2) + "\n")
    energies = driven_path_energies(
        args.output / drive_directory,
        [
            row["window"]
            for row in rows
            if row["window"] in status["completed_windows"]
            and json.loads(
                (
                    args.output / drive_directory / row["window"] / "acceptance.json"
                ).read_text()
            )["native_converged"]
        ],
    )
    for energy in energies:
        energy["accepted"] = next(
            row["accepted"] for row in rows if row["window"] == energy["window"]
        )
    energy_data, _ = write_driven_path_energy_data(energies, bundle)
    for script, arguments in (
        ("plot_reaction_distances.py", [str(distance_path)]),
        (
            "plot_reaction_energies.py",
            [
                str(energy_data),
                str(bundle / "potential_energy.png"),
                "--title",
                "Diagnostic DFTB3/MM PES; first driven window zero; not PMF",
            ],
        ),
    ):
        subprocess.run(
            [str(args.plot_python), str(root / "scripts" / script), *arguments],
            check=True,
        )
    write_checksum_manifest(bundle)
    return bundle
