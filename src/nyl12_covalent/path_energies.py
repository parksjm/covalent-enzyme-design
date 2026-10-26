"""Extract native relay energies and their declared harmonic driving bias."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any


def relay_endpoint_energies(
    run: Path,
    endpoints: list[dict[str, Any]],
    distance_keys: dict[str, str],
) -> list[dict[str, Any]]:
    """Read relative kcal/mol energies for the exact displayed endpoint sequence.

    Args:
        run: Existing Asp272 relay scientific run.
        endpoints: Finalized inspection rows, including explicitly rejected endpoints.
        distance_keys: Explicit acceptance-record keys for attack, donor-H,
            relay-H, and acceptor-H distances.

    Returns:
        Total and driving-bias-removed energies with source digests.
    """
    if not endpoints:
        raise ValueError("No finalized endpoints")
    expected_keys = {"attack", "donor_h", "relay_h", "acceptor_h"}
    if set(distance_keys) != expected_keys:
        raise ValueError("Energy distance-key mapping is incomplete or ambiguous")
    rows = []
    for endpoint in endpoints:
        folder = run / endpoint["window"]
        files = [
            folder / name
            for name in (
                "minimized.mdout",
                "process.json",
                "acceptance.json",
                "restraints.RST",
            )
        ]
        matches = re.findall(r"Etot\s*=\s*([-+\d.EeDd]+)", files[0].read_text())
        if not matches:
            raise ValueError("Missing native DL-FIND energy")
        total = float(matches[-1].replace("D", "E").replace("d", "e"))
        process = json.loads(files[1].read_text())
        if process["stage"] not in (1, 2):
            raise ValueError("Unknown relay stage")
        d = json.loads(files[2].read_text())["distances_A"]
        q = (
            d[distance_keys["donor_h"]] - d[distance_keys["relay_h"]]
            if process["stage"] == 1
            else d[distance_keys["relay_h"]] - d[distance_keys["acceptor_h"]]
        )
        blocks = re.findall(r"&rst\b(.*?)/", files[3].read_text(), flags=re.S | re.I)
        if len(blocks) != 2:
            raise ValueError("Expected two relay restraints")
        bias = 0.0
        for block, value, target in zip(
            blocks,
            (d[distance_keys["attack"]], q),
            (process["attack_target_A"], process["proton_difference_target_A"]),
            strict=True,
        ):
            parameters = {
                k: float(v)
                for k, v in re.findall(r"\b(r[1-4]|rk[23])\s*=\s*([-+\d.eE]+)", block)
            }
            if (
                set(parameters) != {"r1", "r2", "r3", "r4", "rk2", "rk3"}
                or parameters["rk2"] != 100
                or parameters["rk3"] != 100
                or abs(parameters["r2"] - target) > 1e-8
                or abs(parameters["r3"] - target) > 1e-8
                or not parameters["r1"] <= value <= parameters["r4"]
            ):
                raise ValueError("Restraint differs from audited relay harmonic model")
            bias += 100 * (value - target) ** 2
        if not math.isfinite(total + bias):
            raise ValueError("Nonfinite energy")
        rows.append(
            {
                "window": endpoint["window"],
                "stage": process["stage"],
                "accepted": endpoint["accepted"],
                "total_kcal_mol": total,
                "driving_bias_kcal_mol": bias,
                "minus_driving_kcal_mol": total - bias,
                "source_sha256": {
                    str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files
                },
            }
        )
    for row in rows:
        row["relative_total_kcal_mol"] = (
            row["total_kcal_mol"] - rows[0]["total_kcal_mol"]
        )
        row["relative_minus_driving_kcal_mol"] = (
            row["minus_driving_kcal_mol"] - rows[0]["minus_driving_kcal_mol"]
        )
    return rows


def write_energy_artifacts(
    run: Path,
    endpoints: list[dict[str, Any]],
    output: Path,
    plot_python: str,
    plot_script: Path,
    environment: dict[str, str],
    distance_keys: dict[str, str],
) -> None:
    """Write endpoint energy data and its rendered plot into an inspection bundle.

    Args:
        run: Existing reaction-path run directory.
        endpoints: Exact displayed endpoint sequence.
        output: Inspection generation receiving JSON, CSV, and PNG artifacts.
        plot_python: Explicit Python executable used for plotting.
        plot_script: Versioned energy-plot entry point.
        environment: Explicit plotting environment.
        distance_keys: Explicit acceptance-record distance mapping.
    """
    energies = relay_endpoint_energies(run, endpoints, distance_keys)
    json_path = output / "reaction_energies.json"
    json_path.write_text(json.dumps(energies, indent=2) + "\n")
    columns = [key for key in energies[0] if key != "source_sha256"]
    with (output / "reaction_energies.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(energies)
    subprocess.run(
        [
            plot_python,
            str(plot_script),
            str(json_path),
            str(output / "reaction_energies.png"),
        ],
        env=environment,
        check=True,
    )
