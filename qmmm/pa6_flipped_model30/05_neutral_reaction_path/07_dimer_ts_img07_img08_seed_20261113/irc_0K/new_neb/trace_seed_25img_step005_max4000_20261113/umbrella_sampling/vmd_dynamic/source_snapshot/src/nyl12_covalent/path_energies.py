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

from nyl12_covalent.direct_acyl_bias import proton_guard_bias


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
        values = [d[distance_keys["attack"]], q]
        targets = [
            process["attack_target_A"],
            process["proton_difference_target_A"],
        ]
        expected_blocks = 2
        if process["stage"] == 2:
            expected_blocks = 3
            values.append(d[distance_keys["donor_h"]] - d[distance_keys["relay_h"]])
            try:
                targets.append(process["relay_retention_target_A"])
            except KeyError as error:
                raise ValueError(
                    "Stage II lacks the audited relay-retention target"
                ) from error
        if len(blocks) != expected_blocks:
            raise ValueError(
                f"Expected {expected_blocks} restraints for stage {process['stage']}"
            )
        bias = 0.0
        for block, value, target in zip(
            blocks,
            values,
            targets,
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


def driven_path_energies(run: Path, window_names: list[str]) -> list[dict[str, Any]]:
    """Extract reported total and restraint-removed energies for a driven path.

    Args:
        run: Directory containing one subdirectory per coordinate-drive window.
        window_names: Ordered window directory names.

    Returns:
        Auditable energy rows, referenced to the first window.

    Raises:
        ValueError: If a window is missing, unconverged, or lacks one reported
            total-energy/restraint-energy endpoint.
    """
    if not window_names:
        raise ValueError("No driven-path windows supplied")
    rows: list[dict[str, Any]] = []
    for window_name in window_names:
        folder = run / window_name
        files = {
            "mdout": folder / "minimized.mdout",
            "acceptance": folder / "acceptance.json",
            "restraint": folder / "restraints.RST",
        }
        if any(not path.is_file() for path in files.values()):
            raise ValueError(f"Incomplete driven-path window: {window_name}")
        acceptance = json.loads(files["acceptance"].read_text())
        if acceptance.get("native_converged") is not True:
            raise ValueError(f"Unconverged driven-path window: {window_name}")
        output_text = files["mdout"].read_text()
        total_matches = re.findall(r"Etot\s*=\s*([-+\d.EeDd]+)", output_text)
        restraint_matches = re.findall(r"RESTRAINT\s*=\s*([-+\d.EeDd]+)", output_text)
        if not total_matches or not restraint_matches:
            raise ValueError(f"Missing final energy record: {window_name}")
        total = float(total_matches[-1].replace("D", "E").replace("d", "e"))
        restraint = float(restraint_matches[-1].replace("D", "E").replace("d", "e"))
        unbiased = total - restraint
        guard = proton_guard_bias(
            files["restraint"].read_text(),
            acceptance.get("geometry_A", {}).get("substrate_N_HG1_A", float("nan")),
        )
        components = {}
        if guard is not None:
            if not math.isfinite(guard) or guard > restraint + 0.001:
                raise ValueError(f"Invalid substrate guard energy: {window_name}")
            components = {
                "substrate_proton_guard_bias_kcal_mol": guard,
                "composite_bias_kcal_mol": restraint - guard,
            }
        if not all(math.isfinite(value) for value in (total, restraint, unbiased)):
            raise ValueError(f"Nonfinite driven-path energy: {window_name}")
        rows.append(
            {
                "window": window_name,
                "stage": 1,
                "accepted": True,
                "total_kcal_mol": total,
                "driving_bias_kcal_mol": restraint,
                **components,
                "minus_driving_kcal_mol": unbiased,
                "source_sha256": {
                    str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in files.values()
                },
            }
        )
    reference_total = rows[0]["total_kcal_mol"]
    reference_unbiased = rows[0]["minus_driving_kcal_mol"]
    for row in rows:
        row["relative_total_kcal_mol"] = row["total_kcal_mol"] - reference_total
        row["relative_minus_driving_kcal_mol"] = (
            row["minus_driving_kcal_mol"] - reference_unbiased
        )
    return rows


def write_driven_path_energy_data(
    rows: list[dict[str, Any]], output: Path
) -> tuple[Path, Path]:
    """Write driven-path energy JSON and CSV files.

    Args:
        rows: Output from :func:`driven_path_energies`.
        output: Existing output directory.

    Returns:
        Paths to the JSON and CSV files.
    """
    if not rows:
        raise ValueError("No driven-path energy rows supplied")
    json_path = output / "potential_energy.json"
    csv_path = output / "potential_energy.csv"
    json_path.write_text(json.dumps(rows, indent=2) + "\n")
    columns = list(
        dict.fromkeys(key for row in rows for key in row if key != "source_sha256")
    )
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return json_path, csv_path
