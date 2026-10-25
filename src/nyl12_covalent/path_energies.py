"""Extract native relay energies and their declared harmonic driving bias."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any


def relay_endpoint_energies(
    run: Path, endpoints: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Read relative kcal/mol energies for the exact displayed endpoint sequence.

    Args:
        run: Existing Asp272 relay scientific run.
        endpoints: Finalized inspection rows, including explicitly rejected endpoints.

    Returns:
        Total and driving-bias-removed energies with source digests.
    """
    if not endpoints:
        raise ValueError("No finalized endpoints")
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
        q = d["OG_H"] - d["Asp_H"] if process["stage"] == 1 else d["Asp_H"] - d["N_HG"]
        blocks = re.findall(r"&rst\b(.*?)/", files[3].read_text(), flags=re.S | re.I)
        if len(blocks) != 2:
            raise ValueError("Expected two relay restraints")
        bias = 0.0
        for block, value, target in zip(
            blocks,
            (d["attack"], q),
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
