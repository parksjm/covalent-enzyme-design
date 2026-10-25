"""Recompute representative enzyme force evidence from retained worker artifacts."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

_PROBE_LABELS = {
    "Thr_OG1",
    "PA6_C54",
    "PA6_O49",
    "QM_link_parent",
    "MM_link_partner",
    "external_MM_water_O",
}


def _same(left: float, right: float) -> bool:
    return (
        math.isfinite(left)
        and math.isfinite(right)
        and math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-10)
    )


def check_force_evidence(engine: dict, check: dict, directory: Path) -> list[Path]:
    """Verify all six probes at both declared steps against saved energies.

    Args:
        engine: Source-pinned declared engine configuration.
        check: Report claiming accepted representative enzyme forces.
        directory: Preserved force-check directory with isolated worker evidence.

    Returns:
        Evidence paths to hash again before any later engine execution.

    Raises:
        ValueError: Missing coverage, altered workers or inconsistent numerics.
    """
    probes, steps = engine["components"], engine["step_sizes_angstrom"]
    if len(probes) != 6 or {row["label"] for row in probes} != _PROBE_LABELS:
        raise ValueError("Exactly the six declared enzyme probes are required")
    if len({(row["atom"], row["axis"]) for row in probes}) != 6 or any(
        type(row["atom"]) is not int
        or row["atom"] < 1
        or type(row["axis"]) is not int
        or row["axis"] not in (0, 1, 2)
        for row in probes
    ):
        raise ValueError("Probe atom indices and axes must be distinct and valid")
    if (
        len(steps) != 2
        or len(set(steps)) != 2
        or any(
            isinstance(step, bool) or not math.isfinite(step) or step <= 0
            for step in steps
        )
    ):
        raise ValueError("Exactly two explicit distinct finite positive steps required")
    absolute, relative, repeat_tol = [
        engine[key]
        for key in [
            "force_absolute_tolerance",
            "force_relative_tolerance",
            "energy_reproducibility_tolerance",
        ]
    ]
    if (
        not np.isfinite([absolute, relative, repeat_tol]).all()
        or min(absolute, relative, repeat_tol) < 0
    ):
        raise ValueError("Finite nonnegative declared acceptance tolerances required")
    planned = {
        (row["atom"], row["axis"], step): (i, j)
        for i, row in enumerate(probes)
        for j, step in enumerate(steps)
    }
    rows = check["components"]
    keys = [(row["atom"], row["axis"], row["step_angstrom"]) for row in rows]
    if len(keys) != 12 or len(set(keys)) != 12 or set(keys) != set(planned):
        raise ValueError(
            "Force evidence must cover every declared probe and step exactly once"
        )
    paths: list[Path] = []

    def worker(name: str, perturbation: dict | None) -> dict[str, Any]:
        files = [
            directory / name / filename
            for filename in ["request.json", "result.json", "process_status.json"]
        ]
        request, result, process = [json.loads(path.read_text()) for path in files]
        if process["exit_status"] != 0 or type(process["exit_status"]) is not int:
            raise ValueError(
                "Every isolated force worker must have exited successfully"
            )
        if (
            request["config"] != engine
            or request["perturbation"] != perturbation
            or result["perturbation"] != perturbation
        ):
            raise ValueError(
                "Worker configuration or displacement differs from declared proof"
            )
        if (
            result["inventory"] != check["baseline"]["inventory"]
            or result["box"] != check["baseline"]["box"]
        ):
            raise ValueError("Worker inventory or fixed box differs from baseline")
        if not math.isfinite(result["energy_kcal_per_mol"]):
            raise ValueError("Nonfinite saved worker energy")
        paths.extend(files)
        return result

    baseline = worker("baseline", None)
    if baseline != check["baseline"]:
        raise ValueError("Reported baseline differs from preserved worker result")
    repeated = worker("baseline_repeat", None)
    repeat = abs(baseline["energy_kcal_per_mol"] - repeated["energy_kcal_per_mol"])
    if not _same(repeat, check["energy_repeat_difference"]) or repeat > repeat_tol:
        raise ValueError("Baseline energy repeat check failed")
    analytic = baseline["forces"]
    if len(analytic) != 6 or not np.isfinite(analytic).all():
        raise ValueError("Six finite baseline analytic forces required")
    force_path = directory / "baseline/analytic_forces.npy"
    force_array = np.load(force_path, allow_pickle=False)
    if (
        force_array.ndim != 2
        or force_array.shape[1] != 3
        or not np.isfinite(force_array).all()
    ):
        raise ValueError("Baseline full force array is invalid")
    paths.append(force_path)
    for row, key in zip(rows, keys, strict=True):
        i, j = planned[key]
        probe, step = probes[i], steps[j]
        if row["label"] != probe["label"] or row["passed"] is not True:
            raise ValueError("Reported probe identity or acceptance mismatch")
        if probe["atom"] > len(force_array) or not _same(
            force_array[probe["atom"] - 1, probe["axis"]], analytic[i]
        ):
            raise ValueError(
                "Reported analytic force differs from preserved full force array"
            )
        energies = [
            worker(
                f"component_{i}_step_{j}_{name}",
                {"atom": probe["atom"], "axis": probe["axis"], "delta": sign * step},
            )["energy_kcal_per_mol"]
            for name, sign in [("plus", 1), ("minus", -1)]
        ]
        fd = -(energies[0] - energies[1]) / (2 * step)
        error = abs(analytic[i] - fd)
        tolerance = absolute + relative * max(abs(analytic[i]), abs(fd))
        expected = [analytic[i], fd, error, tolerance]
        observed = [
            row[key]
            for key in [
                "analytic_force",
                "finite_difference_force",
                "absolute_error",
                "tolerance",
            ]
        ]
        if (
            not all(
                _same(left, right)
                for left, right in zip(expected, observed, strict=True)
            )
            or error > tolerance
        ):
            raise ValueError(
                "Stored force comparison fails independently recomputed energy differences"
            )
    return paths
