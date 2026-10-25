"""Bounded fresh-process finite-difference checks for the ES DFTB engine."""

from __future__ import annotations
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any
from . import es_dftb_engine as engine


def check_es_dftb_forces(config_path: Path, output_directory: Path) -> dict[str, Any]:
    """Compare analytic forces against central differences in fresh processes.

    Args:
        config_path: Approved JSON inventory, options, environment and tolerances.
        output_directory: Fresh retained evidence directory; must not exist.

    Returns:
        Numeric component comparisons and scientific acceptance, separately from
        each subprocess exit status. No minimization or dynamics is performed.

    Raises:
        ValueError: Invalid configuration, engine output, or source mutation.
        subprocess.CalledProcessError: Any isolated engine process fails.
    """
    config = json.loads(config_path.read_text())
    engine._validate(config)
    output_directory.mkdir(parents=True, exist_ok=False)
    source = output_directory / "source"
    source.mkdir()
    paths = [
        config_path,
        Path(engine.__file__),
        Path(config["python_executable"]),
        Path(config["topology"]),
        Path(config["coordinates"]),
    ]
    paths.append(Path(__file__))
    paths.extend(Path(path) for path in config["parameter_files"])
    records = [engine._record(path) for path in paths]
    engine._write(output_directory / "source_manifest.json", records)
    shutil.copy2(config_path, source / "declared_config.json")
    shutil.copy2(Path(__file__), source / "es_dftb_checks.py")
    worker = source / "es_dftb_engine.py"
    shutil.copy2(Path(engine.__file__), worker)
    env = os.environ.copy()
    env.update(config["environment"])

    def evaluate(name: str, perturbation: dict[str, Any] | None) -> dict[str, Any]:
        directory = output_directory / name
        directory.mkdir()
        request_path = directory / "request.json"
        engine._write(request_path, {"config": config, "perturbation": perturbation})
        command = [
            config["python_executable"],
            str(worker),
            "--worker",
            str(request_path),
        ]
        with (directory / "engine.log").open("w") as log:
            process = subprocess.run(
                command,
                cwd=directory,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        engine._write(
            directory / "process_status.json",
            {"exit_status": process.returncode, "command": command},
        )
        if process.returncode:
            raise subprocess.CalledProcessError(process.returncode, command)
        return json.loads((directory / "result.json").read_text())

    try:
        baseline = evaluate("baseline", None)
        repeated = evaluate("baseline_repeat", None)
        repeat_difference = abs(
            baseline["energy_kcal_per_mol"] - repeated["energy_kcal_per_mol"]
        )
        comparisons = []
        for index, row in enumerate(config["components"]):
            for step_index, step in enumerate(config["step_sizes_angstrom"]):
                energies = [
                    evaluate(
                        f"component_{index}_step_{step_index}_{name}",
                        {
                            "atom": row["atom"],
                            "axis": row["axis"],
                            "delta": sign * step,
                        },
                    )["energy_kcal_per_mol"]
                    for name, sign in (("plus", 1), ("minus", -1))
                ]
                finite_difference = -(energies[0] - energies[1]) / (2 * step)
                analytic = baseline["forces"][index]
                error = abs(analytic - finite_difference)
                tolerance = config["force_absolute_tolerance"] + config[
                    "force_relative_tolerance"
                ] * max(abs(analytic), abs(finite_difference))
                comparisons.append(
                    {
                        **row,
                        "step_angstrom": step,
                        "analytic_force": analytic,
                        "finite_difference_force": finite_difference,
                        "absolute_error": error,
                        "tolerance": tolerance,
                        "passed": math.isfinite(error) and error <= tolerance,
                    }
                )
        unchanged = all(engine._record(Path(row["path"])) == row for row in records)
        accepted = (
            unchanged
            and repeat_difference <= config["energy_reproducibility_tolerance"]
            and all(row["passed"] for row in comparisons)
        )
        result = {
            "accepted": accepted,
            "baseline": baseline,
            "components": comparisons,
            "energy_repeat_difference": repeat_difference,
            "sources_unchanged": unchanged,
            "scope": "Same-snapshot force consistency only; no mechanism or energy-model validation.",
        }
        engine._write(output_directory / "force_check.json", result)
        engine._write(
            output_directory / "scientific_acceptance.json", {"accepted": accepted}
        )
        return result
    except Exception as exc:
        engine._write(
            output_directory / "scientific_acceptance.json",
            {"accepted": False, "error": str(exc)},
        )
        raise
