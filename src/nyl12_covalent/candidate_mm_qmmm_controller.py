"""Sequential, unrestrained MM relaxation and gated QM/MM handoff."""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any, Callable, Sequence

from nyl12_covalent.es_classical_relaxation import prepare_relaxation
from nyl12_covalent.ti1_checks import check_convergence

LOGGER = logging.getLogger(__name__)
DRMS_TARGET = 0.01

EngineRunner = Callable[[Sequence[str]], int]
EndpointGate = Callable[[Path], dict[str, Any]]
QmmmLauncher = Callable[[Path, Path], dict[str, Any]]


class CandidateControllerError(RuntimeError):
    """Raised when the candidate fails a declared process or science gate."""


def run_candidate(
    topology: Path,
    initial_coordinates: Path,
    engine: Path,
    output: Path,
    case_label: str,
    segment_cycles: int,
    total_cycle_ceiling: int,
    endpoint_gate: EndpointGate,
    qmmm_launcher: QmmmLauncher,
    engine_runner: EngineRunner | None = None,
) -> dict[str, Any]:
    """Relax one candidate in immutable MM segments, then conditionally run QM/MM.

    Args:
        topology: Validated Amber topology for the selected candidate.
        initial_coordinates: Candidate restart used for the first segment.
        engine: Path to the intended ``pmemd.cuda`` executable.
        output: New descriptive directory for this candidate attempt.
        case_label: Human-readable candidate identifier.
        segment_cycles: Maximum cycles in each immutable MM attempt.
        total_cycle_ceiling: Sum of segment cycle budgets allowed for this run.
        endpoint_gate: Existing Tier A/B and native Thr233 acceptance check.
        qmmm_launcher: Existing QM/MM preparation and launch callback.
        engine_runner: Optional engine adapter, primarily for tests.

    Returns:
        Controller report with MM segments and the conditional QM/MM result.

    Raises:
        CandidateControllerError: On engine failure, invalid convergence output,
            cycle-ceiling exhaustion, or a failed scientific endpoint gate.
        FileExistsError: If the output directory already exists.
        ValueError: If cycle budgets or the case label are invalid.
    """
    if segment_cycles < 1 or total_cycle_ceiling < 1:
        raise ValueError("Cycle budgets must be positive")
    if not case_label.strip():
        raise ValueError("case_label must be nonempty")
    runner = engine_runner or _run_engine
    output.mkdir(parents=True, exist_ok=False)
    segments_root = output / "mm_segments"
    segments_root.mkdir()
    report: dict[str, Any] = {
        "case_label": case_label,
        "drms_target_kcal_mol_A": DRMS_TARGET,
        "segment_cycle_limit": segment_cycles,
        "total_cycle_ceiling": total_cycle_ceiling,
        "segments": [],
        "status": "mm_in_progress",
        "qmmm_launched": False,
    }
    coordinates = initial_coordinates.resolve(strict=True)
    cycles_budgeted = 0

    while cycles_budgeted < total_cycle_ceiling:
        segment_index = len(report["segments"]) + 1
        cycle_budget = min(segment_cycles, total_cycle_ceiling - cycles_budgeted)
        segment_name = f"segment_{segment_index:03d}_cycles_{cycle_budget:05d}"
        segment_dir = segments_root / segment_name
        prepared = prepare_relaxation(
            topology.resolve(strict=True),
            coordinates,
            engine.resolve(strict=True),
            segment_dir,
            f"{case_label} fully unrestrained MM {segment_name}",
            cycle_budget,
        )
        exit_status = runner(prepared["command"])
        if exit_status != 0:
            report["status"] = "mm_engine_failed"
            report["failed_segment"] = segment_name
            _write_report(output, report)
            raise CandidateControllerError(
                f"MM segment {segment_name} exited with status {exit_status}"
            )

        mdout = segment_dir / "minimized.mdout"
        convergence = check_convergence(mdout, drms=DRMS_TARGET)
        values = convergence["final_energy_rms_gmax"]
        if values is None:
            report["status"] = "mm_convergence_unresolved"
            report["failed_segment"] = segment_name
            _write_report(output, report)
            raise CandidateControllerError(
                f"MM segment {segment_name} has no valid final RMS record"
            )
        segment_report = {
            "name": segment_name,
            "cycle_budget": cycle_budget,
            "input_coordinates": str(coordinates),
            "endpoint_coordinates": str((segment_dir / "minimized.ncrst").resolve()),
            "convergence": convergence,
        }
        report["segments"].append(segment_report)
        cycles_budgeted += cycle_budget
        _write_report(output, report)
        LOGGER.info(
            "MM %s final RMS %.6f kcal mol-1 A-1 (budgeted %d/%d cycles)",
            segment_name,
            values[1],
            cycles_budgeted,
            total_cycle_ceiling,
        )

        if convergence["accepted"]:
            endpoint = segment_dir / "minimized.ncrst"
            gate_report = endpoint_gate(endpoint)
            report["mm_endpoint_gate"] = gate_report
            if not gate_report.get("accepted", False):
                report["status"] = "scientific_gate_failed"
                _write_report(output, report)
                raise CandidateControllerError(
                    "Converged MM endpoint failed the Tier A/B plus native Thr233 gate"
                )
            report["status"] = "qmmm_in_progress"
            _write_report(output, report)
            qmmm_report = qmmm_launcher(segment_dir, output / "qmmm_shells")
            report["qmmm"] = qmmm_report
            report["qmmm_launched"] = True
            report["status"] = "complete"
            _write_report(output, report)
            return report

        if cycles_budgeted == total_cycle_ceiling:
            report["status"] = "mm_cycle_ceiling_failed"
            _write_report(output, report)
            raise CandidateControllerError(
                "MM did not reach DRMS <= 0.01 before the total cycle ceiling"
            )
        coordinates = (segment_dir / "minimized.ncrst").resolve(strict=True)

    raise AssertionError("Sequential controller exited without a terminal state")


def _run_engine(command: Sequence[str]) -> int:
    """Run one Amber segment without changing its declared command."""
    return subprocess.run(list(command), check=False).returncode


def _write_report(output: Path, report: dict[str, Any]) -> None:
    """Write the current controller state as JSON for the outer job log."""
    (output / "controller_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
