"""Audit native DFTB SCC retry history without hiding unrecovered failures."""

from __future__ import annotations

import re

_RECOVERABLE_TRIAL = re.compile(
    r"^\s*QMMM SCC-DFTB:\s*SCC-DFTB FOR STEP\s+(\d+)\s+"
    r"DID NOT CONVERGE AFTER\s+(\d+) cycles\.\s*$",
    re.IGNORECASE,
)
_CONVERGED_TRIAL = re.compile(
    r"^\s*QMMM SCC-DFTB:\s*SCC-DFTB for step\s+(\d+)\s+"
    r"converged in\s+(\d+) cycles\.\s*$",
    re.IGNORECASE,
)
_FATAL_MARKER = re.compile(
    r"(?:\b(?:SCC|SCF)\b.*(?:fail|not\s+converg|will not be accurate)|"
    r"convergence could not be achieved|\bITRMAX\b.*exceed|"
    r"QM calculation failed|\bNaN\b|will not be accurate)",
    re.IGNORECASE,
)
_TEMPERATURE_CHANGE = re.compile(
    r"increasing the electronic temperature|"
    r"energy for this step was calculated with a higher|"
    r"electronic temperature restored",
    re.IGNORECASE,
)
_TELEC = re.compile(
    r"\bTELEC\s*[:=]\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)", re.I
)


def audit_dftb_scc(text: str) -> dict[str, object]:
    """Audit SCC warnings, same-step recoveries, fatal markers, and temperature.

    Only AMBER's exact per-trial SCC warning is recoverable, and only when a
    later native convergence summary for the same step appears in this text.
    Marker-free text is reported as unverified history; callers must retain
    their independent native process, energy, termination, and force gates.

    Args:
        text: Native Sander output text for one calculation.

    Returns:
        A JSON-compatible dictionary describing accepted history and evidence.
    """
    pending: dict[int, list[dict[str, int]]] = {}
    recovered: list[dict[str, int]] = []
    fatal_markers: list[str] = []
    temperature_changes: list[str] = []
    history_observed = False

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        failed = _RECOVERABLE_TRIAL.fullmatch(raw_line)
        if failed:
            history_observed = True
            step, cycles = map(int, failed.groups())
            pending.setdefault(step, []).append({"step": step, "failed_cycles": cycles})
            continue

        converged = _CONVERGED_TRIAL.fullmatch(raw_line)
        if converged:
            history_observed = True
            step, cycles = map(int, converged.groups())
            for trial in pending.pop(step, []):
                recovered.append({**trial, "recovered_cycles": cycles})
            continue

        if _FATAL_MARKER.search(line):
            fatal_markers.append(line[:240])
        if _TEMPERATURE_CHANGE.search(line):
            temperature_changes.append(line[:240])
        for match in _TELEC.finditer(line):
            if float(match.group(1)) != 0.0:
                temperature_changes.append(line[:240])

    unresolved = [trial for trials in pending.values() for trial in trials]
    passed = not fatal_markers and not unresolved and not temperature_changes
    return {
        "passed": passed,
        "history_observed": history_observed,
        "verified_scc_history": history_observed and passed,
        "recovered_trials": recovered,
        "recovered_trial_count": len(recovered),
        "unresolved_trials": unresolved,
        "fatal_markers": fatal_markers,
        "temperature_changes": temperature_changes,
    }
