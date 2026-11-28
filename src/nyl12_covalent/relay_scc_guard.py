"""Bound native zero-temperature SCC recovery without accepting invalid forces."""

from __future__ import annotations

import re

from .dftb_scc_validation import audit_dftb_scc


_ITERATION = re.compile(r"^\s*QMMM SCC-DFTB:\s+(\d+)\s+[-+\d.*]", re.M)
_FAILED = re.compile(r"DID NOT CONVERGE AFTER\s+(\d+) cycles\.", re.I)
_CONVERGED = re.compile(
    r"SCC-DFTB for step\s+\d+\s+converged in\s+(\d+) cycles\.", re.I
)
_RETURNED_GRADIENT = "Energy calculation finished, energy:"


def scc_guard_reason(
    text: str, *, final: bool = False, maximum_iterations: int = 500
) -> str | None:
    """Allow bounded native SCC retries while rejecting invalid returned forces.

    AMBER's70-cycle Broyden reset is a recoverable warning when the same
    electronic evaluation subsequently reaches its original convergence
    criteria. No Hamiltonian, temperature or SCF tolerance changes here.

    Args:
        text: Current native output, possibly ending inside an active SCC trial.
        final: Whether no further same-trial recovery can occur.
        maximum_iterations: Maximum cumulative SCC iterations per evaluation.

    Returns:
        Fatal, exhausted-recovery or unresolved-gradient reason; otherwise None.
    """
    if maximum_iterations < 1:
        raise ValueError("Positive SCC iteration budget required")
    history = audit_dftb_scc(text)
    if history["fatal_markers"] or history["temperature_changes"]:
        return "Native electronic fatal/temperature-change marker"
    iterations = [int(value) for value in _ITERATION.findall(text)]
    converged = [int(value) for value in _CONVERGED.findall(text)]
    failed = [int(value) for value in _FAILED.findall(text)]
    if any(value > maximum_iterations for value in iterations + converged) or any(
        value >= maximum_iterations for value in failed
    ):
        return f"SCC recovery exceeded {maximum_iterations}-iteration evaluation budget"
    # Each force-return boundary must have its own resolved electronic history;
    # repeated native step=1 labels must not let a later evaluation mask failure.
    for completed in text.split(_RETURNED_GRADIENT)[:-1]:
        if not audit_dftb_scc(completed)["passed"]:
            return "Native returned an unresolved SCC electronic gradient"
    if final and not history["passed"]:
        return "Final native output contains unresolved SCC convergence failure"
    return None
