"""Compare printed native NEB replay profiles without certifying convergence."""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any

from .neb_progress import parse_last_complete_neb_report

_START = re.compile(rb"(?m)^[ \t]*NEB Report[ \t]*\r?$")
_NUMBER = rb"([-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?)"
_MAX = re.compile(rb"(?m)^[ \t]*Max grad[ \t]+" + _NUMBER + rb"[^\r\n]*\r?\n")
_RMS = re.compile(rb"(?m)^[ \t]*RMS grad[ \t]+" + _NUMBER + rb"[^\r\n]*\r?\n")


@dataclass(frozen=True)
class ReplayProfile:
    """One complete energy profile and its following max/RMS convergence row."""

    cycle: int
    energies: tuple[float, ...]
    max_gradient: float
    rms_gradient: float
    max_token: str
    rms_token: str
    end_byte: int


def parse_replay_profiles(
    payload: bytes, expected_images: int = 24
) -> dict[int, ReplayProfile]:
    """Parse all completed report/convergence blocks, ignoring an incomplete tail.

    Args:
        payload: Immutable output prefix read once from a native engine.
        expected_images: Number of regular band images.

    Returns:
        Cycle-keyed profiles with source prefix offsets for exact provenance.

    Raises:
        ValueError: If a completed profile has invalid images or duplicate cycles.
    """
    starts = [match.start() for match in _START.finditer(payload)]
    profiles = {}
    for index, start in enumerate(starts):
        stop = starts[index + 1] if index + 1 < len(starts) else len(payload)
        block = payload[start:stop]
        if b"Total path length " not in block:
            continue
        try:
            report = parse_last_complete_neb_report(block, expected_images)
        except ValueError as error:
            if str(error) == "No complete NEB report found":
                continue
            raise
        maximum, rms = _MAX.search(block), _RMS.search(block)
        if report.cycle is None or maximum is None or rms is None:
            continue
        if report.cycle in profiles:
            raise ValueError("Duplicate complete native cycle")
        tokens = [match[1].decode("ascii") for match in (maximum, rms)]
        values = [float(token.replace("D", "E").replace("d", "e")) for token in tokens]
        if not all(math.isfinite(value) for value in values):
            raise ValueError("Nonfinite gradient")
        profiles[report.cycle] = ReplayProfile(
            report.cycle,
            tuple(image.energy_hartree for image in report.images),
            values[0],
            values[1],
            tokens[0],
            tokens[1],
            start + rms.end(),
        )
    return profiles


def printed_quantum(token: str) -> float:
    """Return one last-printed-digit unit for a decimal/scientific number."""
    mantissa, _, exponent = token.lower().replace("d", "e").partition("e")
    decimals = len(mantissa.partition(".")[2])
    return 10.0 ** (int(exponent or "0") - decimals)


def compare_replay(
    pilot: dict[int, ReplayProfile],
    replay: dict[int, ReplayProfile],
) -> dict[str, Any]:
    """Compare same-cycle energies and gradients through the immutable pilot end.

    Args:
        pilot: Complete finished-pilot cycle records.
        replay: Complete replay records parsed from an immutable current prefix.

    Returns:
        Initial-profile, percycle, maximum discrepancy and replay-progress flags.
        Replay agreement is not force convergence, chemical or TS acceptance.

    Raises:
        ValueError: If the pilot has no complete cycle or sequence starts after1.
    """
    if not pilot or sorted(pilot) != list(range(1, max(pilot) + 1)):
        raise ValueError("Pilot must contain complete consecutive cycles starting1")
    rows = []
    for cycle in sorted(set(pilot) & set(replay)):
        old, new = pilot[cycle], replay[cycle]
        if len(old.energies) != len(new.energies):
            raise ValueError("Replay image count changed")
        energy = max(
            abs(a - b) for a, b in zip(old.energies, new.energies, strict=True)
        )
        mtol = (
            0.5 * (printed_quantum(old.max_token) + printed_quantum(new.max_token))
            + 1e-15
        )
        rtol = (
            0.5 * (printed_quantum(old.rms_token) + printed_quantum(new.rms_token))
            + 1e-15
        )
        mdiff = abs(old.max_gradient - new.max_gradient)
        rdiff = abs(old.rms_gradient - new.rms_gradient)
        rows.append(
            dict(
                cycle=cycle,
                energy_max_difference_Eh=energy,
                energy_tolerance_Eh=1.1e-7,
                max_gradient_difference=mdiff,
                max_gradient_tolerance=mtol,
                rms_gradient_difference=rdiff,
                rms_gradient_tolerance=rtol,
                passed=bool(energy <= 1.1e-7 and mdiff <= mtol and rdiff <= rtol),
            )
        )
    failed = [row["cycle"] for row in rows if not row["passed"]]
    matched = len(rows) == len(pilot) and not failed
    passed_pilot = matched and max(replay, default=0) > max(pilot)
    return dict(
        status="diverged"
        if failed
        else (
            "replay_passed_pilot"
            if passed_pilot
            else (
                "matched_pilot_waiting_next_cycle"
                if matched
                else ("replay_matching" if rows else "pending_initial_profile")
            )
        ),
        compared_cycles=len(rows),
        last_compared_cycle=rows[-1]["cycle"] if rows else None,
        pilot_last_cycle=max(pilot),
        replay_last_complete_cycle=max(replay, default=0),
        initial_24image_profile_verified=bool(
            rows and rows[0]["cycle"] == 1 and rows[0]["passed"]
        ),
        replay_matched_pilot=matched,
        replay_passed_pilot=passed_pilot,
        failed_cycles=failed,
        cycles=rows,
        max_energy_difference_Eh=max(
            (row["energy_max_difference_Eh"] for row in rows), default=0.0
        ),
        max_gradient_difference=max(
            (row["max_gradient_difference"] for row in rows), default=0.0
        ),
        max_rms_gradient_difference=max(
            (row["rms_gradient_difference"] for row in rows), default=0.0
        ),
        scientific_success_claim=False,
        stationarity_or_TS_verified=False,
    )
