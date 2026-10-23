#!/usr/bin/env python3
"""Extract Gaussian Cartesian orientation snapshots into multi-frame XYZ."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

LOGGER = logging.getLogger(__name__)
ORIENTATION_RE = re.compile(r"^\s*(Standard|Input) orientation:\s*$")
DASH_RE = re.compile(r"^\s*-{5,}\s*$")
SYMBOLS = [""] + (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn "
    "Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd "
    "Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho "
    "Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra "
    "Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt "
    "Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


@dataclass(frozen=True)
class Frame:
    """One complete Gaussian orientation table."""

    orientation: str
    line: int
    centers: tuple[int, ...]
    atomic_numbers: tuple[int, ...]
    coordinates: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class ParsedLog:
    """Parsed orientation tables and counts for one Gaussian log."""

    frames: tuple[Frame, ...]
    counts: dict[str, int]
    skipped_trailing: dict[str, int]


def _parse_table(lines: Sequence[str], start: int) -> Frame | None:
    """Parse a table, returning None only when EOF truncates that table."""
    match = ORIENTATION_RE.match(lines[start])
    if match is None:
        raise ValueError(f"Expected orientation heading on line {start + 1}")
    orientation = match.group(1)
    cursor = start + 1
    while cursor < len(lines) and not DASH_RE.match(lines[cursor]):
        cursor += 1
    if cursor == len(lines):
        return None
    if cursor + 3 >= len(lines):
        return None
    coordinate_header = lines[cursor + 1].strip()
    column_header = " ".join(lines[cursor + 2].split())
    if not re.search(r"Coordinates\s*\(Angstroms?\)", coordinate_header, re.I):
        raise ValueError(
            f"Orientation table is not in Angstroms near line {cursor + 2}"
        )
    if column_header != "Number Number Type X Y Z":
        raise ValueError(
            f"Unrecognized Gaussian orientation columns near line {cursor + 3}"
        )
    if not DASH_RE.match(lines[cursor + 3]):
        raise ValueError(f"Missing orientation table delimiter near line {cursor + 4}")
    cursor += 4
    centers: list[int] = []
    atomic_numbers: list[int] = []
    coordinates: list[tuple[float, float, float]] = []
    while cursor < len(lines) and not DASH_RE.match(lines[cursor]):
        fields = lines[cursor].split()
        if len(fields) != 6:
            raise ValueError(f"Malformed orientation row on line {cursor + 1}")
        try:
            center = int(fields[0])
            atomic_number = int(fields[1])
            int(fields[2])
            xyz = tuple(float(value.replace("D", "E")) for value in fields[3:])
        except ValueError as error:
            raise ValueError(
                f"Malformed orientation row on line {cursor + 1}"
            ) from error
        if atomic_number < 1 or atomic_number >= len(SYMBOLS):
            raise ValueError(
                f"Unsupported ghost/dummy or atomic number {atomic_number} on line {cursor + 1}"
            )
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Nonfinite coordinate on line {cursor + 1}")
        centers.append(center)
        atomic_numbers.append(atomic_number)
        coordinates.append(xyz)
        cursor += 1
    if cursor == len(lines):
        return None
    if not centers or centers != list(range(1, len(centers) + 1)):
        raise ValueError(
            f"Center numbers must be sequential from 1 near line {start + 1}"
        )
    return Frame(
        orientation,
        start + 1,
        tuple(centers),
        tuple(atomic_numbers),
        tuple(coordinates),
    )


def parse_log(
    text: str, orientation: str = "auto", allow_incomplete: bool = False
) -> ParsedLog:
    """Parse Gaussian orientation snapshots and select the requested stream.

    Args:
        text: Complete decoded Gaussian log contents.
        orientation: ``auto``, ``standard``, or ``input`` orientation stream.
        allow_incomplete: Ignore one final orientation table cut off by EOF.

    Returns:
        Selected frames, complete table counts, and ignored trailing-table count.

    Raises:
        ValueError: For malformed tables, inconsistent atom order, or truncation.
    """
    if orientation not in {"auto", "standard", "input"}:
        raise ValueError("orientation must be auto, standard, or input")
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if ORIENTATION_RE.match(line)]
    parsed: dict[str, list[Frame]] = {"Standard": [], "Input": []}
    skipped = {"Standard": 0, "Input": 0}
    for index, start in enumerate(starts):
        frame = _parse_table(lines, start)
        if frame is None:
            is_trailing = index == len(starts) - 1
            if not allow_incomplete or not is_trailing:
                raise ValueError(
                    f"Truncated {ORIENTATION_RE.match(lines[start]).group(1)} orientation table at line {start + 1}"
                )
            label = ORIENTATION_RE.match(lines[start]).group(1)
            skipped[label] += 1
            continue
        parsed[frame.orientation].append(frame)
    chosen = (
        orientation.title()
        if orientation != "auto"
        else (
            "Standard"
            if starts
            and any(
                "Standard" == ORIENTATION_RE.match(lines[i]).group(1) for i in starts
            )
            else "Input"
        )
    )
    frames = parsed[chosen]
    if not frames:
        raise ValueError(f"No complete {chosen} orientation tables found")
    reference = frames[0].atomic_numbers
    for frame in frames:
        if frame.centers != frames[0].centers or frame.atomic_numbers != reference:
            raise ValueError(
                f"Atom count/order changes in {frame.orientation} orientation at line {frame.line}"
            )
    return ParsedLog(
        tuple(frames), {name: len(parsed[name]) for name in parsed}, skipped
    )


def trajectory_text(parsed: ParsedLog, source_name: str, source_sha256: str) -> str:
    """Serialize selected snapshots as multi-frame XYZ with JSON comments."""
    output: list[str] = []
    for frame_index, frame in enumerate(parsed.frames, 1):
        comment = json.dumps(
            {
                "original_line": frame.line,
                "orientation": frame.orientation,
                "frame_index": frame_index,
                "source_basename": source_name,
                "source_sha256": source_sha256,
            },
            separators=(",", ":"),
        )
        output.append(str(len(frame.atomic_numbers)))
        output.append(comment)
        for atomic_number, xyz in zip(
            frame.atomic_numbers, frame.coordinates, strict=True
        ):
            output.append(
                f"{SYMBOLS[atomic_number]:<2s} {xyz[0]: .10f} {xyz[1]: .10f} {xyz[2]: .10f}"
            )
    return "\n".join(output) + "\n"


def _sha256(data: bytes) -> str:
    """Return a SHA-256 digest for bytes."""
    return hashlib.sha256(data).hexdigest()


def convert(
    log_path: Path, output_path: Path, orientation: str, allow_incomplete: bool
) -> dict[str, object]:
    """Convert a log and create a JSON provenance sidecar without overwrites."""
    sidecar_path = output_path.with_suffix(".json")
    if output_path.resolve() == sidecar_path.resolve():
        raise ValueError("XYZ trajectory and JSON sidecar paths must differ")
    collisions = [path for path in (output_path, sidecar_path) if path.exists()]
    if collisions:
        raise FileExistsError(", ".join(map(str, collisions)))
    source_bytes = log_path.read_bytes()
    source_hash = _sha256(source_bytes)
    parsed = parse_log(
        source_bytes.decode("utf-8", errors="strict"), orientation, allow_incomplete
    )
    xyz = trajectory_text(parsed, log_path.name, source_hash)
    xyz_bytes = xyz.encode("utf-8")
    selected_orientation = parsed.frames[0].orientation
    report: dict[str, object] = {
        "source_basename": log_path.name,
        "source_sha256": source_hash,
        "script_sha256": _sha256(Path(__file__).read_bytes()),
        "requested_orientation": orientation,
        "selected_orientation": selected_orientation,
        "orientation_counts": parsed.counts,
        "ignored_trailing_incomplete_counts": parsed.skipped_trailing,
        "selected_frame_count": len(parsed.frames),
        "atom_count_per_frame": len(parsed.frames[0].atomic_numbers),
        "trajectory_sha256": _sha256(xyz_bytes),
        "frame_selection_note": "Every complete table from the selected orientation stream is retained, including identical snapshots; tables may come from optimization, frequency, ESP, or other Gaussian sections and are not asserted to be optimizer steps.",
    }
    _write_exclusive(output_path, xyz_bytes)
    try:
        _write_exclusive(
            sidecar_path, (json.dumps(report, indent=2) + "\n").encode("utf-8")
        )
    except BaseException:
        output_path.unlink(missing_ok=True)
        raise
    return report


def _write_exclusive(path: Path, data: bytes) -> None:
    """Write bytes exclusively, refusing to overwrite an existing file."""
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(data)


def main() -> None:
    """Run the command-line log converter."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--orientation", choices=("auto", "standard", "input"), default="auto"
    )
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    report = convert(args.log, args.output, args.orientation, args.allow_incomplete)
    LOGGER.info(
        "Wrote %d %s orientation frames to %s",
        report["selected_frame_count"],
        report["selected_orientation"],
        args.output,
    )


if __name__ == "__main__":
    main()
