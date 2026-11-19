"""Inspect native Amber/DL-FIND NEB progress without certifying stationary points."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

HARTREE_TO_KCAL_MOL = 627.5094740631
_CYCLE = re.compile(r"Testing convergence\s+in cycle\s+(\d+)")


@dataclass(frozen=True)
class NebImage:
    """One regular band image or separately reported climbing image."""

    image: int
    energy_hartree: float
    relative_energy_kcal_mol: float
    climbing: bool = False


@dataclass(frozen=True)
class NebReport:
    """Complete profile from an immutable output snapshot.

    The cycle belongs to the following convergence test, when available.
    Report offsets identify exactly which output was consumed. Restart files
    are mutable independent outputs and are never implicitly paired here.
    """

    report_index: int
    cycle: int | None
    report_start_byte: int
    report_end_byte: int
    source_sha256: str
    images: tuple[NebImage, ...]
    climbing_image: NebImage | None
    provisional_minimum_images: tuple[int, ...]


def _image(line: str, climbing: bool, reference: float) -> NebImage:
    # Native format: ('Img ',i4,f15.7,3f10.5,3f8.2,2x,a).
    image = int(line[4:8])
    energy = float(line[8:23])
    if not math.isfinite(energy):
        raise ValueError("NEB image energy is not finite")
    return NebImage(image, energy, (energy - reference) * HARTREE_TO_KCAL_MOL, climbing)


def parse_last_complete_neb_report(payload: bytes, expected_images: int) -> NebReport:
    """Parse the last complete native profile, ignoring a trailing partial report.

    Args:
        payload: Immutable bytes read from Amber output.
        expected_images: Number of regular images, including endpoints but
            excluding DL-FIND's extra climbing-image slot.

    Returns:
        Last complete profile with absolute and reactant-relative energies.
        Strict interior energy dips are provisional candidates only.

    Raises:
        ValueError: If input or a completed profile is invalid, or none exists.
    """
    if expected_images < 2:
        raise ValueError("At least two regular images are required")
    lines = payload.splitlines(keepends=True)
    offset = 0
    start = 0
    report_index = 0
    rows: list[str] | None = None
    completed: NebReport | None = None
    pending_cycle = False
    source_hash = hashlib.sha256(payload).hexdigest()
    for raw_line in lines:
        line = raw_line.decode("ascii", errors="strict").rstrip("\r\n")
        end = offset + len(raw_line)
        if line.strip() == "NEB Report":
            rows = []
            start = offset
            report_index += 1
            pending_cycle = False
        elif rows is not None and line.startswith(("Img ", "Cimg")):
            rows.append(line)
        elif rows is not None and line.startswith("Total path length "):
            # An unterminated final line may still be in the writer's buffer.
            if not raw_line.endswith(b"\n"):
                break
            regular = [row for row in rows if row.startswith("Img ")]
            climbing = [row for row in rows if row.startswith("Cimg")]
            if len(regular) != expected_images or len(climbing) > 1:
                raise ValueError("Completed NEB report has unexpected image count")
            reference = float(regular[0][8:23])
            images = tuple(_image(row, False, reference) for row in regular)
            if [image.image for image in images] != list(range(1, expected_images + 1)):
                raise ValueError("Regular NEB images are duplicated or out of order")
            ci = _image(climbing[0], True, reference) if climbing else None
            if ci is not None and not 1 <= ci.image <= expected_images:
                raise ValueError("Climbing image refers to an invalid regular image")
            minima = tuple(
                images[i].image
                for i in range(1, len(images) - 1)
                if images[i].energy_hartree < images[i - 1].energy_hartree
                and images[i].energy_hartree < images[i + 1].energy_hartree
            )
            completed = NebReport(
                report_index, None, start, end, source_hash, images, ci, minima
            )
            pending_cycle = True
            rows = None
        elif pending_cycle:
            match = _CYCLE.search(line)
            if match and completed is not None:
                completed = NebReport(**{**completed.__dict__, "cycle": int(match[1])})
                pending_cycle = False
        offset = end
    if completed is None:
        raise ValueError("No complete NEB report found")
    return completed


def write_neb_progress(
    source: Path, destination: Path, expected_images: int
) -> NebReport:
    """Write JSON and CSV snapshots of the last complete NEB profile.

    Args:
        source: Amber output file; read once between matching stat checks.
        destination: Output directory, created when absent.
        expected_images: Regular band image count, excluding climbing slot.

    Returns:
        Parsed report. JSON records source time/hash and explicitly states that
        current mutable restart files have not been associated with the report.

    Raises:
        RuntimeError: If the source changes during the snapshot read.
        ValueError: If there is no valid completed profile.
    """
    before = source.stat()
    payload = source.read_bytes()
    after = source.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise RuntimeError("Amber output changed while reading NEB snapshot")
    report = parse_last_complete_neb_report(payload, expected_images)
    destination.mkdir(parents=True, exist_ok=True)
    metadata = {
        **asdict(report),
        "source": str(source.resolve()),
        "source_mtime_ns": after.st_mtime_ns,
        "snapshot_utc": datetime.now(timezone.utc).isoformat(),
        "restart_files_associated": False,
        "scientific_status": "provisional profile; energy dips are not verified minima",
    }
    (destination / "neb_progress.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    with (destination / "neb_progress.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "report_index",
                "cycle",
                "image",
                "climbing",
                "energy_hartree",
                "relative_energy_kcal_mol",
                "provisional_minimum",
            ]
        )
        all_images = (
            *report.images,
            *((report.climbing_image,) if report.climbing_image is not None else ()),
        )
        for image in all_images:
            writer.writerow(
                [
                    report.report_index,
                    report.cycle,
                    image.image,
                    image.climbing,
                    image.energy_hartree,
                    image.relative_energy_kcal_mol,
                    not image.climbing
                    and image.image in report.provisional_minimum_images,
                ]
            )
    return report
