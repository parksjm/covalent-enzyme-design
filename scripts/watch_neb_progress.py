#!/usr/bin/env python3
"""Watch native NEB profiles and export SVG plots without extra dependencies."""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from html import escape
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nyl12_covalent.neb_progress import NebReport, write_neb_progress  # noqa: E402

LOGGER = logging.getLogger(__name__)


def render_profile_svg(report: NebReport) -> str:
    """Render band energies and separate climbing-image marker as SVG.

    Args:
        report: Complete profile with source provenance.

    Returns:
        Standalone SVG with provisional candidates identified.
    """
    energies = [image.relative_energy_kcal_mol for image in report.images]
    if report.climbing_image is not None:
        energies.append(report.climbing_image.relative_energy_kcal_mol)
    low, high = min(energies), max(energies)
    margin = max((high - low) * 0.12, 1.0)
    low, high = low - margin, high + margin

    def position(image: int, energy: float) -> tuple[float, float]:
        return (
            90 + (image - 1) * 690 / (len(report.images) - 1),
            340 - (energy - low) * 260 / (high - low),
        )

    cycle = "pending convergence marker" if report.cycle is None else str(report.cycle)
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="850" height="470" viewBox="0 0 850 470">',
        '<rect width="850" height="470" fill="white"/>',
        '<g font-family="sans-serif" fill="#233248">',
        f'<text x="90" y="28" font-size="19">NEB report {report.report_index}; cycle {escape(cycle)}</text>',
        '<text x="90" y="51" font-size="12">Provisional energy profile; energy dips are not verified minima</text>',
        '<path d="M90 75V340H785" stroke="#556" fill="none"/>',
        '<text x="410" y="385" font-size="14">Regular band image</text>',
        '<text transform="translate(22 300) rotate(-90)" font-size="14">Relative energy (kcal/mol)</text>',
    ]
    for tick in range(6):
        energy = low + (high - low) * tick / 5
        _, y = position(1, energy)
        elements.append(f'<path d="M90 {y:.2f}H780" stroke="#e0e5eb"/>')
        elements.append(
            f'<text x="80" y="{y + 4:.2f}" text-anchor="end" font-size="11">{energy:.2f}</text>'
        )
    points = [
        position(image.image, image.relative_energy_kcal_mol) for image in report.images
    ]
    coordinates = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    elements.append(
        f'<polyline points="{coordinates}" fill="none" stroke="#176b99" stroke-width="2"/>'
    )
    for image, (x, y) in zip(report.images, points, strict=True):
        candidate = image.image in report.provisional_minimum_images
        color = "#be6b12" if candidate else "#176b99"
        elements.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="{color}"/>')
        elements.append(
            f'<text x="{x:.2f}" y="360" text-anchor="middle" font-size="10">{image.image}</text>'
        )
    if report.climbing_image is not None:
        ci = report.climbing_image
        x, y = position(ci.image, ci.relative_energy_kcal_mol)
        elements.append(
            f'<path d="M{x - 5:.2f} {y - 5:.2f}l10 10m0 -10l-10 10" stroke="#b22f4a" stroke-width="2"/>'
        )
    elements.extend(
        [
            '<text x="90" y="410" font-size="12">Blue: regular band; amber: provisional dip; red cross: climbing image</text>',
            f'<text x="90" y="432" font-size="10">Source SHA256: {report.source_sha256}</text>',
            f'<text x="90" y="449" font-size="10">Report bytes [{report.report_start_byte}, {report.report_end_byte}); coordinate association unverified</text>',
            "</g></svg>",
        ]
    )
    return "\n".join(elements) + "\n"


def update_progress_once(
    source: Path, destination: Path, expected_images: int
) -> NebReport:
    """Export a coherent energy profile, archive it, and log candidate dips.

    Args:
        source: Growing native Amber output.
        destination: Monitoring output directory.
        expected_images: Regular image count excluding the extra CI slot.

    Returns:
        Parsed profile. No mutable restart or XYZ coordinates are associated.

    Raises:
        ValueError: If a complete valid profile is unavailable.
        RuntimeError: If the source changes during the snapshot read.
    """
    report = write_neb_progress(source, destination, expected_images)
    svg = render_profile_svg(report)
    temporary = destination / ".neb_profile.svg.tmp"
    temporary.write_text(svg, encoding="utf-8")
    temporary.replace(destination / "neb_profile.svg")
    cycle = "pending" if report.cycle is None else f"cycle_{report.cycle:06d}"
    archive = destination / "reports" / f"report_{report.report_index:06d}_{cycle}"
    if not archive.exists():
        archive.mkdir(parents=True)
        for name in ("neb_progress.json", "neb_progress.csv", "neb_profile.svg"):
            (archive / name).write_bytes((destination / name).read_bytes())
    candidates = destination / "provisional_candidates.jsonl"
    identity = [report.report_index, report.report_end_byte]
    recorded = set()
    if candidates.exists():
        for line in candidates.read_text(encoding="utf-8").splitlines():
            record = json.loads(line)
            recorded.add(tuple(record["report_identity"]))
    if tuple(identity) not in recorded and report.provisional_minimum_images:
        metadata = json.loads((destination / "neb_progress.json").read_text())
        record = {
            "report_identity": identity,
            "cycle": report.cycle,
            "images": report.provisional_minimum_images,
            "source_sha256": report.source_sha256,
            "snapshot_utc": metadata["snapshot_utc"],
            "source_mtime_ns": metadata["source_mtime_ns"],
            "coordinate_association": "unverified; no release launched",
            "status": "provisional energy dips; stationary points not established",
        }
        with candidates.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
        LOGGER.info(
            "Report %d provisional dips: %s",
            report.report_index,
            report.provisional_minimum_images,
        )
    return report


def main(argv: list[str] | None = None) -> int:
    """Run one export or periodically watch a growing native NEB output.

    Args:
        argv: Optional command-line arguments.

    Returns:
        Zero on completion. Unexpected parsing or filesystem errors propagate.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--expected-images", type=int, required=True)
    parser.add_argument("--interval-seconds", type=float, default=60.0)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--stop-file", type=Path)
    args = parser.parse_args(argv)
    if args.interval_seconds < 1:
        parser.error("--interval-seconds must be at least 1")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    seen_source = False
    while args.stop_file is None or not args.stop_file.exists():
        if not args.source.exists() and not seen_source and not args.once:
            LOGGER.info("Waiting for Amber output %s", args.source)
        else:
            seen_source = True
            try:
                report = update_progress_once(
                    args.source, args.destination, args.expected_images
                )
                LOGGER.info(
                    "Exported complete report %d (cycle %s)",
                    report.report_index,
                    report.cycle,
                )
            except ValueError as error:
                if args.once or str(error) != "No complete NEB report found":
                    raise
                LOGGER.info("Waiting for first complete NEB report")
            except RuntimeError as error:
                if (
                    args.once
                    or str(error) != "Amber output changed while reading NEB snapshot"
                ):
                    raise
                LOGGER.info("Source advanced during snapshot; retrying next interval")
        if args.once:
            return 0
        time.sleep(args.interval_seconds)
    LOGGER.info("Watcher stopped by stop file %s", args.stop_file)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
