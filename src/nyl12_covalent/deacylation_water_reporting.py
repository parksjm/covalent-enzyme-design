"""Output writers for the PA6 deacylation hydration reconnaissance."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any

from nyl12_covalent.deacylation_water_geometry import _plot, sha256

LOGGER = logging.getLogger(__name__)


def write_reports(
    output: Path,
    summary: dict[str, Any],
    water_rows: list[dict[str, Any]],
    nearest_rows: list[tuple[Any, ...]],
    minima: dict[int, dict[str, Any]],
    frame_count: int,
    code_paths: list[Path],
    source_paths: list[Path],
) -> None:
    """Write JSON, CSV, Markdown, plots, and hash manifest for one analysis."""
    nearest_path = output / "nearest_water_per_frame.csv"
    with nearest_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "frame_index",
                "stage",
                "time_ps",
                "water_residue_index_1based",
                "oxygen_atom_index_1based",
                "O_C_A",
                "attack_angle_deg",
            ]
        )
        writer.writerows(nearest_rows)
    summary["files"] = {
        "per_frame_metrics": "per_frame_metrics.csv",
        "nearest_water_per_frame": nearest_path.name,
        "per_water_summary": "per_water_summary.json",
        "report_markdown": "report.md",
    }
    (output / "per_water_summary.json").write_text(
        json.dumps(
            {"schema_version": 1, "waters": water_rows}, indent=2, allow_nan=False
        )
        + "\n"
    )
    (output / "report.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    _plot(output, nearest_rows)
    total_near = sum(item["near_attack_frame_count"] for item in minima.values())
    stage_text = ", ".join(
        f"{key}: {value}" for key, value in summary["frame_count_by_stage"].items()
    )
    report_lines = [
        "# Deacylation hydration geometry reconnaissance",
        "",
        f"Frames analyzed: {frame_count} ({stage_text}).",
        f"All {summary['water_count']:,} water oxygens were evaluated in every frame with triclinic minimum-image distances.",
        "The exact 0 ps endpoint is included before the saved heating and 300 K equilibration frames.",
        "",
        f"Pending, unaccepted extension segments excluded: "
        f"{', '.join(item['segment'] for item in summary['pending_extension_segments']) or 'none'}.",
        "",
        "Near-attack annotation: water O-C ≤3.5 Å and O-C-Ocarbonyl 95–125°. "
        "This is exploratory geometry only and does not establish reactivity or proton transfer.",
        "",
        f"Water-frame observations meeting the near-attack gate: {total_near}.",
        "",
        "Top five waters ranked by minimum sampled O-C distance:",
        "",
        "| Water residue (1-based) | O atom (1-based) | H atoms (1-based) | Minimum O-C (Å) | Angle at minimum (°) | Near-attack frames |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    report_lines.extend(
        f"| {w['water_residue_index_1based']} | {w['oxygen_atom_index_1based']} | "
        f"{', '.join(map(str, w['hydrogen_atom_indices_1based']))} | "
        f"{w['min_attack_distance_A']:.3f} | {w['min_attack_angle_deg']:.1f} | "
        f"{w['near_attack_frame_count']} |"
        for w in water_rows[:5]
    )
    endpoint = nearest_rows[0]
    final = nearest_rows[-1]
    report_lines.extend(
        [
            "",
            f"At t=0 ps the nearest water is residue {endpoint[3]}/O {endpoint[4]} "
            f"at {endpoint[5]:.3f} Å and {endpoint[6]:.1f}°. At the final {final[2]:.2f} ps "
            f"frame it is {final[5]:.3f} Å and {final[6]:.1f}°.",
            "Water-H···Nalpha values are descriptive only because Nalpha is protonated. "
            "Asp OD1 and OD2 are checked independently; outputs include selected acceptor atom, "
            "O···A distance, H···A distance, and O-H···A angle.",
            "Sampling is short and not converged; these geometries do not indicate a reaction mechanism.",
            "",
        ]
    )
    (output / "report.md").write_text("\n".join(report_lines))
    manifest = {str(path): sha256(path) for path in [*source_paths, *code_paths]}
    output_paths = [
        output / name
        for name in (
            "per_frame_metrics.csv",
            "nearest_water_per_frame.csv",
            "per_water_summary.json",
            "report.json",
            "report.md",
            "nearest_water_geometry.png",
        )
    ]
    manifest.update({str(path): sha256(path) for path in output_paths})
    (output / "source_manifest.json").write_text(
        json.dumps({"sha256": manifest}, indent=2, allow_nan=False) + "\n"
    )
    LOGGER.info("Wrote exploratory hydration analysis for %d frames", frame_count)
