"""Plot all complete native NEB reports from one immutable output snapshot."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from pathlib import Path
import re
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = next(
    p for p in Path(__file__).resolve().parents if (p / "src/nyl12_covalent").is_dir()
)
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_progress import (  # noqa: E402
    HARTREE_TO_KCAL_MOL,
    parse_last_complete_neb_report,
)

LOGGER = logging.getLogger(__name__)
REPORT = re.compile(rb"(?m)^[ \t]*NEB Report[ \t]*\r?$")
GRAD = re.compile(
    rb"(?m)^\s*(Max|RMS) grad\s+([-+\d.EeDd]+)\s+Target:\s+([-+\d.EeDd]+)[^\n]*\n"
)


def parse_reports(payload: bytes, images: int) -> tuple[list[dict[str, Any]], bool]:
    """Parse every complete report; reject invalid interior or completed records.

    Args:
        payload: Exact frozen native bytes, read once from the live output.
        images: Declared regular band size, excluding the climbing slot.

    Returns:
        Ordered reports with optional complete convergence evidence, and whether
        one final incomplete report was ignored.
    """
    starts = [match.start() for match in REPORT.finditer(payload)]
    rows, ignored = [], False
    for index, start in enumerate(starts):
        end = starts[index + 1] if index + 1 < len(starts) else len(payload)
        block = payload[start:end]
        complete = block[: block.rfind(b"\n") + 1]
        try:
            report = parse_last_complete_neb_report(complete, images)
        except ValueError as error:
            if (
                index == len(starts) - 1
                and str(error) == "No complete NEB report found"
            ):
                ignored = True
                continue
            raise
        values = [image.relative_energy_kcal_mol for image in report.images]
        maximum = max(report.images, key=lambda image: image.relative_energy_kcal_mol)
        row = {
            "report_index": index + 1,
            "cycle": report.cycle,
            "maximum_kcal_mol": maximum.relative_energy_kcal_mol,
            "maximum_image": maximum.image,
            "endpoint_kcal_mol": values[-1],
            "image_relative_kcal_mol": values,
            "image_energy_hartree": [image.energy_hartree for image in report.images],
            "source_start_byte": start,
            "source_report_end_byte": start + report.report_end_byte,
            **{
                key: None
                for key in (
                    "max_grad",
                    "max_target",
                    "rms_grad",
                    "rms_target",
                    "max_ratio",
                    "rms_ratio",
                )
            },
        }
        following = complete[report.report_end_byte :]
        matched = list(GRAD.finditer(following)) if report.cycle is not None else []
        if len(matched) >= 2:
            if [match[1] for match in matched[:2]] != [b"Max", b"RMS"]:
                raise ValueError(
                    "Convergence gradient rows are duplicated or out of order"
                )
            for label, match in zip(("max", "rms"), matched[:2], strict=True):
                value, target = [
                    float(item.replace(b"D", b"E").replace(b"d", b"e"))
                    for item in (match[2], match[3])
                ]
                if not np.isfinite([value, target]).all() or value <= 0 or target <= 0:
                    raise ValueError(
                        "Native gradient and target must be finite positive"
                    )
                row.update(
                    {
                        f"{label}_grad": value,
                        f"{label}_target": target,
                        f"{label}_ratio": value / target,
                    }
                )
        rows.append(row)
    if not rows:
        raise ValueError("No complete native NEB profiles in frozen snapshot")
    cycles = [row["cycle"] for row in rows if row["cycle"] is not None]
    if any(right <= left for left, right in zip(cycles, cycles[1:])):
        raise ValueError("Native NEB convergence cycles are not strictly increasing")
    return rows, ignored


def main() -> None:
    """Freeze live output once and save three-panel evolution plus full JSON/CSV."""
    output = Path(__file__).resolve().parent
    source = output.parent / "calculation/minimized.mdout"
    snapshot = output / "minimized_snapshot.mdout"
    if snapshot.exists():
        raise FileExistsError("Frozen snapshot already exists; use a new plot attempt")
    payload = source.read_bytes()  # The only read of the live native output.
    with snapshot.open("xb") as stream:
        stream.write(payload)
    rows, ignored = parse_reports(payload, 24)
    selected = sorted(
        {
            0,
            len(rows) - 1,
            *[
                i
                for i, row in enumerate(rows)
                if row["cycle"] in (25, 50, 75, 100, 125)
            ],
        }
    )
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), constrained_layout=True)
    for index in selected:
        row = rows[index]
        axes[0].plot(
            range(1, 25),
            row["image_relative_kcal_mol"],
            marker=".",
            label=f"Cycle {row['cycle']}"
            if row["cycle"] is not None
            else f"Report {row['report_index']} (cycle pending)",
        )
    axes[0].set(
        xlabel="Regular image",
        ylabel="Same-band image1 ΔE (kcal/mol)",
        title="Selected native NEB profiles (chronological)",
    )
    axes[0].legend(ncol=4, fontsize=9)
    x = [
        row["cycle"] if row["cycle"] is not None else row["report_index"]
        for row in rows
    ]
    axes[1].plot(
        x, [row["maximum_kcal_mol"] for row in rows], color="#9b3a23", linewidth=1.5
    )
    axes[1].scatter(
        [x[0], x[-1]],
        [rows[0]["maximum_kcal_mol"], rows[-1]["maximum_kcal_mol"]],
        color="#222222",
        zorder=3,
    )
    axes[1].set(
        xlabel="Native convergence cycle (report index if pending)",
        ylabel="Profile maximum ΔE (kcal/mol)",
        title="Every complete report: nonmonotonic profile evolution",
    )
    force_rows = [row for row in rows if row["max_ratio"] is not None]
    if not force_rows:
        raise ValueError("No complete paired native Max/RMS gradient rows")
    for label, color in (("max", "#305a9d"), ("rms", "#25866a")):
        axes[2].semilogy(
            [r["cycle"] for r in force_rows],
            [r[f"{label}_ratio"] for r in force_rows],
            label=f"{label.upper()} gradient / native target",
            color=color,
        )
    axes[2].axhline(
        1, color="black", linestyle="--", linewidth=1, label="Native target"
    )
    axes[2].set(
        xlabel="Native convergence cycle",
        ylabel="Gradient / native target (log scale)",
        title="Native convergence evidence; gradient units retained as printed",
    )
    axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.suptitle(
        "Relay-acyl NEB evolution: diagnostic PES, not a verified TS barrier or PMF",
        fontsize=12,
    )
    for suffix in ("png", "svg"):
        fig.savefig(output / f"neb_evolution.{suffix}", dpi=180)
    plt.close(fig)
    last = rows[-20:]
    maxima = np.array([row["maximum_kcal_mol"] for row in last])
    summary = {
        "first": rows[0],
        "latest": rows[-1],
        "last20_count": len(last),
        "last20_maximum_delta_kcal_mol": float(maxima[-1] - maxima[0]),
        "last20_range_kcal_mol": [float(maxima.min()), float(maxima.max())],
        "last20_linear_slope_kcal_mol_per_report": float(
            np.polyfit(range(len(last)), maxima, 1)[0]
        )
        if len(last) > 1
        else None,
        "latest_complete_force_record": force_rows[-1],
    }
    record = {
        "snapshot": str(snapshot),
        "source": str(source),
        "snapshot_sha256": hashlib.sha256(payload).hexdigest(),
        "snapshot_bytes": len(payload),
        "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
        "regular_images": 24,
        "complete_reports": len(rows),
        "trailing_incomplete_report_ignored": ignored,
        "selected_report_indices": [rows[i]["report_index"] for i in selected],
        "energy_reference": "Each report referenced to its own regular image1",
        "reports": rows,
        "summary": summary,
        "source_hashes": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                Path(__file__).resolve(),
                ROOT / "src/nyl12_covalent/neb_progress.py",
            )
        },
        "interpretation": "Native NEB image energies and convergence diagnostics; neither a verified transition-state barrier nor a PMF.",
    }
    (output / "neb_evolution.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n"
    )
    with (output / "neb_evolution.csv").open("w", newline="") as stream:
        flat = [
            {key: value for key, value in row.items() if not isinstance(value, list)}
            for row in rows
        ]
        writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    LOGGER.info(
        "Reports=%d firstcycle=%s peak=%.3f latestcycle=%s peak=%.3f last20delta=%.3f",
        len(rows),
        rows[0]["cycle"],
        rows[0]["maximum_kcal_mol"],
        rows[-1]["cycle"],
        rows[-1]["maximum_kcal_mol"],
        summary["last20_maximum_delta_kcal_mol"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
