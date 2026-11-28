#!/usr/bin/env python3
"""Validate MCCE titration occupancies and plot selected beta-chain sites."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
from pathlib import Path
from typing import Any, Sequence

LOGGER = logging.getLogger(__name__)
TARGETS = (
    ("NTR", 1, "Thr233 N-terminus"),
    ("ASP", 40, "Asp272"),
    ("ASP", 42, "Asp274"),
)
CHAINS = "EFGH"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_head3(path: Path) -> dict[str, tuple[tuple[str, int, str], float]]:
    conformers: dict[str, tuple[tuple[str, int, str], float]] = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].split()[:5] == [
        "iConf",
        "CONFORMER",
        "FL",
        "occ",
        "crg",
    ]:
        raise ValueError("head3.lst has an unexpected or missing header")
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split()
        if not fields:
            continue
        if len(fields) < 5:
            raise ValueError(
                f"head3.lst line {line_number}: expected at least 5 fields"
            )
        name = fields[1]
        if len(name) < 10:
            raise ValueError(
                f"head3.lst line {line_number}: malformed conformer id {name!r}"
            )
        try:
            charge = float(fields[4])
            resnum = int(name[6:10])
        except ValueError as exc:
            raise ValueError(
                f"head3.lst line {line_number}: invalid charge or residue number"
            ) from exc
        if not math.isfinite(charge):
            raise ValueError(f"head3.lst line {line_number}: nonfinite charge")
        key = (name[:3], resnum, name[5])
        if name in conformers:
            raise ValueError(
                f"head3.lst line {line_number}: duplicate conformer id {name}"
            )
        conformers[name] = (key, charge)
    return conformers


def _read_fort38(path: Path) -> tuple[list[float], dict[str, list[float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError("fort.38 is empty")
    header = lines[0].split()
    if not header or header[0] != "ph":
        raise ValueError("fort.38 must start with a 'ph' header")
    try:
        ph = [float(value) for value in header[1:]]
    except ValueError as exc:
        raise ValueError("fort.38 has a nonnumeric pH value") from exc
    if not ph or any(not math.isfinite(value) for value in ph):
        raise ValueError("fort.38 pH grid must be nonempty and finite")
    if any(right <= left for left, right in zip(ph, ph[1:])):
        raise ValueError("fort.38 pH values must be strictly increasing")
    rows: dict[str, list[float]] = {}
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split()
        if not fields:
            continue
        name = fields[0]
        if name in rows:
            raise ValueError(
                f"fort.38 line {line_number}: duplicate conformer id {name}"
            )
        if len(fields) != len(ph) + 1:
            raise ValueError(
                f"fort.38 line {line_number}: expected {len(ph)} occupancies"
            )
        try:
            values = [float(value) for value in fields[1:]]
        except ValueError as exc:
            raise ValueError(
                f"fort.38 line {line_number}: nonnumeric occupancy"
            ) from exc
        if any(
            not math.isfinite(value) or value < 0.0 or value > 1.0 for value in values
        ):
            raise ValueError(
                f"fort.38 line {line_number}: occupancies must be finite and in [0, 1]"
            )
        rows[name] = values
    return ph, rows


def load_curves(native: Path) -> tuple[list[float], dict[str, Any]]:
    """Validate native MCCE files and sum charged-conformer occupancies.

    Args:
        native: Directory containing head3.lst and fort.38.

    Returns:
        The pH values and curve records grouped by target and chain.
    """
    conformers = _read_head3(native / "head3.lst")
    ph, rows = _read_fort38(native / "fort.38")
    if set(conformers) != set(rows):
        missing = sorted(set(conformers) - set(rows))
        extra = sorted(set(rows) - set(conformers))
        raise ValueError(
            f"head3.lst/fort.38 conformer IDs differ; missing={missing[:5]}, extra={extra[:5]}"
        )
    totals: dict[tuple[str, int, str], list[float]] = {}
    for name, (key, _charge) in conformers.items():
        values = rows[name]
        running = totals.setdefault(key, [0.0] * len(ph))
        for index, occupancy in enumerate(values):
            running[index] += occupancy
    for key, values in totals.items():
        for index, total in enumerate(values):
            if abs(total - 1.0) > 0.005:
                raise ValueError(
                    f"occupancy sum for {key} at pH {ph[index]} is {total:.6f}, outside 1±0.005"
                )

    curves: dict[str, Any] = {}
    for resname, resnum, label in TARGETS:
        target: dict[str, Any] = {}
        for chain in CHAINS:
            selected = [
                (name, charge)
                for name, (key, charge) in conformers.items()
                if key == (resname, resnum, chain)
            ]
            if not selected:
                raise ValueError(f"missing target {resname}{resnum} chain {chain}")
            protonated = [0.0] * len(ph)
            for name, charge in selected:
                include = charge > 0.5 if resname == "NTR" else abs(charge) < 0.5
                if include:
                    for index, occupancy in enumerate(rows[name]):
                        protonated[index] += occupancy
            if any(value < -1e-10 or value > 1.005 for value in protonated):
                raise ValueError(
                    f"protonated occupancy outside valid range for {label} chain {chain}"
                )
            target[chain] = {
                "protonated_fraction": protonated,
                "conformers": [name for name, _ in selected],
            }
        curves[f"{resname}{resnum}"] = {"label": label, "chains": target}
    return ph, curves


def _plot(ph: list[float], curves: dict[str, Any], output: Path) -> None:
    """Render the three selected pilot curves to a PNG file."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(zip(CHAINS, ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")))
    figure, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
    for axis, record in zip(axes, curves.values()):
        for chain, series in record["chains"].items():
            axis.plot(
                ph,
                series["protonated_fraction"],
                marker="o",
                markersize=3,
                color=colors[chain],
                label=f"Chain {chain}",
            )
        axis.axhline(0.5, color="0.35", linewidth=1, linestyle="--")
        axis.set_title(record["label"])
        axis.set_xlim(0, 14)
        axis.set_ylim(0, 1)
        axis.set_xlabel("pH")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Protonated fraction")
    axes[-1].legend(frameon=False, fontsize=8)
    figure.suptitle("MCCE titration curves — pilot; not validated")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main(argv: Sequence[str] | None = None) -> int:
    """Validate native titration files and write a JSON report and plot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument(
        "--output", required=True, type=Path, help="new or existing output directory"
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    ph, curves = load_curves(args.native)
    args.output.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).resolve()
    report = {
        "status": "pilot; not validated",
        "interpretation": "No fitted or validated pKa values are reported.",
        "pH": ph,
        "curves": curves,
        "source_sha256": {
            "head3.lst": _sha256(args.native / "head3.lst"),
            "fort.38": _sha256(args.native / "fort.38"),
            "script": _sha256(script),
        },
    }
    with (args.output / "titration_curves.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    _plot(ph, curves, args.output / "titration_curves.png")
    LOGGER.info("Wrote pilot curves to %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
