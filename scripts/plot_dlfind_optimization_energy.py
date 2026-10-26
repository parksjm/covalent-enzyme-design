#!/usr/bin/env python3
"""Extract DL-FIND energies from an Amber output and plot optimization progress."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LOGGER = logging.getLogger(__name__)
HARTREE_TO_KCAL_MOL = 627.5094740631
ENERGY_PATTERN = re.compile(
    rb"Energy calculation finished, energy:\s*([+-]?\d+\.\d+E[+-]\d+)"
)
CYCLE_PATTERN = re.compile(rb"Testing convergence\s+in cycle\s+(\d+)")


@dataclass(frozen=True)
class EnergyRecord:
    """One completed DL-FIND cycle and its reported total energy."""

    cycle: int
    energy_hartree: float


def parse_energy_records(payload: bytes) -> list[EnergyRecord]:
    """Pair each reported energy with its following DL-FIND cycle.

    Args:
        payload: Immutable snapshot of an Amber/DL-FIND output file.

    Returns:
        Ordered, gap-free energy records.

    Raises:
        ValueError: If records are missing, duplicated, or nonsequential.
    """
    records: list[EnergyRecord] = []
    pending_energy: float | None = None
    for line in payload.splitlines():
        energy_match = ENERGY_PATTERN.search(line)
        if energy_match:
            if pending_energy is not None:
                raise ValueError("Found a second energy before a cycle record")
            pending_energy = float(energy_match.group(1))
            continue
        cycle_match = CYCLE_PATTERN.search(line)
        if cycle_match and pending_energy is not None:
            records.append(
                EnergyRecord(
                    cycle=int(cycle_match.group(1)), energy_hartree=pending_energy
                )
            )
            pending_energy = None
    if not records:
        raise ValueError("No complete DL-FIND energy/cycle pairs found")
    cycles = [record.cycle for record in records]
    expected = list(range(cycles[0], cycles[-1] + 1))
    if cycles != expected:
        raise ValueError("DL-FIND cycle records are duplicated or nonsequential")
    return records


def write_csv(records: list[EnergyRecord], destination: Path) -> None:
    """Write absolute and relative energies as a CSV table.

    Args:
        records: Ordered optimization records.
        destination: CSV output path.
    """
    reference = records[0].energy_hartree
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cycle", "energy_hartree", "relative_energy_kcal_mol"])
        for record in records:
            writer.writerow(
                [
                    record.cycle,
                    f"{record.energy_hartree:.10f}",
                    f"{(record.energy_hartree - reference) * HARTREE_TO_KCAL_MOL:.8f}",
                ]
            )


def plot_energy(records: list[EnergyRecord], destination: Path) -> None:
    """Plot absolute and relative energy across DL-FIND cycles.

    Args:
        records: Ordered optimization records.
        destination: PNG output path.
    """
    cycles = [record.cycle for record in records]
    energies = [record.energy_hartree for record in records]
    reference = energies[0]
    relative = [(energy - reference) * HARTREE_TO_KCAL_MOL for energy in energies]
    figure, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(cycles, energies, color="#1f4e79", linewidth=1.2)
    axes[0].set_ylabel("Energy (Hartree)")
    axes[0].grid(alpha=0.25)
    axes[0].set_title("Unrestrained E-S DL-FIND minimization")
    axes[1].plot(cycles, relative, color="#a33b20", linewidth=1.2)
    axes[1].axhline(0.0, color="black", linewidth=0.7, alpha=0.5)
    axes[1].set_xlabel("DL-FIND optimization cycle")
    axes[1].set_ylabel("Relative energy (kcal/mol)\nfrom cycle 0")
    axes[1].grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def analyze(source: Path, output_dir: Path) -> dict[str, object]:
    """Extract, tabulate, and plot a snapshot of optimization energies.

    Args:
        source: Active Amber/DL-FIND output file.
        output_dir: Analysis destination directory.

    Returns:
        Analysis metadata.
    """
    payload = source.read_bytes()
    records = parse_energy_records(payload)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "optimization_energy.csv"
    plot_path = output_dir / "optimization_energy.png"
    write_csv(records, csv_path)
    plot_energy(records, plot_path)
    metadata: dict[str, object] = {
        "source": source.resolve().as_posix(),
        "source_snapshot_bytes": len(payload),
        "source_snapshot_sha256": hashlib.sha256(payload).hexdigest(),
        "first_cycle": records[0].cycle,
        "last_complete_cycle": records[-1].cycle,
        "complete_records": len(records),
        "reference_energy_hartree": records[0].energy_hartree,
        "final_energy_hartree": records[-1].energy_hartree,
        "net_change_kcal_mol": (records[-1].energy_hartree - records[0].energy_hartree)
        * HARTREE_TO_KCAL_MOL,
        "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
        "csv": csv_path.resolve().as_posix(),
        "plot": plot_path.resolve().as_posix(),
    }
    (output_dir / "optimization_energy.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    LOGGER.info(
        "Plotted %d complete cycles through cycle %d",
        len(records),
        records[-1].cycle,
    )
    return metadata


def main() -> int:
    """Parse command-line arguments and run the energy analysis.

    Returns:
        Process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    analyze(args.source, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
