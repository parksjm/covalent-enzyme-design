"""Plot retained forward minimization endpoints; this is no reaction barrier."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
import shutil
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

COMPONENTS = (
    "BOND",
    "ANGLE",
    "DIHED",
    "VDWAALS",
    "EEL",
    "HBOND",
    "1-4 VDW",
    "1-4 EEL",
    "RESTRAINT",
    "DFTBESCF",
    "CMAP",
)
NUMBER = r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)"


def main() -> None:
    """Validate endpoints, write auditable data, and plot energy and geometry."""
    source, output = (Path(value).resolve() for value in sys.argv[1:])
    status = json.loads((source / "status.json").read_text())
    assert status["state"] == "stopped"
    assert len(status["completed_windows"]) == 18
    assert not list(source.glob("reverse_*")), "Unexpected reverse artifacts"
    hashes = {}
    rows = []
    for index, name in enumerate(status["completed_windows"]):
        folder = source / name
        for filename in (
            "process.json",
            "acceptance.json",
            "minimized.mdout",
            "restraints.RST",
            "minimize.mdin",
        ):
            path = folder / filename
            hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        process = json.loads((folder / "process.json").read_text())
        accept = json.loads((folder / "acceptance.json").read_text())
        text = (folder / "minimized.mdout").read_text()
        block = text[text.rfind("   NSTEP") :]
        summary = re.search(
            r"NUMBER\s*\n\s*(\d+)\s+" + NUMBER + r"\s+" + NUMBER + r"\s+" + NUMBER,
            block,
        )
        assert summary, name
        row = {
            "window": index,
            "name": name,
            "final_cycle": int(summary[1]),
            "printed_energy_low_precision": float(summary[2]),
            "RMS": float(summary[3]),
            "GMAX": float(summary[4]),
        }
        for label in COMPONENTS:
            found = re.search(
                r"(?<![\w-])" + re.escape(label) + r"\s*=\s*" + NUMBER, block
            )
            assert found, (name, label)
            row[label] = float(found[1])
        row.update(accept["distances_A"])
        row["target_attack_A"] = process["target_attack_A"]
        row["target_proton_A"] = process["target_proton_A"]
        row["actual_proton_A"] = row["OG_H"] - row["N_H"]
        da = row["attack"] - row["target_attack_A"]
        dp = row["actual_proton_A"] - row["target_proton_A"]
        assert abs(da) <= 2 and abs(dp) <= 2
        row["attack_bias_kcal_mol"] = 25 * da**2
        row["proton_bias_kcal_mol"] = 25 * dp**2
        row["driving_bias_kcal_mol"] = 25 * (da**2 + dp**2)
        bond = re.search(r"NMR restraints: Bond\s*=\s*" + NUMBER, block)
        general = re.search(r"Gen\. Dist\. Coord\.\s*=\s*" + NUMBER, block)
        assert bond and general
        row["printed_nmr_bias_kcal_mol"] = float(bond[1]) + float(general[1])
        row["bias_crosscheck_error_kcal_mol"] = (
            row["driving_bias_kcal_mol"] - row["printed_nmr_bias_kcal_mol"]
        )
        assert abs(row["bias_crosscheck_error_kcal_mol"]) <= 0.002, name
        row["full_energy_kcal_mol"] = sum(row[key] for key in COMPONENTS)
        row["minus_driving_kcal_mol"] = (
            row["full_energy_kcal_mol"] - row["driving_bias_kcal_mol"]
        )
        rows.append(row)
    for row in rows:
        row["relative_full_kcal_mol"] = (
            row["full_energy_kcal_mol"] - rows[0]["full_energy_kcal_mol"]
        )
        row["relative_minus_driving_kcal_mol"] = (
            row["minus_driving_kcal_mol"] - rows[0]["minus_driving_kcal_mol"]
        )
    with (output / "forward_endpoint_energies.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    x = np.arange(len(rows))

    def values(key: str) -> list[float]:
        return [row[key] for row in rows]

    axes[0].plot(
        x,
        values("relative_full_kcal_mol"),
        "o-",
        label="Full system incl. driving bias",
    )
    axes[0].plot(
        x,
        values("relative_minus_driving_kcal_mol"),
        "s-",
        label="Minus driving bias (MM positional restraints retained)",
    )
    axes[0].set_ylabel("Relative energy (kcal/mol)\nreference: forward_00")
    axes[1].plot(x, values("driving_bias_kcal_mol"), "o-", label="Total driving bias")
    axes[1].plot(x, values("attack_bias_kcal_mol"), "--", label="Attack restraint")
    axes[1].plot(x, values("proton_bias_kcal_mol"), ":", label="Proton restraint")
    axes[1].set_ylabel("Driving energy (kcal/mol)")
    for actual, target, label, color in (
        ("attack", "target_attack_A", "OG1–C attack", "C0"),
        ("actual_proton_A", "target_proton_A", "OG1–H minus N–H", "C1"),
    ):
        axes[2].plot(x, values(actual), "o-", color=color, label=f"{label}: actual")
        axes[2].plot(x, values(target), "--", color=color, label=f"{label}: target")
    axes[2].set_ylabel("Coordinate (Å)")
    axes[2].set_xlabel("Forward window (sequential minimized endpoints)")
    axes[2].set_xticks(x)
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)
    fig.suptitle(
        "FAILED TI1 construction — forward endpoints only\n"
        "Not a barrier or free energy; reverse NOT RUN",
        color="darkred",
    )
    fig.text(
        0.5,
        0.015,
        "Endpoint: attack 2.608 Å vs 1.550 Å target; RMS 0.04279, "
        "GMAX 6.8292, 7 cycles.\nEnergies include the full solvated system; "
        "subtracting driving bias does not validate a reaction path.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 0.95))
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"forward_endpoint_energies.{suffix}", dpi=180)
    plt.close(fig)
    script = Path(__file__).resolve()
    hashes[str(script)] = hashlib.sha256(script.read_bytes()).hexdigest()
    hashes[str(source / "status.json")] = hashlib.sha256(
        (source / "status.json").read_bytes()
    ).hexdigest()
    shutil.copy2(script, output / "plotting_job" / "plot_energies.snapshot.py")
    report = {
        "scientific_status": "FAILED TI1 construction; reverse NOT RUN",
        "process_acceptance": "18 endpoints parsed; bias crosschecks passed",
        "scientific_acceptance": False,
        "source_status": status,
        "component_sum": list(COMPONENTS),
        "rows": rows,
        "input_sha256": hashes,
        "python": sys.executable,
        "matplotlib": matplotlib.__version__,
        "numpy": np.__version__,
        "caveat": "Minus driving retains MM positional restraints; not barrier/free energy",
    }
    (output / "forward_endpoint_energies.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    logging.info(
        "Wrote plots and data; max bias residual %.6g kcal/mol",
        max(abs(row["bias_crosscheck_error_kcal_mol"]) for row in rows),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
