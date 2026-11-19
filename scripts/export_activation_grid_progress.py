"""Export an incremental relaxed grid and ordered, per-row VMD trajectories."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
GRID = ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison/08_activation_attack_grid"
sys.path.insert(0, str(GRID.parent / "07_water_to_asp272_ti1_like"))
from export_water_path import export  # noqa: E402

LOGGER = logging.getLogger(__name__)


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write a table with stable fields shared by all rows."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def correct_vmd_bundle(folder: Path) -> None:
    """Add N–HG1 dynamics and remove static bonds involving transferring protons."""
    manifest = json.loads((folder / "manifest.json").read_text())
    mapping = {
        old: new
        for new, old in enumerate(manifest["stripped_to_original_atom_indices_0based"])
    }
    n, og, hg, wo, wh = [mapping[i - 1] for i in (13905, 13916, 13917, 20317, 20319)]
    path = folder / "view_path.vmd"
    script = path.read_text()
    cleanup = (
        "set sel [atomselect $molid all]\nset bonds [$sel getbonds]\n"
        f"foreach pair {{{{{og} {hg}}} {{{n} {hg}}} {{{wo} {hg}}} {{{wo} {wh}}}}} {{\n"
        " lassign $pair a b\n"
        " foreach {i j} [list $a $b $b $a] {\n"
        "  set keep {}; foreach k [lindex $bonds $i] {if {$k != $j} {lappend keep $k}}\n"
        "  lset bonds $i $keep\n }\n}\n$sel setbonds $bonds\n$sel delete\n"
    )
    script = script.replace("mol delrep 0 $molid\n", "mol delrep 0 $molid\n" + cleanup)
    for a, b in ((n, hg), (wo, wh)):
        script += (
            f"\nmol representation DynamicBonds 1.4 0.1 1\n"
            f'mol selection "index {a} {b}"\nmol addrep $molid\n'
            f"label add Bonds $molid/{a} $molid/{b}\n"
        )
    path.write_text(
        script.replace(
            "Alternative Asp272-COOH / neutral-Thr-N path",
            "Independent activation/attack grid row",
        )
    )
    manifest["state"] = (
        "independent q_N activation/attack grid row; inspect proton inventory"
    )
    manifest["grid_exporter_sha256"] = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    checksum = folder / "MANIFEST.sha256"
    names = [line.split("  ", 1)[1] for line in checksum.read_text().splitlines()]
    checksum.write_text(
        "".join(
            f"{hashlib.sha256((folder / name).read_bytes()).hexdigest()}  {name}\n"
            for name in names
        )
    )


def refresh() -> bool:
    """Refresh grid plots and newly completed row frames; return terminal status."""
    protocol = json.loads((GRID / "protocol.json").read_text())
    data = json.loads((GRID / "grid_points.json").read_text())["points"]
    output = GRID / "progress"
    output.mkdir(exist_ok=True)
    completed = [
        (key, value) for key, value in data.items() if value["status"] == "completed"
    ]
    if completed:
        write_csv(
            output / "grid_points.csv",
            [{"point": key, **value} for key, value in completed],
        )
        baseline = min(value["physical_energy_kcal_mol"] for _, value in completed)
        fig, ax = plt.subplots(figsize=(7, 5))
        sc = ax.scatter(
            [v["attack_actual_A"] for _, v in completed],
            [v["q_actual_A"] for _, v in completed],
            c=[v["physical_energy_kcal_mol"] - baseline for _, v in completed],
            s=140,
            cmap="viridis",
        )
        ax.set(
            xlabel="Measured OG1–carbonyl C (Å)",
            ylabel="Measured q_N (Å)",
            title=f"Relaxed PES grid: {len(completed)}/25 points",
        )
        fig.colorbar(sc, ax=ax, label="Energy above lowest completed point (kcal/mol)")
        fig.tight_layout()
        fig.savefig(output / "grid_heatmap.png", dpi=180)
        plt.close(fig)
    topology = GRID.parent / "00_system/ash274_solvated.prmtop"
    terminal = []
    for row in range(1, 6):
        row_dir = GRID / f"row_{row:03d}"
        plan = json.loads((row_dir / "row.json").read_text())["stages"]
        records = []
        restarts = []
        hashes = set()
        for stage in plan:
            folder = GRID / stage["directory"]
            acceptance = folder / "acceptance.json"
            if not acceptance.exists():
                break
            record = json.loads(acceptance.read_text())
            if not record.get("accepted"):
                break
            digest = record["restart_sha256"]
            if digest not in hashes:
                records.append((folder.name, record))
                restarts.append(folder / "optimized.ncrst")
                hashes.add(digest)
        terminal.append(
            (row_dir / plan[-1]["stage"] / "acceptance.json").exists()
            or any(
                (GRID / s["directory"] / "acceptance.json").exists()
                and not json.loads(
                    (GRID / s["directory"] / "acceptance.json").read_text()
                ).get("accepted")
                for s in plan
            )
        )
        if not records:
            continue
        bundle = output / f"row_{row:03d}"
        bundle.mkdir(exist_ok=True)
        signature = "\n".join(record["restart_sha256"] for _, record in records)
        marker = bundle / "exported_frames.txt"
        if marker.exists() and marker.read_text() == signature:
            continue
        base = records[0][1]["energies"]["physical_energy_kcal_mol"]
        write_csv(
            bundle / "energies.csv",
            [
                {
                    "window": name,
                    "energy_kcal_mol": rec["energies"]["physical_energy_kcal_mol"]
                    - base,
                    "energy_bias_inclusive_kcal_mol": rec["energies"][
                        "bias_inclusive_energy_kcal_mol"
                    ]
                    - base,
                }
                for name, rec in records
            ],
        )
        write_csv(
            bundle / "distances.csv",
            [{"window": name, **rec["geometry"]} for name, rec in records],
        )
        export(
            topology,
            restarts,
            bundle / "energies.csv",
            bundle / "distances.csv",
            bundle,
        )
        correct_vmd_bundle(bundle)
        marker.write_text(signature)
    LOGGER.info(
        "Exported %d grid points; terminal rows %d/5; protocol %s",
        len(completed),
        sum(terminal),
        protocol["status"],
    )
    return all(terminal)


def main() -> None:
    """Refresh once or monitor with a bounded deadline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    deadline = time.monotonic() + 48 * 3600
    while True:
        if refresh() or not args.watch:
            return
        if time.monotonic() >= deadline:
            raise TimeoutError("Grid export monitor exceeded 48-hour bound")
        time.sleep(120)


if __name__ == "__main__":
    main()
