#!/usr/bin/env python3
"""Rebuild a portable VMD bundle from accepted water-assisted path windows."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import parmed
from parmed.amber.mask import AmberMask
from water_path_vmd_support import (
    ang,
    dist,
    plot_svg,
    read_energy,
    sha256,
    render_manifest,
)

CPPTRAJ = Path("/data/jpc/code/ambertools26/venv/bin/cpptraj")
LOGGER = logging.getLogger(__name__)


def export(
    topology: Path,
    restarts: list[Path],
    energy_csv: Path,
    distance_csv: Path,
    output: Path,
) -> None:
    """Rebuild the cumulative VMD trajectory and provenance after accepted windows.

    Args:
        topology: Full-system ASH274 Amber topology.
        restarts: Accepted window restarts in path order.
        energy_csv: CSV with window and energy_kcal_mol columns.
        distance_csv: CSV with window and path distance/angle columns.
        output: Incrementally rebuilt portable VMD bundle directory.
    """
    if not restarts:
        raise ValueError("At least one accepted restart is required")
    if len({sha256(path) for path in restarts}) != len(restarts):
        raise ValueError("Duplicate restart content would duplicate a trajectory frame")
    if not CPPTRAJ.is_file():
        raise FileNotFoundError(CPPTRAJ)
    structure = parmed.load_file(str(topology))
    prep = (
        Path(__file__).parent.parent
        / "06_qmmm_optimization/replica_1_tier_A_seed/water_qm_reoptimization/preparation.json"
    )
    qm_mask = json.loads(prep.read_text(encoding="utf-8"))["qm_region"]["mask"]
    qm = [
        i for i, value in enumerate(AmberMask(structure, qm_mask).Selection()) if value
    ]
    keep = [
        i
        for i, value in enumerate(AmberMask(structure, ":1-1316,1576").Selection())
        if value
    ]
    if not qm or not set(qm) <= set(keep):
        raise ValueError("QM selection is empty or absent from stripped topology")
    energy = read_energy(energy_csv)
    with distance_csv.open(newline="", encoding="utf-8") as stream:
        distance_rows = list(csv.DictReader(stream))
    distances = {row["window"]: row for row in distance_rows}
    if len(distances) != len(distance_rows):
        raise ValueError("Duplicate window in distance ledger")
    names = [restart.parent.name for restart in restarts]
    if any(name not in energy or name not in distances for name in names):
        raise ValueError(
            "Energy and distance ledgers must contain each accepted window"
        )
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="water_path_export_") as tmp_name:
        tmp = Path(tmp_name)
        stripped = parmed.load_file(str(topology))
        stripped.strip("! :1-1316,1576")
        stripped.save(str(tmp / "path.parm7"), overwrite=True)
        lines = [
            f"parm {topology}",
            *(f"trajin {p}" for p in restarts),
            "strip !:1-1316,1576",
            f"trajout {tmp / 'path.dcd'} dcd",
            "run",
            "quit",
        ]
        control = tmp / "export.cpptraj"
        control.write_text("\n".join(lines) + "\n", encoding="utf-8")
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = "/data/jpc/code/ambertools26/venv/lib"
        with (tmp / "cpptraj.log").open("w", encoding="utf-8") as log:
            subprocess.run(
                [str(CPPTRAJ), "-i", str(control)],
                check=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
            )
        rows: list[dict[str, Any]] = []
        for name, restart in zip(names, restarts, strict=True):
            endpoint = parmed.load_file(str(topology), xyz=str(restart))
            xyz = np.asarray(endpoint.coordinates, dtype=float)
            if xyz.shape != (len(structure.atoms), 3) or not np.isfinite(xyz).all():
                raise ValueError(f"Invalid coordinates in {restart}")
            row: dict[str, Any] = {
                "window": name,
                **energy[name],
                "Thr_OG1_HG1_A": dist(xyz, endpoint, (929, "OG1"), (929, "HG1")),
                "Thr_HG1_water_O_A": dist(xyz, endpoint, (929, "HG1"), (1576, "O")),
                "water_H2_Asp272_OD1_A": dist(
                    xyz, endpoint, (1576, "H2"), (968, "OD1")
                ),
                "OG1_water_O_A": dist(xyz, endpoint, (929, "OG1"), (1576, "O")),
                "water_O_H2_A": dist(xyz, endpoint, (1576, "O"), (1576, "H2")),
                "Thr_OG1_HG1_water_O_deg": ang(
                    xyz, endpoint, (929, "OG1"), (929, "HG1"), (1576, "O")
                ),
                "water_O_H2_Asp272_OD1_deg": ang(
                    xyz, endpoint, (1576, "O"), (1576, "H2"), (968, "OD1")
                ),
            }
            row.update(
                {
                    key: value
                    for key, value in distances[name].items()
                    if key != "window" and value
                }
            )
            rows.append(row)
        with (tmp / "path.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        plot_svg(
            tmp / "energy.svg",
            "PES scan energy (not a PMF or barrier)",
            names,
            {
                "Physical Δ(EPtot−RESTRAINT)": [r["energy_kcal_mol"] for r in rows],
                "Bias-inclusive ΔEPtot": [
                    r["energy_bias_inclusive_kcal_mol"] for r in rows
                ],
            },
            "ΔE (kcal/mol)",
        )
        keys = [
            key for key in rows[0] if key.endswith("_A") and key != "energy_kcal_mol"
        ]
        plot_svg(
            tmp / "distances.svg",
            "Reactive and catalytic distances",
            names,
            {key: [float(r[key]) for r in rows] for key in keys},
            "Distance (Å)",
        )
        angle_keys = [key for key in rows[0] if key.endswith("_deg")]
        if angle_keys:
            plot_svg(
                tmp / "angles.svg",
                "Attack and relay angles",
                names,
                {key: [float(row[key]) for row in rows] for key in angle_keys},
                "Angle (degrees)",
            )
        index_map = {old: new for new, old in enumerate(keep)}
        mapped_qm = [index_map[i] for i in qm]
        manifest = {
            "state": "alternative Asp272-COOH / neutral-Thr-N TI1-like path",
            "source_code": {
                str(source.resolve()): sha256(source)
                for source in (
                    Path(__file__),
                    Path(__file__).with_name("water_path_vmd_support.py"),
                )
            },
            "topology": {"path": str(topology), "sha256": sha256(topology)},
            "accepted_restarts": [
                {"path": str(p), "sha256": sha256(p)} for p in restarts
            ],
            "energy_ledger": {"path": str(energy_csv), "sha256": sha256(energy_csv)},
            "distance_ledger": {
                "path": str(distance_csv),
                "sha256": sha256(distance_csv),
            },
            "energy_definition": {
                "energy_kcal_mol": "relative EPtot minus RESTRAINT (physical PES energy)",
                "energy_bias_inclusive_kcal_mol": "relative EPtot including path bias",
                "interpretation": "PES scan; not a PMF or validated barrier",
            },
            "frames": len(restarts),
            "stripped_atom_count": len(keep),
            "stripped_to_original_atom_indices_0based": keep,
            "qm_indices_stripped_0based": mapped_qm,
            "retained_water_residue_1based": 1576,
        }
        (tmp / "view_path.vmd").write_text(
            render_manifest(structure, manifest), encoding="utf-8"
        )
        (tmp / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        payload = [
            "path.parm7",
            "path.dcd",
            "path.csv",
            "energy.svg",
            "distances.svg",
            "view_path.vmd",
            "manifest.json",
            "cpptraj.log",
            "export.cpptraj",
        ]
        if angle_keys:
            payload.append("angles.svg")
        (tmp / "MANIFEST.sha256").write_text(
            "".join(f"{sha256(tmp / name)}  {name}\n" for name in payload),
            encoding="utf-8",
        )
        payload.append("MANIFEST.sha256")
        for name in payload:
            shutil.copy2(tmp / name, output / name)
    LOGGER.info("Rebuilt %d cumulative accepted frames in %s", len(restarts), output)


def main() -> None:
    """Parse arguments and create the VMD bundle."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument(
        "--accepted-restarts", type=Path, action="append", required=True
    )
    parser.add_argument("--energy-ledger", type=Path, required=True)
    parser.add_argument("--distance-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    export(
        args.topology.resolve(),
        [p.resolve() for p in args.accepted_restarts],
        args.energy_ledger.resolve(),
        args.distance_ledger.resolve(),
        args.output.resolve(),
    )


if __name__ == "__main__":
    main()
