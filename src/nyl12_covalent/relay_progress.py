"""Accepted-only relay campaign inspection helpers; never evaluate energies."""

from __future__ import annotations

import csv
import importlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
import parmed
from parmed.amber.mask import AmberMask

from .relay_progress_vmd import vmd_script

ROOT = Path(__file__).resolve().parents[2]
COMPARISON = ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison"
EXPORT_SOURCE = COMPARISON / "07_water_to_asp272_ti1_like"


def legacy_exporter() -> Any:
    """Load the retained exporter from its original directory."""
    sys.path.insert(0, str(EXPORT_SOURCE))
    return importlib.import_module("export_water_path")


def accepted_rows(campaign: Path) -> list[dict[str, Any]]:
    """Read complete accepted JSONL records in stable stage/iteration order."""
    rows = []
    for path in sorted(campaign.rglob("accepted_iterations.jsonl")):
        if "progress" in path.relative_to(campaign).parts:
            continue
        previous = -1
        for line in path.read_text().splitlines(keepends=True):
            if not line.endswith("\n"):
                continue
            row = json.loads(line)
            iteration = row["iteration"]
            if type(iteration) is not int or iteration <= previous:
                raise ValueError(f"Nonmonotonic accepted iterations: {path}")
            previous = iteration
            energies = [
                row[key]
                for key in (
                    "physical_energy_kcal_mol",
                    "bias_energy_kcal_mol",
                    "total_energy_kcal_mol",
                )
            ]
            if not np.isfinite(energies).all():
                raise ValueError(f"Nonfinite accepted energy: {path}")
            rows.append({**row, "stage": str(path.parent.relative_to(campaign))})
    return rows


def snapshot_frames(
    campaign: Path, rows: list[dict[str, Any]]
) -> list[tuple[Path, dict[str, Any]]]:
    """Select stage initial/every25 and completed final accepted geometries."""
    stages: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        stages.setdefault(row["stage"], []).append(row)
    frames = []
    for stage, records in stages.items():
        directory = campaign / stage
        for row in records:
            iteration = row["iteration"]
            if iteration == 0 or iteration % 25 == 0:
                path = directory / f"accepted_{iteration:05d}.ncrst"
                if path.is_file():
                    frames.append((path, row))
        final = directory / "optimized.ncrst"
        if (directory / "process.json").is_file() and final.is_file():
            frames.append((final, records[-1]))
    return frames


def write_table(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a scalar CSV inspection view while JSON retains nested evidence."""
    flat = [
        {
            **{k: v for k, v in row.items() if not isinstance(v, (dict, list))},
            **{k: v for k, v in row["state"].items() if isinstance(v, (int, float))},
        }
        for row in rows
    ]
    fields = list(dict.fromkeys(key for row in flat for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flat)


def plot_progress(output: Path, rows: list[dict[str, Any]]) -> None:
    """Plot stage-labelled physical/bias energies, reaction distances and CNs."""
    plot = legacy_exporter().plot_svg
    labels = [f"{r['stage']}:{r['iteration']}" for r in rows]
    baseline = rows[0]["physical_energy_kcal_mol"]
    plot(
        output / "energy.svg",
        "Accepted states; diagnostic PES, not PMF/barrier",
        labels,
        {
            "Physical ΔE": [r["physical_energy_kcal_mol"] - baseline for r in rows],
            "Total ΔE": [r["total_energy_kcal_mol"] - baseline for r in rows],
            "Bias E": [r["bias_energy_kcal_mol"] for r in rows],
        },
        "kcal/mol",
    )
    for filename, keys, title, unit in (
        (
            "reaction_distances.svg",
            ("rOG1C", "rCN", "rCO", "rNwater", "rwaterLeave"),
            "Accepted reactive distances by stage",
            "Å",
        ),
        (
            "coordination.svg",
            ("CN_N", "CN_OG1", "CN_water", "CN_leave"),
            "All-real-QM-H smooth coordination by stage",
            "Smooth count",
        ),
    ):
        plot(
            output / filename,
            title,
            labels,
            {key: [r["state"][key] for r in rows] for key in keys},
            unit,
        )


def verify_mapping(case: dict[str, Any]) -> Any:
    """Require the exporter QM102/54H/retained-water mapping to equal the case."""
    topology = parmed.load_file(case["topology"])
    prep = (
        COMPARISON
        / "06_qmmm_optimization/replica_1_tier_A_seed/water_qm_reoptimization/preparation.json"
    )
    mask = json.loads(prep.read_text())["qm_region"]["mask"]
    qm = list(AmberMask(topology, mask).Selected())
    declared = [int(i) - 1 for i in case["qm_atom_ids"]]
    h = [i for i in qm if topology.atoms[i].atomic_number == 1]
    if qm != declared or len(qm) != 102 or len(h) != 54:
        raise ValueError("Case differs from retained exporter QM102/54H mapping")
    water = [a.idx + 1 for a in topology.residues[1575].atoms]
    if water != [20317, 20318, 20319] or not set(water) <= set(case["qm_atom_ids"]):
        raise ValueError("Expected unchanged retained water1576 atom IDs20317-20319")
    return topology


def export_movie(
    output: Path, frames: list[tuple[Path, dict[str, Any]]], case: dict[str, Any]
) -> None:
    """Reuse validated coordinate stripping after strict QM/water map checks."""
    if not frames:
        return
    exporter = legacy_exporter()
    topology = verify_mapping(case)
    links, energies, distances, records = [], [], [], []
    seen: dict[str, int] = {}
    baseline = frames[0][1]["physical_energy_kcal_mol"]
    for path, row in frames:
        identity = exporter.sha256(path)
        if identity in seen:
            records.append(
                {
                    "stage": row["stage"],
                    "iteration": row["iteration"],
                    "source_path": str(path),
                    "movie_frame_0based": seen[identity],
                    "sha256": identity,
                }
            )
            continue
        seen[identity] = len(links)
        records.append(
            {
                "stage": row["stage"],
                "iteration": row["iteration"],
                "source_path": str(path),
                "movie_frame_0based": len(links),
                "sha256": identity,
            }
        )
        name = f"frame_{len(links):05d}"
        directory = output / "frame_links" / name
        directory.mkdir(parents=True)
        link = directory / "accepted.ncrst"
        link.symlink_to(path.resolve())
        links.append(link)
        energies.append(
            {
                "window": name,
                "energy_kcal_mol": row["physical_energy_kcal_mol"] - baseline,
                "energy_bias_inclusive_kcal_mol": row["total_energy_kcal_mol"]
                - baseline,
            }
        )
        distances.append(
            {
                "window": name,
                **{
                    k: v for k, v in row["state"].items() if isinstance(v, (int, float))
                },
            }
        )
    for filename, rows in (
        ("snapshot_energies.csv", energies),
        ("snapshot_coordinates.csv", distances),
    ):
        with (output / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    bundle = output / "vmd"
    exporter.export(
        Path(case["topology"]),
        links,
        output / "snapshot_energies.csv",
        output / "snapshot_coordinates.csv",
        bundle,
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.update(
        {
            "state": "Accepted diagnostic states; primary water screen negative; endpoint test independent from accepted ES; no validated relay/barrier claim",
            "accepted_stage_frame_records": records,
            "accepted_restarts": [
                {"path": record["source_path"], "sha256": record["sha256"]}
                for record in records
            ],
            "coordination_definition": "All54 real QM hydrogens; no fixed proton identities",
        }
    )
    (bundle / "view_path.vmd").write_text(vmd_script(manifest, case, topology))
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    checksum = bundle / "MANIFEST.sha256"
    names = [line.split("  ", 1)[1] for line in checksum.read_text().splitlines()]
    checksum.write_text(
        "".join(f"{exporter.sha256(bundle / name)}  {name}\n" for name in names)
    )
