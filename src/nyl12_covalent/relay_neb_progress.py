"""Publish coherent relay NEB image geometry and compact inspection bundles."""

from __future__ import annotations

from dataclasses import asdict
import gzip
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
from typing import Any

import numpy as np
from parmed.amber import Rst7

from .active_site_geometry import angle_degrees
from .neb_progress import HARTREE_TO_KCAL_MOL, NebReport
from .relay_case import load_case, state_record
from .relay_progress import EXPORT_SOURCE, export_movie, legacy_exporter, write_table

ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    """Return the SHA-256 of a retained monitor source or compact artifact."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    """Atomically replace a small JSON record in its destination filesystem."""
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def source_paths() -> list[Path]:
    """Return the precise monitor, snapshot and movie sources retained at milestones."""
    package = Path(__file__).parent
    return [
        ROOT / "scripts/watch_relay_acyl_neb.py",
        ROOT / "scripts/run_relay_acyl_campaign.py",
        Path(__file__).resolve(),
        *[
            package / name
            for name in (
                "neb_xyz_snapshot.py",
                "neb_snapshot.py",
                "neb_progress.py",
                "relay_case.py",
                "relay_progress.py",
                "relay_progress_vmd.py",
                "proton_coordination.py",
            )
        ],
        EXPORT_SOURCE / "export_water_path.py",
        EXPORT_SOURCE / "water_path_vmd_support.py",
    ]


def retain_sources(output: Path, identity: dict[str, Any], paths: list[Path]) -> Path:
    """Retain one immutable source milestone per distinct source/input identity."""
    root = output / "source_milestones"
    root.mkdir(exist_ok=True)
    for manifest in sorted(root.glob("milestone_*/manifest.json")):
        if json.loads(manifest.read_text()) == identity:
            return manifest.parent
    folder = root / f"milestone_{len(list(root.glob('milestone_*'))) + 1:04d}"
    folder.mkdir()
    for index, path in enumerate(paths):
        shutil.copy2(path, folder / f"{index:02d}_{path.name}")
    atomic_json(folder / "manifest.json", identity)
    return folder


def image_geometry(atom_map: dict[str, int], xyz: np.ndarray) -> dict[str, float]:
    """Measure mapped transfer-H, oxyanion and attack geometry for one NEB image.

    Args:
        atom_map: Explicit case mapping with one-based native atom identifiers.
        xyz: Full-system coordinates in zero-based topology order, in angstrom.

    Returns:
        Scalar distances and BD angle; water_h1 is distinct from original HG1.
    """
    indices = {key: index - 1 for key, index in atom_map.items()}
    pairs = {
        "rWaterH1_OG1": ("water_h1", "thr_og"),
        "rWaterH1_Nleave": ("water_h1", "leaving_n"),
        "rO_PheN": ("oxygen", "phe_n"),
        "rO_AsnND2": ("oxygen", "asn_nd2"),
    }
    result = {
        label: float(np.linalg.norm(xyz[indices[a]] - xyz[indices[b]]))
        for label, (a, b) in pairs.items()
    }
    result["BD_angle_deg"] = angle_degrees(
        xyz[indices["thr_og"]], xyz[indices["carbon"]], xyz[indices["oxygen"]]
    )
    return result


def publish_band(
    config: dict[str, Any],
    band: Path,
    report: NebReport,
    output: Path,
    identity: dict[str, Any],
    terminal: bool,
    milestone: Path,
) -> dict[str, Any]:
    """Publish image-matched energies, geometry, stripped movie and provenance.

    Args:
        config: Explicit NEB configuration including the original case_path.
        band: Validated temporary full-system XYZ reconstruction directory.
        report: Native physical image energies proven associated with this band.
        output: NEB progress directory; existing generations remain immutable.
        identity: Complete report/config/source identity for deduplication.
        terminal: Whether the calculation launcher has ended.
        milestone: Retained monitor source/input snapshot directory.

    Returns:
        Compact progress record after the latest symlink switches atomically.
    """
    case = load_case(Path(config["case_path"]))
    number = len(list(output.glob("update_[0-9]*"))) + 1
    name = f"update_{number:06d}"
    pending = output / f".building_{number:06d}"
    pending.mkdir()
    rows, frames = [], []
    for image in report.images:
        path = band / f"optimized.{image.image:04d}.ncrst"
        xyz = np.asarray(Rst7.open(str(path)).coordinates).reshape(-1, 3)
        energy = image.energy_hartree * HARTREE_TO_KCAL_MOL
        row = {
            "stage": f"NEB_cycle_{report.cycle}",
            "iteration": image.image,
            "image": image.image,
            "cycle": report.cycle,
            "physical_energy_kcal_mol": energy,
            "total_energy_kcal_mol": energy,
            "bias_energy_kcal_mol": 0.0,
            "relative_energy_kcal_mol": image.relative_energy_kcal_mol,
            "state": {**state_record(case, xyz), **image_geometry(case.atom_map, xyz)},
        }
        rows.append(row)
        frames.append((path, row))
    atomic_json(pending / "images.json", rows)
    write_table(pending / "images.csv", rows)
    plot = legacy_exporter().plot_svg
    labels = [str(row["image"]) for row in rows]
    plot(
        pending / "energy.svg",
        "NEB path candidate; not a verified TS or PMF",
        labels,
        {"Physical image ΔE": [r["relative_energy_kcal_mol"] for r in rows]},
        "kcal/mol",
    )
    for filename, keys, title, unit in [
        (
            "reaction_distances.svg",
            (
                "rOG1C",
                "rCN",
                "rCO",
                "rNwater",
                "rwaterLeave",
                "rWaterH1_OG1",
                "rWaterH1_Nleave",
                "rO_PheN",
                "rO_AsnND2",
            ),
            "NEB path candidate: reactive image distances",
            "Å",
        ),
        (
            "attack_angle.svg",
            ("BD_angle_deg",),
            "NEB path candidate: OG1-C-O Burgi-Dunitz angle",
            "Degrees",
        ),
        (
            "coordination.svg",
            ("CN_N", "CN_OG1", "CN_water", "CN_leave"),
            "NEB path candidate: all-real-QM-H coordination",
            "Smooth count",
        ),
    ]:
        plot(
            pending / filename,
            title,
            labels,
            {key: [r["state"][key] for r in rows] for key in keys},
            unit,
        )
    export_movie(pending, frames, case.config)
    movie_manifest = pending / "vmd/manifest.json"
    movie_record = json.loads(movie_manifest.read_text())
    movie_record["state"] = (
        "Coherent NEB path candidate; not a verified transition state or PMF"
    )
    movie_record["image_energy_coordinate_association"] = (
        "Ordered native XYZ save loop with fixed-reference reconstruction"
    )
    atomic_json(movie_manifest, movie_record)
    checksum = pending / "vmd/MANIFEST.sha256"
    names = [line.split("  ", 1)[1] for line in checksum.read_text().splitlines()]
    checksum.write_text(
        "".join(f"{digest(checksum.parent / name)}  {name}\n" for name in names)
    )
    metadata = json.loads((band / "band_snapshot.json").read_text())
    atomic_json(pending / "band_snapshot.json", metadata)
    with gzip.open(pending / "source_nebpath.xyz.gz", "wb") as stream:
        stream.write((band / "source_nebpath.xyz").read_bytes())
    output_bytes = (band / "source_output_snapshot.mdout").read_bytes()
    (pending / "native_report.mdout").write_bytes(
        output_bytes[report.report_start_byte :]
    )
    record = {
        "identity": identity,
        "terminal": terminal,
        "latest_update": name,
        "report_index": report.report_index,
        "cycle": report.cycle,
        "regular_images": len(rows),
        "source_milestone": str(milestone),
        "report": asdict(report),
        "interpretation": "Path candidate, not a verified transition state, barrier or PMF.",
        "coordinate_source": "Proven ordered native XYZ save loop; numbered native restart files were not read.",
        "energy_reference": "Physical image energy relative to regular image1; NEB springs excluded.",
        "reconstructed_restarts": "Temporary derived NetCDFs; hashes in band_snapshot.json. Regenerate from retained native XYZ and fixed reference.",
        "retained_native_xyz_sha256": metadata["source_xyz_sha256"],
        "local_native_evidence_sha256": {
            name: digest(pending / name)
            for name in ("native_report.mdout", "source_nebpath.xyz.gz")
        },
        "geometry_atom_indices_0based": {
            key: case.atom_map[key] - 1
            for key in (
                "water_h1",
                "proton",
                "thr_og",
                "leaving_n",
                "carbon",
                "oxygen",
                "phe_n",
                "asn_nd2",
            )
        },
        "transfer_hydrogen_definition": "water_h1 (native20318 in the primary case) tracks actual OG1-to-leaving-N transfer; original ThrHG1 is separate atom_map.proton (native13917), transferred onto water. All-H ownership remains audited independently.",
    }
    if terminal:
        from scripts.run_relay_acyl_campaign import classify

        status = json.loads((output.parent / "job/status.json").read_text())
        native_text = output_bytes.decode("ascii")
        errors = re.findall(
            r"(?:SCC|SCF)[^\n]*(?:failed|not converg|failure)[^\n]*", native_text, re.I
        )
        endpoint = classify(case, xyz)
        endpoint.pop("state")  # Full ownership/state evidence remains in images.json.
        maximum = max(report.images, key=lambda image: image.relative_energy_kcal_mol)
        acceptance = {
            "process": status,
            "native_optimisation_converged": "Optimisation converged" in native_text,
            "scc_error_count": len(errors),
            "scc_error_lines": errors[-10:],
            "endpoint_image": len(rows),
            "provisional_endpoint_chemistry": endpoint,
            "curvature_verified": False,
            "transition_state_verified": False,
            "profile_maximum_image": maximum.image,
            "profile_maximum_relative_kcal_mol": maximum.relative_energy_kcal_mol,
            "profile_maximum_is_ts_barrier": False,
            "automatic_ts_launch": False,
            "interpretation": "Native path candidate only. Process, NEB convergence, endpoint chemistry and curvature verification are distinct; profile maximum is not a validated TS barrier or free energy.",
        }
        atomic_json(pending / "scientific_acceptance.json", acceptance)
        atomic_json(output.parent / "scientific_acceptance.json", acceptance)
    atomic_json(pending / "progress.json", record)
    # Temporary reconstructed-frame links are intentionally absent from the laptop payload.
    payload = [
        p
        for p in pending.rglob("*")
        if p.is_file()
        and not p.is_symlink()
        and "frame_links" not in p.relative_to(pending).parts
        and p.name not in {"native_report.mdout", "source_nebpath.xyz.gz"}
    ]
    (pending / "MANIFEST.sha256").write_text(
        "".join(
            f"{digest(path)}  {path.relative_to(pending)}\n" for path in sorted(payload)
        )
    )
    pending.rename(output / name)
    pointer = output / f".latest_{os.getpid()}"
    pointer.symlink_to(name, target_is_directory=True)
    pointer.replace(output / "latest")
    atomic_json(output / "progress.json", record)
    return record


def record_terminal_failure(config_path: Path, error: Exception) -> Path | None:
    """Record a terminal export failure without claiming unavailable band chemistry.

    Args:
        config_path: NEB configuration whose launcher status must prove completion.
        error: Unexpected exporter or native-band validation failure.

    Returns:
        Written stage acceptance path, or None while the native run remains live.
    """
    status_path = config_path.parent / "job/status.json"
    if not status_path.is_file():
        return None
    status = json.loads(status_path.read_text())
    if status["ended_utc"] is None or type(status["exit_status"]) is not int:
        return None
    mdout = config_path.parent / "calculation/minimized.mdout"
    text = mdout.read_text() if mdout.is_file() else ""
    errors = re.findall(
        r"(?:SCC|SCF)[^\n]*(?:failed|not converg|failure)[^\n]*", text, re.I
    )
    record = {
        "process": status,
        "native_optimisation_converged": "Optimisation converged" in text,
        "scc_error_count": len(errors),
        "scc_error_lines": errors[-10:],
        "coherent_final_band_exported": False,
        "provisional_endpoint_chemistry": None,
        "curvature_verified": False,
        "transition_state_verified": False,
        "profile_maximum_is_ts_barrier": False,
        "automatic_ts_launch": False,
        "export_error": f"{type(error).__name__}: {error}",
    }
    path = config_path.parent / "scientific_acceptance.json"
    atomic_json(path, record)
    return path
