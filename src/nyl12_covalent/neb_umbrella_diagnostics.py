"""Two-dimensional sample/bias diagnostics and compact endpoint-only umbrella views."""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7
from scipy.io import netcdf_file

from .neb_umbrella_bias import HarmonicWindow, bias_energy
from .neb_umbrella_campaign import digest, write_json


def trajectory_cvs(path: Path, roles: dict[str, int]) -> np.ndarray:
    """Read both CVs from native coordinate frames, never assume a DUMPAVE layout."""
    ids = [roles[k] - 1 for k in ("carbon", "leaving_n", "thr_og", "proton")]
    with netcdf_file(path, "r", mmap=False) as stream:
        xyz = np.array(stream.variables["coordinates"][:, ids, :], dtype=float)
    if xyz.ndim != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("Missing/nonfinite trajectory coordinate samples")
    c, n, o, h = (xyz[:, i, :] for i in range(4))

    def distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.linalg.norm(a - b, axis=1)

    return np.column_stack(
        (distance(c, n) - distance(c, o), distance(o, h) - distance(n, h))
    )


def block_diagnostics(
    cfg: dict[str, Any], row: dict[str, Any], stage: Path, topology: Any
) -> None:
    """Keep every valid frame and evaluate its exact bias at every target window."""
    import re
    from parmed.amber.mask import AmberMask
    from .neb_umbrella_observables import extract_observables

    qm_mask = re.findall(
        r"\bqmmask\s*=\s*'([^']+)'", (stage / "run.mdin").read_text(), re.I
    )[0]
    qm_ids = [i + 1 for i in AmberMask(topology, qm_mask).Selected()]
    extract_observables(
        Path(cfg["topology"]),
        stage / "trajectory.nc",
        cfg["roles_1based"],
        qm_ids,
        stage / "observables",
    )
    with (stage / "observables/observables.csv").open() as stream:
        cv = np.array(
            [[float(r["q_acyl"]), float(r["q_h"])] for r in csv.DictReader(stream)]
        )
    windows = [HarmonicWindow(**r["bias"]) for r in cfg["windows"]]
    energy = np.array(
        [[bias_energy(tuple(sample), window) for sample in cv] for window in windows]
    )
    beta = 1 / (0.00198720425864083 * cfg["temperature_K"])
    pending = stage / "samples_and_cross_bias.pending.npz"
    np.savez_compressed(
        pending, cv_values=cv, bias_kcal_mol=energy, u_kn_bias_only=beta * energy
    )
    pending.replace(stage / "samples_and_cross_bias.npz")
    write_json(
        stage / "analysis_request.json",
        {
            "samples": len(cv),
            "cv_names": ["C-Nleave minus C-OG1", "OG1-HG1 minus Nleave-HG1"],
            "bias_matrix_shape": list(energy.shape),
            "sample_source_sha256": digest(stage / "trajectory.nc"),
            "bias_convention": cfg["bias_convention"],
            "sample_filtering": "NONE: protonbranches/shared/unassignedHretained",
            "2D_MBAR": "Pendingdedicatedanalysis/correlation/overlapreview; no1DPMFparserused",
            "common_physical_energy": "SameHamiltonian; commonper-framephysicaltermcancelsinMBARweights",
            "all_window_bias_definitions": [r["bias"] for r in cfg["windows"]],
        },
    )
    # Only scalar diagnostic data are read across windows; each writer owns its file.
    profiles = []
    for other in cfg["windows"]:
        paths = sorted(
            stage.parents[1].glob(
                f"{other['label']}/03_production*/samples_and_cross_bias.npz"
            )
        )
        if paths:
            samples = []
            for p in paths:
                with np.load(p) as data:
                    samples.append(data["cv_values"].copy())
            profiles.append((other["label"], np.concatenate(samples)))
    bounds = np.array(
        [
            [
                min(r["q_acyl_center"] for r in cfg["windows"]) - 0.5,
                max(r["q_acyl_center"] for r in cfg["windows"]) + 0.5,
            ],
            [
                min(r["q_h_center"] for r in cfg["windows"]) - 0.5,
                max(r["q_h_center"] for r in cfg["windows"]) + 0.5,
            ],
        ]
    )
    all_samples = np.concatenate([values for _, values in profiles])
    bounds[:, 0] = np.minimum(bounds[:, 0], all_samples.min(axis=0) - 0.01)
    bounds[:, 1] = np.maximum(bounds[:, 1], all_samples.max(axis=0) + 0.01)
    histograms = [
        np.histogram2d(values[:, 0], values[:, 1], bins=40, range=bounds)[0]
        for _, values in profiles
    ]
    overlaps = []
    for i in range(len(profiles) - 1):
        a, b = histograms[i], histograms[i + 1]
        overlaps.append(
            {
                "first": profiles[i][0],
                "second": profiles[i + 1][0],
                "histogram_intersection": float(
                    np.minimum(a / a.sum(), b / b.sum()).sum()
                ),
            }
        )
    write_json(
        stage / "overlap_checkpoint.json",
        {
            "available_windows": len(profiles),
            "overlaps": overlaps,
            "interpretation": "Diagnosticfinite-sample2Dhistogramintersection; notsamplingconvergence",
        },
    )
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    fig, axis = plt.subplots(figsize=(6, 5))
    axis.scatter(cv[:, 0], cv[:, 1], s=8, alpha=0.5)
    axis.scatter(
        [r["q_acyl_center"] for r in cfg["windows"]],
        [r["q_h_center"] for r in cfg["windows"]],
        s=9,
        marker="x",
        color="black",
    )
    axis.set(
        xlabel="qAcyl (Å)",
        ylabel="qH (Å)",
        title=f"{row['label']}: all production-block samples",
    )
    fig.tight_layout()
    fig.savefig(stage / "cv_coverage.png", dpi=160)
    plt.close(fig)


def export_view(
    cfg: dict[str, Any], row: dict[str, Any], folder: Path, topology: Any
) -> None:
    """Publish immutable portable PDB/XYZ views of production block endpoints only."""
    stages = sorted(
        p for p in folder.glob("03_production*") if (p / "acceptance.json").is_file()
    )
    if not stages:
        return
    roles = cfg["roles_1based"]
    seed = Rst7.open(row["seed"])
    topology.coordinates = np.asarray(seed.coordinates).reshape(-1, 3)
    # Source QM98 plus the complete P6 substrate; explicit atom indices, no residue offsets.
    import re

    text = (folder / row["stages"][0]["label"] / "run.mdin").read_text()
    qm_mask = re.findall(r"\bqmmask\s*=\s*'([^']+)'", text, re.I)[0]
    from parmed.amber.mask import AmberMask

    qm = set(AmberMask(topology, qm_mask).Selected())
    from .relay_progress_vmd import substrate_indices

    substrate = substrate_indices(topology, roles["carbon"] - 1, roles["leaving_n"] - 1)
    keep = sorted(qm | substrate)
    mapping = {old: new for new, old in enumerate(keep)}
    selection = "@" + ",".join(str(i + 1) for i in keep)
    compact = topology[selection]
    frames = [
        np.asarray(Rst7.open(str(p / "restart.ncrst")).coordinates).reshape(-1, 3)[keep]
        for p in stages
    ]
    output = folder / "progress"
    output.mkdir(exist_ok=True)
    generation = output / f"block_{len(stages):02d}"
    generation.mkdir()
    compact.coordinates = frames[0]
    compact.save(str(generation / "view.pdb"), overwrite=False)
    if len(frames) > 1:
        with (generation / "endpoints.xyz").open("w") as stream:
            for xyz in frames[1:]:
                stream.write(f"{len(keep)}\nProductionblockendpoint:notMDtrajectory\n")
                for atom, point in zip(compact.atoms, xyz, strict=True):
                    stream.write(
                        f"{atom.element_name} {point[0]:.7f} {point[1]:.7f} {point[2]:.7f}\n"
                    )
    hydrogens = {i for i in qm if topology.atoms[i].atomic_number == 1}
    hetero = {i for i in qm if topology.atoms[i].atomic_number in (7, 8)}
    c, n, o = (roles[k] - 1 for k in ("carbon", "leaving_n", "thr_og"))
    bonds = []
    for original in keep:
        neighbors = []
        for atom in topology.atoms[original].bond_partners:
            j = atom.idx
            if j not in mapping:
                continue
            if (
                (original in hydrogens and j in hetero)
                or (j in hydrogens and original in hetero)
                or {original, j} in ({c, n}, {c, o})
            ):
                continue
            neighbors.append(mapping[j])
        bonds.append("{" + " ".join(map(str, neighbors)) + "}")
    lines = [
        "set here [file dirname [file normalize [info script]]]",
        "set molid [mol new [file join $here view.pdb] type pdb waitfor all]",
    ]
    if len(frames) > 1:
        lines.append(
            "mol addfile [file join $here endpoints.xyz] type xyz waitfor all molid $molid"
        )
    lines += [
        "mol delrep 0 $molid",
        "set sel [atomselect $molid all]",
        "$sel setbonds {" + " ".join(bonds) + "}",
        "$sel delete",
        "color Element H white",
        "color Element N blue",
        "color Element O red",
        "mol representation Licorice 0.18 12 12",
        "mol color Element",
        "mol selection all",
        "mol addrep $molid",
    ]
    for first, second in ((c, o), (c, n)):
        lines += [
            "mol representation DynamicBonds 1.85 0.1 12",
            f'mol selection "index {mapping[first]} {mapping[second]}"',
            "mol addrep $molid",
        ]
    for h in sorted(hydrogens):
        lines += [
            "mol representation DynamicBonds 1.4 0.1 12",
            'mol selection "index '
            + " ".join(map(str, [mapping[h], *[mapping[a] for a in sorted(hetero)]]))
            + '"',
            "mol addrep $molid",
        ]
    (generation / "view_umbrella.vmd").write_text("\n".join(lines) + "\n")
    write_json(
        generation / "manifest.json",
        {
            "frames": len(frames),
            "source_stages": [str(p) for p in stages],
            "original_atom_ids_1based": [i + 1 for i in keep],
            "interpretation": "Productionblockendpoints only;notcontinuousMDmovie",
            "window": row["label"],
            "bias": row["bias"],
        },
    )
    (generation / "MANIFEST.sha256").write_text(
        "".join(
            f"{digest(p)}  {p.name}\n"
            for p in sorted(generation.iterdir())
            if p.is_file()
        )
    )
    pointer = output / f".latest_{os.getpid()}"
    pointer.symlink_to(generation.name)
    pointer.replace(output / "latest")
