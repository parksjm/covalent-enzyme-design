"""Audit all Asp272/Asp274 oxygen pairs without changing protonation.

Amber numbering is verified against the four beta-subunit Thr233 starts.
Trajectory coordinates are streamed one frame at a time from Amber NetCDF.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

LOGGER = logging.getLogger(__name__)
THR_STARTS = (929, 1025, 1121, 1217)
PAIR_NAMES = ("OD1_OD1", "OD1_OD2", "OD2_OD1", "OD2_OD2")


def file_hash(path: Path) -> str:
    """Return a SHA-256 digest using bounded reads."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pair_distances(coordinates: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    """Return four oxygen distances in Angstrom for each site.

    Args:
        coordinates: Array of shape (sites, 4, 3), ordered as Asp272
            OD1/OD2 followed by Asp274 OD1/OD2.
        box: Optional lengths followed by alpha, beta and gamma angles.

    Returns:
        Array of shape (sites, 4), in PAIR_NAMES order.

    Raises:
        ValueError: Coordinates or box parameters are invalid.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 3 or coordinates.shape[1:] != (4, 3):
        raise ValueError("Expected oxygen coordinates with shape (sites, 4, 3)")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("Nonfinite coordinates")
    delta = coordinates[:, :2, None, :] - coordinates[:, None, 2:, :]
    if box is not None:
        box = np.asarray(box, dtype=float)
        if box.shape != (6,) or not np.all(np.isfinite(box)):
            raise ValueError("Expected three box lengths and three angles")
        if np.any(box[:3] <= 0) or np.any(box[3:] <= 0) or np.any(box[3:] >= 180):
            raise ValueError("Invalid periodic box lengths or angles")
        alpha, beta, gamma = np.deg2rad(box[3:])
        ca, cb, cg, sg = np.cos(alpha), np.cos(beta), np.cos(gamma), np.sin(gamma)
        cy = (ca - cb * cg) / sg
        cz2 = 1 - cb * cb - cy * cy
        if cz2 <= 0:
            raise ValueError("Degenerate periodic box")
        cell = np.array([[box[0], 0, 0], [box[1] * cg, box[1] * sg, 0],
                         [box[2] * cb, box[2] * cy, box[2] * np.sqrt(cz2)]])
        # Conventional Amber orthogonal/truncated-octahedral cells are reduced;
        # 27 neighbors correct fractional rounding for nonorthogonal cells.
        if np.linalg.cond(cell) > 3:
            raise ValueError("Box too skewed for the bounded nearest-image search")
        fractional = delta @ np.linalg.inv(cell)
        centered = fractional - np.rint(fractional)
        offsets = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
        candidates = (centered[..., None, :] + offsets) @ cell
        return np.linalg.norm(candidates, axis=-1).min(axis=-1).reshape(
            coordinates.shape[0], 4
        )
    return np.linalg.norm(delta, axis=-1).reshape(coordinates.shape[0], 4)


def atom_mapping(topology: Any) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Verify native residue mapping and return zero-based oxygen indices."""
    indices, records = [], []
    for site, start in enumerate(THR_STARTS, 1):
        if topology.residues[start - 1].name != "THR":
            raise ValueError(f"Site {site}: Amber residue {start} is not THR233")
        site_indices, residues = [], []
        for native in (272, 274):
            amber = start + native - 233
            residue = topology.residues[amber - 1]
            if residue.name != "ASP":
                raise ValueError(f"Native Asp{native}, Amber {amber}: {residue.name}")
            for name in ("OD1", "OD2"):
                atoms = [atom for atom in residue.atoms if atom.name == name]
                if len(atoms) != 1:
                    raise ValueError(f"Amber {amber}: expected unique {name}")
                site_indices.append(atoms[0].idx)
            residues.append({"native": native, "amber": amber})
        indices.append(site_indices)
        records.append({"site": site, "thr233_amber": start,
                        "aspartates": residues, "atom_indices_0based": site_indices})
    return np.asarray(indices), records


def named_paths(values: list[str]) -> list[tuple[str, Path]]:
    """Parse distinct LABEL=PATH input arguments."""
    result = []
    for value in values:
        label, separator, source = value.partition("=")
        if not separator or not label or Path(label).name != label:
            raise ValueError("Inputs must use a simple LABEL=PATH")
        path = Path(source).resolve(strict=True)
        result.append((label, path))
    if len({label for label, _ in result}) != len(result):
        raise ValueError("Duplicate input labels")
    return result


def write_summary(label: str, distances: np.ndarray) -> list[dict[str, Any]]:
    """Summarize each pair and the naming-invariant minimum per site."""
    records = []
    for site in range(distances.shape[1]):
        series = np.column_stack((distances[:, site], distances[:, site].min(axis=1)))
        for column, pair in enumerate((*PAIR_NAMES, "minimum")):
            values = series[:, column]
            records.append({"input": label, "site": site + 1, "pair": pair,
                            "frames": len(values), "mean_A": float(values.mean()),
                            "std_A": float(values.std()), "min_A": float(values.min()),
                            "p05_A": float(np.quantile(values, 0.05)),
                            "median_A": float(np.median(values)),
                            "p95_A": float(np.quantile(values, 0.95)),
                            "max_A": float(values.max()),
                            "fraction_lt_3A": float(np.mean(values < 3.0)),
                            "fraction_lt_3p5A": float(np.mean(values < 3.5))})
    return records


def analyze_trajectory(label: str, path: Path, indices: np.ndarray,
                       atom_count: int, output: Path
                       ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Stream every Amber trajectory frame and save distances and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.io import netcdf_file

    distances, times = [], []
    with netcdf_file(path, "r", mmap=True) as trajectory:
        variables = trajectory.variables
        shape = variables["coordinates"].shape
        if shape[1:] != (atom_count, 3):
            raise ValueError(f"Unexpected trajectory shape: {shape}")
        coordinate_units = getattr(variables["coordinates"], "units", b"")
        if coordinate_units not in (b"angstrom", "angstrom"):
            raise ValueError(f"Unrecognized coordinate units: {coordinate_units!r}")
        if "cell_lengths" not in variables or "cell_angles" not in variables:
            raise ValueError("Periodic trajectory lacks box lengths/angles")
        units = getattr(variables["time"], "units", b"")
        if units not in (b"picosecond", b"picoseconds", "picosecond", "picoseconds"):
            raise ValueError(f"Unrecognized trajectory time units: {units!r}")
        for frame in range(variables["coordinates"].shape[0]):
            oxygen = np.array(variables["coordinates"][frame][indices], dtype=float)
            box = np.concatenate((variables["cell_lengths"][frame],
                                  variables["cell_angles"][frame]))
            distances.append(pair_distances(oxygen, box))
            times.append(float(variables["time"][frame]))
        del variables
    values, time = np.asarray(distances), np.asarray(times)
    if not len(time) or not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0):
        raise ValueError("Trajectory time must be finite, nonempty and increasing")
    with (output / f"{label}_distances.csv").open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(("frame", "time_ps", "site", *PAIR_NAMES, "minimum_A"))
        for frame, timestamp in enumerate(time):
            for site in range(4):
                row = values[frame, site]
                writer.writerow((frame, timestamp, site + 1, *row, row.min()))
    figure, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    for site, axis in enumerate(axes):
        for pair, name in enumerate(PAIR_NAMES):
            axis.plot(time / 1000, values[:, site, pair], lw=0.5, alpha=0.45,
                      label=name)
        axis.plot(time / 1000, values[:, site].min(axis=1), color="black",
                  lw=0.8, label="minimum")
        axis.set_ylabel(f"Site {site + 1}\nO–O (Å)")
    axes[0].legend(ncol=5, fontsize=8)
    axes[0].set_title(label)
    axes[-1].set_xlabel("Stored trajectory time (ns)")
    figure.tight_layout()
    figure.savefig(output / f"{label}_distances.png", dpi=160)
    plt.close(figure)
    return write_summary(label, values), {"frames": len(time),
                                         "first_time_ps": float(time[0]),
                                         "last_time_ps": float(time[-1])}


def main() -> None:
    """Run the explicit topology/restart/trajectory distance audit."""
    import parmed

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--restart", action="append", default=[])
    parser.add_argument("--trajectory", action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    restarts, trajectories = named_paths(args.restart), named_paths(args.trajectory)
    if not restarts and not trajectories:
        raise ValueError("At least one input is required")
    if len({x[0] for x in restarts + trajectories}) != len(restarts + trajectories):
        raise ValueError("Input labels must be globally unique")
    args.output.mkdir(parents=True, exist_ok=False)
    topology = parmed.load_file(str(args.topology))
    indices, mapping = atom_mapping(topology)
    result: dict[str, Any] = {"mapping": mapping, "distance_units": "angstrom",
                              "periodicity": "triclinic minimum image, 27 neighbor search",
                              "thresholds": "descriptive, not acceptance gates",
                              "inputs": [], "summary": []}
    for label, path in restarts:
        structure = parmed.load_file(str(args.topology), xyz=str(path))
        if structure.box is None:
            raise ValueError(f"Periodic restart lacks box: {path}")
        values = pair_distances(np.asarray(structure.coordinates)[indices],
                                np.asarray(structure.box))
        result["summary"].extend(write_summary(label, values[None]))
        result["inputs"].append({"label": label, "path": str(path),
                                  "sha256": file_hash(path), "frames": 1})
    for label, path in trajectories:
        LOGGER.info("Streaming %s", path)
        records, metadata = analyze_trajectory(label, path, indices,
                                                len(topology.atoms), args.output)
        result["summary"].extend(records)
        result["inputs"].append({"label": label, "path": str(path),
                                  "sha256": file_hash(path), **metadata})
    result["topology"] = {"path": str(args.topology.resolve()),
                            "sha256": file_hash(args.topology)}
    result["script_sha256"] = file_hash(Path(__file__))
    (args.output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    with (args.output / "summary.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result["summary"][0]))
        writer.writeheader()
        writer.writerows(result["summary"])
    (args.output / "scientific_acceptance.json").write_text(json.dumps({
        "mapping_verified": True, "all_saved_frames_analyzed": True,
        "all_four_pairs_measured": True, "protonation_changed": False,
        "sampling_convergence_assessed": False,
        "limitation": "Fixed classical protonation; O-O distances alone do not establish proton sharing."
    }, indent=2) + "\n")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
