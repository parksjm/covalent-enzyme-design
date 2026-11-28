"""Full-precision, uncensored umbrella CVs and proton-label-independent chemistry."""

from __future__ import annotations

import csv
import hashlib
from itertools import product
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from parmed import unit
from parmed.amber import AmberFormat
from parmed.geometry import box_lengths_and_angles_to_vectors
from scipy.io import netcdf_file

_SHIFTS = np.asarray(list(product((-1, 0, 1), repeat=3)), dtype=float)


def _distances(
    xyz: np.ndarray, first: Sequence[int], second: Sequence[int], box: np.ndarray
) -> np.ndarray:
    cell = np.asarray(
        box_lengths_and_angles_to_vectors(*box).value_in_unit(unit.angstroms)
    )
    if not np.isfinite(cell).all() or np.linalg.det(cell) <= 0:
        raise ValueError("Invalid periodic cell")
    # Neighbor search around fractional rounding handles Amber's reduced
    # orthorhombic and truncated-octahedron cells, including wrapped molecules.
    delta = xyz[np.asarray(first) - 1, None, :] - xyz[np.asarray(second) - 1, :]
    frac = delta @ np.linalg.inv(cell)
    base = frac - np.rint(frac)
    vectors = (base[..., None, :] - _SHIFTS) @ cell
    return np.sqrt(np.min(np.sum(vectors * vectors, axis=-1), axis=-1))


def frame_observables(
    xyz: np.ndarray,
    box: np.ndarray,
    roles: Mapping[str, int],
    hydrogen_pool: Sequence[int],
    owner_atom_ids: Sequence[int],
    initial_owners: Mapping[int, int],
) -> dict:
    """Measure one frame; preserve chemistry changes as flags, never exclusions.

    Args:
        xyz: Full system coordinates in Angstrom, shape (atoms,3).
        box: Amber lengths/angles, Angstrom/degrees, shape (6,).
        roles: Old mapped carbon, leaving_n, og1, hg1, nalpha, nalpha_h1/h2 IDs.
        hydrogen_pool: One-based real QM H initially bonded to N/O, no link H.
        owner_atom_ids: One-based QM N/O IDs and initial N/O parents.
        initial_owners: Topology-based initial N/O parent of each pool hydrogen.

    Returns:
        Raw distances, labeled CVs, coordination CVs and uncensored owner flags.
    """
    xyz = np.asarray(xyz, dtype=float)
    box = np.asarray(box, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError("Nonfinite or invalid coordinates")
    if box.shape != (6,) or not np.isfinite(box).all() or not hydrogen_pool:
        raise ValueError("Need finite periodic box and a nonempty proton pool")
    c, n, o, h, na = [roles[k] for k in ("carbon", "leaving_n", "og1", "hg1", "nalpha")]
    if min([*roles.values(), *hydrogen_pool, *owner_atom_ids]) < 1 or max(
        [*roles.values(), *hydrogen_pool, *owner_atom_ids]
    ) > len(xyz):
        raise ValueError("Atom role/pool ID outside coordinate array")

    def d(a: int, b: int) -> float:
        return float(_distances(xyz, [a], [b], box)[0, 0])

    cn, co, oh, nh = d(c, n), d(c, o), d(o, h), d(n, h)
    pool_dist = _distances(xyz, [n, o, na], hydrogen_pool, box)
    coordination = (1 / (1 + (pool_dist / 1.25) ** 6)).sum(axis=1)
    distances = _distances(xyz, hydrogen_pool, owner_atom_ids, box)
    nearest = distances.argmin(axis=1)
    owners = {
        int(atom): {
            "owner": int(owner_atom_ids[int(j)]),
            "distance_A": float(distances[i, j]),
        }
        for i, (atom, j) in enumerate(zip(hydrogen_pool, nearest, strict=True))
    }
    changed = [
        int(hid) for hid in hydrogen_pool if owners[hid]["owner"] != initial_owners[hid]
    ]
    owner_distances = {
        int(hid): {
            int(owner): float(distances[i, j]) for j, owner in enumerate(owner_atom_ids)
        }
        for i, hid in enumerate(hydrogen_pool)
    }
    ambiguous = [
        int(hid)
        for hid in hydrogen_pool
        if sum(distance <= 1.30 for distance in owner_distances[int(hid)].values()) > 1
    ]
    return {
        "r_c_n": cn,
        "r_c_og1": co,
        "r_og1_hg1": oh,
        "r_nleave_hg1": nh,
        "q_acyl": cn - co,
        "q_h": oh - nh,
        "cn_nleave_hpool": float(coordination[0]),
        "cn_og1_hpool": float(coordination[1]),
        "cn_progress": float(coordination[0] - coordination[1]),
        "cn_nalpha_hpool": float(coordination[2]),
        "r_nalpha_h1": d(na, roles["nalpha_h1"]),
        "r_nalpha_h2": d(na, roles["nalpha_h2"]),
        "hg1_owner": owners[h]["owner"],
        "nalpha_h1_owner": owners[roles["nalpha_h1"]]["owner"],
        "nalpha_h2_owner": owners[roles["nalpha_h2"]]["owner"],
        "owner_changed": bool(changed),
        "changed_hydrogen_ids": changed,
        "hydrogen_owners": owners,
        "hydrogen_owner_distances_A": owner_distances,
        "ambiguous_hydrogens_within_1_30_A": ambiguous,
    }


def topology_proton_pool(
    topology: Path, qm_atom_ids: Sequence[int]
) -> tuple[list[int], dict[int, int], list[int]]:
    """Return real QM N/O-bound hydrogens and ownership atoms from topology.

    Args:
        topology: Amber topology containing the full atom inventory.
        qm_atom_ids: One-based real QM atom IDs, excluding link atoms.

    Returns:
        Sorted hydrogen pool, initial H-to-parent map, and sorted N/O owners.
    """
    data = AmberFormat(str(topology)).parm_data
    natom = int(data["POINTERS"][0])
    z = data["ATOMIC_NUMBER"]
    qm = set(int(atom) for atom in qm_atom_ids)
    if not qm or min(qm) < 1 or max(qm) > natom:
        raise ValueError("QM IDs outside topology")
    owners: dict[int, int] = {}
    for flag in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
        bonds = data[flag]
        for i in range(0, len(bonds), 3):
            a, b = int(bonds[i]) // 3 + 1, int(bonds[i + 1]) // 3 + 1
            for hydrogen, parent in ((a, b), (b, a)):
                if hydrogen in qm and z[hydrogen - 1] == 1 and z[parent - 1] in (7, 8):
                    if hydrogen in owners:
                        raise ValueError("Hydrogen has multiple initial N/O parents")
                    owners[hydrogen] = parent
    pool = sorted(owners)
    owner_ids = sorted(
        {atom for atom in qm if z[atom - 1] in (7, 8)} | set(owners.values())
    )
    if not pool:
        raise ValueError("QM topology contains no N/O-bound real hydrogen")
    return pool, owners, owner_ids


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_observables(
    topology: Path,
    trajectory: Path,
    roles: Mapping[str, int],
    qm_atom_ids: Sequence[int],
    output_dir: Path,
) -> dict:
    """Stream every full-precision frame into CSV plus a source-pinned receipt.

    Uses memory-mapped Amber NetCDF and only one full coordinate frame at once.
    Numeric invalidity is an error; valid chemical branches remain in the CSV.
    """
    output_dir.mkdir(parents=True, exist_ok=False)
    pool, owners, hetero = topology_proton_pool(topology, qm_atom_ids)
    qm = set(qm_atom_ids)
    natom = int(AmberFormat(str(topology)).parm_data["POINTERS"][0])
    tracked = [roles[k] for k in ("hg1", "nalpha_h1", "nalpha_h2")]
    if not set(tracked) <= set(pool):
        raise ValueError("Tracked protons must belong to real QM N/O proton pool")
    records = 0
    flagged = 0
    csv_path = output_dir / "observables.csv"
    with (
        netcdf_file(str(trajectory), "r", mmap=True) as nc,
        csv_path.open("w") as output,
    ):
        coords = nc.variables["coordinates"]
        if coords.shape[1:] != (natom, 3):
            raise ValueError("Trajectory/topology atom count mismatch")
        writer = None
        for index in range(coords.shape[0]):
            box = np.r_[
                nc.variables["cell_lengths"][index], nc.variables["cell_angles"][index]
            ]
            row = {
                "frame": index + 1,
                "time_ps": float(nc.variables["time"][index]),
                **frame_observables(coords[index], box, roles, pool, hetero, owners),
            }
            flagged += row["owner_changed"]
            row["changed_hydrogen_ids"] = json.dumps(row["changed_hydrogen_ids"])
            row["hydrogen_owners"] = json.dumps(row["hydrogen_owners"], sort_keys=True)
            if writer is None:
                writer = csv.DictWriter(output, fieldnames=list(row))
                writer.writeheader()
            writer.writerow(row)
            records += 1
        del coords
    if records == 0:
        raise ValueError("Empty trajectory")
    receipt = {
        "frames": records,
        "owner_changed_frames": flagged,
        "frames_censored": 0,
        "source_sha256": {
            str(p): _hash(p) for p in (topology, trajectory, Path(__file__))
        },
        "csv": str(csv_path),
        "csv_sha256": _hash(csv_path),
        "roles_1based": dict(roles),
        "qm_atom_ids": sorted(qm),
        "hydrogen_pool_1based": pool,
        "initial_owners_1based": owners,
        "coordination_r0_A": 1.25,
        "coordination_exponent": 6,
        "coordination_function": "1/(1+(r/r0)^6)",
        "policy": "All valid frames retained; proton ownership changes flag only; no additional bias",
        "precision": "Native NetCDF precision, never three-decimal DUMPAVE",
        "pbc": "Minimum image over 27 neighboring reduced Amber cell images",
    }
    (output_dir / "observables.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt
