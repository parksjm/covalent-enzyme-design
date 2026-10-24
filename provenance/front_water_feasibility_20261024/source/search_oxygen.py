"""Deterministic fixed-pocket front-side relay oxygen feasibility search."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
SOURCE = (
    ROOT
    / "output_nyl12_local_repair_20261021/pa6_l4_alpha_beta/attempt_05_tip3p_preconditioned/calculation"
)
SEED = ROOT / "output_nyl12_oxyanion_cluster_20261023/seed_01"
USER = ROOT / "output_nyl12_user_water_seed_20261023/native_seed_jmp.xyz"
# Bondi vdW radii, angstrom. Severe steric overlap means >0.40 A overlap.
RADII = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 16: 1.80}
LOGGER = logging.getLogger(__name__)


def front_frame(n: np.ndarray, og: np.ndarray, preferred: np.ndarray) -> tuple:
    """Define an atom-based front hemisphere from the preferred oxygen.

    Args:
        n: Terminal nitrogen coordinate.
        og: Thr OG1 coordinate.
        preferred: User-supplied oxygen coordinate in the same frame.

    Returns:
        Axis N-to-OG, front unit vector perpendicular to it, and right vector.
    """
    axis = og - n
    axis /= np.linalg.norm(axis)
    front = preferred - n - np.dot(preferred - n, axis) * axis
    if np.linalg.norm(front) < 1e-8:
        raise ValueError("Preferred oxygen does not define a front direction")
    front /= np.linalg.norm(front)
    return axis, front, np.cross(axis, front)


def oxygen_grid(n: np.ndarray, og: np.ndarray, preferred: np.ndarray) -> tuple:
    """Generate the explicit 0.025 A distance/0.5 degree front-side grid.

    Args:
        n: Terminal nitrogen coordinate.
        og: Thr oxygen coordinate.
        preferred: User oxygen defining front.

    Returns:
        Coordinates and normalized front projections, ordered deterministically.
    """
    axis, front, right = front_frame(n, og, preferred)
    distance = np.linalg.norm(og - n)
    rn, ro, phi = np.meshgrid(
        np.linspace(2.6, 3.2, 25),
        np.linspace(2.6, 3.2, 25),
        np.deg2rad(np.arange(-89.75, 90.0, 0.5)),
        indexing="ij",
    )
    x = (rn**2 - ro**2 + distance**2) / (2 * distance)
    rho2 = rn**2 - x**2
    valid = rho2 > 0
    phi, x, rho = phi[valid], x[valid], np.sqrt(rho2[valid])
    xyz = (
        n
        + x[:, None] * axis
        + rho[:, None] * (np.cos(phi)[:, None] * front + np.sin(phi)[:, None] * right)
    )
    return xyz, np.cos(phi)


def load_pocket() -> tuple:
    """Load immutable source geometry and extra constructed cluster caps.

    Returns:
        ParmEd model, precise seed coordinates, provenance and full environment.
    """
    parm = AmberParm(str(SOURCE / "model.prmtop"), str(SOURCE / "prepared.rst7"))
    lines = (SEED / "native_seed.xyz").read_text().splitlines()[2:]
    seed = np.array([[float(v) for v in line.split()[1:]] for line in lines])
    provenance = json.loads((SEED / "coordinate_provenance.json").read_text())
    for item in provenance:
        if not item["constructed_cap_H"] and not item["origin_id"].startswith("WAT:"):
            source_idx = item["source_atom_index_1based"] - 1
            if not np.allclose(
                seed[item["index"]], parm.coordinates[source_idx], atol=1e-6, rtol=0
            ):
                raise ValueError(f"Source frame mismatch: {item['origin_id']}")
    nh_ids = {
        item["source_atom_index_1based"] - 1
        for item in provenance
        if item["origin_id"] in {"THR:H1", "THR:H2", "THR:HG1"}
    }
    keep = [a.idx for a in parm.atoms if a.residue.idx != 329 and a.idx not in nh_ids]
    coords = [parm.coordinates[i] for i in keep]
    radii = [RADII[parm.atoms[i].atomic_number] for i in keep]
    names = [
        f"{parm.atoms[i].residue.name}:{parm.atoms[i].residue.idx + 1}:{parm.atoms[i].name}"
        for i in keep
    ]
    donor_parent = []
    for i in keep:
        a = parm.atoms[i]
        parent = [b for b in a.bond_partners if b.atomic_number in {7, 8}]
        donor_parent.append(
            parm.coordinates[parent[0].idx]
            if a.atomic_number == 1 and parent
            else [np.nan] * 3
        )
    for item in provenance:
        if item["constructed_cap_H"]:
            coords.append(seed[item["index"]])
            radii.append(RADII[1])
            names.append(item["origin_id"])
            donor_parent.append([np.nan] * 3)
    return (
        parm,
        seed,
        provenance,
        (np.array(coords), np.array(radii), names, np.array(donor_parent)),
    )


def screen_oxygen(
    points: np.ndarray,
    environment: tuple,
    n: np.ndarray,
    max_n_distance: float = 3.2,
) -> tuple:
    """Calculate worst steric clearance at every oxygen, allowing real H bonds.

    Args:
        points: Oxygen grid coordinates.
        environment: Coordinates, radii, names, and N/O-bound hydrogen parents.
        n: Nitrogen center bounding every oxygen within 3.2 A.
        max_n_distance: Explicit largest candidate N-O distance in angstrom.

    Returns:
        Minimum clearance, worst collider index and retained environment indices.
    """
    coords, radii, names, parents = environment
    # Triangle inequality: omitted atoms cannot overlap O at any candidate.
    if np.max(np.linalg.norm(points - n, axis=1)) > max_n_distance + 1e-8:
        raise ValueError("Candidate lies outside declared screening bound")
    local = np.flatnonzero(
        np.linalg.norm(coords - n, axis=1)
        <= max_n_distance + 1.52 + max(RADII.values())
    )
    coords, radii, parents = coords[local], radii[local], parents[local]
    minima, indices = [], []
    threshold = 1.52 + radii - 0.40
    intended = np.array(
        [
            names[i] in {"T1P:233:N", "T1P:233:OG1"}
            or (names[i].split(":")[1:] in [["233", "N"], ["233", "OG1"]])
            for i in local
        ]
    )
    threshold[intended] = 2.6
    for chunk in np.array_split(points, max(1, int(np.ceil(len(points) / 2048)))):
        distances = np.linalg.norm(chunk[:, None, :] - coords[None, :, :], axis=2)
        margins = distances - threshold
        # Only genuine donor-H...O geometry can permit short O-H approach.
        h_to_parent = parents - coords
        h_to_o = chunk[:, None, :] - coords
        cosine = np.sum(h_to_parent * h_to_o, axis=2) / (
            np.linalg.norm(h_to_parent, axis=1)[None, :] * distances
        )
        donor_distance = np.linalg.norm(chunk[:, None, :] - parents[None, :, :], axis=2)
        hbond = (
            (cosine <= np.cos(np.deg2rad(140)))
            & (donor_distance >= 2.6)
            & (donor_distance <= 3.3)
        )
        margins[hbond] = distances[hbond] - 1.50
        indices.append(local[np.argmin(margins, axis=1)])
        minima.append(np.min(margins, axis=1))
    return np.concatenate(minima), np.concatenate(indices), local


def main() -> None:
    """Search fixed-pocket front oxygen positions and retain ranked evidence."""
    logging.basicConfig(level=logging.INFO)
    parm, seed, provenance, environment = load_pocket()
    preferred = np.array(
        [float(v) for v in USER.read_text().splitlines()[39].split()[1:]]
    )
    expected = np.array([-3.833305, 8.999226, -1.874442])
    if not np.allclose(preferred, expected, atol=1e-6):
        raise ValueError("User oxygen index/coordinate mismatch")
    n, og = seed[0], seed[3]
    points, front = oxygen_grid(n, og, preferred)
    margins, worst, local = screen_oxygen(points, environment, n)
    passed = margins >= -1e-8
    rank = np.lexsort((np.linalg.norm(points - preferred, axis=1), -margins))
    best = rank[:20]
    np.savez(
        OUT / "oxygen_grid.npz",
        coordinates=points,
        front=front,
        margin=margins,
        collider=worst,
    )
    report = {
        "scientific_status": "oxygen_positions_exist"
        if passed.any()
        else "no_oxygen_grid_position_clears_fixed_pocket",
        "not_a_continuous_impossibility_proof": True,
        "grid_count": len(points),
        "oxygen_pass_count": int(passed.sum()),
        "distance_step_A": 0.025,
        "azimuth_step_degrees": 0.5,
        "front_definition": "positive dot of O projection perpendicular to N-OG with user O projection; normalized radial projection cos(phi)",
        "front_frame": [v.tolist() for v in front_frame(n, og, preferred)],
        "screen": {
            "radii_A": RADII,
            "allowed_vdw_overlap_A": 0.4,
            "intended_N_OG_minimum_A": 2.6,
            "genuine_N/O-H_to_waterO_minimum_A": 1.5,
            "hbond_angle_minimum_degrees": 140,
            "donor_O_distance_A": [2.6, 3.3],
        },
        "environment_count": len(environment[0]),
        "local_count": len(local),
        "local_selection_proof": "N-distance <=3.2+1.52+max(radius)=6.52A; omitted atoms cannot clash by triangle inequality",
        "unchanged_nonwater_heavy_coordinates": True,
        "front_sensitivity": {
            str(cut): int(np.sum(passed & (front >= cut)))
            for cut in [0.0, 0.25, 0.5, 0.75]
        },
        "best_clearance_positions": [
            {
                "oxygen_A": points[i].tolist(),
                "margin_A": float(margins[i]),
                "worst_collider": environment[2][worst[i]],
                "front_projection": float(front[i]),
                "distance_from_user_O_A": float(np.linalg.norm(points[i] - preferred)),
            }
            for i in best
        ],
        "source_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                SOURCE / "model.prmtop",
                SOURCE / "prepared.rst7",
                SEED / "native_seed.xyz",
                SEED / "native_seed.sdf",
                SEED / "atom_map.json",
                SEED / "coordinate_provenance.json",
                USER,
                Path(__file__),
            ]
        },
    }
    (OUT / "oxygen_feasibility.json").write_text(json.dumps(report, indent=2) + "\n")
    LOGGER.info(
        "Oxygen points %d; sterically clear %d; best margin %.4f A",
        len(points),
        passed.sum(),
        margins[rank[0]],
    )


if __name__ == "__main__":
    main()
