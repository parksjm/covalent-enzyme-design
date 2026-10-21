"""Bounded rigid-geometry repair search for the TI1 relay water."""

from __future__ import annotations

from typing import Any

import numpy as np

from .relay_feasibility import _validate_inventory
from .relay_repair_geometry import (
    axis_local_environment_indices,
    contact_limits_for_indices,
    donor_geometry,
    moved_atom_contact_screen,
    oxygen_grid,
    rotate_points,
    water_geometry,
)

TOLERANCE_A = 0.10


def _water_screen(
    topology: Any,
    water_xyz: np.ndarray,
    coordinates: np.ndarray,
    environment_indices: list[int],
) -> dict[str, Any]:
    """Apply existing relay-water oxygen/hydrogen limits to an environment."""
    oxygen_limits, hydrogen_limits = contact_limits_for_indices(
        topology, environment_indices
    )
    environment = coordinates[np.asarray(environment_indices, dtype=int)]
    oxygen_distances = np.linalg.norm(environment - water_xyz[0], axis=1)
    hydrogen_distances = np.linalg.norm(
        water_xyz[1:, None, :] - environment[None, :, :], axis=2
    )
    shortfalls = np.concatenate(
        (
            np.maximum(oxygen_limits - oxygen_distances, 0.0),
            np.maximum(hydrogen_limits[None, :] - hydrogen_distances, 0.0).ravel(),
        )
    )
    maximum = float(shortfalls.max(initial=0.0))
    return {
        "accepted": maximum <= TOLERANCE_A,
        "maximum_shortfall_A": maximum,
        "violation_count": int(np.count_nonzero(shortfalls > TOLERANCE_A)),
    }


def _site_screen(
    point: np.ndarray, coordinates: np.ndarray, indices: list[int], limits: np.ndarray
) -> bool:
    if not indices:
        return True
    distances = np.linalg.norm(coordinates[np.asarray(indices)] - point, axis=1)
    return bool(np.max(np.maximum(limits - distances, 0.0), initial=0.0) <= TOLERANCE_A)


def _rmsd(coordinates: np.ndarray, reference: np.ndarray, indices: list[int]) -> float:
    if not indices:
        return 0.0
    delta = coordinates[np.asarray(indices)] - reference[np.asarray(indices)]
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _inventory(topology: Any, coordinates: np.ndarray, target: str) -> dict[str, Any]:
    inventory, _, water_indices = _validate_inventory(topology, coordinates)
    asn = {atom.name: atom for atom in topology.residues[185].atoms}
    thr = {atom.name: atom for atom in topology.residues[928].atoms}
    ligand = {atom.name: atom for atom in topology.residues[1312].atoms}
    if target not in ligand:
        raise ValueError(
            f"Target oxyanion atom {target!r} is absent from ligand residue 1313"
        )
    if not {"CB", "CG", "CA", "ND2"}.issubset(asn):
        raise ValueError("ASN186 requires CA, CB, CG, and ND2")
    if asn["CG"] not in asn["CB"].bond_partners:
        raise ValueError("ASN186 CG must be bonded to CB")
    if thr["CA"] not in thr["N"].bond_partners:
        raise ValueError("Thr929 CA must be bonded to N")
    if set(water_indices) != {a.idx for a in topology.residues[1313].atoms}:
        raise ValueError("Water inventory is inconsistent")
    inventory.update(
        {
            "target_legacy_atom_name": target,
            "asn186_rotatable_axes": ["CA-CB", "CB-CG"],
            "thr929_nh3_axis": "N-CA",
        }
    )
    return inventory


def _pair_screen_between(
    topology: Any,
    moved: list[int],
    other: list[int],
    coordinates: np.ndarray,
) -> dict[str, Any]:
    """Screen cross-group contacts using the explicit atom-pair thresholds."""
    violations = 0
    maximum = 0.0
    from .relay_repair_geometry import pair_contact_limit

    for i in moved:
        for j in other:
            a, b = topology.atoms[i], topology.atoms[j]
            distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            shortfall = max(0.0, pair_contact_limit(a, b) - distance)
            maximum = max(maximum, shortfall)
            violations += shortfall > TOLERANCE_A
    return {
        "accepted": maximum <= TOLERANCE_A,
        "maximum_shortfall_A": maximum,
        "violation_count": int(violations),
    }


def search_repair(
    topology: Any,
    coords: np.ndarray,
    target: str,
    step_degrees: int = 15,
) -> tuple[dict[str, Any], np.ndarray | None]:
    """Search bounded ASN186/NH3+ rotations and rigid relay-water geometries.

    Args:
        topology: Hydrogenated TI1 topology with the validated 1314-residue
            inventory and bond connectivity.
        coords: Finite coordinates in topology order, in Angstrom.
        target: Archival substrate identifier pa6_l2 or poly_xama_maxa.
        step_degrees: Angular sampling increment dividing 360 degrees.

    Returns:
        A JSON-serializable audit report and accepted copied coordinates, or
        ``None`` coordinates when every candidate fails an explicit gate.

    Raises:
        ValueError: For invalid inventory, connectivity, target, or angular grid.
    """
    coordinates = np.asarray(coords, dtype=float)
    if step_degrees < 1 or 360 % step_degrees:
        raise ValueError("step_degrees must be a positive divisor of 360")
    legacy_target = target
    target = {"pa6_l2": "O49", "poly_xama_maxa": "O034"}[target]
    inventory = _inventory(topology, coordinates, target)
    asn = {a.name: a for a in topology.residues[185].atoms}
    thr = {a.name: a for a in topology.residues[928].atoms}
    wat = {a.name: a for a in topology.residues[1313].atoms}
    asn_side = sorted(set(_component(asn["CB"], asn["CA"])))
    chi2_side = sorted(set(_component(asn["CG"], asn["CB"])))
    beta_h = [thr[name].idx for name in ("H1", "H2", "H3")]
    water_ids = [wat[name].idx for name in ("O", "H1", "H2")]
    all_nonwater = [a.idx for a in topology.atoms if a.idx not in water_ids]
    fixed = [i for i in all_nonwater if i not in set(asn_side + beta_h)]
    n_idx, ca_idx, og_idx = thr["N"].idx, thr["CA"].idx, thr["OG1"].idx
    oxyanion_idx = next(
        a.idx for a in topology.residues[1312].atoms if a.name == target
    )
    step_values = list(range(0, 360, step_degrees))
    sites = oxygen_grid(coordinates[n_idx], coordinates[og_idx])
    fixed_o_limits, _ = contact_limits_for_indices(topology, fixed)
    fixed_pass = [
        site
        for site in sites
        if _site_screen(site[0], coordinates, fixed, fixed_o_limits)
    ]
    asn_o_limits, _ = contact_limits_for_indices(topology, asn_side)
    beta_o_limits, _ = contact_limits_for_indices(topology, beta_h)

    # Baselines are retained even when they fail; candidate gates are never relaxed.
    chi1_origin = coordinates[asn["CA"].idx]
    chi1_axis = coordinates[asn["CB"].idx] - chi1_origin
    side_radius = max(
        float(np.linalg.norm(coordinates[i] - coordinates[asn["CB"].idx]))
        for i in asn_side
    )
    side_local = axis_local_environment_indices(
        coordinates, all_nonwater, chi1_origin, chi1_axis, side_radius
    )
    side_baseline = moved_atom_contact_screen(
        topology, coordinates, asn_side, side_local
    )
    n_axis = coordinates[n_idx] - coordinates[ca_idx]
    beta_local = axis_local_environment_indices(
        coordinates,
        all_nonwater,
        coordinates[n_idx],
        n_axis,
        max(float(np.linalg.norm(coordinates[i] - coordinates[n_idx])) for i in beta_h),
    )
    beta_baseline = moved_atom_contact_screen(topology, coordinates, beta_h, beta_local)
    original_water_xyz = coordinates[np.asarray(water_ids)]
    original_water_screen = _water_screen(
        topology, original_water_xyz, coordinates, all_nonwater
    )
    original_donors = []
    for atom_name in ("H1", "H2", "H3"):
        hidx = thr[atom_name].idx
        distance, angle = donor_geometry(
            coordinates[n_idx], coordinates[hidx], coordinates[water_ids[0]]
        )
        original_donors.append(
            {"hydrogen_name": atom_name, "H_to_O_A": distance, "N_H_O_angle_deg": angle}
        )
    original_water_acceptors = []
    for atom_name in ("H1", "H2"):
        hidx = wat[atom_name].idx
        distance, angle = donor_geometry(
            coordinates[water_ids[0]], coordinates[hidx], coordinates[og_idx]
        )
        original_water_acceptors.append(
            {
                "hydrogen_name": atom_name,
                "H_to_OG1_A": distance,
                "O_H_OG1_angle_deg": angle,
            }
        )

    candidate_rows: list[
        tuple[tuple[float, float, float], np.ndarray, dict[str, Any]]
    ] = []
    rotation_audit: list[dict[str, Any]] = []
    failure_counts = {
        key: 0
        for key in (
            "asn_contact",
            "asn_oxyanion_distance",
            "oxygen_site",
            "donor_geometry",
            "water_orientation",
            "moved_contact",
            "water_full_environment",
            "accepted",
        )
    }
    accepted_count = 0
    best_failure = "no oxygen site passed optimistic fixed-environment screen"
    # Precompute fixed-environment water orientations for each optimistic O site.
    orientation_cache: dict[tuple[int, int], tuple[np.ndarray, dict[str, Any]]] = {}
    for site_id, site in enumerate(fixed_pass):
        for roll in range(0, 360, 10):
            water_xyz = water_geometry(site[0], coordinates[og_idx], roll)
            screen = _water_screen(topology, water_xyz, coordinates, fixed)
            orientation_cache[(site_id, roll)] = (water_xyz, screen)

    nh3_variants = []
    for hdelta in step_values:
        nh3 = rotate_points(coordinates, beta_h, coordinates[n_idx], n_axis, hdelta)
        beta_contact = moved_atom_contact_screen(
            topology, nh3, beta_h, [i for i in beta_local if i not in asn_side]
        )
        nh3_variants.append((hdelta, nh3, beta_contact))

    for chi1_delta in step_values:
        chi1_coords = rotate_points(
            coordinates,
            asn_side,
            coordinates[asn["CA"].idx],
            coordinates[asn["CB"].idx] - coordinates[asn["CA"].idx],
            chi1_delta,
        )
        for chi2_delta in step_values:
            rotamer = rotate_points(
                chi1_coords,
                chi2_side,
                chi1_coords[asn["CB"].idx],
                chi1_coords[asn["CG"].idx] - chi1_coords[asn["CB"].idx],
                chi2_delta,
            )
            asn_record = {
                "chi1_delta_deg": chi1_delta,
                "chi2_delta_deg": chi2_delta,
                "asn_contact_failures": 0,
                "oxygen_site_failures": 0,
                "nh3_rotamers_tested": 0,
                "donor_geometry_failures": 0,
                "water_orientation_failures": 0,
                "moved_contact_failures": 0,
                "water_full_environment_failures": 0,
                "accepted_candidates": 0,
            }
            side_local_rot = axis_local_environment_indices(
                rotamer,
                all_nonwater,
                rotamer[asn["CA"].idx],
                rotamer[asn["CB"].idx] - rotamer[asn["CA"].idx],
                side_radius,
            )
            d3 = float(np.linalg.norm(rotamer[asn["ND2"].idx] - rotamer[oxyanion_idx]))
            asn_record["asn_ND2_oxyanion_distance_A"] = d3
            if not 2.4 <= d3 <= 3.5:
                failure_counts["asn_oxyanion_distance"] += 1
                rotation_audit.append(asn_record)
                continue
            asn_contact = moved_atom_contact_screen(
                topology,
                rotamer,
                asn_side,
                [i for i in side_local_rot if i not in beta_h],
            )
            if not asn_contact["accepted"]:
                asn_record["asn_contact_failures"] = 1
                failure_counts["asn_contact"] += 1
                rotation_audit.append(asn_record)
                continue
            compatible_sites = [
                (site_id, site)
                for site_id, site in enumerate(fixed_pass)
                if _site_screen(site[0], rotamer, asn_side, asn_o_limits)
            ]
            asn_record["oxygen_site_failures"] = len(fixed_pass) - len(compatible_sites)
            failure_counts["oxygen_site"] += len(fixed_pass) - len(compatible_sites)
            for hdelta, nh3, beta_contact in nh3_variants:
                asn_record["nh3_rotamers_tested"] += 1
                if not beta_contact["accepted"]:
                    asn_record["moved_contact_failures"] += 1
                    failure_counts["moved_contact"] += 1
                    continue
                combined = rotamer.copy()
                combined[beta_h] = nh3[beta_h]
                cross_contact = _pair_screen_between(
                    topology, asn_side, beta_h, combined
                )
                if not cross_contact["accepted"]:
                    asn_record["moved_contact_failures"] += 1
                    failure_counts["moved_contact"] += 1
                    continue
                for site_id, site in compatible_sites:
                    if not _site_screen(site[0], nh3, beta_h, beta_o_limits):
                        asn_record["oxygen_site_failures"] += 1
                        failure_counts["oxygen_site"] += 1
                        continue
                    donor_options = []
                    for hname in ("H1", "H2", "H3"):
                        hidx = thr[hname].idx
                        d, angle = donor_geometry(nh3[n_idx], nh3[hidx], site[0])
                        if 1.4 <= d <= 2.5 and angle >= 140.0:
                            donor_options.append((hname, hidx, d, angle))
                    if not donor_options:
                        asn_record["donor_geometry_failures"] += 1
                        failure_counts["donor_geometry"] += 1
                        continue
                    for roll in range(0, 360, 10):
                        water_xyz, fixed_water_gate = orientation_cache[(site_id, roll)]
                        if not fixed_water_gate["accepted"]:
                            asn_record["water_orientation_failures"] += 1
                            failure_counts["water_orientation"] += 1
                            continue
                        # Test the three water atoms against both moved groups.
                        moved_water_gate = _water_screen(
                            topology, water_xyz, rotamer, asn_side
                        )
                        beta_water_gate = _water_screen(
                            topology, water_xyz, nh3, beta_h
                        )
                        if (
                            not moved_water_gate["accepted"]
                            or not beta_water_gate["accepted"]
                        ):
                            asn_record["moved_contact_failures"] += 1
                            failure_counts["moved_contact"] += 1
                            continue
                        candidate = coordinates.copy()
                        candidate[asn_side] = rotamer[asn_side]
                        candidate[beta_h] = nh3[beta_h]
                        candidate[water_ids] = water_xyz
                        # Bounded group and cross-contact screens cover all pairs.
                        # Independently verify the selected candidate globally.
                        moved_gate = {
                            "accepted": True,
                            "screen": "bounded groups and cross contacts",
                        }
                        water_gate = {
                            "accepted": True,
                            "screen": "fixed and moved environment partition",
                        }
                        accepted_count += 1
                        asn_record["accepted_candidates"] += 1
                        failure_counts["accepted"] += 1
                        heavy_asn = [
                            i for i in asn_side if topology.atoms[i].atomic_number > 1
                        ]
                        asn_rmsd = _rmsd(candidate, coordinates, heavy_asn)
                        beta_rmsd = _rmsd(candidate, coordinates, beta_h)
                        penalty = abs(site[1] - 2.85) + abs(site[2] - 2.85)
                        rank = (asn_rmsd, beta_rmsd, 0.05 * penalty)
                        selected_donor = min(donor_options, key=lambda row: row[2])
                        _, oh_og_angle = donor_geometry(
                            water_xyz[0], water_xyz[1], coordinates[og_idx]
                        )
                        details = {
                            "accepted": True,
                            "status": "geometrically_accepted_seed",
                            "target_legacy_atom_name": target,
                            "asn_ND2_oxyanion_distance_A": d3,
                            "chi1_delta_deg": chi1_delta,
                            "chi2_delta_deg": chi2_delta,
                            "nh3_rotation_delta_deg": hdelta,
                            "water_h2_roll_deg": roll,
                            "donor_hydrogen_name": selected_donor[0],
                            "nh3_donor_H_to_O_A": selected_donor[2],
                            "N_H_O_angle_deg": selected_donor[3],
                            "water_H1_to_OG1_A": float(
                                np.linalg.norm(water_xyz[1] - coordinates[og_idx])
                            ),
                            "O_H1_OG1_angle_deg": oh_og_angle,
                            "N_to_water_O_A": site[1],
                            "OG1_to_water_O_A": site[2],
                            "oxyanion_to_water_O_A": float(
                                np.linalg.norm(water_xyz[0] - coordinates[oxyanion_idx])
                            ),
                            "asn186_heavy_RMSD_A": asn_rmsd,
                            "thr929_NH3_hydrogen_RMSD_A": beta_rmsd,
                            "water_distance_penalty_A": penalty,
                            "moved_contact_gate": moved_gate,
                            "water_full_environment_gate": water_gate,
                            "changed_atom_indices_zero_based": sorted(
                                i
                                for i in asn_side + beta_h + water_ids
                                if np.linalg.norm(candidate[i] - coordinates[i])
                                > 1.0e-10
                            ),
                        }
                        if not candidate_rows or rank < candidate_rows[0][0]:
                            candidate_rows[:] = [(rank, candidate, details)]
            rotation_audit.append(asn_record)
    report: dict[str, Any] = {
        "method": "bounded rigid-geometry relay repair; finite-grid rejection screens",
        "accepted": bool(candidate_rows),
        "status": "geometrically_accepted_seed"
        if candidate_rows
        else "rejected_no_accepted_candidate",
        "inventory": inventory,
        "archival_target_alias": legacy_target,
        "target_legacy_atom_name": target,
        "grid": {
            "oxygen_grid": "9x9x72 exact sphere-intersection grid",
            "n_and_og1_distance_range_A": [2.6, 3.4],
            "oxygen_grid_site_count": len(sites),
            "optimistic_fixed_screen_site_count": len(fixed_pass),
            "optimistic_fixed_screen_sites": [
                {
                    "oxygen_xyz_A": s[0].tolist(),
                    "n_to_oxygen_A": s[1],
                    "og1_to_oxygen_A": s[2],
                    "distance_grid_indices": [s[3], s[4]],
                    "azimuth_index": s[5],
                }
                for s in fixed_pass
            ],
            "nh3_rotation_values_deg": step_values,
            "asn_chi1_chi2_values_deg": step_values,
            "water_h2_roll_values_deg": list(range(0, 360, 10)),
            "tolerance_A": TOLERANCE_A,
        },
        "original_geometry": {
            "water_contact_gate": original_water_screen,
            "thr929_N_H_O_to_original_water": original_donors,
            "original_water_O_H_OG1": original_water_acceptors,
            "asn186_moved_group_contact_gate": side_baseline,
            "thr929_NH3_moved_group_contact_gate": beta_baseline,
        },
        "stage_failure_counts": failure_counts,
        "per_chi_rotation_audit": rotation_audit,
        "accepted_candidate_count": accepted_count,
        "selected_geometry": None,
        "best_failure_reason": None if candidate_rows else best_failure,
        "changed_atom_indices_zero_based": [],
        "ranking": "lexicographic ASN186 sidechain heavy RMSD, Thr929 NH3 hydrogen RMSD, then 0.05 times summed deviation from 2.85 A N/OG1-water distances; geometric preference only, not energy",
        "limitations": (
            "Finite angular and distance grids screen geometry only. This does not calculate energy, "
            "validate force-field parameters, prove stereochemical optimization, or infer a proton-transfer pathway."
        ),
    }
    if not candidate_rows:
        # Distinguish empty optimistic oxygen feasibility from later rejection stages.
        if fixed_pass:
            best_failure = "all ASN/NH3/water orientations failed explicit geometry or contact screens"
            report["best_failure_reason"] = best_failure
        return report, None
    candidate_rows.sort(key=lambda row: row[0])
    _, selected_coords, selected_details = candidate_rows[0]
    final_moved = moved_atom_contact_screen(
        topology, selected_coords, asn_side + beta_h, all_nonwater
    )
    final_water = _water_screen(
        topology, selected_coords[water_ids], selected_coords, all_nonwater
    )
    if not final_moved["accepted"] or not final_water["accepted"]:
        raise ValueError(
            "Conservative screen disagrees with full-environment verification"
        )
    selected_details["moved_contact_gate"] = final_moved
    selected_details["water_full_environment_gate"] = final_water
    report["selected_geometry"] = selected_details
    report["changed_atom_indices_zero_based"] = selected_details[
        "changed_atom_indices_zero_based"
    ]
    report["accepted_candidate_count"] = accepted_count
    return report, selected_coords


def _component(atom: Any, cut_partner: Any) -> list[int]:
    """Get the component on the selected side of an ASN bond cut."""
    from .relay_repair_geometry import component_after_cut

    return component_after_cut(atom, cut_partner)
