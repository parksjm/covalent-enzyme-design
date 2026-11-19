"""Validate and project saved endpoint D gradients on all real-QM hydrogens."""

from __future__ import annotations

from typing import Any

import numpy as np

from run_reference import HERE, ROOT, TOPOLOGY, inspect_source
from analyze import BOHR, HARTREE, deck, distance, dual_directions, engrad, fchk_array
from reference_common import digest, write_json
from nyl12_covalent.proton_coordination import coordination_transfer


def select() -> dict[str, Any]:
    """Return matched real-QM gradients and occupancy/attack projections.

    Directions hold boundary QM atoms and external charges fixed, so no moving
    link-parent chain rule is required. These are electronic QM derivatives.
    """
    system, _, mapping = inspect_source()
    ids = mapping["qm_real_atom_ids"]
    index = {v: i for i, v in enumerate(ids)}
    source = HERE / "extern_export/old.gau_job.inp"
    fchk = HERE / "gaussian_direct/direct_force.fchk"
    xtb = HERE / "gfn2_xtb/gfn2.qm_link.engrad"
    _, x = deck(source)
    dft_x = fchk_array(fchk, "Current cartesian coordinates").reshape(-1, 3) * BOHR
    dft_z = fchk_array(fchk, "Atomic numbers").astype(int)
    xtb_z, xtb_x, xtb_g = engrad(xtb)
    dft_g = fchk_array(fchk, "Cartesian Gradient").reshape(-1, 3)
    expected_z = np.array([system.atoms[i - 1].atomic_number for i in ids] + [1] * 8)
    if not np.array_equal(dft_z, expected_z) or not np.array_equal(xtb_z, expected_z):
        raise ValueError("QM/link gradient atomic order differs")
    if dft_g.shape != (110, 3) or xtb_g.shape != (110, 3):
        raise ValueError("Expected 110-atom analytic gradients")
    if not np.isfinite(dft_g).all() or not np.isfinite(xtb_g).all():
        raise ValueError("Nonfinite gradient")
    errors = {
        "DFT_export_max_A": float(abs(dft_x - x).max()),
        "xTB_export_max_A": float(abs(xtb_x - x).max()),
        "restart_export_max_A": float(
            abs(system.coordinates[np.array(ids) - 1] - x[:102]).max()
        ),
    }
    if max(errors.values()) > 1e-6:
        raise ValueError(f"Geometry mismatch: {errors}")
    h = np.flatnonzero(expected_z[:102] == 1)
    if len(h) != 54:
        raise ValueError("Expected all54 real hydrogen atoms")
    real_x = x[:102]
    n, o, tagged_h, c = [index[i] for i in (13905, 13916, 13917, 19920)]
    occupancy = coordination_transfer(real_x, o, n, h)
    attack, attack_j = distance(real_x, o, c)
    oh, oh_j = distance(real_x, o, tagged_h)
    nh, nh_j = distance(real_x, n, tagged_h)
    cut_qm = [
        i
        for i in ids
        if any(b.idx + 1 not in index for b in system.atoms[i - 1].bond_partners)
    ]
    arrays = {
        "real_QM_amber_ids": np.array(ids),
        "coordinates_A": x,
        "DFT_gradient_Eh_bohr": dft_g,
        "GFN2_gradient_Eh_bohr": xtb_g,
    }
    rows = []
    for cv, cv_j in (("q_occ", occupancy.gradient), ("r_OH", oh_j), ("r_NH", nh_j)):
        directions = dual_directions(np.array([cv_j, attack_j]))
        if np.any(directions[:, [index[i] for i in cut_qm]] != 0):
            raise ValueError("Moving cut parent requires explicit link chain rule")
        arrays[f"{cv}_directions_A_per_CV"] = directions
        for axis, name in enumerate((cv, f"attack_hold_{cv}")):
            for method, gradient in (
                ("DFT_fixed_PC_QM", dft_g),
                ("GFN2_fixed_PC_QM", xtb_g),
            ):
                rows.append(
                    {
                        "method": method,
                        "direction": name,
                        "derivative": float(
                            np.sum(gradient[:102] * HARTREE / BOHR * directions[axis])
                        ),
                        "units": "kcal/mol per dimensionless q"
                        if name == "q_occ"
                        else "kcal/mol/A",
                    }
                )
    # Three-CV duals preserve occupancy and attack to first order.
    for label, donor in (("Phe113_N_O", 1638), ("Asn186_ND2_O", 2732)):
        separation, contact_j = distance(real_x, index[donor], index[19921])
        jacobian = np.array([contact_j, occupancy.gradient, attack_j])
        direction = -dual_directions(jacobian)[0]
        if np.any(direction[[index[i] for i in cut_qm]] != 0):
            raise ValueError("Oxyanion direction moves a boundary parent")
        response = jacobian.reshape(3, -1) @ direction.ravel()
        if not np.allclose(response, [-1.0, 0.0, 0.0], atol=1e-10, rtol=0):
            raise ValueError("Oxyanion contact dual constraint failure")
        arrays[f"{label}_shorten_directions_A_per_A"] = direction
        for method, gradient in (
            ("DFT_fixed_PC_QM", dft_g),
            ("GFN2_fixed_PC_QM", xtb_g),
        ):
            rows.append(
                {
                    "method": method,
                    "direction": f"shorten_{label}_hold_q_attack",
                    "distance_A": separation,
                    "derivative": float(
                        np.sum(gradient[:102] * HARTREE / BOHR * direction)
                    ),
                    "units": "kcal/mol/A of shortening",
                    "first_order_CV_response": response.tolist(),
                }
            )
    np.savez_compressed(HERE / "gradient_arrays.npz", **arrays)
    result = {
        "verification": errors,
        "real_H": 54,
        "cut_QM_amber_ids": cut_qm,
        "q_occ": occupancy.value,
        "attack_A": attack,
        "r_OH_A": oh,
        "r_NH_A": nh,
        "projections": rows,
        "scope": "Finite fixed-PC electronic QM derivatives",
        "sources_sha256": {
            str(p): digest(p)
            for p in (
                source,
                fchk,
                xtb,
                TOPOLOGY,
                ROOT / "src/nyl12_covalent/proton_coordination.py",
            )
        },
    }
    write_json(HERE / "gradient_projections.json", result)
    return result
