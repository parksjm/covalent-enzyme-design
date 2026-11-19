"""Prepare trace-derived ASP272/ASH274 native NEB inputs; never launch engines."""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from nyl12_covalent.native_sander_cli import source_sha256  # noqa: E402
from nyl12_covalent.relay_case import write_json  # noqa: E402
from nyl12_covalent.relay_neb_native import native_initial_band, native_neb_input  # noqa: E402
from nyl12_covalent.relay_neb_trace import load_trace, proton_audit, resample_trace  # noqa: E402
from scripts.prepare_dimer_continuation import write_restart  # noqa: E402

LOGGER = logging.getLogger(__name__)
NATIVE = Path(
    "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src"
)


def prepare(case_path: Path, output: Path) -> dict[str, Any]:
    """Prepare13 controls and audit the installed24regular-image initializer.

    Args:
        case_path: Completed latest relay campaign fixed case JSON.
        output: New preparation directory, conventionally campaign/04_neb.

    Returns:
        Finite preflight evidence; config.json defines an unlaunched native job.

    Raises:
        ValueError: If source, atom/proton, exterior or predicted band gates fail.
        FileExistsError: If output already exists; nothing is overwritten.
    """
    case_path, output = case_path.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    trace = load_trace(case_path)
    model = trace.model
    active_ids = np.flatnonzero(model.active)
    local = {int(atom): i for i, atom in enumerate(active_ids)}
    qm_positions = np.array([local[int(atom)] for atom in model.qm])
    controls, brackets, arcs = resample_trace(
        trace.active_coordinates, qm_positions, 13
    )
    if not np.array_equal(
        controls[0], trace.active_coordinates[0]
    ) or not np.array_equal(controls[-1], trace.active_coordinates[-1]):
        raise ValueError(
            "Control endpoints differ from exact accepted source endpoints"
        )
    band, prediction = native_initial_band(controls, qm_positions)

    def audit_frames(frames: np.ndarray) -> list[dict[str, Any]]:
        """Audit real H labels and reactive geometry in reconstructed full images."""
        rows = []
        ids = {key: value - 1 for key, value in model.atom_map.items()}
        for i, active_xyz in enumerate(frames):
            xyz = model.coordinates.copy()
            xyz[model.active] = active_xyz
            row = {"image": i, **proton_audit(model, xyz)}
            row["distances_A"] = {
                label: float(np.linalg.norm(xyz[ids[a]] - xyz[ids[b]]))
                for label, a, b in (
                    ("OG1_C", "thr_og", "carbon"),
                    ("C_Nleaving", "carbon", "leaving_n"),
                    ("C_O", "carbon", "oxygen"),
                    ("water_HG1", "water_o", "proton"),
                    ("OG1_H20318", "thr_og", "water_h1"),
                    ("Nleave_H20318", "leaving_n", "water_h1"),
                )
            }
            rows.append(row)
        return rows

    boundary = [
        (bond.atom1.idx, bond.atom2.idx)
        for bond in model.topology.bonds
        if model.active[bond.atom1.idx] != model.active[bond.atom2.idx]
    ]
    boundary_changes = []
    for endpoint in (0, -1):
        original, aligned = model.coordinates.copy(), model.coordinates.copy()
        original[model.active], aligned[model.active] = (
            controls[endpoint],
            band[endpoint],
        )
        boundary_changes.extend(
            abs(
                float(np.linalg.norm(aligned[a] - aligned[b]))
                - float(np.linalg.norm(original[a] - original[b]))
            )
            for a, b in boundary
        )
    boundary_max = max(boundary_changes, default=0.0)
    prediction["alignment_boundary_bond_max_change_A"] = boundary_max
    prediction["alignment_boundary_warning"] = bool(
        boundary_max > 0.05 or prediction["max_active_alignment_shift_A"] > 0.1
    )
    control_audits = audit_frames(controls)
    predicted_audits = audit_frames(band)
    midpoint_audits = audit_frames((band[:-1] + band[1:]) / 2)
    output.mkdir(parents=True, exist_ok=False)
    images = output / "inputs/images"
    images.mkdir(parents=True)
    paths = []
    for i, active_xyz in enumerate(controls):
        xyz = model.coordinates.copy()
        xyz[model.active] = active_xyz
        path = images / f"control_{i:02d}.ncrst"
        write_restart(path, xyz, model.restart)
        paths.append(path)
    predicted = output / "predicted_native_initialization"
    predicted.mkdir()
    for i, active_xyz in enumerate(band):
        xyz = model.coordinates.copy()
        xyz[model.active] = active_xyz
        write_restart(predicted / f"regular_{i + 1:02d}.ncrst", xyz, model.restart)
    calculation = output / "calculation"
    calculation.mkdir()
    mdin = native_neb_input(model.config["mdin"], paths)
    (calculation / "minimize.mdin").write_text(mdin)
    source_copies = output / "source_snapshot"
    source_copies.mkdir()
    copies, seen = {}, {}
    files = [
        case_path,
        Path(model.config["topology"]),
        Path(model.config["calibration"]),
        *trace.paths,
        Path(__file__),
        ROOT / "src/nyl12_covalent/relay_neb_trace.py",
        ROOT / "src/nyl12_covalent/relay_neb_native.py",
        ROOT / "src/nyl12_covalent/relay_case.py",
        NATIVE / "sander/dlfind_module.F90",
        NATIVE / "libdlfind/src/dlf_neb.F90",
        NATIVE / "libdlfind/src/dlf_coords.F90",
    ]
    hashes = dict(trace.sha256)
    for i, path in enumerate(files):
        digest = source_sha256(path)
        hashes[str(path)] = digest
        if digest not in seen:
            copied = source_copies / f"source_{i:04d}_{path.parent.name}_{path.name}"
            shutil.copy2(path, copied)
            if source_sha256(copied) != digest:
                raise ValueError("Retained NEB source snapshot mismatch")
            seen[digest] = str(copied)
        copies[str(path)] = seen[digest]
    reactive_steps = {
        key: np.diff([row["distances_A"][key] for row in predicted_audits]).tolist()
        for key in predicted_audits[0]["distances_A"]
    }
    preflight = dict(
        status="prepared_not_launched; actual native startup verification required",
        direction="unbiased_ES_return_to_unbiased_acylenzyme",
        source_frame_count=len(trace.paths),
        exact_ES_source=str(trace.paths[0]),
        exact_acyl_source=str(trace.paths[-1]),
        control_endpoints_exact_source_coordinates=True,
        source_paths=[str(p) for p in trace.paths],
        source_sha256=hashes,
        source_copies=copies,
        source_audits=trace.audits,
        source_arc=arcs,
        control_brackets=brackets,
        control_audits=control_audits,
        predicted_regular_audits=predicted_audits,
        predicted_midpoint_audits=midpoint_audits,
        predicted_native_initialization=prediction,
        adjacent_reactive_distance_changes_A=reactive_steps,
        transition_window025_checkpoint_count=sum(
            p.parent.name == "window_025" for p in trace.paths
        ),
        atom_order_and_all54H_labels_preserved=True,
        frozen_exterior_exact=True,
        qm_atoms=102,
        real_qm_H=54,
        qm_charge=-1,
        qm_link_H=8,
        active_atoms=4662,
        electronic_temperature_K=0,
        seed_controls=13,
        regular_images=24,
        neb_nimage=25,
        maxcycle=16000,
        native_energy_budget=112000,
        amber_maxcyc=100000,
        maxstep_bohr=0.05,
        endpoint_policy="FROZEN",
        climbing_image_enabled=False,
        wtmask_role="QM102 weights affect image distances/dimer orientation; all4662 active QM/MM forces retained",
        external_reaction_or_positional_restraints=False,
        endpoint_hessian_verified=False,
        native_band_force_converged=False,
        barrier_or_PMF_claim=False,
    )
    write_json(output / "preflight.json", preflight)
    active_masks = re.findall(r"\bactive\s*=\s*'([^']+)'", model.config["mdin"], re.I)
    if len(active_masks) != 1:
        raise ValueError(
            "Expected one exact current-case active mask for monitor contract"
        )
    config = dict(
        case_path=str(case_path),
        reactant=str(paths[0]),
        active_mask=active_masks[0],
        active_atoms=4662,
        regular_images=24,
        topology=model.config["topology"],
        sander=str(Path(model.config["amber_home"]) / "bin/sander"),
        start_restart=str(paths[0]),
        base_input=str(calculation / "minimize.mdin"),
        calculation=str(calculation),
        source_case=str(case_path),
        images=[str(p) for p in paths],
        expected_regular_images=24,
        native_neb_nimage=25,
        method="DFTB3/3OB 0K periodic QM/MM",
        qm_atom_ids=model.config["qm_atom_ids"],
        active_atom_ids=model.config["active_atom_ids"],
        atom_map=model.atom_map,
        qm_charge=-1,
        endpoint_policy="FROZEN",
        trajectory_authority="nebpath.xyz; numbered native restarts unusable",
        startup_verification_required=True,
        source_sha256=hashes,
    )
    write_json(output / "config.json", config)
    LOGGER.info(
        "Prepared13 controls for24 regular images; source=%d, max active alignment shift %.6g A",
        len(trace.paths),
        prediction["max_active_alignment_shift_A"],
    )
    return preflight


def main() -> None:
    """Prepare a new audited band without launching a native engine."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    prepare(args.case, args.output)


if __name__ == "__main__":
    main()
