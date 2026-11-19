"""Prepare or run exactly one endpoint-D finite-embedding reference pair."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import shutil
import sys
from typing import Any

import numpy as np
import parmed

HERE = Path(__file__).resolve().parent
GRID = HERE.parent.parent
REFERENCE = GRID / "reference_checks"
ROOT = HERE.parents[5]
sys.path.insert(0, str(REFERENCE))
sys.path.insert(0, str(REFERENCE / "gradient_analysis"))
sys.path.insert(0, str(ROOT / "src"))
from reference_common import BASELINE, TOPOLOGY, digest, journal, write_json  # noqa: E402
from reference_calculators import capture_gaussian_export, run_xtb  # noqa: E402
from reference_embedding import (  # noqa: E402
    expected_mm_ids,
    parse_external_input,
    validate_pc_membership,
    verify_baseline_pc_ids,
)
from reference_gaussian import run_direct_gaussian  # noqa: E402
from nyl12_covalent.proton_coordination import coordination_transfer  # noqa: E402

SOURCE = HERE.parent / "pilot/point_D"
RESTART = SOURCE / "accepted_00131_continuation_final.ncrst"
ENVIRONMENT = {
    "AMBERHOME": "/data/jpc/code/ambertools26/venv",
    "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
    "GAUSS_EXEDIR": "/usr/local/g16",
    "g16root": "/usr/local",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
}


def inspect_source() -> tuple[Any, dict[str, Any], dict[str, Any]]:
    """Validate the restart, fixed atom map and nearest-heavy proton inventory."""
    acceptance = json.loads((SOURCE / "scientific_acceptance.json").read_text())
    process = json.loads((SOURCE / "process.json").read_text())
    if (
        not acceptance["biased_force_converged"]
        or process["termination"] != "force_converged"
    ):
        raise ValueError("D pilot force convergence is not recorded")
    system = parmed.load_file(str(TOPOLOGY), xyz=str(RESTART))
    membership = json.loads((GRID / "raw_mm_membership.json").read_text())
    pc = json.loads((REFERENCE / "pc_manifest.json").read_text())
    ids = membership["qm_real_atom_ids"]
    if len(ids) != 102 or expected_mm_ids(membership) != pc["ordered_amber_atom_ids"]:
        raise ValueError("Fixed QM/MM atom membership changed")
    baseline = BASELINE / "old.gau_job.inp"
    validate_pc_membership(baseline, membership, pc["charge_order_e"])
    verify_baseline_pc_ids(baseline, membership, pc["ordered_amber_atom_ids"])
    x = system.coordinates[np.array(ids) - 1]
    z = np.array([system.atoms[i - 1].atomic_number for i in ids])
    h = np.flatnonzero(z == 1)
    heavy = np.flatnonzero(z != 1)
    if len(h) != 54:
        raise ValueError("Expected all 54 real QM hydrogens")
    owners = [ids[heavy[np.argmin(np.linalg.norm(x[heavy] - x[i], axis=1))]] for i in h]
    counts = {str(i): owners.count(i) for i in (13905, 13916, 20317, 14578, 14579)}
    if counts != {"13905": 3, "13916": 0, "20317": 2, "14578": 0, "14579": 1}:
        raise ValueError(f"Endpoint D proton inventory changed: {counts}")
    index = {atom_id: i for i, atom_id in enumerate(ids)}
    occupancy = coordination_transfer(x, index[13916], index[13905], h)
    attack = float(np.linalg.norm(x[index[13916]] - x[index[19920]]))
    if abs(occupancy.value - acceptance["endpoint"]["q"]) > 1e-10:
        raise ValueError("Restart q differs from accepted endpoint")
    if abs(attack - acceptance["endpoint"]["attack_A"]) > 1e-10:
        raise ValueError("Restart attack distance differs from accepted endpoint")
    from analyze import distance, dual_directions

    _, attack_j = distance(x, index[13916], index[19920])
    cut = [
        index[i]
        for i in ids
        if any(b.idx + 1 not in index for b in system.atoms[i - 1].bond_partners)
    ]
    contact_checks = {}
    for name, donor in (("Phe113_N_O", 1638), ("Asn186_ND2_O", 2732)):
        separation, contact_j = distance(x, index[donor], index[19921])
        jacobian = np.array([contact_j, occupancy.gradient, attack_j])
        direction = -dual_directions(jacobian)[0]
        response = jacobian.reshape(3, -1) @ direction.ravel()
        if np.any(direction[cut] != 0) or not np.allclose(
            response, [-1, 0, 0], atol=1e-10, rtol=0
        ):
            raise ValueError("Invalid fixed-q/attack oxyanion shortening direction")
        contact_checks[name] = {
            "distance_A": separation,
            "first_order_CV_response": response.tolist(),
            "boundary_motion": False,
            "Jacobian_rank": 3,
        }
    return (
        system,
        pc,
        {
            "qm_real_atom_ids": ids,
            "real_H": len(h),
            "nearest_heavy_H_counts": counts,
            "q": occupancy.value,
            "attack_A": attack,
            "point_charges": len(pc["charge_order_e"]),
            "oxyanion_shortening_preflight": contact_checks,
        },
    )


def prepare() -> None:
    """Write immutable hashes and preflight evidence without invoking any engine."""
    if (HERE / "preflight.json").exists():
        raise FileExistsError("Preflight already prepared; preserve this attempt")
    journal("intent", "Prepare endpoint D reference; no electronic-structure calls")
    _, _, mapping = inspect_source()
    files = [
        RESTART,
        SOURCE / "scientific_acceptance.json",
        SOURCE / "process.json",
        TOPOLOGY,
        GRID / "raw_mm_membership.json",
        REFERENCE / "pc_manifest.json",
        BASELINE / "minimize.mdin",
        BASELINE / "gau_job.tpl",
        BASELINE / "old.gau_job.inp",
        Path(__file__),
        HERE / "select_gradients.py",
        REFERENCE / "gradient_analysis/analyze.py",
        ROOT / "src/nyl12_covalent/proton_coordination.py",
    ]
    files.extend(
        REFERENCE / name
        for name in (
            "reference_common.py",
            "reference_calculators.py",
            "reference_embedding.py",
            "reference_gaussian.py",
        )
    )
    snapshot = HERE / "source_snapshot"
    snapshot.mkdir()
    for path in files:
        if path.suffix == ".py":
            shutil.copy2(path, snapshot / path.name)
    shutil.copy2(RESTART, HERE / "input.ncrst")
    hashes = {str(p): digest(p) for p in files}
    hashes[str(HERE / "input.ncrst")] = digest(HERE / "input.ncrst")
    write_json(
        HERE / "preflight.json",
        {
            "status": "prepared_no_engine_invoked",
            "sources_sha256": hashes,
            "mapping": mapping,
            "qm_charge": -1,
            "multiplicity": 1,
            "method": "wB97XD/def2SVP analytic Force; GFN2-xTB gradient",
            "resources": {
                "Gaussian_CPUs": 8,
                "Gaussian_memory": "32GB",
                "xTB_threads": 1,
                "sequential": True,
            },
            "environment_overrides": ENVIRONMENT,
            "command": [sys.executable, str(Path(__file__)), "--run"],
            "cwd": str(HERE),
            "scope": "Finite nonperiodic full-MM-PC electronic gradients; not complete DFT/MM gradients",
            "acceptance": "Normal exits;110 QM/link atoms;330 gradient values; exact restart/export mapping; H54 projections",
        },
    )
    journal("outcome", "Prepare endpoint D reference", "Engine-free preflight complete")


def run() -> None:
    """Export D and execute the reviewed two-method finite-embedding benchmark."""
    preflight = json.loads((HERE / "preflight.json").read_text())
    for name, expected in preflight["sources_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Preflight input/source changed: {name}")
    journal("intent", "Run bounded endpoint D matched finite-embedding pair")
    system, pc, mapping = inspect_source()
    environment = os.environ.copy() | ENVIRONMENT
    environment["PATH"] = "/usr/local/g16:" + os.environ["PATH"]
    exported = capture_gaussian_export(
        HERE / "extern_export",
        HERE / "input.ncrst",
        (BASELINE / "minimize.mdin").read_text(),
        environment,
    )
    embedding = parse_external_input(exported)
    membership = json.loads((GRID / "raw_mm_membership.json").read_text())
    validate_pc_membership(exported, membership, pc["charge_order_e"])
    qm_x = np.array([row[1:] for row in embedding["atoms"][:102]])
    pc_x = np.array([row[:3] for row in embedding["point_charges"]])
    qm_ids = np.array(mapping["qm_real_atom_ids"]) - 1
    pc_ids = np.array(pc["ordered_amber_atom_ids"]) - 1
    if not np.allclose(qm_x, system.coordinates[qm_ids], atol=1e-5, rtol=0):
        raise ValueError("D export real-QM coordinates/order changed")
    if not np.allclose(pc_x, system.coordinates[pc_ids], atol=1e-5, rtol=0):
        raise ValueError("D export MM-PC coordinates/order changed")
    gaussian = run_direct_gaussian(
        exported,
        HERE / "gaussian_direct",
        HERE / "extern_export/single_point.mdout",
        environment,
    )
    xtb_folder = HERE / "gfn2_xtb"
    xtb_folder.mkdir()
    journal("intent", "Run endpoint D GFN2-xTB gradient")
    xtb_energy = run_xtb(xtb_folder, exported, environment, pc["charge_order_e"])
    write_json(
        xtb_folder / "process.json",
        {
            "status": "finished",
            "energy_kcal_mol": xtb_energy,
            "outputs_sha256": {
                p.name: digest(p) for p in xtb_folder.iterdir() if p.is_file()
            },
        },
    )
    journal("outcome", "Run endpoint D GFN2-xTB gradient", "Normal termination")
    from select_gradients import select

    projections = select()
    write_json(
        HERE / "results.json",
        {
            "status": "calculations_complete",
            "source_restart_sha256": digest(RESTART),
            "gaussian": gaussian,
            "GFN2_energy_kcal_mol": xtb_energy,
            "projections": projections,
            "interpretation": "Finite PC QM electronic gradients only; no periodic DFT/MM or barrier claim",
        },
    )
    journal(
        "outcome", "Run bounded endpoint D matched finite-embedding pair", "Complete"
    )


def main() -> None:
    """Require an explicit run flag; otherwise prepare without engine execution."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    if args.run:
        run()
    else:
        prepare()


if __name__ == "__main__":
    main()
