"""Continue only the two unfinished grid references, preserving completed work."""

from concurrent.futures import ThreadPoolExecutor
import json
import logging
import os
import shutil

from reference_common import BASELINE, GRID, HERE, digest, journal, write_json
from reference_calculators import capture_gaussian_export, run_xtb
from reference_embedding import parse_external_input, validate_pc_membership
from reference_gaussian import run_direct_gaussian


def evaluate(point: dict) -> dict:
    """Reuse the verified DFTB3 result and finish the embedded DFT/xTB checks."""
    folder = (
        HERE
        / "evaluations"
        / ("02_activated" if point["class"] == "activated" else "03_attack")
    )
    restart = folder / "input.ncrst"
    if digest(restart) != point["restart_sha256"]:
        raise ValueError(f"Restart changed: {restart}")
    dftb = json.loads((folder / "dftb3/process.json").read_text())
    if dftb["status"] != "finished" or dftb["exit_status"] != 0:
        raise ValueError(f"DFTB3 reference not complete: {folder}")
    for name, expected in dftb["outputs_sha256"].items():
        if digest(folder / "dftb3" / name) != expected:
            raise ValueError(f"DFTB3 output changed: {name}")
    failed = folder / "extern_export"
    if (folder / "extern_export_failed_pathlength").exists():
        raise FileExistsError("Continuation already attempted")
    shutil.move(failed, folder / "extern_export_failed_pathlength")
    environment = os.environ.copy() | {
        "AMBERHOME": "/data/jpc/code/ambertools26/venv",
        "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
        "GAUSS_EXEDIR": "/usr/local/g16",
        "g16root": "/usr/local",
        "PATH": "/usr/local/g16:" + os.environ["PATH"],
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    template = (BASELINE / "minimize.mdin").read_text()
    exported = capture_gaussian_export(
        folder / "extern_export", restart, template, environment
    )
    membership = json.loads((GRID / "raw_mm_membership.json").read_text())
    baseline_pc = json.loads((HERE / "pc_manifest.json").read_text())
    charges = validate_pc_membership(
        exported, membership, baseline_pc["charge_order_e"]
    )
    expected_ids = baseline_pc["ordered_amber_atom_ids"]
    # Match every exported charge position to this state's original Amber atom ID.
    import parmed
    from reference_common import TOPOLOGY
    import numpy as np

    system = parmed.load_file(str(TOPOLOGY), xyz=str(restart))
    coordinates = np.array(
        [p[:3] for p in parse_external_input(exported)["point_charges"]]
    )
    if not np.allclose(
        coordinates, system.coordinates[np.array(expected_ids) - 1], atol=1e-5, rtol=0
    ):
        raise ValueError("Selected-state external charge atom ordering changed")
    gaussian = run_direct_gaussian(
        exported,
        folder / "gaussian_direct",
        folder / "extern_export/single_point.mdout",
        environment,
    )
    xtb = folder / "gfn2_xtb"
    xtb.mkdir()
    xtb_energy = run_xtb(xtb, exported, environment, baseline_pc["charge_order_e"])
    write_json(
        xtb / "process.json",
        {
            "status": "finished",
            "energy_kcal_mol": xtb_energy,
            "output_sha256": digest(xtb / "gfn2.out"),
        },
    )
    return {
        "state": point["key"],
        "class": point["class"],
        "source_restart_sha256": point["restart_sha256"],
        "energies_kcal_mol": {
            "dftb3_sander": dftb["energy_kcal_mol"],
            "wb97xd_pc_self_corrected": gaussian["corrected_gaussian_energy_kcal_mol"],
            "gfn2_xtb": xtb_energy,
        },
        "gaussian": gaussian,
        "point_charge_count": len(charges),
        "gradient_analysis": "pending",
    }


def main() -> None:
    """Resume a bounded, explicitly failed selection; leave scientific review pending."""
    logging.basicConfig(level=logging.INFO)
    journal("intent", "Continue the two selected references after export-path repair")
    selection = json.loads((HERE / "selection.json").read_text())
    points = selection["selected_grid_states"]
    if len(points) != 2 or {p["class"] for p in points} != {"activated", "attack"}:
        raise ValueError("Continuation is pinned to the two inspected selected states")
    old = json.loads((HERE / "results.json").read_text())
    if len(old["states"]) != 1 or old["states"][0]["state"] != "ES":
        raise ValueError("Unexpected prior reference results")
    shutil.copy2(HERE / "queue_status.json", HERE / "queue_status_export_failure.json")
    write_json(
        HERE / "queue_status.json",
        {
            "schema_version": 1,
            "status": "running",
            "reason": "continuing selected references; ES retained",
        },
    )
    records = old["states"]
    try:
        with ThreadPoolExecutor(max_workers=2) as workers:
            for result in workers.map(evaluate, points):
                records.append(result)
                write_json(
                    HERE / "results.json", {"status": "running", "states": records}
                )
    except Exception as error:
        write_json(
            HERE / "queue_status.json",
            {"schema_version": 1, "status": "failed", "reason": str(error)},
        )
        journal("outcome", "Reference continuation failed", str(error))
        raise
    write_json(
        HERE / "results.json",
        {
            "status": "calculations_complete",
            "states": records,
            "relative_energies_kcal_mol": [
                {
                    "state": r["state"],
                    **{
                        m: e - records[0]["energies_kcal_mol"][m]
                        for m, e in r["energies_kcal_mol"].items()
                    },
                }
                for r in records
            ],
            "interpretation": "Matched finite embedding; periodic PES and gradient projection remain separate",
        },
    )
    write_json(
        HERE / "queue_status.json",
        {
            "schema_version": 1,
            "status": "complete",
            "method_comparison": {
                "status": "inconclusive",
                "reviewed_by": None,
                "qualitative_agreement": None,
                "raw_results": "reference_checks/results.json",
            },
            "required_state_count": {"expected": 3, "completed": 3},
            "source_grid_protocol_sha256": digest(GRID / "protocol.json"),
            "fixed_pc_manifest_path": "reference_checks/pc_manifest.json",
            "fixed_pc_manifest_sha256": digest(HERE / "pc_manifest.json"),
            "process_records": [
                str(p.relative_to(GRID))
                for p in (HERE / "evaluations").rglob("process.json")
            ],
        },
    )
    journal(
        "outcome",
        "Selected-state references complete",
        "Three total states; scientific comparison pending",
    )


if __name__ == "__main__":
    main()
