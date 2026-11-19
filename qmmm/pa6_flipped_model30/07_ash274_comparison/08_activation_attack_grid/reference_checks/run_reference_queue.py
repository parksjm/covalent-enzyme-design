"""Run the bounded matched-embedding comparison queue with at most two workers."""

from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import logging
import os
from pathlib import Path
import shutil
import time
from typing import Any
from queue_consumer import accepted_points, all_rows_terminal, read_json, select_points
from reference_common import (
    HERE,
    GRID,
    BASELINE,
    BASELINE_JOB,
    MAX_WAIT_SECONDS,
    digest,
    write_json,
    journal,
)
from reference_embedding import (
    parse_external_input,
    expected_mm_ids,
    verify_baseline_pc_ids,
    validate_pc_membership,
)
from reference_gaussian import run_direct_gaussian
from reference_calculators import run_sander, run_xtb, capture_gaussian_export

LOGGER = logging.getLogger(__name__)


def external_input(folder: Path) -> Path:
    """Return the single successful Sander-exported Gaussian input."""
    path = folder / "old.gau_job.inp"
    if not path.is_file():
        raise FileNotFoundError(
            f"Sander did not retain expected external-charge input: {path}"
        )
    return path


def prepare_and_run() -> None:
    """Wait for baseline, run ES checks immediately, then queue up to three grid states."""
    started = time.monotonic()
    while True:
        if (
            BASELINE_JOB.is_file()
            and read_json(BASELINE_JOB).get("exit_status") is not None
        ):
            break
        if time.monotonic() - started > MAX_WAIT_SECONDS:
            raise TimeoutError("ES baseline did not finish within 48 hours")
        time.sleep(300)
    baseline_job = read_json(BASELINE_JOB)
    if baseline_job.get("exit_status") != 0:
        raise RuntimeError(
            "ES Sander export process failed; preserve its status record"
        )
    queue = read_json(HERE / "queue.json")
    baseline_inp = external_input(BASELINE)
    baseline_embedding = parse_external_input(baseline_inp)
    membership = read_json(GRID / "raw_mm_membership.json")
    ordered_ids = expected_mm_ids(membership)
    expected_count = len(ordered_ids)
    if len(baseline_embedding["point_charges"]) != expected_count:
        raise ValueError(
            "ES baseline PC count does not equal all non-QM atoms minus eight link parents"
        )
    if baseline_embedding["charge"] != -1 or len(baseline_embedding["atoms"]) != 110:
        raise ValueError(
            "ES baseline Gaussian model is not 102 real atoms + 8 links, charge −1"
        )
    verify_baseline_pc_ids(baseline_inp, membership, ordered_ids)
    baseline_charges = [charge[3] for charge in baseline_embedding["point_charges"]]
    pc_manifest = {
        "schema_version": 1,
        "source": str(baseline_inp),
        "source_sha256": digest(baseline_inp),
        "point_charge_count": expected_count,
        "charge_order_e": baseline_charges,
        "ordered_amber_atom_ids": ordered_ids,
        "atom_membership_source": "all non-QM topology atoms excluding the eight bonded QM/MM link parents; input-coordinate mapping verified at ES",
        "native_adjusted_charges": True,
        "qm_charge": baseline_embedding["charge"],
        "qm_link_atom_count": len(baseline_embedding["atoms"]),
    }
    write_json(HERE / "pc_manifest.json", pc_manifest)
    base_template = (BASELINE / "minimize.mdin").read_text(encoding="utf-8")
    environment = os.environ.copy() | {
        "AMBERHOME": "/data/jpc/code/ambertools26/venv",
        "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
        "GAUSS_EXEDIR": "/usr/local/g16",
        "g16root": "/usr/local",
        "PATH": "/usr/local/g16:" + os.environ.get("PATH", ""),
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    evaluations = HERE / "evaluations"
    evaluations.mkdir(exist_ok=False)
    records: list[dict[str, Any]] = []
    es_dir = evaluations / "01_es"
    es_dir.mkdir()
    dftb_es = run_sander(
        es_dir / "dftb3", BASELINE / "input.ncrst", "DFTB3", base_template, environment
    )
    gaussian_es = run_direct_gaussian(
        baseline_inp,
        es_dir / "gaussian_direct",
        BASELINE / "single_point.mdout",
        environment,
    )
    xtb_es_dir = es_dir / "gfn2_xtb"
    xtb_es_dir.mkdir()
    xtb_es = run_xtb(xtb_es_dir, baseline_inp, environment, baseline_charges)
    records.append(
        {
            "state": "ES",
            "energies_kcal_mol": {
                "dftb3_sander": dftb_es,
                "wb97xd_pc_self_corrected": gaussian_es[
                    "corrected_gaussian_energy_kcal_mol"
                ],
                "gfn2_xtb": xtb_es,
            },
            "gaussian": gaussian_es,
            "gradient_analysis": "raw FCHK/xTB artifacts retained; DFTB3-versus-DFT reaction projection pending",
        }
    )
    write_json(
        HERE / "results.json", {"status": "ES_complete_waiting_grid", "states": records}
    )
    while not all_rows_terminal():
        if time.monotonic() - started > MAX_WAIT_SECONDS:
            write_json(
                HERE / "selection.json",
                {"status": "timed_out", "selected_grid_states": []},
            )
            raise TimeoutError(
                "Grid rows did not reach terminal status within 48 hours"
            )
        LOGGER.info(
            "ES baseline references complete; waiting for five grid row statuses"
        )
        time.sleep(300)
    selected = select_points(accepted_points(queue))
    write_json(
        HERE / "selection.json",
        {
            "status": "selected",
            "selected_grid_states": selected,
            "maximum_selected_grid_states": 3,
            "grid_protocol_sha256": digest(GRID / "protocol.json"),
            "grid_points_sha256": digest(GRID / "grid_points.json"),
        },
    )

    def evaluate_selected(item: tuple[int, dict[str, Any]]) -> dict[str, Any]:
        index, point = item
        records: list[dict[str, Any]] = []
        state_dir = evaluations / f"{index + 1:02d}_{point['class']}"
        state_dir.mkdir()
        source_restart = Path(point["restart"])
        if digest(source_restart) != point["restart_sha256"]:
            raise ValueError(f"Selected restart changed: {source_restart}")
        acceptance = Path(point["acceptance_record"])
        shutil.copy2(acceptance, state_dir / "source_acceptance.json")
        shutil.copy2(source_restart, state_dir / "input.ncrst")
        dftb_energy = run_sander(
            state_dir / "dftb3",
            state_dir / "input.ncrst",
            "DFTB3",
            base_template,
            environment,
        )
        ext_inp = capture_gaussian_export(
            state_dir / "extern_export",
            state_dir / "input.ncrst",
            base_template,
            environment,
        )
        ext_embedding = parse_external_input(ext_inp)
        charges = validate_pc_membership(ext_inp, membership, baseline_charges)
        if len(charges) != len(ordered_ids):
            raise ValueError(
                "Selected-state point-charge membership differs from frozen Amber atom IDs"
            )
        if len(ext_embedding["atoms"]) != 110 or ext_embedding["charge"] != -1:
            raise ValueError(f"QM/link identity changed at {point['key']}")
        gaussian = run_direct_gaussian(
            ext_inp,
            state_dir / "gaussian_direct",
            state_dir / "extern_export/single_point.mdout",
            environment,
        )
        xtb_dir = state_dir / "gfn2_xtb"
        xtb_dir.mkdir()
        xtb_energy = run_xtb(xtb_dir, ext_inp, environment, baseline_charges)
        records.append(
            {
                "state": point["key"],
                "class": point["class"],
                "source_restart_sha256": point["restart_sha256"],
                "energies_kcal_mol": {
                    "dftb3_sander": dftb_energy,
                    "wb97xd_pc_self_corrected": gaussian[
                        "corrected_gaussian_energy_kcal_mol"
                    ],
                    "gfn2_xtb": xtb_energy,
                },
                "gaussian": gaussian,
                "point_charge_count": len(charges),
                "gradient_analysis": "raw FCHK/xTB artifacts retained; DFTB3-versus-DFT reaction projection pending",
            }
        )
        return records[0]

    with ThreadPoolExecutor(max_workers=2) as pool:
        for result in pool.map(evaluate_selected, enumerate(selected, 1)):
            records.append(result)
            write_json(HERE / "results.json", {"status": "running", "states": records})
    write_json(
        HERE / "results.json",
        {
            "status": "calculations_complete",
            "states": records,
            "energy_note": "Gaussian direct Force energy is SCF Done minus its reported point-charge self energy; do not use the field-query EXTERNESCF from the failed baseline.",
            "absolute_energy_comparison_to_periodic_native_qmmm": "prohibited",
            "reaction_coordinate_gradient_analysis": "pending",
        },
    )
    manifest_path = HERE / "pc_manifest.json"
    write_json(
        HERE / "queue_status.json",
        {
            "schema_version": 1,
            "status": "complete",
            "method_comparison": {
                "status": "inconclusive",
                "reviewed_by": None,
                "qualitative_agreement": None,
                "basis": "Raw energies and analytic-gradient artifacts retained; reaction-coordinate gradient analysis and lead review pending.",
                "raw_results": "reference_checks/results.json",
            },
            "required_state_count": {
                "expected": len(records),
                "completed": len(records),
            },
            "source_grid_protocol_sha256": digest(GRID / "protocol.json"),
            "fixed_pc_manifest_path": str(manifest_path.relative_to(GRID)),
            "fixed_pc_manifest_sha256": digest(manifest_path),
            "process_records": [
                str(path.relative_to(GRID))
                for path in sorted(evaluations.rglob("process.json"))
            ],
            "baseline_job_record": str(BASELINE_JOB.relative_to(GRID)),
            "baseline_job_record_sha256": digest(BASELINE_JOB),
            "queue_status_note": "Candidate-release queue must remain pending unless method_comparison.status is changed to pass after review.",
        },
    )


def main() -> None:
    """Execute the reviewed bounded reference queue."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    journal(
        "intent",
        "Start the bounded three-method reference queue; wait at most 48 hours for baseline and grid completion",
    )
    try:
        prepare_and_run()
    except Exception as error:
        write_json(HERE / "queue_status.json", {"schema_version": 1, "status": "failed", "reason": str(error)})
        journal("outcome", "Reference queue failed; outputs preserved", str(error))
        raise
    journal(
        "outcome",
        "Complete the three-method reference queue",
        "all prepared ES/grid state calculations completed",
    )


if __name__ == "__main__":
    main()
