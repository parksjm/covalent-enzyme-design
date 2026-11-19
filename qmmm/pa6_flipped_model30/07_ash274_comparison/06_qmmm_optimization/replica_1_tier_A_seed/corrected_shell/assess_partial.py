"""Classify preserved outputs after abnormal ASH274 run termination."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import parmed

RUN = Path(__file__).resolve().parent
CASE = RUN.parent
TOPOLOGY = CASE.parent.parent / "00_system/ash274_solvated.prmtop"


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_coordinate(structure: Any, residue_number: int, name: str) -> np.ndarray:
    """Return one exact residue atom coordinate."""
    matches = [
        atom
        for atom in structure.residues[residue_number - 1].atoms
        if atom.name.strip() == name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one atom {residue_number}:{name}, found {len(matches)}"
        )
    return np.asarray(structure.coordinates[matches[0].idx], dtype=float)


def main() -> None:
    """Write independent process, force, chemistry, tier, and water findings."""
    acceptance_path = RUN / "acceptance.json"
    if acceptance_path.exists():
        raise FileExistsError(acceptance_path)
    mdout = (RUN / "minimized.mdout").read_text(errors="replace")
    normal_markers = (
        ("FINAL RESULTS" in mdout or "Run   done at" in mdout)
        and "FATAL" not in mdout
        and "NaN" not in mdout
    )
    spec = importlib.util.spec_from_file_location(
        "corrected_ash274_runner", RUN / "prepare_and_run_ash274_qmmm.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load corrected runner for endpoint assessment")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.CASE = RUN
    scientific = runner.assess()
    initial = parmed.load_file(str(TOPOLOGY), xyz=str(RUN / "starting.ncrst"))
    endpoint = parmed.load_file(str(TOPOLOGY), xyz=str(RUN / "optimized.ncrst"))
    water = endpoint.residues[1575]
    names = {atom.name.strip() for atom in water.atoms}
    oxygen_name = "O" if "O" in names else "OW"
    hydrogen_names = [name for name in ("H1", "H2") if name in names]
    if len(hydrogen_names) != 2:
        raise ValueError(f"Unexpected water 1576 atom names: {sorted(names)}")
    initial_o = atom_coordinate(initial, 1576, oxygen_name)
    endpoint_o = atom_coordinate(endpoint, 1576, oxygen_name)
    og1 = atom_coordinate(endpoint, 929, "OG1")
    h1 = atom_coordinate(endpoint, 1576, hydrogen_names[0])
    h2 = atom_coordinate(endpoint, 1576, hydrogen_names[1])
    water_report = {
        "residue": 1576,
        "included_in_active_shell": True,
        "atom_names": sorted(names),
        "initial_O_to_Thr929_OG1_A": float(
            np.linalg.norm(initial_o - atom_coordinate(initial, 929, "OG1"))
        ),
        "saved_restart_O_to_Thr929_OG1_A": float(np.linalg.norm(endpoint_o - og1)),
        "saved_restart_O_displacement_from_start_A": float(
            np.linalg.norm(endpoint_o - initial_o)
        ),
        "saved_restart_OH1_A": float(np.linalg.norm(endpoint_o - h1)),
        "saved_restart_OH2_A": float(np.linalg.norm(endpoint_o - h2)),
        "saved_restart_H1_to_Thr929_OG1_A": float(np.linalg.norm(h1 - og1)),
        "saved_restart_H2_to_Thr929_OG1_A": float(np.linalg.norm(h2 - og1)),
        "interpretation": "Geometry is from the last written partial restart, not a converged endpoint.",
    }
    process = {
        "status": "terminated_abnormally_before_normal_completion",
        "exit_status": None,
        "exit_status_known": False,
        "normal_engine_markers": normal_markers,
        "force_convergence_marker": scientific["force_convergence_marker"],
        "last_saved_restart": "optimized.ncrst",
        "last_saved_restart_sha256": sha256(RUN / "optimized.ncrst"),
        "observed_utc": datetime.now(timezone.utc).isoformat(),
        "termination_cause": "Unknown; wrapper and Sander process disappeared before process.json finalization.",
    }
    report = {
        "process": process,
        "force_convergence": {
            "accepted": scientific["force_convergence_marker"],
            "criteria": "DL-FIND reports Convergence reached and Optimisation converged before maxcycle=2000",
            "mdout_marker_counts": {
                marker: mdout.count(marker)
                for marker in (
                    "Convergence reached",
                    "Optimisation converged",
                    "Run   done at",
                    "FINAL RESULTS",
                )
            },
            "mdout_tail": "\n".join(mdout.splitlines()[-12:]),
        },
        "chemistry": {
            "proton_retained_on_ash274": scientific["gates"][
                "proton_retained_on_ash274"
            ],
            "neutral_thr929_n_retained": scientific["gates"][
                "neutral_thr233_n_retained"
            ],
            "scissile_cn_intact": scientific["gates"]["scissile_cn_intact"],
            "proton_distances_A": scientific["proton_distances"],
        },
        "tier_geometry": {
            "tier": scientific["tier"],
            "tier_a_or_b_positive_phi": scientific["gates"]["tier_a_or_b_positive_phi"],
            "geometry": scientific["geometry"],
            "interpretation": "Provisional partial-restart geometry; force convergence and normal process completion failed.",
        },
        "water_1576_geometry": water_report,
        "overall_accepted": False,
        "calculation_count": 1,
    }
    acceptance_path.write_text(json.dumps(report, indent=2) + "\n")
    process_path = RUN / "process.json"
    process_record = json.loads(process_path.read_text())
    process_record.update(process)
    process_record["output_sha256"] = {
        name: sha256(RUN / name)
        for name in (
            "minimize.mdin",
            "starting.ncrst",
            "minimized.mdout",
            "optimized.ncrst",
            "optimization.nc",
            "engine.log",
        )
        if (RUN / name).is_file()
    }
    process_record["acceptance_sha256"] = sha256(acceptance_path)
    process_path.write_text(json.dumps(process_record, indent=2) + "\n")
    (RUN / "status.json").write_text(
        json.dumps(
            {
                "exit_status": None,
                "process_status": process["status"],
                "scientific_acceptance": report,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
