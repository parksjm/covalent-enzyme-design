"""Run one prepared grid row and record stage and point provenance."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any

import numpy as np
import parmed

from prepare_grid import ATOM, HERE, ROOT, TOPOLOGY, q_n

LOGGER = logging.getLogger(__name__)
SANDER = Path("/data/jpc/code/ambertools26/venv/bin/sander")
ENERGY_RE = {
    "bias_inclusive_energy_kcal_mol": re.compile(r"EPtot\s*=\s*([-+0-9.EeDd]+)"),
    "restraint_energy_kcal_mol": re.compile(r"RESTRAINT\s*=\s*([-+0-9.EeDd]+)"),
}


def digest(path: Path) -> str:
    """Return a file SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def last_energy(path: Path, pattern: re.Pattern[str]) -> float:
    """Read one finite energy from a completed Amber output."""
    text = path.read_text(errors="replace")
    vals = pattern.findall(text)
    if not vals or "Energy calculation finished" not in text:
        raise ValueError(f"Missing completed energy field {pattern.pattern}: {path}")
    value = float(vals[-1].replace("D", "E").replace("d", "e"))
    if not np.isfinite(value):
        raise ValueError(f"Nonfinite energy in {path}")
    return value


def measured(restart: Path) -> dict[str, float]:
    """Measure key reaction coordinates and full monitored proton distances."""
    system = parmed.load_file(str(TOPOLOGY), xyz=str(restart))
    xyz = np.asarray(system.coordinates, dtype=float)

    def dist(a: int, b: int) -> float:
        return float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))

    result = {"q_N_A": q_n(xyz), "attack_A": dist(ATOM["thr_og1"], ATOM["substrate_c"])}
    pairs = {
        "thr_og1_hg1_A": ("thr_og1", "thr_hg1"),
        "thr_n_hg1_A": ("thr_n", "thr_hg1"),
        "thr_n_h1_A": ("thr_n", "thr_h1"),
        "thr_n_h2_A": ("thr_n", "thr_h2"),
        "water_o_h1_A": ("water_o", "water_h1"),
        "water_o_h2_A": ("water_o", "water_h2"),
        "water_o_hg1_A": ("water_o", "thr_hg1"),
        "asp272_od1_h2_A": ("asp272_od1", "water_h2"),
        "ash274_od2_hd2_A": ("ash274_od2", "ash274_hd2"),
        "ash274_hd2_water_o_A": ("ash274_hd2", "water_o"),
        "substrate_c_o_A": ("substrate_c", "substrate_o"),
        "substrate_c_n_A": ("substrate_c", "leaving_n"),
        "leaving_n_hn_A": ("leaving_n", "leaving_hn"),
        "thr_n_water_h1_A": ("thr_n", "water_h1"),
        "thr_n_water_h2_A": ("thr_n", "water_h2"),
        "leaving_n_hg1_A": ("leaving_n", "thr_hg1"),
        "leaving_n_water_h1_A": ("leaving_n", "water_h1"),
        "leaving_n_water_h2_A": ("leaving_n", "water_h2"),
        "substrate_o_phe113_n_A": ("substrate_o", "phe113_n"),
        "substrate_o_asn186_nd2_A": ("substrate_o", "asn186_nd2"),
    }
    result.update({key: dist(ATOM[a], ATOM[b]) for key, (a, b) in pairs.items()})
    return result


def journal(phase: str, action: str, result: str | None = None) -> None:
    """Write an append-only action record before or after an engine stage."""
    cmd = [
        sys.executable,
        str(ROOT / "scripts/record_agent_action.py"),
        "--log",
        str(ROOT / "output_agent_actions/20261021.jsonl"),
        "--agent",
        "activation_grid",
        "--phase",
        phase,
        "--permission-mode",
        "escalated",
        "--reason",
        "Run prepared activation grid row",
        "--action",
        action,
        "--rollback",
        "Preserve failed stages; retry only in a new attempt directory.",
    ]
    if result is not None:
        cmd.extend(("--result", result))
    subprocess.run(cmd, check=True, cwd=ROOT)


def update_points(records: list[dict[str, Any]]) -> None:
    """Atomically update completed point fields under an advisory file lock."""
    path = HERE / "grid_points.json"
    with (HERE / ".grid_points.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        data = json.loads(path.read_text())
        for key, values in records:
            data["points"][key].update(values)
        temp = path.with_suffix(".json.tmp")
        temp.write_text(json.dumps(data, indent=2) + "\n")
        os.replace(temp, path)


def run_row(row_number: int) -> None:
    """Run all sequential preparation and attack stages for one q row."""
    row_dir = HERE / f"row_{row_number:03d}"
    manifest = json.loads((row_dir / "row.json").read_text())
    prior_restart = row_dir / "starting.ncrst"
    point_records = []
    protocol_path = HERE / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    if protocol.get("status") != "prepared_not_run":
        raise ValueError("Protocol is not in prepared_not_run state")
    for source_path, expected_hash in protocol["source_sha256"].items():
        if digest(Path(source_path)) != expected_hash:
            raise ValueError(f"Prepared source hash mismatch: {source_path}")
    if not SANDER.is_file() or not (HERE / "grid_points.json").is_file():
        raise FileNotFoundError("Sander executable or grid progress manifest missing")
    row_manifest_hash = digest(row_dir / "row.json")
    row_status_path = row_dir / "row_status.json"
    if row_status_path.exists():
        raise FileExistsError(
            f"Row already started; refusing resume: {row_status_path}"
        )
    row_status = {
        "row": row_number,
        "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_sha256": digest(protocol_path),
        "row_manifest_sha256": row_manifest_hash,
    }
    row_status_path.write_text(json.dumps(row_status, indent=2) + "\n")
    for item in manifest["stages"]:
        stage_dir = HERE / item["directory"]
        if (stage_dir / "process.json").exists():
            raise FileExistsError(f"Refusing to overwrite stage {stage_dir}")
        shutil.copy2(prior_restart, stage_dir / "starting.ncrst")
        before = {
            name: digest(stage_dir / name)
            for name in ("drive.mdin", "restraints.RST", "starting.ncrst")
        }
        for name in ("drive.mdin", "restraints.RST"):
            if before[name] != item["input_sha256"][name]:
                raise ValueError(f"Prepared stage input changed: {stage_dir / name}")
        command = [
            str(SANDER),
            "-O",
            "-i",
            "drive.mdin",
            "-o",
            "minimized.mdout",
            "-p",
            str(TOPOLOGY),
            "-c",
            "starting.ncrst",
            "-r",
            "optimized.ncrst",
            "-x",
            "optimization.nc",
            "-inf",
            "mdinfo",
        ]
        env = os.environ.copy() | {
            "AMBERHOME": "/data/jpc/code/ambertools26/venv",
            "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
        }
        journal("intent", f"Run row {row_number} {item['stage']}: {' '.join(command)}")
        start = datetime.now(timezone.utc)
        proc = {
            "status": "running",
            "command": command,
            "row": row_number,
            "stage": item,
            "started_utc": start.isoformat(),
            "permission_mode": "escalated",
            "source_protocol_sha256": digest(protocol_path),
            "runner_sha256": digest(Path(__file__)),
            "row_manifest_sha256": row_manifest_hash,
            "environment": {
                key: env[key]
                for key in (
                    "AMBERHOME",
                    "LD_LIBRARY_PATH",
                    "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                )
            },
            "inputs_sha256": before,
        }
        (stage_dir / "process.json").write_text(json.dumps(proc, indent=2) + "\n")
        with (stage_dir / "engine.log").open("w") as log:
            process = subprocess.run(
                command,
                cwd=stage_dir,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        text = (
            (stage_dir / "minimized.mdout").read_text(errors="replace")
            if (stage_dir / "minimized.mdout").exists()
            else ""
        )
        converged = "Convergence reached" in text and "Optimisation converged" in text
        scc_ok = not any(
            term in text.lower()
            for term in ("scc not converged", "scf not converged", "scf failed")
        )
        completed = (
            process.returncode == 0
            and converged
            and scc_ok
            and (stage_dir / "optimized.ncrst").is_file()
        )
        proc = {
            "status": "finished" if completed else "failed",
            "exit_status": process.returncode,
            "force_converged": converged,
            "scc_converged": scc_ok,
            "started_utc": start.isoformat(),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "inputs_sha256": before,
            "command": command,
            "row": row_number,
            "stage": item,
            "permission_mode": "escalated",
            "source_protocol_sha256": digest(protocol_path),
            "runner_sha256": digest(Path(__file__)),
            "row_manifest_sha256": row_manifest_hash,
            "environment": {
                key: env[key]
                for key in (
                    "AMBERHOME",
                    "LD_LIBRARY_PATH",
                    "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                )
            },
            "outputs_sha256": {
                p.name: digest(p)
                for p in stage_dir.iterdir()
                if p.is_file() and p.name not in {"process.json"}
            },
        }
        (stage_dir / "process.json").write_text(json.dumps(proc, indent=2) + "\n")
        if not completed:
            (stage_dir / "acceptance.json").write_text(
                json.dumps(
                    {"accepted": False, "reason": "process/SCC/force failure"}, indent=2
                )
                + "\n"
            )
            if item["grid_key"]:
                update_points(
                    [
                        (
                            item["grid_key"],
                            {
                                "status": "failed",
                                "chemistry_class": "unclassified",
                                "failure_reason": "process/SCC/force failure",
                            },
                        )
                    ]
                )
            row_status.update(
                {
                    "status": "failed",
                    "failed_stage": item["stage"],
                    "finished_utc": datetime.now(timezone.utc).isoformat(),
                }
            )
            row_status_path.write_text(json.dumps(row_status, indent=2) + "\n")
            journal(
                "outcome",
                f"Run row {row_number} {item['stage']}",
                "failed; row stopped, artifacts preserved",
            )
            raise RuntimeError(
                f"Calculation failed explicit convergence gates: {stage_dir}"
            )
        energies = {
            key: last_energy(stage_dir / "minimized.mdout", pat)
            for key, pat in ENERGY_RE.items()
        }
        energies["physical_energy_kcal_mol"] = (
            energies["bias_inclusive_energy_kcal_mol"]
            - energies["restraint_energy_kcal_mol"]
        )
        geometry = measured(stage_dir / "optimized.ncrst")
        out_hash = digest(stage_dir / "optimized.ncrst")
        chem = (
            "retained_candidate"
            if geometry["substrate_c_n_A"] <= 1.55
            and geometry["thr_n_h1_A"] <= 1.25
            and geometry["thr_n_h2_A"] <= 1.25
            and geometry["ash274_od2_hd2_A"] <= 1.25
            and geometry["leaving_n_hn_A"] <= 1.25
            and geometry["water_o_h1_A"] <= 1.25
            else "alternate_chemistry"
        )
        accept = {
            "accepted": chem == "retained_candidate",
            "process_accepted": True,
            "force_converged": True,
            "scc_converged": True,
            "chemistry_class": chem,
            "geometry": geometry,
            "energies": energies,
            "restart_sha256": out_hash,
            "restraint_targets": {
                k: item[k] for k in ("q_target_A", "attack_target_A")
            },
        }
        (stage_dir / "acceptance.json").write_text(json.dumps(accept, indent=2) + "\n")
        if item["grid_key"]:
            point_records.append(
                (
                    item["grid_key"],
                    {
                        "status": "completed",
                        "q_actual_A": geometry["q_N_A"],
                        "attack_actual_A": geometry["attack_A"],
                        "physical_energy_kcal_mol": energies[
                            "physical_energy_kcal_mol"
                        ],
                        "bias_inclusive_energy_kcal_mol": energies[
                            "bias_inclusive_energy_kcal_mol"
                        ],
                        "restraint_energy_kcal_mol": energies[
                            "restraint_energy_kcal_mol"
                        ],
                        "restart_sha256": out_hash,
                        "chemistry_class": chem,
                        "accepted": chem == "retained_candidate",
                        "geometry": geometry,
                    },
                )
            )
            update_points(point_records)
            point_records.clear()
        journal(
            "outcome",
            f"Run row {row_number} {item['stage']}",
            f"completed; chemistry={chem}; restart_sha256={out_hash}",
        )
        prior_restart = stage_dir / "optimized.ncrst"
        row_status.update(
            {
                "current_stage": item["stage"],
                "last_completed_stage": item["stage"],
                "updated_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        row_status_path.write_text(json.dumps(row_status, indent=2) + "\n")
    row_status.update(
        {"status": "completed", "finished_utc": datetime.now(timezone.utc).isoformat()}
    )
    row_status_path.write_text(json.dumps(row_status, indent=2) + "\n")


def main() -> None:
    """Run one row; preparation is required and no resume/overwrite is allowed."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--row", type=int, required=True, choices=range(1, 6))
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    run_row(args.row)


if __name__ == "__main__":
    main()
