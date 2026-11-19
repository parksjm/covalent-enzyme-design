"""Prepare and sequentially run validated unbiased release candidates."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Any

from release_queue_support import (
    GRID,
    PROTOCOL,
    QUEUE_GATE,
    ROOT,
    STATUS_PATH,
    TOPOLOGY,
    digest,
    geometry,
    load_json,
)

SOURCE_CASE = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/06_qmmm_optimization/replica_1_tier_A_seed/water_qm_reoptimization"
)
SOURCE_MDIN = SOURCE_CASE / "minimize.mdin"
SANDER = Path("/data/jpc/code/ambertools26/venv/bin/sander")
ENERGIES = {
    "ep_tot_kcal_mol": re.compile(r"EPtot\s*=\s*([-+0-9.EeDd]+)"),
    "restraint_kcal_mol": re.compile(r"RESTRAINT\s*=\s*([-+0-9.EeDd]+)"),
}


def journal(phase: str, action: str, result: str | None = None) -> None:
    """Append a required action record for each unbiased release."""
    command = [
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
        "Run gated candidate release",
        "--action",
        action,
        "--rollback",
        "Preserve terminal release evidence; retries require a new attempt directory.",
    ]
    if result is not None:
        command.extend(("--result", result))
    subprocess.run(command, check=True, cwd=ROOT)


def write_status(status: str, reason: str, **extra: Any) -> None:
    """Write an append-safe single terminal queue status document."""
    record = {
        "schema_version": 1,
        "status": status,
        "reason": reason,
        "updated_utc": datetime.now(timezone.utc).isoformat(),
    } | extra
    STATUS_PATH.write_text(json.dumps(record, indent=2) + "\n")


def energy(mdout: Path) -> dict[str, float]:
    """Extract total and restraint energy from a completed Amber output."""
    text = mdout.read_text(errors="replace")
    output = {}
    for name, pattern in ENERGIES.items():
        values = pattern.findall(text)
        if not values or "Energy calculation finished" not in text:
            raise ValueError(f"Missing completed energy field {name} in {mdout}")
        output[name] = float(values[-1].replace("D", "E").replace("d", "e"))
        if not math.isfinite(output[name]):
            raise ValueError(f"Nonfinite {name} in {mdout}")
    output["physical_kcal_mol"] = (
        output["ep_tot_kcal_mol"] - output["restraint_kcal_mol"]
    )
    return output


def run_release(folder: Path, candidate: dict[str, Any]) -> dict[str, Any]:
    """Run and independently classify one sequential unbiased release."""
    if (folder / "process.json").exists():
        raise FileExistsError(f"Release already has process provenance: {folder}")
    release_manifest = load_json(folder / "release.json")
    if digest(folder / "starting.ncrst") != release_manifest["starting_restart_sha256"]:
        raise ValueError(f"Release starting restart hash mismatch: {folder}")
    if (
        digest(folder / "unrestrained_relax.mdin")
        != release_manifest["prepared_unrestrained_mdin_sha256"]
    ):
        raise ValueError(f"Release mdin hash mismatch: {folder}")
    if digest(QUEUE_GATE) != release_manifest["reference_gate_sha256"]:
        raise ValueError("Reference comparison gate changed after release preparation")
    for source, expected in release_manifest["runtime_source_sha256"].items():
        if digest(Path(source)) != expected:
            raise ValueError(
                f"Release runner source changed after preparation: {source}"
            )
    snapshot = GRID / release_manifest["runtime_source_snapshot"]
    for source, expected in release_manifest["runtime_source_sha256"].items():
        if digest(snapshot / Path(source).name) != expected:
            raise ValueError(f"Release source snapshot mismatch: {source}")
    env = os.environ.copy() | {
        "AMBERHOME": "/data/jpc/code/ambertools26/venv",
        "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    command = [
        str(SANDER),
        "-O",
        "-i",
        "unrestrained_relax.mdin",
        "-p",
        str(TOPOLOGY),
        "-c",
        "starting.ncrst",
        "-o",
        "release.mdout",
        "-r",
        "release.ncrst",
        "-x",
        "release.nc",
        "-inf",
        "release.mdinfo",
    ]
    journal(
        "intent",
        f"Run {candidate['class_name']} release from {candidate['key']}: {' '.join(command)}",
    )
    process = {
        "status": "running",
        "command": command,
        "candidate_class": candidate["class_name"],
        "grid_point": candidate["key"],
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "permission_mode": "escalated",
        "environment": {
            k: env[k]
            for k in (
                "AMBERHOME",
                "LD_LIBRARY_PATH",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
            )
        },
        "source_protocol_sha256": digest(PROTOCOL),
        "reference_gate_sha256": digest(QUEUE_GATE),
        "candidate_restart_sha256": digest(folder / "starting.ncrst"),
        "unrestrained_mdin_sha256": digest(folder / "unrestrained_relax.mdin"),
        "topology_sha256": digest(TOPOLOGY),
        "runner_sha256": digest(Path(__file__)),
        "runtime_source_sha256": load_json(folder / "release.json")[
            "runtime_source_sha256"
        ],
    }
    (folder / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    with (folder / "engine.log").open("w") as log:
        result = subprocess.run(
            command,
            cwd=folder,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    output = (
        (folder / "release.mdout").read_text(errors="replace")
        if (folder / "release.mdout").is_file()
        else ""
    )
    converged = (
        result.returncode == 0
        and "Convergence reached" in output
        and "Optimisation converged" in output
    )
    error_pattern = re.compile(
        r"(?:SCC|SCF).{0,100}(?:fail|not converg|did not converge|maximum iteration)|(?:fail|not converg).{0,100}(?:SCC|SCF)",
        re.I,
    )
    scc_errors = error_pattern.findall(output)
    scc_ok = "DFTBESCF=" in output and not scc_errors
    process.update(
        {
            "status": "finished" if result.returncode == 0 else "failed",
            "exit_status": result.returncode,
            "force_converged": converged,
            "scc_converged": scc_ok,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "output_sha256": {
                p.name: digest(p)
                for p in folder.iterdir()
                if p.is_file() and p.name != "process.json"
            },
        }
    )
    (folder / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    accepted = False
    measured = None
    energies = None
    chemistry_gates: dict[str, bool] = {}
    reason = "process/SCC/force failure"
    if (
        process["status"] == "finished"
        and converged
        and scc_ok
        and (folder / "release.ncrst").is_file()
    ):
        measured = geometry(folder / "release.ncrst")
        energies = energy(folder / "release.mdout")
        chemistry_gates = {
            "thr_n_h1_h2_retained": measured["thr_n_h1_A"] <= 1.25
            and measured["thr_n_h2_A"] <= 1.25,
            "ash274_hd2_retained": measured["ash274_od2_hd2_A"] <= 1.25
            and measured["ash274_hd2_water_o_A"] >= 1.30,
            "water_h1_retained": measured["water_o_h1_A"] <= 1.25,
            "substrate_leaving_hn_retained": measured["leaving_n_hn_A"] <= 1.25,
            "scissile_cn_retained": measured["substrate_c_n_A"] <= 1.55,
        }
        if candidate["class_name"] == "activated":
            chemistry_gates["activated_thr_n_hg1"] = measured["thr_n_hg1_A"] <= 1.25
        else:
            chemistry_gates.update(
                {
                    "ti_like_thr_n_hg1": measured["thr_n_hg1_A"] <= 1.25,
                    "attack_geometry_retained": measured["attack_A"] <= 1.75,
                    "tetrahedral_co": measured["substrate_c_o_A"] >= 1.30,
                }
            )
        accepted = all(chemistry_gates.values())
        reason = (
            "all release gates passed"
            if accepted
            else "unrestrained release changed required geometry or proton inventory"
        )
    scientific = {
        "process_status": process["status"],
        "exit_status": process["exit_status"],
        "force_converged": converged,
        "scc_converged": scc_ok,
        "chemistry_gates": chemistry_gates,
        "chemistry_accepted": accepted,
        "accepted": accepted,
        "reason": reason,
        "geometry": measured,
        "energies": energies,
    }
    (folder / "scientific_acceptance.json").write_text(
        json.dumps(scientific, indent=2) + "\n"
    )
    journal(
        "outcome",
        f"Run {candidate['class_name']} release from {candidate['key']}",
        f"process={process['status']}; force={converged}; SCC={scc_ok}; chemistry_accepted={accepted}: {reason}",
    )
    return scientific


def prepare_releases(
    candidates: list[dict[str, Any]], decisions: dict[str, Any], gate: dict[str, Any]
) -> list[Path]:
    """Write at most two immutable release inputs and source manifests."""
    root = GRID / "candidate_releases"
    if root.exists():
        raise FileExistsError(f"Release preparation exists; refusing overwrite: {root}")
    root.mkdir()
    source_files = [
        Path(__file__),
        Path(__file__).with_name("release_queue_support.py"),
        Path(__file__).with_name("queue_candidate_releases.py"),
    ]
    snapshot = root / "release_queue_source_snapshot"
    snapshot.mkdir()
    source_hashes = {str(path): digest(path) for path in source_files}
    for source in source_files:
        shutil.copy2(source, snapshot / source.name)
        if digest(snapshot / source.name) != source_hashes[str(source)]:
            raise ValueError(f"Release runtime source snapshot mismatch: {source}")
    release_dirs = []
    for index, candidate in enumerate(candidates, 1):
        folder = root / f"release_{index:02d}_{candidate['class_name']}"
        folder.mkdir()
        shutil.copy2(candidate["restart_path"], folder / "starting.ncrst")
        shutil.copy2(SOURCE_MDIN, folder / "unrestrained_relax.mdin")
        manifest = {
            "candidate_class": candidate["class_name"],
            "grid_point": candidate["key"],
            "grid_restart_sha256": candidate["point"]["restart_sha256"],
            "starting_restart_sha256": digest(folder / "starting.ncrst"),
            "source_unrestrained_mdin_sha256": digest(SOURCE_MDIN),
            "prepared_unrestrained_mdin_sha256": digest(
                folder / "unrestrained_relax.mdin"
            ),
            "physical_grid_energy_kcal_mol": candidate["energy"],
            "geometry": candidate["geometry"],
            "source_acceptance_record": str(
                candidate["acceptance_path"].relative_to(GRID)
            ),
            "reference_gate_sha256": digest(QUEUE_GATE),
            "runtime_source_sha256": source_hashes,
            "runtime_source_snapshot": str(snapshot.relative_to(GRID)),
            "reference_gate_source_sha256": gate["source_grid_protocol_sha256"],
            "grid_protocol_sha256": digest(PROTOCOL),
            "topology_sha256": digest(TOPOLOGY),
            "status": "prepared_not_run",
        }
        (folder / "release.json").write_text(json.dumps(manifest, indent=2) + "\n")
        release_dirs.append(folder)
    write_status(
        "prepared_not_run",
        "reference gate passed; candidates selected",
        selected=decisions,
        reference_gate_sha256=digest(QUEUE_GATE),
        grid_protocol_sha256=digest(PROTOCOL),
        releases=[str(p.relative_to(GRID)) for p in release_dirs],
    )
    return release_dirs
