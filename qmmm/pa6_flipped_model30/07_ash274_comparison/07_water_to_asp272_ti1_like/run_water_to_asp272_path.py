"""Execute prepared water-to-Asp272 QM/MM windows inside tmux."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Any

import numpy as np
import parmed

from prepare_water_to_asp272_path import (
    CASE,
    DISTANCE_KEYS,
    JOURNAL,
    ROOT,
    SANDER,
    SOURCE_CASE,
    TOPOLOGY,
    WINDOWS,
    geometry,
    sha256,
)

LOGGER = logging.getLogger(__name__)
RECORD_ACTION = ROOT / "scripts/record_agent_action.py"
EXPORTER = CASE / "export_water_path.py"
ENERGY_PATTERNS = {
    "EPtot_bias_inclusive_kcal_mol": re.compile(r"EPtot\s*=\s*([-+0-9.EeDd]+)"),
    "restraint_kcal_mol": re.compile(r"RESTRAINT\s*=\s*([-+0-9.EeDd]+)"),
}
REACTIVE_CONTINUITY_KEYS = (
    "thr_og1_hg1_A",
    "water_o_hg1_A",
    "water_o_h2_A",
    "asp272_od1_h2_A",
    "og1_substrate_c_A",
    "substrate_c_o_A",
    "substrate_c_n_A",
)
CONTINUITY_LIMIT_A = 0.35
RESUME_WINDOW = 5


def journal(phase: str, action: str, result: str | None = None) -> None:
    """Append a required pre/post action record to the active journal."""
    command = [
        sys.executable,
        str(RECORD_ACTION),
        "--log",
        str(JOURNAL),
        "--agent",
        "water_ti1_path_script",
        "--phase",
        phase,
        "--permission-mode",
        "tool",
        "--reason",
        "Execute reviewed water-to-Asp272 path window",
        "--action",
        action,
        "--rollback",
        "Preserve all completed windows; any retry uses a new attempt directory.",
    ]
    if result is not None:
        command.extend(("--result", result))
    subprocess.run(command, check=True, cwd=ROOT)


def last_energy(path: Path, pattern: re.Pattern[str]) -> float:
    """Read the last finite energy matching a verified Amber output field."""
    text = path.read_text(errors="replace")
    matches = pattern.findall(text)
    if not matches or "Energy calculation finished" not in text:
        raise ValueError(f"Missing energy field or completion marker in {path}")
    value = float(matches[-1].replace("D", "E").replace("d", "e"))
    if not np.isfinite(value):
        raise ValueError(f"Nonfinite energy in {path}")
    return value


def energy_record(mdout: Path) -> dict[str, float]:
    """Return restraint-inclusive and physical hybrid energy from one Sander block."""
    values = {
        name: last_energy(mdout, pattern) for name, pattern in ENERGY_PATTERNS.items()
    }
    values["physical_kcal_mol"] = (
        values["EPtot_bias_inclusive_kcal_mol"] - values["restraint_kcal_mol"]
    )
    return values


def interim_gates(values: dict[str, float], thr_h3_present: bool) -> dict[str, bool]:
    """Require the intact scissile amide and unaltered N-H inventories."""
    return {
        "scissile_cn_intact": values["substrate_c_n_A"] <= 1.55,
        "thr_n_remains_neutral_nh2": not thr_h3_present
        and values["thr_n_h1_A"] <= 1.25
        and values["thr_n_h2_A"] <= 1.25,
        "substrate_leaving_n_keeps_hn": values["leaving_n_hn_A"] <= 1.25,
        "water_h1_stays_on_water": values["water_o_h1_A"] <= 1.25,
        "ash274_hd2_stays_on_od2": values["ash274_od2_hd2_A"] <= 1.25
        and values["ash274_hd2_water_o_A"] >= 1.30,
        "no_unintended_thr_n_protonation": min(
            values["thr_n_hg1_A"],
            values["thr_n_water_h1_A"],
            values["thr_n_water_h2_A"],
        )
        >= 1.30,
        "no_unintended_leaving_n_protonation": min(
            values["leaving_n_hg1_A"],
            values["leaving_n_water_h1_A"],
            values["leaving_n_water_h2_A"],
        )
        >= 1.30,
    }


def endpoint_gates(values: dict[str, float], thr_h3_present: bool) -> dict[str, bool]:
    """Classify the alternative Asp272-protonated TI1-like endpoint."""
    return interim_gates(values, thr_h3_present) | {
        "thr_proton_transferred_to_water": values["thr_og1_hg1_A"] >= 1.30
        and values["water_o_hg1_A"] <= 1.20,
        "water_h2_transferred_to_asp272": values["water_o_h2_A"] >= 1.30
        and values["asp272_od1_h2_A"] <= 1.20,
        "tetrahedral_attack_geometry": values["og1_substrate_c_A"] <= 1.75
        and values["substrate_c_o_A"] >= 1.30,
        "scissile_amide_remains_intact": values["substrate_c_n_A"] <= 1.55,
        "oxyanion_hole_retained": values["substrate_o_phe113_n_A"] <= 3.5
        and values["substrate_o_asn186_nd2_A"] <= 3.5,
    }


def adjacent_changes(
    previous: dict[str, float], current: dict[str, float]
) -> dict[str, float]:
    """Calculate coordinate changes between neighboring accepted windows."""
    return {key: current[key] - previous[key] for key in DISTANCE_KEYS}


def continuity_passes(changes: dict[str, float]) -> bool:
    """Apply the adjacent-window gate to declared reactive distances only."""
    return all(
        abs(changes[key]) <= CONTINUITY_LIMIT_A for key in REACTIVE_CONTINUITY_KEYS
    )


def verify_output_hashes(window: Path, acceptance: dict[str, Any]) -> None:
    """Require every previously recorded output digest to match disk."""
    for name, expected in acceptance["output_sha256"].items():
        path = window / name
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"Window output hash mismatch: {path}")
    restart = window / "optimized.ncrst"
    if sha256(restart) != acceptance["restart_sha256"]:
        raise ValueError(f"Window restart hash mismatch: {restart}")


def resume_history(
    protocol: dict[str, Any], initial_energy: dict[str, float]
) -> list[dict[str, Any]]:
    """Validate windows 001–004 and re-evaluate frozen window 005 for resume."""
    if protocol.get("status") != "path_running":
        raise ValueError("Resume requires protocol status path_running")
    accepted_window = protocol.get("last_accepted_window")
    if accepted_window not in {"window_004", "window_005"}:
        raise ValueError(
            "Resume requires last_accepted_window window_004 or window_005"
        )

    records: list[dict[str, Any]] = [
        {
            "window": "window_000",
            "geometry": protocol["coordinates"]["starting_geometry"],
            "energy": initial_energy,
            "relative_physical_kcal_mol": 0.0,
            "relative_EPtot_bias_inclusive_kcal_mol": 0.0,
            "accepted": True,
        }
    ]
    for index in range(1, RESUME_WINDOW):
        window = CASE / f"window_{index:03d}"
        acceptance = json.loads((window / "acceptance.json").read_text())
        if not acceptance.get("accepted"):
            raise ValueError(f"Prior window is not accepted: {window.name}")
        verify_output_hashes(window, acceptance)
        process = json.loads((window / "process.json").read_text())
        previous_restart = CASE / f"window_{index - 1:03d}" / "optimized.ncrst"
        if process.get("starting_restart_sha256") != sha256(previous_restart):
            raise ValueError(
                f"Prior starting restart provenance mismatch: {window.name}"
            )
        records.append(acceptance)

    window = CASE / f"window_{RESUME_WINDOW:03d}"
    revision = window / "acceptance_revision_1.json"
    if accepted_window == "window_005" and not revision.is_file():
        raise FileNotFoundError(f"Resumed acceptance backup is missing: {revision}")
    process = json.loads((window / "process.json").read_text())
    previous_restart = CASE / "window_004" / "optimized.ncrst"
    if process.get("status") != "finished" or process.get("exit_status") != 0:
        raise ValueError("window_005 process record is not successfully finished")
    if process.get("starting_restart_sha256") != sha256(previous_restart):
        raise ValueError("window_005 starting restart provenance mismatch")
    if sha256(window / "starting.ncrst") != process["starting_restart_sha256"]:
        raise ValueError("window_005 starting restart file does not match provenance")
    if process.get("drive_mdin_sha256") != sha256(window / "drive.mdin"):
        raise ValueError("window_005 drive input hash mismatch")
    if process.get("restraints_sha256") != sha256(window / "restraints.RST"):
        raise ValueError("window_005 restraint input hash mismatch")
    for name in ("minimized.mdout", "optimized.ncrst", "optimization.nc"):
        if not (window / name).is_file():
            raise FileNotFoundError(
                f"Required completed output missing: {window / name}"
            )
    output_text = (window / "minimized.mdout").read_text(errors="replace")
    if (
        "Convergence reached" not in output_text
        or "Optimisation converged" not in output_text
    ):
        raise ValueError("window_005 output lacks force-convergence markers")
    old_acceptance = json.loads((window / "acceptance.json").read_text())
    verify_output_hashes(window, old_acceptance)
    if revision.exists():
        original_acceptance = json.loads(revision.read_text())
        verify_output_hashes(window, original_acceptance)
        if (
            not old_acceptance.get("accepted")
            or old_acceptance.get("re_evaluated_under_protocol_revision") != 5
        ):
            raise ValueError("Existing acceptance revision is incomplete or invalid")
    elif old_acceptance.get("re_evaluated_under_protocol_revision") == 5:
        raise ValueError("Revision-5 acceptance exists without its preserved original")
    if not (window / "window.json").is_file():
        raise FileNotFoundError("window_005 prepared input manifest is missing")
    if revision.exists():
        records.append(old_acceptance)
    return records


def write_ledgers(records: list[dict[str, Any]]) -> None:
    """Atomically refresh compact energy and full raw-distance CSV ledgers."""
    temp = CASE / "energy_ledger.csv.tmp"
    with temp.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("window", "energy_kcal_mol", "energy_bias_inclusive_kcal_mol"),
        )
        writer.writeheader()
        writer.writerows(
            {
                "window": row["window"],
                "energy_kcal_mol": row["relative_physical_kcal_mol"],
                "energy_bias_inclusive_kcal_mol": row[
                    "relative_EPtot_bias_inclusive_kcal_mol"
                ],
            }
            for row in records
        )
    temp.replace(CASE / "energy_ledger.csv")
    temp = CASE / "distance_ledger.csv.tmp"
    fields = ["window", *DISTANCE_KEYS]
    with temp.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {"window": row["window"], **row["geometry"]} for row in records
        )
    temp.replace(CASE / "distance_ledger.csv")


def export_progress(records: list[dict[str, Any]]) -> None:
    """Rebuild cumulative VMD progress and plots after an accepted window."""
    if not EXPORTER.exists():
        raise FileNotFoundError(f"Required live path exporter is absent: {EXPORTER}")
    command = [sys.executable, str(EXPORTER), "--topology", str(TOPOLOGY)]
    for row in records:
        command.extend(
            ("--accepted-restarts", str(CASE / row["window"] / "optimized.ncrst"))
        )
    command.extend(("--energy-ledger", str(CASE / "energy_ledger.csv")))
    command.extend(("--distance-ledger", str(CASE / "distance_ledger.csv")))
    command.extend(("--output", str(CASE / "vmd_progress")))
    journal(
        "intent",
        "Refresh cumulative VMD/energy/distance progress after accepted window",
    )
    result = subprocess.run(
        command, cwd=CASE, check=False, capture_output=True, text=True
    )
    if result.returncode:
        journal(
            "outcome",
            "Refresh cumulative VMD/energy/distance progress",
            f"failed rc={result.returncode}: {result.stderr[-1000:]}",
        )
        raise RuntimeError(
            "Live path export failed; accepted scientific window is preserved"
        )
    journal("outcome", "Refresh cumulative VMD/energy/distance progress", "completed")


def run_release(protocol: dict[str, Any]) -> dict[str, Any]:
    """Run a separate unrestrained release only after the driven endpoint passes."""
    release = CASE / "release_pending"
    if protocol.get("status") != "endpoint_accepted":
        raise ValueError("Unrestrained release requires an accepted driven endpoint")
    if (release / "release.mdout").exists():
        raise FileExistsError("Release output exists; refusing to overwrite it")
    source_window = CASE / str(protocol["last_accepted_window"])
    source_restart = source_window / "optimized.ncrst"
    shutil.copy2(source_restart, release / "starting.ncrst")
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
        "intent", f"Run separate unrestrained endpoint release: {' '.join(command)}"
    )
    process = {
        "phase": "unrestrained_release",
        "command": command,
        "starting_restart": str(source_restart.resolve()),
        "starting_restart_sha256": sha256(source_restart),
        "input_sha256": sha256(release / "unrestrained_relax.mdin"),
        "runner_sha256": sha256(Path(__file__)),
        "exporter_sha256": sha256(EXPORTER),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
    }
    (release / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    environment = os.environ.copy() | {
        "AMBERHOME": "/data/jpc/code/ambertools26/venv",
        "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
        "OMP_NUM_THREADS": "1",
    }
    with (release / "engine.log").open("w") as log:
        result = subprocess.run(
            command,
            cwd=release,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    output = (release / "release.mdout").read_text(errors="replace")
    converged = (
        result.returncode == 0
        and "Convergence reached" in output
        and "Optimisation converged" in output
    )
    process.update(
        {
            "exit_status": result.returncode,
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "status": "finished" if result.returncode == 0 else "failed",
        }
    )
    (release / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    if not converged:
        acceptance = {
            "accepted": False,
            "force_converged": False,
            "exit_status": result.returncode,
        }
    else:
        endpoint = parmed.load_file(str(TOPOLOGY), xyz=str(release / "release.ncrst"))
        values = geometry(np.asarray(endpoint.coordinates, dtype=float))
        final = endpoint_gates(values, False)
        acceptance = {
            "accepted": all(final.values()),
            "force_converged": True,
            "gates": final,
            "geometry": values,
            "energy": energy_record(release / "release.mdout"),
            "restart_sha256": sha256(release / "release.ncrst"),
        }
    (release / "acceptance.json").write_text(json.dumps(acceptance, indent=2) + "\n")
    protocol["release"] = acceptance
    protocol["status"] = (
        "endpoint_accepted_release_accepted"
        if acceptance["accepted"]
        else "endpoint_accepted_release_failed"
    )
    (CASE / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    journal(
        "outcome",
        "Run separate unrestrained endpoint release",
        json.dumps(acceptance, separators=(",", ":")),
    )
    return protocol


def preflight(resume: bool = False) -> dict[str, Any]:
    """Validate source, restrained-energy parsing, windows, and exporter without launching."""
    protocol = json.loads((CASE / "protocol.json").read_text())
    expected_status = "path_running" if resume else "prepared_not_run"
    if protocol.get("status") != expected_status:
        raise ValueError(f"Preflight expects protocol status {expected_status}")
    if not EXPORTER.is_file():
        raise FileNotFoundError(f"Required live exporter not found: {EXPORTER}")
    for source_name, expected_hash in protocol["sources_sha256"].items():
        path = Path(source_name)
        if not path.is_file() or expected_hash != sha256(path):
            raise ValueError(
                f"Prepared source hash changed since input generation: {path}"
            )
    source_energy = energy_record(SOURCE_CASE / "minimized.mdout")
    restraint_example = energy_record(
        ROOT
        / "qmmm/pa6_flipped_model30/05_neutral_reaction_path/04_direct_acylenzyme"
        / "work/01_coupled_coordinate_drive/asynchronous_05/minimized.mdout"
    )
    if source_energy["restraint_kcal_mol"] != 0.0:
        raise ValueError(
            "Source E-S baseline unexpectedly has a nonzero restraint energy"
        )
    if restraint_example["restraint_kcal_mol"] == 0.0:
        raise ValueError(
            "Restraint parser did not recover the validated nonzero restraint energy"
        )
    if not all(
        (CASE / f"window_{index:03d}" / name).is_file()
        for index in range(1, WINDOWS + 1)
        for name in ("drive.mdin", "restraints.RST", "window.json")
    ):
        raise FileNotFoundError("Prepared path window inputs are incomplete")
    if resume:
        if protocol.get("protocol_revision") != 5:
            raise ValueError("Resume requires frozen protocol revision 5")
        accepted_window = protocol.get("last_accepted_window")
        if accepted_window not in {"window_004", "window_005"}:
            raise ValueError(
                "Resume requires last_accepted_window window_004 or window_005"
            )
        if any(
            (CASE / f"window_{index:03d}" / "minimized.mdout").exists()
            for index in range(RESUME_WINDOW + 1, WINDOWS + 1)
        ):
            raise FileExistsError("Later calculation output exists; refusing resume")
        validated_resume_windows = len(resume_history(protocol, source_energy))
    else:
        validated_resume_windows = 0
    source, _, start = __import__("prepare_water_to_asp272_path").verify_source()
    if protocol["qm"]["real_atoms"] != 102 or protocol["qm"]["charge"] != -1:
        raise ValueError("Prepared QM inventory changed")
    if protocol["active"]["atoms"] != 4662:
        raise ValueError("Prepared active mask changed")
    if (
        abs(
            start["relay_Q_A"]
            - protocol["coordinates"]["starting_geometry"]["relay_Q_A"]
        )
        > 1e-8
    ):
        raise ValueError("Prepared start geometry differs from accepted source")
    previous_q = start["relay_Q_A"]
    previous_attack = start["og1_substrate_c_A"]
    q_steps = []
    attack_steps = []
    for window in protocol["windows"]:
        q_target = float(window["relay_Q_target_A"])
        attack_target = float(window["attack_target_A"])
        q_steps.append(abs(q_target - previous_q))
        attack_steps.append(abs(attack_target - previous_attack))
        previous_q = q_target
        previous_attack = attack_target
    max_q_step = max(q_steps)
    max_attack_step = max(attack_steps)
    if max_q_step > 0.15 or max_attack_step > 0.08:
        raise ValueError("Prepared target increments exceed the reviewed step bounds")
    return {
        "passed": True,
        "status": "preflight_only_no_calculation_launched",
        "source_energy": source_energy,
        "restraint_syntax_energy_smoke": restraint_example,
        "source_atoms": len(source.atoms),
        "window_count": WINDOWS,
        "resume": resume,
        "validated_resume_windows": validated_resume_windows,
        "maximum_target_steps_A": {
            "Q": max_q_step,
            "attack": max_attack_step,
        },
        "runner_sha256": sha256(Path(__file__)),
        "exporter_sha256": sha256(EXPORTER),
    }


def execute(resume: bool = False) -> dict[str, Any]:
    """Run windows sequentially, gate each endpoint, and preserve all outcomes."""
    protocol_path = CASE / "protocol.json"
    if not protocol_path.is_file():
        raise FileNotFoundError(
            "Run --prepare and review protocol.json before execution"
        )
    protocol = json.loads(protocol_path.read_text())
    preflight(resume=resume)
    if resume:
        history = resume_history(
            protocol, energy_record(SOURCE_CASE / "minimized.mdout")
        )
    else:
        if protocol.get("status") != "prepared_not_run":
            raise ValueError("Path protocol is not in its prepared state")
        existing = [
            CASE / f"window_{index:03d}" / "minimized.mdout"
            for index in range(1, WINDOWS + 1)
            if (CASE / f"window_{index:03d}" / "minimized.mdout").exists()
        ]
        if existing:
            raise FileExistsError(
                f"Existing calculation output prevents overwrite: {existing[0]}"
            )
        history = None

    initial_energy = energy_record(SOURCE_CASE / "minimized.mdout")
    source_topology = parmed.load_file(str(TOPOLOGY))
    thr_h3_present = "H3" in {
        atom.name.strip() for atom in source_topology.residues[928].atoms
    }
    if thr_h3_present:
        raise ValueError("Canonical source has unexpected Thr929 H3 proton")
    exporter_hash = sha256(EXPORTER)
    runner_hash = sha256(Path(__file__))
    protocol["execution_sources_sha256"] = {
        "runner": runner_hash,
        "exporter": exporter_hash,
        "sander": sha256(SANDER),
    }
    protocol_path.write_text(json.dumps(protocol, indent=2) + "\n")
    records: list[dict[str, Any]] = [
        {
            "window": "window_000",
            "geometry": protocol["coordinates"]["starting_geometry"],
            "energy": initial_energy,
            "relative_physical_kcal_mol": 0.0,
            "relative_EPtot_bias_inclusive_kcal_mol": 0.0,
            "accepted": True,
        }
    ]
    if resume:
        records = history
        if records[-1]["window"] != "window_005":
            previous_values = records[-1]["geometry"]
            window = CASE / "window_005"
            system = parmed.load_file(
                str(TOPOLOGY), xyz=str(window / "optimized.ncrst")
            )
            values = geometry(np.asarray(system.coordinates, dtype=float))
            energy = energy_record(window / "minimized.mdout")
            changes = adjacent_changes(previous_values, values)
            continuity = continuity_passes(changes)
            gates = interim_gates(values, thr_h3_present)
            final_gates = endpoint_gates(values, thr_h3_present)
            acceptance = {
                "window": window.name,
                "geometry": values,
                "energy": energy,
                "relative_physical_kcal_mol": energy["physical_kcal_mol"]
                - initial_energy["physical_kcal_mol"],
                "relative_EPtot_bias_inclusive_kcal_mol": energy[
                    "EPtot_bias_inclusive_kcal_mol"
                ]
                - initial_energy["EPtot_bias_inclusive_kcal_mol"],
                "energy_change_from_previous_kcal_mol": energy["physical_kcal_mol"]
                - records[-1]["energy"]["physical_kcal_mol"],
                "geometry_change_from_previous": changes,
                "continuity_keys": list(REACTIVE_CONTINUITY_KEYS),
                "continuity_limit_A": CONTINUITY_LIMIT_A,
                "continuity": continuity,
                "gates": interim_gates(values, thr_h3_present),
                "accepted": all(gates.values()) and continuity,
                "endpoint_gates": final_gates,
                "endpoint_chemistry_accepted": all(final_gates.values()),
                "restart_sha256": sha256(window / "optimized.ncrst"),
                "output_sha256": {
                    name: sha256(window / name)
                    for name in (
                        "minimized.mdout",
                        "optimized.ncrst",
                        "optimization.nc",
                        "mdinfo",
                        "restraint_values.dat",
                    )
                    if (window / name).is_file()
                },
                "re_evaluated_under_protocol_revision": 5,
            }
            if not acceptance["accepted"]:
                raise ValueError(
                    "Existing window_005 outputs fail revision-5 acceptance"
                )
            revision = window / "acceptance_revision_1.json"
            with (window / "acceptance.json").open("rb") as source:
                with revision.open("xb") as destination:
                    shutil.copyfileobj(source, destination)
            (window / "acceptance.json").write_text(
                json.dumps(acceptance, indent=2) + "\n"
            )
            records.append(acceptance)
        protocol["last_accepted_window"] = "window_005"
        protocol["status"] = "path_running"
        protocol_path.write_text(json.dumps(protocol, indent=2) + "\n")
        write_ledgers(records)
        export_progress(records)
    else:
        write_ledgers(records)
    first_index = RESUME_WINDOW + 1 if resume else 1
    for index in range(first_index, WINDOWS + 1):
        window = CASE / f"window_{index:03d}"
        previous = CASE / f"window_{index - 1:03d}" / "optimized.ncrst"
        if not previous.is_file():
            raise FileNotFoundError(f"Missing prior accepted restart: {previous}")
        shutil.copy2(previous, window / "starting.ncrst")
        command = [
            str(SANDER),
            "-O",
            "-i",
            "drive.mdin",
            "-p",
            str(TOPOLOGY),
            "-c",
            "starting.ncrst",
            "-o",
            "minimized.mdout",
            "-r",
            "optimized.ncrst",
            "-x",
            "optimization.nc",
            "-inf",
            "mdinfo",
        ]
        journal(
            "intent", f"Run {window.name} DFTB3/MM path window: {' '.join(command)}"
        )
        process = {
            "window": window.name,
            "command": command,
            "cwd": str(window.resolve()),
            "starting_restart": str(previous.resolve()),
            "starting_restart_sha256": sha256(previous),
            "drive_mdin_sha256": sha256(window / "drive.mdin"),
            "restraints_sha256": sha256(window / "restraints.RST"),
            "status": "running",
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "runner_sha256": runner_hash,
            "exporter_sha256": exporter_hash,
        }
        (window / "process.json").write_text(json.dumps(process, indent=2) + "\n")
        environment = os.environ.copy() | {
            "AMBERHOME": "/data/jpc/code/ambertools26/venv",
            "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
            "OMP_NUM_THREADS": "1",
        }
        with (window / "engine.log").open("w") as log:
            result = subprocess.run(
                command,
                cwd=window,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        process.update(
            {
                "exit_status": result.returncode,
                "status": "finished" if result.returncode == 0 else "failed",
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        (window / "process.json").write_text(json.dumps(process, indent=2) + "\n")
        mdout = window / "minimized.mdout"
        output = mdout.read_text(errors="replace") if mdout.is_file() else ""
        converged = (
            result.returncode == 0
            and "Convergence reached" in output
            and "Optimisation converged" in output
        )
        if not converged:
            journal(
                "outcome",
                f"Run {window.name} DFTB3/MM path window",
                f"not accepted; rc={result.returncode}; DL-FIND convergence marker absent",
            )
            (window / "acceptance.json").write_text(
                json.dumps(
                    {
                        "accepted": False,
                        "reason": "process or force convergence gate failed",
                        "exit_status": result.returncode,
                    },
                    indent=2,
                )
                + "\n"
            )
            raise RuntimeError(f"{window.name} failed; outputs are preserved")
        system = parmed.load_file(str(TOPOLOGY), xyz=str(window / "optimized.ncrst"))
        values = geometry(np.asarray(system.coordinates, dtype=float))
        energy = energy_record(mdout)
        previous_values = records[-1]["geometry"]
        changes = adjacent_changes(previous_values, values)
        continuity = continuity_passes(changes)
        gates = interim_gates(values, thr_h3_present)
        final_gates = endpoint_gates(values, thr_h3_present)
        accepted = all(gates.values()) and continuity
        row = {
            "window": window.name,
            "geometry": values,
            "energy": energy,
            "relative_physical_kcal_mol": energy["physical_kcal_mol"]
            - initial_energy["physical_kcal_mol"],
            "relative_EPtot_bias_inclusive_kcal_mol": energy[
                "EPtot_bias_inclusive_kcal_mol"
            ]
            - initial_energy["EPtot_bias_inclusive_kcal_mol"],
            "energy_change_from_previous_kcal_mol": energy["physical_kcal_mol"]
            - records[-1]["energy"]["physical_kcal_mol"],
            "geometry_change_from_previous": changes,
            "continuity_keys": list(REACTIVE_CONTINUITY_KEYS),
            "continuity_limit_A": CONTINUITY_LIMIT_A,
            "continuity": continuity,
            "gates": gates,
            "accepted": accepted,
            "endpoint_gates": final_gates,
            "endpoint_chemistry_accepted": all(final_gates.values()),
        }
        row["restart_sha256"] = sha256(window / "optimized.ncrst")
        row["output_sha256"] = {
            name: sha256(window / name)
            for name in (
                "minimized.mdout",
                "optimized.ncrst",
                "optimization.nc",
                "mdinfo",
                "restraint_values.dat",
            )
            if (window / name).is_file()
        }
        (window / "acceptance.json").write_text(json.dumps(row, indent=2) + "\n")
        journal(
            "outcome",
            f"Run {window.name} DFTB3/MM path window",
            json.dumps(
                {
                    "accepted_for_continuation": accepted,
                    "geometry": values,
                    "energy_kcal_mol": energy,
                },
                separators=(",", ":"),
            ),
        )
        if not accepted:
            raise RuntimeError(
                f"{window.name} failed proton-inventory or continuity gate; outputs are preserved"
            )
        records.append(row)
        write_ledgers(records)
        export_progress(records)
        protocol["last_accepted_window"] = window.name
        protocol["status"] = "path_running"
        protocol_path.write_text(json.dumps(protocol, indent=2) + "\n")
        if row["endpoint_chemistry_accepted"]:
            protocol["status"] = "endpoint_accepted"
            protocol["endpoint_chemistry"] = row["endpoint_gates"]
            protocol["last_accepted_window"] = window.name
            protocol_path.write_text(json.dumps(protocol, indent=2) + "\n")
            break
    else:
        protocol["status"] = "endpoint_negative_or_incomplete"
        protocol["endpoint_chemistry"] = records[-1].get("endpoint_gates")
        protocol["last_accepted_window"] = records[-1]["window"]
        protocol_path.write_text(json.dumps(protocol, indent=2) + "\n")
        return protocol
    return run_release(protocol)


def main() -> None:
    """Run no-launch preflight or execute the path inside an approved tmux session."""
    import argparse

    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--execute", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume only from validated window_005 outputs under protocol revision 5",
    )
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.preflight:
        LOGGER.info(
            "No-launch preflight: %s",
            json.dumps(preflight(resume=args.resume), sort_keys=True),
        )
        return
    if not os.environ.get("TMUX"):
        raise RuntimeError(
            "Scientific path execution requires an existing tmux session"
        )
    result = execute(resume=args.resume)
    LOGGER.info("Path finished with status %s", result["status"])


if __name__ == "__main__":
    main()
