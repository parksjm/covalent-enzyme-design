"""One-off continuation of reviewed Nyl12 analytic relaxation attempt2."""

from __future__ import annotations
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask


def sha(path: Path) -> str:
    """Return the file digest without modification."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path: Path, data: object) -> None:
    """Preserve structured process or scientific evidence."""
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def derive_mdin(original: str) -> tuple[str, dict]:
    """Change only cycle budget and output cadences, with reversible proof."""
    changes = {"maxcyc": (25, 500), "ntpr": (1, 10), "ntwx": (1, 5), "ntwr": (1, 25)}
    derived = original
    for key, (before, after) in changes.items():
        pattern = rf"(\b{key}\s*=\s*){before}(?=\s*[,/\n])"
        derived, count = re.subn(
            pattern, lambda match: match[1] + str(after), derived, flags=re.I
        )
        if count != 1:
            raise ValueError(f"Expected exactly one original {key}={before}")
    reversed_text = derived
    for key, (before, after) in changes.items():
        pattern = rf"(\b{key}\s*=\s*){after}(?=\s*[,/\n])"
        reversed_text, count = re.subn(
            pattern, lambda match: match[1] + str(before), reversed_text, flags=re.I
        )
        if count != 1:
            raise ValueError("Cannot verify unique reversible namelist change")
    if reversed_text != original:
        raise ValueError("Unexpected Hamiltonian changes")
    return derived, {
        key: {"before": before, "after": after}
        for key, (before, after) in changes.items()
    }


def prepare(source: Path, output: Path) -> dict:
    """Prepare a fresh continuation with the original restraint reference."""
    source = source.resolve()
    metadata_path, process_path = source / "preparation.json", source / "process.json"
    metadata, previous = (
        json.loads(metadata_path.read_text()),
        json.loads(process_path.read_text()),
    )
    if previous["exit_status"] != 0:
        raise ValueError("Source attempt must have exited successfully")
    old_command = metadata["command"]
    topology = Path(old_command[old_command.index("-p") + 1])
    original_mdin = Path(old_command[old_command.index("-i") + 1])
    reference = Path(old_command[old_command.index("-ref") + 1])
    restart = source / "optimized.ncrst"
    if reference.resolve() != (source / "starting.ncrst").resolve():
        raise ValueError("Original attempt2 restraint reference must be retained")
    if metadata["restraint"]["force_constant_kcal_mol_A2"] != 10.0 or metadata[
        "restraint"
    ]["reference"] != str(reference):
        raise ValueError("Retain original10kcal/MM restraints and reference")
    if (
        metadata["region_audit"]["real_atom_count"] != 98
        or metadata["region_audit"]["qmcharge"] != -2
    ):
        raise ValueError("Expected primary98-atom charge-2 core")
    for record in metadata["sources"]:
        if sha(Path(record["path"])) != record["sha256"]:
            raise ValueError(f"Original source changed: {record['path']}")
    original_text = original_mdin.read_text()
    settings = dict(re.findall(r"(\w+)\s*=\s*('[^']*'|[^,\n/]+)", original_text))
    numeric = {
        "imin": 1,
        "ntmin": 1,
        "ncyc": 10,
        "dx0": 0.02,
        "drms": 0.05,
        "ntb": 1,
        "igb": 0,
        "ntc": 1,
        "ntf": 1,
        "cut": 10,
        "ifqnt": 1,
        "ntr": 1,
        "restraint_wt": 10,
        "qmcharge": -2,
        "adjust_q": 2,
        "qm_ewald": 1,
        "qm_pme": 1,
        "qmmm_int": 1,
        "qmcut": 10,
        "scfconv": 1e-10,
        "lnk_dis": 1.09,
        "lnk_atomic_no": 1,
        "lnk_method": 1,
    }
    if any(
        float(settings[key].replace("d", "e")) != value
        for key, value in numeric.items()
    ):
        raise ValueError("Original controls differ from approved case")
    if settings["qm_theory"].strip("'") != "DFTB3":
        raise ValueError("DFTB3 required")
    if any(name in settings for name in ["nmropt", "jar", "icfe", "ievb"]):
        raise ValueError("Unexpected reaction-coordinate/alternative-potential control")
    parm = AmberParm(str(topology), xyz=str(restart))
    ref = AmberParm(str(topology), xyz=str(reference))
    if (
        len(parm.atoms) != 103251
        or not np.isfinite(parm.coordinates).all()
        or not np.isfinite(ref.coordinates).all()
    ):
        raise ValueError("Expected finite full103251-atom endpoint and reference")
    if parm.get_coordinates().shape[0] != 1 or ref.get_coordinates().shape[0] != 1:
        raise ValueError("Exactly one endpoint and reference frame required")
    expected_box = np.asarray(metadata["source_box"])
    if not all(
        np.allclose(item.box, expected_box, rtol=0, atol=1e-10) for item in [parm, ref]
    ):
        raise ValueError("Endpoint/reference fixed box changed")
    qm = AmberMask(parm, settings["qmmask"].strip("'")).Selection()
    restrained = AmberMask(parm, settings["restraintmask"].strip("'")).Selection()
    if (
        sum(qm) != 98
        or sum(restrained) != 103153
        or any(q + r != 1 for q, r in zip(qm, restrained, strict=True))
    ):
        raise ValueError("MM restraints must complement98 free QM atoms")
    derived, changes = derive_mdin(original_text)
    output.mkdir(exist_ok=False)
    output = output.resolve()
    new_mdin = output / "continuation.mdin"
    new_mdin.write_text(derived)
    replacements = {
        "-i": new_mdin,
        "-c": restart,
        "-ref": reference,
        "-o": output / "minimized.mdout",
        "-r": output / "optimized.ncrst",
        "-x": output / "optimization.nc",
        "-inf": output / "relaxation.mdinfo",
    }
    command = old_command.copy()
    for flag, path in replacements.items():
        command[command.index(flag) + 1] = str(path)
    source_paths = list(
        dict.fromkeys(
            [
                *[Path(row["path"]) for row in metadata["sources"]],
                metadata_path,
                process_path,
                original_mdin,
                restart,
                reference,
                source / "minimized.mdout",
                Path(__file__).resolve(),
            ]
        )
    )
    frozen = output / "source"
    frozen.mkdir()
    frozen_records = []
    for index, path in enumerate(source_paths):
        target = frozen / f"{index:02d}_{path.name}"
        digest = sha(path)
        shutil.copy2(path, target)
        if sha(target) != digest:
            raise ValueError("Frozen source copy mismatch")
        frozen_records.append(
            {
                "path": str(path),
                "sha256": digest,
                "frozen_path": str(target),
                "frozen_sha256": digest,
            }
        )
    manifest = {
        "command": command,
        "environment": metadata["environment"],
        "sources": frozen_records,
        "derived_mdin_sha256": sha(new_mdin),
        "original_mdin_sha256": sha(original_mdin),
        "namelist_changes": changes,
        "source_attempt": str(source),
        "source_box": expected_box.tolist(),
        "start_coordinates": str(restart),
        "start_coordinates_sha256": sha(restart),
        "restraint_reference": str(reference),
        "restraint_reference_sha256": sha(reference),
        "atom_count": len(parm.atoms),
        "qm_real_atoms": sum(qm),
        "restrained_MM_atoms": sum(restrained),
        "scope": "Analytic-gradient continuation; failed FD unresolved; no MD/TS or new chemical model.",
        "trajectory_policy": "ntwx5 gives force-evaluation snapshots, not MD times or guaranteed cycle frames.",
        "endpoint_and_runner_review_required_before_execute": True,
    }
    write(output / "preparation.json", manifest)
    write(
        output / "scientific_acceptance.json",
        {
            "process_completed": False,
            "optimization_converged": False,
            "geometry_accepted": False,
            "status": "Prepared only; root endpoint and runner review required",
        },
    )
    return manifest


def execute(output: Path) -> None:
    """Execute only after root review of endpoint chemistry and this runner."""
    output = output.resolve()
    metadata = json.loads((output / "preparation.json").read_text())
    if (output / "process.json").exists():
        raise ValueError("Never rerun an existing engine attempt")
    before = [
        {"path": row["path"], "sha256": sha(Path(row["path"]))}
        for row in metadata["sources"]
    ]
    if any(
        item["sha256"] != record["sha256"]
        for item, record in zip(before, metadata["sources"], strict=True)
    ):
        raise ValueError("A source changed after preparation")
    if sha(output / "continuation.mdin") != metadata["derived_mdin_sha256"]:
        raise ValueError("Derived minimization input changed after review")
    write(output / "source_hashes_before.json", before)
    status = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "exit_status": None,
        "command": metadata["command"],
        "environment": metadata["environment"],
    }
    write(output / "process.json", status)
    environment = os.environ.copy()
    environment.update(metadata["environment"])
    with (output / "engine.log").open("w") as log:
        process = subprocess.run(
            metadata["command"],
            cwd=output,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    status.update(
        exit_status=process.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    write(output / "process.json", status)
    after = [
        {"path": row["path"], "sha256": sha(Path(row["path"]))}
        for row in metadata["sources"]
    ]
    write(output / "source_hashes_after.json", after)
    if before != after:
        raise ValueError("Original sources changed during continuation")
    if process.returncode:
        raise RuntimeError(f"Sander continuation exited{process.returncode}")
    command = metadata["command"]
    topology = command[command.index("-p") + 1]
    final = AmberParm(topology, xyz=str(output / "optimized.ncrst"))
    box_ok = np.allclose(final.box, metadata["source_box"], rtol=0, atol=1e-10)
    trajectory = AmberParm(topology, xyz=str(output / "optimization.nc"))
    frames = trajectory.get_coordinates()
    if (
        len(final.atoms) != 103251
        or not np.isfinite(final.coordinates).all()
        or not box_ok
        or frames.shape[1:] != (103251, 3)
        or not np.isfinite(frames).all()
    ):
        raise ValueError(
            "Invalid final coordinates, fixed box or optimization snapshots"
        )
    rows = re.findall(
        r"^\s*(\d+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+(\S+)\s+(\d+)\s*$",
        (output / "minimized.mdout").read_text(),
        re.M,
    )
    if not rows:
        raise ValueError("No final minimization gradient table found")
    last = rows[-1]
    write(
        output / "scientific_acceptance.json",
        {
            "process_completed": True,
            "optimization_converged": float(last[2]) <= 0.05,
            "final_reported_cycle": int(last[0]),
            "final_energy": float(last[1]),
            "final_RMS_gradient": float(last[2]),
            "target_RMS_gradient": 0.05,
            "fixed_box_preserved": bool(box_ok),
            "trajectory_force_evaluation_frames": len(frames),
            "geometry_accepted": False,
            "geometry_status": "Separate endpoint/proton-state review required",
            "convergence_scope": "ReportedRMS only; no minimum, ES/TI/TS state or chemical-model validation established",
        },
    )
    write(
        output / "output_hashes.json",
        [
            {"path": str(path), "sha256": sha(path)}
            for path in [
                output / name
                for name in [
                    "minimized.mdout",
                    "optimized.ncrst",
                    "optimization.nc",
                    "engine.log",
                    "process.json",
                    "scientific_acceptance.json",
                ]
            ]
        ],
    )


def main() -> None:
    """Prepare first; root-authorized execute is a separate explicit action."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--source-attempt", type=Path, required=True)
    prep.add_argument("--output", type=Path, required=True)
    run = sub.add_parser("execute")
    run.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.source_attempt, args.output)
    else:
        execute(args.output)


if __name__ == "__main__":
    main()
