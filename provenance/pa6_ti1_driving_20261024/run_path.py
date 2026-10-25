"""Sequential, fail-closed forward/reverse TI1 construction pilot."""

from __future__ import annotations
import argparse
import hashlib
import json
import logging
from pathlib import Path
import re
import shutil
import subprocess
from datetime import datetime, timezone
import numpy as np
from parmed.amber import AmberParm, NetCDFTraj

ROOT = Path("/data/jpc/covalent-enzyme-design")
PILOT = ROOT / "output_nyl12_pa6_dftb_pilot_20261024"
BASE = PILOT / "analytic_relaxation_continuation_attempt1"
CONFIG = ROOT / "configs/pa6_ti1_coordinate_driving_20261024.json"


def digest(path: Path) -> str:
    """Hash one retained input."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path: Path, value: object) -> None:
    """Write current status alongside immutable attempts."""
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def targets(fraction: float) -> tuple[float, float]:
    """Return separately controlled attack and proton targets in Angstrom."""
    return 3.141038 + fraction * (1.55 - 3.141038), -1.584863 + fraction * (
        0.5 + 1.584863
    )


def prepare(output: Path) -> None:
    """Freeze the reviewed baseline and source before execution."""
    output.mkdir(exist_ok=False)
    source = output / "source"
    source.mkdir()
    for path in [
        Path(__file__),
        CONFIG,
        BASE / "preparation.json",
        BASE / "continuation.mdin",
    ]:
        shutil.copy2(path, source / path.name)
    metadata = json.loads((BASE / "preparation.json").read_text())
    command = metadata["command"]
    topology = Path(command[command.index("-p") + 1])
    reference = Path(command[command.index("-ref") + 1])
    initial = BASE / "optimized.ncrst"
    parm = AmberParm(str(topology), xyz=str(initial))
    expected = {
        13905: "N",
        13916: "OG1",
        13917: "HG1",
        19903: "C54",
        19904: "O49",
        19905: "N59",
    }
    if len(parm.atoms) != 103251 or any(
        parm.atoms[i - 1].name != name for i, name in expected.items()
    ):
        raise ValueError("Reactive atom identities or full-system inventory differ")
    write(
        output / "inputs.json",
        {
            "command": command,
            "environment": metadata["environment"],
            "topology": str(topology),
            "reference": str(reference),
            "initial": str(initial),
            "box": parm.box.tolist(),
            "sha256": {
                str(p): digest(p)
                for p in [topology, reference, initial, CONFIG, Path(__file__)]
            },
        },
    )
    write(
        output / "status.json",
        {"state": "prepared", "started_utc": None, "completed_windows": []},
    )


def one_window(
    output: Path, label: str, fraction: float, restart: Path, inputs: dict
) -> Path:
    """Minimize one static two-coordinate window and retain its evidence."""
    import os

    folder = output / label
    folder.mkdir(exist_ok=False)
    attack, proton = targets(fraction)
    rst = folder / "restraints.RST"
    blocks = []
    for atoms, weights, target in [
        ("13916,19903,0,0,0,0,0,0", "1,0,0,0", attack),
        ("13916,13917,13905,13917,0,0,0,0", "1,-1,0,0", proton),
    ]:
        blocks.append(
            f"&rst iat={atoms}, rstwt={weights}, r1={target - 2:.9f}, r2={target:.9f}, r3={target:.9f}, r4={target + 2:.9f}, rk2=25, rk3=25, /\n"
        )
    rst.write_text("".join(blocks))
    text = (output / "source/continuation.mdin").read_text()
    for key, value in {"maxcyc": 1000, "ntpr": 10, "ntwx": 10, "ntwr": 50}.items():
        text, count = re.subn(
            rf"(\b{key}\s*=\s*)[0-9]+", lambda m: m[1] + str(value), text, count=1
        )
        if count != 1:
            raise ValueError("Missing minimization control " + key)
    if "nmropt" in text.lower():
        raise ValueError("Baseline already contains reaction restraints")
    text = text.replace("&cntrl", "&cntrl\n nmropt=1,", 1)
    text += "\n&wt type='END' /\nDISANG=restraints.RST\nDUMPAVE=restraint_values.dat\n"
    mdin = folder / "minimize.mdin"
    mdin.write_text(text)
    command = inputs["command"].copy()
    for flag, path in {
        "-i": mdin,
        "-c": restart,
        "-ref": Path(inputs["reference"]),
        "-o": folder / "minimized.mdout",
        "-r": folder / "optimized.ncrst",
        "-x": folder / "optimization.nc",
        "-inf": folder / "mdinfo",
    }.items():
        command[command.index(flag) + 1] = str(path)
    record = {
        "command": command,
        "environment": inputs["environment"],
        "target_attack_A": attack,
        "target_proton_A": proton,
        "input_hashes": {
            str(p): digest(p) for p in [mdin, rst, restart, Path(inputs["reference"])]
        },
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    write(folder / "process.json", record)
    environment = os.environ.copy()
    environment.update(inputs["environment"])
    logging.info("Starting %s: attack %.4f, proton %.4f", label, attack, proton)
    with (folder / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=folder,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    record.update(
        exit_status=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    write(folder / "process.json", record)
    if result.returncode:
        raise RuntimeError("Engine failed " + label)
    rows = re.findall(
        r"^\s*(\d+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+(\S+)\s+(\d+)\s*$",
        (folder / "minimized.mdout").read_text(),
        re.M,
    )
    if not rows or float(rows[-1][2]) > 0.05:
        raise ValueError("Window did not reach RMS target: " + label)
    final = folder / "optimized.ncrst"
    parm = AmberParm(inputs["topology"], xyz=str(final))
    if not np.isfinite(parm.coordinates).all() or not np.allclose(
        parm.box, inputs["box"], rtol=0, atol=1e-10
    ):
        raise ValueError("Invalid coordinates or changed box")
    xyz = parm.coordinates
    pairs = {
        "attack": (13916, 19903),
        "OG_H": (13916, 13917),
        "N_H": (13905, 13917),
        "C_O": (19903, 19904),
        "C_N": (19903, 19905),
    }
    distances = {
        name: float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))
        for name, (a, b) in pairs.items()
    }
    write(
        folder / "acceptance.json",
        {
            "RMS": float(rows[-1][2]),
            "distances_A": distances,
            "converged": True,
            "chemical_validation": False,
            "restart_sha256": digest(final),
            "mdout_sha256": digest(folder / "minimized.mdout"),
        },
    )
    trajectory = NetCDFTraj.open_old(str(folder / "optimization.nc"))
    saved_frames = len(trajectory.coordinates)
    trajectory.close()
    write(
        folder / "display_frames.json",
        {
            "saved_optimizer_frames": saved_frames,
            "exact_endpoint_added": True,
            "zero_frame_policy": "A window converging before ntwx has only its exact endpoint",
        },
    )
    trajectory_input = f"trajin {folder / 'optimization.nc'}\n" if saved_frames else ""
    box = inputs["box"]
    cpp = folder / "export.cpptraj"
    cpp.write_text(
        f"parm {inputs['topology']}\n{trajectory_input}trajin {final}\nbox x {box[0]} y {box[1]} z {box[2]} alpha {box[3]} beta {box[4]} gamma {box[5]}\nstrip !:1-1313,3328,10325,20216\ntrajout {folder / 'optimization_stripped.dcd'} dcd namdcell\nrun\n"
    )
    with (folder / "export.log").open("w") as handle:
        subprocess.run(
            [
                str(Path(inputs["environment"]["AMBERHOME"]) / "bin/cpptraj"),
                "-i",
                str(cpp),
            ],
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return final


def run(output: Path) -> None:
    """Run one forward path and conditionally one reverse path, sequentially."""
    inputs = json.loads((output / "inputs.json").read_text())
    state = json.loads((output / "status.json").read_text())
    if state["state"] != "prepared":
        raise ValueError("Never rerun an existing attempt")
    for path, expected in inputs["sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError("Source changed before execution: " + path)
    state.update(state="running", started_utc=datetime.now(timezone.utc).isoformat())
    write(output / "status.json", state)
    restart = Path(inputs["initial"])
    try:
        for direction, indices in [
            ("forward", range(18)),
            ("reverse", range(16, -1, -1)),
        ]:
            if direction == "reverse":
                d = json.loads((output / "forward_17/acceptance.json").read_text())[
                    "distances_A"
                ]
                if not (
                    1.3 < d["attack"] < 1.9
                    and 1.2 < d["C_O"] < 1.6
                    and 1.2 < d["C_N"] < 1.7
                    and d["N_H"] <= 1.25
                    and d["OG_H"] >= 1.3
                ):
                    raise ValueError(
                        "Forward endpoint fails candidate TI1 geometry gate; reverse not launched"
                    )
            for index in indices:
                label = f"{direction}_{index:02d}"
                restart = one_window(output, label, index / 17, restart, inputs)
                state["completed_windows"].append(label)
                write(output / "status.json", state)
        state["state"] = "completed"
    except Exception as error:
        state.update(state="stopped", reason=str(error))
        raise
    finally:
        state["updated_utc"] = datetime.now(timezone.utc).isoformat()
        write(output / "status.json", state)


def main() -> None:
    """Prepare or launch an explicit fresh sequential path attempt."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.run:
        run(args.output.resolve())
    else:
        prepare(args.output.resolve())


if __name__ == "__main__":
    main()
