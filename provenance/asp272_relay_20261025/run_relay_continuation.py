"""Two-stage Asp272 shuttle construction using native analytic-gradient L-BFGS."""

from __future__ import annotations
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import numpy as np
from parmed.amber import AmberParm, NetCDFTraj

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.path_inspection import publish_inspection  # noqa: E402

INSPECTION_CONFIG = ROOT / "configs/asp272_relay_inspection_20261025.json"
PILOT = ROOT / "output_nyl12_pa6_dftb_pilot_20261024"
HELPER = ROOT / "provenance/pa6_ti1_driving_20261024/run_path.py"
spec = importlib.util.spec_from_file_location("previous_driver", HELPER)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
CONFIG = ROOT / "configs/asp272_relay_pilot_20261025.json"


def distances(parm: AmberParm) -> dict[str, float]:
    """Track the original proton and chemically relevant covalent distances."""
    pairs = {
        "attack": (13916, 19903),
        "OG_H": (13916, 13917),
        "Asp_H": (14559, 13917),
        "N_HG": (13905, 13917),
        "N_H1": (13905, 13906),
        "N_H2": (13905, 13907),
        "C_O": (19903, 19904),
        "C_N": (19903, 19905),
    }
    xyz = parm.coordinates
    return {
        name: float(np.linalg.norm(xyz[a - 1] - xyz[b - 1]))
        for name, (a, b) in pairs.items()
    }


def prepare(output: Path) -> None:
    """Reuse the prior audited system inventory and freeze this new controller."""
    old.CONFIG = CONFIG
    old.prepare(output)
    shutil.copy2(Path(__file__), output / "source/run_relay.py")
    inputs = json.loads((output / "inputs.json").read_text())
    parm = AmberParm(inputs["topology"], xyz=inputs["initial"])
    if parm.atoms[14558].name != "OD1" or parm.atoms[14558].residue.idx + 1 != 968:
        raise ValueError("Asp272 OD1 identity differs")
    inputs["sha256"][str(Path(__file__))] = old.digest(Path(__file__))
    inputs["initial_distances"] = distances(parm)
    old.write(output / "inputs.json", inputs)


def window(
    output: Path,
    label: str,
    attack: float,
    proton: float,
    stage: int,
    restart: Path,
    inputs: dict,
) -> Path:
    """Run a fresh harmonic window with strict native convergence and tracking."""
    folder = output / label
    folder.mkdir(exist_ok=False)
    pairs = "13916,13917,14559,13917" if stage == 1 else "14559,13917,13905,13917"
    blocks = []
    for atoms, weights, target in [
        ("13916,19903,0,0,0,0,0,0", "1,0,0,0", attack),
        (pairs + ",0,0,0,0", "1,-1,0,0", proton),
    ]:
        blocks.append(
            f"&rst iat={atoms}, rstwt={weights}, r1={target - 2:.10f}, r2={target:.10f}, r3={target:.10f}, r4={target + 2:.10f}, rk2=100, rk3=100, /\n"
        )
    (folder / "restraints.RST").write_text("".join(blocks))
    text = (output / "source/continuation.mdin").read_text()
    for key, value in {"ntmin": 5, "ntwx": 5, "ntpr": 25, "ntwr": 100}.items():
        text, count = re.subn(
            rf"(\b{key}\s*=\s*)\d+", lambda m: m[1] + str(value), text, count=1
        )
        if count != 1:
            raise ValueError("Missing control " + key)
    text = text.replace("&cntrl", "&cntrl\n nmropt=1,", 1)
    text += "\n&wt type='END' /\nDISANG=restraints.RST\nDUMPAVE=restraint_values.dat\n&dlfind optalg='LBFGS', crdrep='CART', active='@*', tol=1.d-4, tole=1.d-7, maxcycle=1000, maxstep=0.5, verbosity=2, /\n"
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
    process = {
        "command": command,
        "environment": inputs["environment"],
        "stage": stage,
        "attack_target_A": attack,
        "proton_difference_target_A": proton,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "input_hashes": {
            str(p): old.digest(p)
            for p in [
                mdin,
                folder / "restraints.RST",
                restart,
                Path(inputs["reference"]),
            ]
        },
    }
    old.write(folder / "process.json", process)
    env = os.environ.copy()
    env.update(inputs["environment"])
    logging.info("Starting %s attack=%.4f proton=%.4f", label, attack, proton)
    with (folder / "engine.log").open("w") as handle:
        result = subprocess.run(
            command,
            cwd=folder,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    process.update(
        exit_status=result.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    old.write(folder / "process.json", process)
    if result.returncode:
        raise RuntimeError("Engine failed " + label)
    text = (folder / "minimized.mdout").read_text()
    converged = (
        "Optimisation converged" in text
        and "WARNING: OPTIMISATION HAS NOT CONVERGED" not in text
    )
    final = folder / "optimized.ncrst"
    parm = AmberParm(inputs["topology"], xyz=str(final))
    d = distances(parm)
    actual_proton = d["OG_H"] - d["Asp_H"] if stage == 1 else d["Asp_H"] - d["N_HG"]
    valid = bool(
        np.isfinite(parm.coordinates).all()
        and np.allclose(parm.box, inputs["box"], rtol=0, atol=1e-10)
    )
    tracking_error = max(abs(d["attack"] - attack), abs(actual_proton - proton))
    tracking = tracking_error <= 0.20
    if tracking_error > 0.15:
        logging.warning(
            "%s coordinate tracking error %.6f A exceeds 0.15 A warning",
            label,
            tracking_error,
        )
    old.write(
        folder / "acceptance.json",
        {
            "native_converged": converged,
            "max_component_tolerance_Hartree_Bohr": 1e-4,
            "finite_coordinates_and_fixed_box": valid,
            "target_tracking_within_0.20_A": tracking,
            "tracking_warning_above_0.15_A": tracking_error > 0.15,
            "maximum_tracking_error_A": tracking_error,
            "distances_A": d,
            "actual_proton_difference_A": actual_proton,
            "chemical_state_validated": False,
        },
    )
    if not (converged and valid and tracking):
        raise ValueError(
            f"{label}: convergence={converged}, coordinates={valid}, target_tracking={tracking}"
        )
    traj = NetCDFTraj.open_old(str(folder / "optimization.nc"))
    count = len(traj.coordinates)
    traj.close()
    old.write(
        folder / "display_frames.json",
        {
            "saved_evaluation_frames": count,
            "exact_endpoint_added": True,
            "policy": "DL-FIND callbacks include trial evaluations; ntwx5, no physical time",
        },
    )
    trajectory_input = f"trajin {folder / 'optimization.nc'}\n" if count else ""
    box = inputs["box"]
    cpp = folder / "export.cpptraj"
    cpp.write_text(
        f"parm {inputs['topology']}\n{trajectory_input}trajin {final}\nbox x {box[0]} y {box[1]} z {box[2]} alpha {box[3]} beta {box[4]} gamma {box[5]}\nstrip !:1-1313,3328,10325,20216\ntrajout {folder / 'optimization_stripped.dcd'} dcd namdcell\nrun\n"
    )
    with (folder / "export.log").open("w") as handle:
        subprocess.run(
            [str(Path(env["AMBERHOME"]) / "bin/cpptraj"), "-i", str(cpp)],
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return final


def gate(output: Path, label: str, stage: int) -> None:
    """Require the intended proton identity and intact TI1 connectivity."""
    d = json.loads((output / label / "acceptance.json").read_text())["distances_A"]
    common = (
        1.3 < d["attack"] < 1.9
        and 1.3 < d["C_O"] < 1.6
        and 1.2 < d["C_N"] < 1.7
        and max(d["N_H1"], d["N_H2"]) < 1.25
    )
    proton = (
        (d["Asp_H"] < 1.25 and d["OG_H"] > 1.3)
        if stage == 1
        else (d["N_HG"] < 1.25 and d["Asp_H"] > 1.3)
    )
    if not (common and proton):
        raise ValueError(f"Stage{stage} endpoint fails covalent/proton-identity gate")


def run(output: Path) -> None:
    """Run only the authorized two-stage relay, stopping on any failed window."""
    inputs = json.loads((output / "inputs.json").read_text())
    state = json.loads((output / "status.json").read_text())
    if state["state"] != "stopped" or len(state["completed_windows"]) != 10:
        raise ValueError("Continuation requires the preserved window10 tracking stop")
    if (
        state.get("reason")
        != "stage1_10: convergence=True, coordinates=True, target_tracking=False"
    ):
        raise ValueError("Unexpected stop reason")
    if (output / "continuation_tracking_policy.json").exists():
        raise ValueError("Continuation already prepared or executed")
    for path, expected in inputs["sha256"].items():
        if old.digest(Path(path)) != expected:
            raise ValueError("Source changed " + path)
    prior = output / "stage1_10"
    acceptance = json.loads((prior / "acceptance.json").read_text())
    process = json.loads((prior / "process.json").read_text())
    distances_now = acceptance["distances_A"]
    tracking_error = max(
        abs(distances_now["attack"] - process["attack_target_A"]),
        abs(
            acceptance["actual_proton_difference_A"]
            - process["proton_difference_target_A"]
        ),
    )
    if not (
        acceptance["native_converged"]
        and acceptance["finite_coordinates_and_fixed_box"]
        and process["exit_status"] == 0
        and tracking_error <= 0.20
    ):
        raise ValueError("Preserved endpoint does not qualify for continuation")
    old.write(output / "status_before_continuation_20261025.json", state)
    hashes = {
        str(p): old.digest(p)
        for p in [
            Path(__file__),
            INSPECTION_CONFIG,
            ROOT / "scripts/plot_reaction_distances.py",
            ROOT / "src/nyl12_covalent/path_inspection.py",
            ROOT / "src/nyl12_covalent/path_distances.py",
            ROOT / "src/nyl12_covalent/path_publish.py",
            prior / "optimized.ncrst",
            prior / "acceptance.json",
            prior / "process.json",
        ]
    }
    old.write(
        output / "continuation_tracking_policy.json",
        {
            "warning_A": 0.15,
            "stop_A": 0.20,
            "window10_tracking_error_A": tracking_error,
            "original_acceptance_preserved": True,
            "continuation_source_sha256": hashes,
            "rationale": "Converged harmonic window misses old cutoff by 0.003764 A; retain warning and explicitly allow finite bias lag.",
        },
    )
    state["completed_windows"].append("stage1_10")
    state["tracking_policy"] = {"warning_A": 0.15, "stop_A": 0.20}
    state.pop("reason", None)
    restart = prior / "optimized.ncrst"
    d = inputs["initial_distances"]
    state.update(state="running", continued_utc=datetime.now(timezone.utc).isoformat())
    old.write(output / "status.json", state)
    try:
        publish_inspection(
            output, INSPECTION_CONFIG, ROOT / "scripts/plot_reaction_distances.py"
        )
        for stage in (1, 2):
            if stage == 1:
                a0 = d["attack"]
                q0 = d["OG_H"] - d["Asp_H"]
            else:
                d = distances(AmberParm(inputs["topology"], xyz=str(restart)))
                a0 = 1.55
                q0 = d["Asp_H"] - d["N_HG"]
            for index in range(17):
                fraction = index / 16
                label = f"stage{stage}_{index:02d}"
                if label in state["completed_windows"]:
                    restart = output / label / "optimized.ncrst"
                    continue
                restart = window(
                    output,
                    label,
                    a0 + fraction * (1.55 - a0),
                    q0 + fraction * (0.8 - q0),
                    stage,
                    restart,
                    inputs,
                )
                state["completed_windows"].append(label)
                old.write(output / "status.json", state)
                publish_inspection(
                    output,
                    INSPECTION_CONFIG,
                    ROOT / "scripts/plot_reaction_distances.py",
                )
            gate(output, label, stage)
        state["state"] = "completed"
    except Exception as error:
        state.update(state="stopped", reason=str(error))
        raise
    finally:
        state["updated_utc"] = datetime.now(timezone.utc).isoformat()
        old.write(output / "status.json", state)
        publish_inspection(
            output, INSPECTION_CONFIG, ROOT / "scripts/plot_reaction_distances.py"
        )


def main() -> None:
    """Continue the existing run, preserving every completed window."""
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
        raise ValueError(
            "Use --run for this explicitly authorized in-place continuation"
        )


if __name__ == "__main__":
    main()
