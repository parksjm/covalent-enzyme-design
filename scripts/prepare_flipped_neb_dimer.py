"""Prepare native DL-FIND dimer inputs from the final flipped NEB band."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import shlex
import shutil
import sys
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from nyl12_covalent.neb_xyz_snapshot import reconstruct_xyz_band


def digest(path: Path) -> str:
    """Return a source-file SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def option(text: str, name: str) -> str:
    """Extract exactly one quoted namelist option."""
    values = re.findall(rf"\b{name}\s*=\s*'([^']*)'", text, re.I)
    if len(values) != 1:
        raise ValueError(f"Expected unique {name}")
    return values[0]


def main() -> None:
    """Prepare immutable TS and displaced-mode inputs without launching."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neb-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interface-source", type=Path, required=True)
    parser.add_argument("--session", default="model30_true_dimer_ts_20261113")
    args = parser.parse_args()
    neb, output = args.neb_dir.resolve(strict=True), args.output.resolve()
    project = Path(__file__).resolve().parents[1]
    source = args.interface_source.resolve(strict=True)
    interface = source.read_text()
    if not all(
        x in interface for x in ("options%dimer", "options%coords2(1)", "iopt = 3")
    ):
        raise ValueError("Native dimer interface not verified")
    calculation = neb / "calculation"
    mdin_path, mdout_path = (
        calculation / "minimize.mdin",
        calculation / "minimized.mdout",
    )
    xyz_path, top_path = calculation / "nebpath.xyz", neb / "inputs/system.prmtop"
    mdin, mdout = mdin_path.read_text(), mdout_path.read_text()
    active_mask, weight_mask = option(mdin, "active"), option(mdin, "wtmask")
    reference_path = Path(option(mdin, "refcrd"))
    sources = [
        mdin_path,
        mdout_path,
        xyz_path,
        top_path,
        reference_path,
        source,
        Path(__file__).resolve(),
        project / "src/nyl12_covalent/neb_xyz_snapshot.py",
    ]
    hashes = {str(path): digest(path) for path in sources}
    position = mdout.rindex("NEB Report\n")
    end = mdout.index("Total path length", position)
    cycles = re.findall(r"No climbing image used in cycle\s+(\d+)", mdout[:position])
    rows = [
        (int(i), float(e))
        for i, e in re.findall(
            r"^Img\s+(\d+)\s+([-+\d.Ee]+)", mdout[position:end], re.M
        )
    ]
    if (
        not cycles
        or int(cycles[-1]) != 499
        or [i for i, _ in rows] != list(range(1, 17))
    ):
        raise ValueError("Expected complete native cycle499 with16 regular images")
    tail = mdout[end:]
    closure = re.search(r"Testing convergence\s+in cycle\s+499\b", tail)
    if closure is None:
        raise ValueError("Final native save loop not closed at cycle499")
    saves = re.findall(r"Saving img\s+(\d+) to restart:", tail[: closure.start()])
    if [int(value) for value in saves] != list(range(1, 17)):
        raise ValueError("Final native band lacks ordered saves1 through16")
    peak = max(range(len(rows)), key=lambda index: rows[index][1])
    if peak != 7:
        raise ValueError("Expected native image8 maximum")
    topology, reference = AmberParm(str(top_path)), Rst7.open(str(reference_path))
    selected = np.array(
        [
            index
            for index, flag in enumerate(AmberMask(topology, active_mask).Selection())
            if flag and topology.atoms[index].mass > 0.1
        ],
        dtype=int,
    )
    frames, active = reconstruct_xyz_band(
        xyz_path.read_bytes(),
        topology,
        reference,
        active_mask,
        expected_images=16,
        expected_active=len(selected),
    )
    if list(active) != list(selected):
        raise ValueError("Active mapping differs from topology mask")
    weights = np.flatnonzero(AmberMask(topology, weight_mask).Selection())
    if not len(weights) or not set(weights).issubset(set(selected)):
        raise ValueError("Rotation mask must lie within active region")
    start = np.asarray(frames[peak], dtype=float)
    tangent = np.asarray(frames[peak + 1]) - np.asarray(frames[peak - 1])
    frozen = np.ones(len(topology.atoms), dtype=bool)
    frozen[selected] = False
    norm = float(np.linalg.norm(tangent[weights]))
    if not np.isfinite(norm) or norm <= 1e-10 or np.max(np.abs(tangent[frozen])) > 1e-7:
        raise ValueError("Invalid tangent or changed frozen atoms")
    mode = start + 0.01 * tangent / norm
    if hashes != {str(path): digest(path) for path in sources}:
        raise RuntimeError("Native source changed during capture")
    output.mkdir(parents=True, exist_ok=False)
    inputs, work = output / "inputs", output / "calculation"
    inputs.mkdir()
    work.mkdir()
    for index, path in enumerate(sources):
        shutil.copy2(path, inputs / f"source_{index:02d}_{path.name}")
    shutil.copy2(top_path, inputs / "system.prmtop")
    shutil.copy2(reference_path, inputs / "reference.ncrst")
    for name, coordinates in (("ts_guess", start), ("ts_mode_guess", mode)):
        restart = Rst7(natom=len(coordinates))
        restart.coordinates, restart.box = coordinates, reference.box
        restart.write(str(inputs / f"{name}.ncrst"), netcdf=True)
    dimer = (
        "&dlfind\n dimer=.true., dimer_mode=1, dimer_tolrot=2., dimer_maxrot=20,\n"
        " optalg='LBFGS', crdrep='CART',\n"
        f" refcrd='{inputs / 'reference.ncrst'}',\n"
        f" coords2(1)='{inputs / 'ts_mode_guess.ncrst'}',\n"
        f" wtmask='{weight_mask}', active='{active_mask}',\n"
        " tol=1.d-4, tole=1.d-5, maxcycle=1000, maxstep=0.1, verbosity=2,\n/\n"
    )
    (work / "dimer.mdin").write_text(mdin[: mdin.lower().index("&dlfind")] + dimer)
    old = json.loads((neb / "job/command.json").read_text())
    command = [
        old["command"][0],
        "-O",
        "-i",
        "dimer.mdin",
        "-p",
        str(inputs / "system.prmtop"),
        "-c",
        str(inputs / "ts_guess.ncrst"),
        "-o",
        "dimer.mdout",
        "-r",
        "optimized.ncrst",
        "-x",
        "optimization.nc",
        "-inf",
        "mdinfo",
    ]
    record = {
        "source_sha256": hashes,
        "native_cycle": 499,
        "native_image_1based": 8,
        "native_energy_hartree": rows[peak][1],
        "algorithm": "Native DL-FIND dimer mode1",
        "same_cntrl_qmmm_topology_protonation": True,
        "native_xyz_cycle_energy_match_verified": False,
        "native_cycle499_ordered_save_loop_verified": True,
        "tangent": "central difference of final native images7 and9; norm0.01A on wtmask",
        "scientific_status": "Prepared; requires source energy match then saddle/Hessian validation",
        "prepared_sha256": {
            str(p.relative_to(output)): digest(p)
            for p in [*inputs.iterdir(), work / "dimer.mdin"]
        },
    }
    (output / "preparation.json").write_text(json.dumps(record, indent=2) + "\n")
    launch = [
        sys.executable,
        str(project / "scripts/run_nyl12_job.py"),
        "launch",
        "--session",
        args.session,
        "--output-dir",
        str(output / "job"),
        "--cwd",
        str(work),
    ]
    for path in [*inputs.iterdir(), work / "dimer.mdin", output / "preparation.json"]:
        launch += ["--input", str(path)]
    launch += ["--input", command[0]]
    for key, value in old["env_overrides"].items():
        launch += ["--env", f"{key}={value}"]
    launch += ["--", *command]
    guard = (
        "import json,sys; "
        "assert json.load(open(sys.argv[1]))['native_xyz_cycle_energy_match_verified'], "
        "'Validate reconstructed image energy before launch'"
    )
    (output / "launch_after_validation.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        + shlex.join([sys.executable, "-c", guard, str(output / "preparation.json")])
        + "\n"
        + shlex.join(launch)
        + "\n"
    )
    logging.info("Prepared %s; native image-energy association still required", output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
