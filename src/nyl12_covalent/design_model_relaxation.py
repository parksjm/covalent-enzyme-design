"""Prepare and minimize one four-ligand PA6 design model with Amber."""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import parmed

from .design_pa6_mapping import convert_four_ligands, export_canonical_heavy
from .design_relaxation_topology import leap_input, validate_four_ligand_topology

LOGGER = logging.getLogger(__name__)
ION_PAIR_FACTOR = 0.150 * 6.02214076e-4
DRMS_TARGET = 0.05


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _paths(config: dict[str, Any]) -> dict[str, Path]:
    names = ("release", "mol2", "mapping", "tleap", "engine", "amberhome")
    missing = set((*names, "gpu")) - config.keys()
    if missing:
        raise ValueError(f"Missing config keys: {sorted(missing)}")
    result = {name: Path(config[name]) for name in names}
    for name, path in result.items():
        if not path.is_absolute():
            raise ValueError(f"{name} must be absolute: {path}")
        if not path.exists() or (name in {"tleap", "engine"} and not path.is_file()):
            raise FileNotFoundError(path)
    if not str(config["gpu"]).strip():
        raise ValueError("gpu must be an explicit device selector")
    return result


def _run(
    argv: list[str],
    cwd: Path,
    env: dict[str, str],
    log: Path,
    command_log: Path,
    inputs: list[Path],
    events: list[dict[str, Any]],
) -> subprocess.CompletedProcess[bytes]:
    with command_log.open("a") as handle:
        handle.write(
            json.dumps(
                {
                    "event": "planned",
                    "argv": argv,
                    "cwd": str(cwd),
                    "input_sha256": {str(path): _hash(path) for path in inputs},
                },
                sort_keys=True,
            )
            + "\n"
        )
    with log.open("w") as handle:
        result = subprocess.run(
            argv, cwd=cwd, env=env, stdout=handle, stderr=subprocess.STDOUT, check=False
        )
    events.append({"argv": argv, "returncode": result.returncode})
    return result


def _environment(paths: dict[str, Path], gpu: str) -> dict[str, str]:
    env = os.environ.copy()
    env["AMBERHOME"] = str(paths["amberhome"])
    libs = (paths["amberhome"] / "lib", paths["engine"].parent.parent / "lib")
    env["LD_LIBRARY_PATH"] = ":".join(
        [
            *(str(path) for path in libs),
            *([env["LD_LIBRARY_PATH"]] if env.get("LD_LIBRARY_PATH") else []),
        ]
    )
    env.update(CUDA_VISIBLE_DEVICES=gpu, OMP_NUM_THREADS="1")
    return env


def _leap_ok(result: subprocess.CompletedProcess[bytes], log: Path) -> None:
    text = log.read_text(errors="replace")
    if result.returncode:
        raise RuntimeError(f"LEaP exited {result.returncode}: {log}")
    if "Errors = 0" not in text or "FATAL" in text.upper():
        raise ValueError(f"LEaP log is not clean: {log}")


def _volume(topology: Path, coordinates: Path) -> float:
    box = parmed.load_file(str(topology), xyz=str(coordinates)).box
    if box is None or len(box) < 6:
        raise ValueError("No six-value periodic box found")
    a, b, c, alpha, beta, gamma = map(float, box[:6])
    if not all(map(math.isfinite, (a, b, c, alpha, beta, gamma))):
        raise ValueError("Periodic box contains non-finite values")
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    factor = 1 + 2 * ca * cb * cg - ca**2 - cb**2 - cg**2
    if min(a, b, c) <= 0 or factor <= 0:
        raise ValueError("Periodic box dimensions or angles are invalid")
    return a * b * c * math.sqrt(factor)


def _mdin(atom_count: int) -> str:
    if atom_count <= 0:
        raise ValueError("Atom count must be positive")
    return (
        "PA6 unrestrained MM minimization\n&cntrl\n"
        f" imin=1, maxcyc=10000, ncyc=1000, drms={DRMS_TARGET},\n"
        " ntb=1, cut=10.0, ntr=0, nmropt=0, ntc=1, ntf=1,\n"
        f" ntpr=100, ntwx=100, ntwr=1000, ioutfm=1, ntxo=2, ntwprt={atom_count},\n/\n"
    )


def _final_results(text: str) -> tuple[float, float]:
    if "FINAL RESULTS" not in text:
        raise ValueError("Amber output is missing FINAL RESULTS")
    rows = re.findall(
        r"^\s*\d+\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+[-+0-9.EeDd]+\s+\S+\s+\d+",
        text,
        re.MULTILINE,
    )
    if not rows:
        raise ValueError("Cannot parse final energy and RMS gradient")
    energy, rms = (float(x.replace("D", "E").replace("d", "e")) for x in rows[-1])
    if not math.isfinite(energy) or not math.isfinite(rms):
        raise ValueError("Final energy or RMS gradient is non-finite")
    return energy, rms


def _run_leap(
    system: Path,
    release: Path,
    ions: tuple[int, int] | None,
    paths: dict[str, Path],
    env: dict[str, str],
    base: list[Path],
    command_log: Path,
    events: list[dict[str, Any]],
    stem: str,
) -> None:
    input_path, log = system / f"{stem}.in", system / f"{stem}.log"
    input_path.write_text(leap_input(release, ions))
    result = _run(
        [str(paths["tleap"]), "-I", str(release), "-f", input_path.name],
        system,
        env,
        log,
        command_log,
        [input_path, *base, paths["tleap"]],
        events,
    )
    _leap_ok(result, log)


def _process_record(
    output: Path,
    status: str,
    stage: str,
    events: list[dict[str, Any]],
    error: str | None = None,
) -> None:
    data: dict[str, Any] = {"status": status, "stage": stage, "events": events}
    data.update({"error": error} if error else {})
    (output / "process_acceptance.json").write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n"
    )


def relax_model(source: Path, output: Path, config: dict[str, Any]) -> dict[str, Any]:
    """Prepare and minimize a four-ligand PA6 model.

    Args:
        source: Canonical eight-protein/four-ligand PDB.
        output: New artifact directory.
        config: Required absolute Amber, release and mapping paths plus ``gpu``.

    Returns:
        Scientific acceptance; process status is written separately.

    Raises:
        FileExistsError: If output exists. RuntimeError or ValueError on failed
            preparation or execution; process failures are recorded.
    """
    if output.exists():
        raise FileExistsError(output)
    if not source.is_absolute() or not source.is_file():
        raise FileNotFoundError(source)
    paths = _paths(config)
    output.mkdir(parents=True, exist_ok=False)
    events: list[dict[str, Any]] = []
    stage = "input_preservation"
    _process_record(output, "running", stage, events)
    try:
        inputs = output / "inputs"
        inputs.mkdir()
        copies = (
            ("source_model.pdb", source),
            ("whole_l4_source.mol2", paths["mol2"]),
            ("whole_l4_context_mapping.json", paths["mapping"]),
        )
        for name, path in copies:
            shutil.copy2(path, inputs / name)
        release = inputs / "parameters"
        shutil.copytree(paths["release"], release)
        system = output / "system"
        convert_four_ligands(
            inputs / copies[0][0],
            inputs / copies[1][0],
            inputs / copies[2][0],
            release,
            system,
        )
        conversion = system / "conversion.json"
        base = [
            *(inputs / name for name, _ in copies),
            conversion,
            system / "complex_heavy_modular.pdb",
            *(path for path in release.iterdir() if path.is_file()),
        ]
        env = _environment(paths, str(config["gpu"]))
        command_log = output / "command_log.jsonl"

        stage = "dry_leap"
        _run_leap(system, release, None, paths, env, base, command_log, events, "count")
        dry = validate_four_ligand_topology(system / "dry.prmtop", release, False)
        volume = _volume(system / "noions.prmtop", system / "noions.rst7")
        pairs = round(ION_PAIR_FACTOR * volume)
        sodium, chloride = pairs + 8, pairs
        stage = "salted_leap"
        _run_leap(
            system,
            release,
            (sodium, chloride),
            paths,
            env,
            base,
            command_log,
            events,
            "final",
        )
        solvated = validate_four_ligand_topology(
            system / "complex.prmtop", release, True
        )

        stage = "minimization"
        minimization = output / "minimization"
        minimization.mkdir()
        mdin = minimization / "min.in"
        mdin.write_text(_mdin(int(solvated["atom_count"])))
        out, rst = minimization / "min.out", minimization / "min.rst7"
        argv = [
            str(paths["engine"]),
            "-O",
            "-i",
            str(mdin),
            "-p",
            str(system / "complex.prmtop"),
            "-c",
            str(system / "complex.rst7"),
            "-o",
            str(out),
            "-r",
            str(rst),
            "-x",
            str(minimization / "min.nc"),
            "-inf",
            str(minimization / "min.mdinfo"),
        ]
        result = _run(
            argv,
            minimization,
            env,
            minimization / "engine.log",
            command_log,
            [mdin, system / "complex.prmtop", system / "complex.rst7", paths["engine"]],
            events,
        )
        _process_record(
            output, "completed" if result.returncode == 0 else "failed", stage, events
        )
        if result.returncode:
            raise RuntimeError(f"Amber minimization exited {result.returncode}")
        energy, rms = _final_results(out.read_text(errors="replace"))
        endpoint = parmed.load_file(str(system / "complex.prmtop"), xyz=str(rst))
        endpoint.strip(":WAT,Na+,Cl-")
        top, coords = (
            minimization / "endpoint_complex.prmtop",
            minimization / "endpoint_complex.rst7",
        )
        endpoint.save(str(top), overwrite=False)
        endpoint.save(str(coords), format="ncrst", overwrite=False)
        pdb = output / "relaxed.pdb"
        export = export_canonical_heavy(top, coords, conversion, pdb)
        acceptance = {
            "status": "minimization_complete",
            "source": str(inputs / copies[0][0]),
            "endpoint": str(pdb),
            "parameter_release_id": dry["release_id"],
            "dry_topology": dry,
            "solvated_topology": solvated,
            "salt": dict(
                concentration_molar=0.150,
                box_volume_angstrom3=volume,
                neutralizing_plus_salt_pairs=pairs,
                sodium_count=sodium,
                chloride_count=chloride,
            ),
            "method": dict(
                protein_force_field="ff19SB",
                ligand_force_field="released modular PA6 RESP/GAFF2",
                water_model="TIP3P",
                solvent_padding_angstrom=12.0,
                protonation="Neutral beta Thr; ASP272/ASP274 anionic at all sites.",
                restraints="none; ntr=0, nmropt=0, ntc=1, ntf=1",
                max_cycles=10000,
                drms_target_kcal_mol_angstrom=DRMS_TARGET,
            ),
            "final_energy_kcal_mol": energy,
            "final_rms_gradient": rms,
            "force_converged": rms <= DRMS_TARGET,
            "process_exit_code": result.returncode,
            "heavy_export": export,
            "interpretation": "Process completion and force convergence are independent.",
        }
        (output / "acceptance.json").write_text(
            json.dumps(acceptance, indent=2, sort_keys=True) + "\n"
        )
        _process_record(output, "completed", stage, events)
        LOGGER.info("Finished PA6 design minimization; RMS gradient %.6g", rms)
        return acceptance
    except Exception as exc:
        _process_record(output, "failed", stage, events, repr(exc))
        raise
