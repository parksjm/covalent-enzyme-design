"""Separate, provenance-preserving DFTB3/MM repair of a crowded TI1 seed."""

from __future__ import annotations

import json
import logging
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from parmed.amber import AmberMask, AmberParm

from .local_repair_checks import audit_local_model, local_restraint_selection
from .local_repair_engine import run_checked_engine
from .local_repair_model import QM_MASK, export_pdb, prepare_model
from .local_water_premin import preminimize_water
from .ti1_attempt import _sha, _snapshot
from .ti1_checks import check_convergence

LOGGER = logging.getLogger(__name__)


def _write_json(path: Path, payload: object) -> None:
    with path.open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")


def write_stage_input(
    directory: Path,
    mask: str,
    pairs: list[list[int]],
    cycles: int,
    drms: float,
    bridge_force: float,
) -> Path:
    """Write one explicit restrained repair problem.

    Args:
        directory: Fresh existing stage directory.
        mask: Outer heavy atoms restrained to the original prepared reference.
        pairs: Two nitrogen/OG1 to water-O pairs in Amber atom numbering.
        cycles: Positive maximum minimization cycle budget.
        drms: Positive gradient tolerance.
        bridge_force: Nonnegative Amber distance restraint force constant.

    Returns:
        The Sander minimization input path.
    """
    if cycles < 1 or drms <= 0 or bridge_force < 0 or len(pairs) != 2:
        raise ValueError("Invalid repair stage specification")
    mdin = directory / "minimize.mdin"
    mdin.write_text(
        "PA6 L4 TI1 local DFTB3/MM repair; explicit QM water\n"
        " &cntrl\n"
        f" imin=1, ntmin=1, ncyc={min(20, cycles)}, maxcyc={cycles}, drms={drms},\n"
        " ntb=0, igb=8, cut=999.0, ntc=1, ntf=1, ntpr=10, ntwr=50, ntwx=1, ioutfm=0,\n"
        f" ntr=1, restraint_wt=5.0, restraintmask='{mask}',\n"
        " ifqnt=1, nmropt=1,\n /\n &qmmm\n"
        f" qmmask='{QM_MASK}', qmcharge=0, qm_theory='DFTB3',\n"
        " scfconv=1.0d-8, dftb_maxiter=200, itrmax=2000, verbosity=1,\n /\n"
        " &wt type='END' /\nDISANG=bridge.RST\n"
    )
    (directory / "bridge.RST").write_text(
        "".join(
            f" &rst iat={first},{second}, r1=2.2, r2=2.6, r3=3.2, r4=3.6, rk2={bridge_force}, rk3={bridge_force}, /\n"
            for first, second in pairs
        )
    )
    return mdin


def run_repair(
    topology: Path,
    coordinates: Path,
    cif: Path,
    amber_bin: Path,
    output: Path,
    prepare_only: bool = False,
) -> None:
    """Prepare and run three separately recorded local-repair stages in tmux.

    Args:
        topology: Product-gated source octamer topology.
        coordinates: Matching original topology coordinates.
        cif: Two-contact AF3 alpha/beta water model.
        amber_bin: Explicit AmberTools executable directory.
        output: Fresh calculation directory.
        prepare_only: Stop after preparation and inventory/mask verification.
    """
    executable = amber_bin / "sander"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(executable)
    output.mkdir(parents=True, exist_ok=False)
    inputs = {}
    for name, path in (
        ("source.prmtop", topology),
        ("source.rst7", coordinates),
        ("af3.cif", cif),
    ):
        shutil.copy2(path, output / name)
        inputs[name] = {"path": str(path.absolute()), "sha256": _sha(output / name)}
    source_hashes = _snapshot(output)
    slko = amber_bin.parent / "dat/slko"
    data = sorted(p for p in slko.rglob("*") if p.is_file())
    if not data:
        raise FileNotFoundError(f"DFTB parameters missing: {slko}")
    shutil.copytree(slko, output / "engine_data_snapshot")
    data_hashes = {str(p): _sha(p) for p in data}
    _write_json(
        output / "manifest.json",
        {
            "inputs": inputs,
            "source_sha256": source_hashes,
            "sander": str(executable.absolute()),
            "sander_sha256": _sha(executable),
            "engine_data_sha256": data_hashes,
            "system_scope": "single alpha/beta diagnostic; octamer validation pending",
            "positional_reference": "prepared.rst7 for every stage",
            "method": "DFTB3/MM, GB8, QM charge0; inherited ff19SB/GAFF/TIP3P parameters",
            "parameter_quality": "reactive classical parameters provisional; no MD/energy ranking validation",
        },
    )
    prepare_model(
        output / "source.prmtop", output / "source.rst7", output / "af3.cif", output
    )
    model = AmberParm(str(output / "model.prmtop"), xyz=str(output / "prepared.rst7"))
    initial = audit_local_model(model)
    _write_json(output / "initial_audit.json", initial)
    if not initial["bond_screen_passed"]:
        raise ValueError("Prepared reactive bond/proton geometry fails before repair")
    selection = local_restraint_selection(model)
    _write_json(output / "selection.json", selection)
    residues = ",".join(str(i) for i in selection["local_residues_1based"])
    masks = [":1-329&!@H=", f"!(:{residues})&!@H=", f"!(:{residues})&!@H="]
    expected = [
        selection["initial_heavy_atom_indices_1based"],
        selection["outer_heavy_atom_indices_1based"],
        selection["outer_heavy_atom_indices_1based"],
    ]
    for mask, indices in zip(masks, expected, strict=True):
        if set(AmberMask(model, mask).Selected()) != {i - 1 for i in indices}:
            raise ValueError("Positional mask does not match explicit atom selection")
    if prepare_only:
        LOGGER.info(
            "Preparation and chemistry/mask preflight complete; engine not started"
        )
        return
    environment = os.environ.copy()
    overrides = {
        "AMBERHOME": str(amber_bin.parent.absolute()),
        "OMP_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "1",
        "LD_LIBRARY_PATH": str(amber_bin.parent.absolute() / "lib"),
    }
    if environment.get("LD_LIBRARY_PATH"):
        overrides["LD_LIBRARY_PATH"] += ":" + environment["LD_LIBRARY_PATH"]
    environment.update(overrides)
    current = preminimize_water(
        output, executable, environment, selection["bridge_atoms_1based"]
    )
    stages = [
        ("01_hydrogens_water", 100, 0.5, 10.0),
        ("02_local_pocket", 300, 0.05, 5.0),
        ("03_bridge_release", 300, 0.05, 0.0),
    ]
    results = []
    for (name, cycles, drms, force), mask in zip(stages, masks, strict=True):
        stage = output / name
        stage.mkdir()
        shutil.copy2(current, stage / "start.rst7")
        mdin = write_stage_input(
            stage, mask, selection["bridge_atoms_1based"], cycles, drms, force
        )
        restart = stage / "optimized.rst7"
        mdout = stage / "minimize.mdout"
        command = [
            str(executable.absolute()),
            "-i",
            str(mdin),
            "-o",
            str(mdout),
            "-p",
            str(output / "model.prmtop"),
            "-c",
            str(stage / "start.rst7"),
            "-ref",
            str(output / "prepared.rst7"),
            "-r",
            str(restart),
            "-inf",
            str(stage / "minimize.mdinfo"),
            "-x",
            str(stage / "optimization.mdcrd"),
        ]
        _write_json(
            stage / "command.json",
            {
                "argv": command,
                "environment_overrides": overrides,
                "starting_coordinates_sha256": _sha(stage / "start.rst7"),
                "topology_sha256": _sha(output / "model.prmtop"),
                "reference_sha256": _sha(output / "prepared.rst7"),
                "input_sha256": _sha(mdin),
                "bridge_sha256": _sha(stage / "bridge.RST"),
                "drms": drms,
                "bridge_force": force,
            },
        )
        process = {
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "ended_utc": None,
            "exit_status": None,
        }
        LOGGER.info("Starting %s; %d cycles, bridge force %s", name, cycles, force)
        completed = run_checked_engine(
            command,
            cwd=stage,
            env=environment,
            log_path=stage / "sander.log",
            mdout=mdout,
        )
        process.update(
            ended_utc=datetime.now(timezone.utc).isoformat(),
            **completed,
        )
        if restart.is_file():
            process["restart_sha256"] = _sha(restart)
        _write_json(stage / "process.json", process)
        if (
            completed["exit_status"]
            or completed["stopped_on_scf_failure"]
            or not restart.is_file()
        ):
            _write_json(
                stage / "acceptance.json",
                {
                    "accepted_restrained_geometry": False,
                    "reason": "Engine failure, missing restart, or electronically unconverged forces",
                    "stopped_on_scf_failure": completed["stopped_on_scf_failure"],
                },
            )
            raise RuntimeError(f"Sander stage {name} failed; inspect retained logs")
        model = AmberParm(str(output / "model.prmtop"), xyz=str(restart))
        audit = audit_local_model(model)
        chirality_pass = all(
            initial["signed_volumes_A3"][key] * value > 0
            for key, value in audit["signed_volumes_A3"].items()
        )
        audit["chirality_retained"] = chirality_pass
        audit["convergence"] = check_convergence(mdout, drms)
        audit["accepted_restrained_geometry"] = (
            audit["geometry_screen_passed"]
            and chirality_pass
            and audit["convergence"]["accepted"]
        )
        _write_json(stage / "acceptance.json", audit)
        export_pdb(model, stage / "optimized.pdb")
        results.append(
            {
                "stage": name,
                "accepted_restrained_geometry": audit["accepted_restrained_geometry"],
                "bridge_force": force,
            }
        )
        _write_json(stage / "progress.json", {"stages_completed": results})
        LOGGER.info(
            "%s: bridge %s, maximum contact shortfall %.3f A",
            name,
            audit["bridge_distances_A"],
            audit["water_maximum_contact_shortfall_A"],
        )
        if (
            not audit["bond_screen_passed"]
            or not chirality_pass
            or audit["convergence"]["engine_failure_lines"]
        ):
            raise ValueError(
                f"Chemistry/engine screen failed after {name}; later stages withheld"
            )
        current = restart
    _write_json(
        output / "acceptance.json",
        {
            "stages": results,
            "accepted_local_geometry": results[-1]["accepted_restrained_geometry"],
            "stationary_point_validated": False,
            "octamer_compatibility_validated": False,
        },
    )
