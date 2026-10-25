"""Source-checked preparation and execution of a bounded ES relaxation pilot."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm
from nyl12_covalent.es_qm_region import audit_es_qm_region
from nyl12_covalent.es_dftb_relaxation_gate import _force_gate, _sha


def _write(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _mdin(audit: dict) -> str:
    mask = audit["qmmask"]
    return (
        "Bounded ES pilot: no guarantee of convergence or retained ES\n"
        "&cntrl\n imin=1, maxcyc=25, ntmin=1, ncyc=10, dx0=0.02, drms=0.05,\n"
        " ntb=1, igb=0, ntc=1, ntf=1, cut=10.0, ifqnt=1, jfastw=4,\n"
        f" ntr=1, restraint_wt=10.0, restraintmask='!({mask})',\n"
        " ntpr=1, ntwx=1, ntwr=1, ioutfm=1, ntxo=2,\n/\n"
        "&qmmm\n"
        f" qmmask='{mask}', qmcharge=-2, qm_theory='DFTB3',\n"
        " qmmm_int=1, adjust_q=2, qm_ewald=1, qm_pme=1, qmcut=10.0,\n"
        " scfconv=1.0d-10, lnk_dis=1.09, lnk_atomic_no=1, lnk_method=1,\n/\n"
    )


def prepare_es_relaxation(
    engine_config: Path,
    force_check: Path,
    region_config: Path,
    sander: Path,
    output: Path,
    pilot_decision: Path | None = None,
) -> dict[str, Any]:
    """Prepare a source-pinned, force-gated 25-cycle restrained-MM ES pilot.

    Args:
        engine_config: Exact live configuration accepted by the enzyme FD check.
        force_check: Accepted enzyme force check and its source evidence.
        region_config: Approved explicit 98-atom primary region configuration.
        sander: Explicit installed Sander executable.
        output: Fresh directory; existing paths are rejected.
        pilot_decision: Explicit bounded analytic-pilot decision; absent means strict FD.

    Returns:
        Retained command, explicit environment and audited preparation metadata.
        Preparation alone never launches an engine.

    Raises:
        ValueError: Source identity, force acceptance, physics or restart fails.
        FileExistsError: The output directory already exists.
    """
    declared = json.loads(region_config.read_text())
    variants = [
        row for row in declared["variants"] if row["name"] == "primary_106_nuclei"
    ]
    if len(variants) != 1:
        raise ValueError("Exactly one approved primary ES core required")
    topology, coordinates = Path(declared["topology"]), Path(declared["coordinates"])
    audit = audit_es_qm_region(topology, coordinates, variants[0])
    if (
        audit["real_atom_count"],
        audit["link_H_count"],
        audit["electron_count_with_links"],
    ) != (98, 8, 368):
        raise ValueError("Relaxation requires the approved 98-real/8-link primary core")
    engine = _force_gate(
        engine_config, force_check, audit, region_config, pilot_decision
    )
    if (
        Path(engine["topology"]).resolve() != topology.resolve()
        or Path(engine["coordinates"]).resolve() != coordinates.resolve()
    ):
        raise ValueError("Enzyme force check and region audit use different sources")
    if not sander.is_file() or not os.access(sander, os.X_OK):
        raise ValueError("Explicit executable Sander prerequisite missing")
    if (
        sander.resolve()
        != (Path(engine["environment"]["AMBERHOME"]) / "bin/sander").resolve()
    ):
        raise ValueError("Sander must belong to the force-checked Amber installation")
    parm = AmberParm(str(topology), xyz=str(coordinates))
    xyz, box = parm.coordinates.copy(), np.asarray(parm.box, dtype=float)
    if (
        parm.get_coordinates().shape[0] != 1
        or box.shape != (6,)
        or not np.isfinite(box).all()
    ):
        raise ValueError("Exactly one source frame with a finite periodic box required")
    if engine["_baseline_force_atom_count"] != len(parm.atoms) or not np.array_equal(
        box, engine["_baseline_box"]
    ):
        raise ValueError(
            "Analytic baseline force size or fixed box does not match source frame"
        )
    selected = {row["index"] for row in audit["qm_atoms"]}
    output.mkdir(exist_ok=False)
    output = output.resolve()
    restart = output / "starting.ncrst"
    parm.save(str(restart), format="ncrst")
    reloaded = AmberParm(str(topology), xyz=str(restart))
    if not np.array_equal(xyz, reloaded.coordinates) or not np.array_equal(
        box, reloaded.box
    ):
        raise ValueError("NetCDF restart changed source coordinates or periodic box")
    mdin = output / "relaxation.mdin"
    mdin.write_text(_mdin(audit))
    command = [
        str(sander.resolve()),
        "-O",
        "-i",
        str(mdin),
        "-p",
        str(topology.resolve()),
        "-c",
        str(restart),
        "-ref",
        str(restart),
        "-o",
        str(output / "minimized.mdout"),
        "-r",
        str(output / "optimized.ncrst"),
        "-x",
        str(output / "optimization.nc"),
        "-inf",
        str(output / "relaxation.mdinfo"),
    ]
    implementation_sources = [
        Path(__file__).with_name(name)
        for name in [
            "es_qm_region.py",
            "es_dftb_relaxation_gate.py",
            "es_dftb_analytic_pilot.py",
            "es_dftb_force_evidence.py",
        ]
    ]
    implementation_sources.append(
        Path(__file__).resolve().parents[2] / "scripts/run_nyl12_es_dftb_relaxation.py"
    )
    paths = [
        *implementation_sources,
        engine_config,
        force_check,
        region_config,
        *([pilot_decision] if pilot_decision is not None else []),
        topology,
        coordinates,
        sander,
        Path(__file__),
        restart,
        mdin,
        force_check.parent / "source_manifest.json",
        force_check.parent / "source/declared_config.json",
        *[Path(path) for path in engine["_force_source_paths"]],
    ]
    paths = list(dict.fromkeys(path.resolve() for path in paths))
    manifest = [{"path": str(path), "sha256": _sha(path)} for path in paths]
    source = output / "source"
    source.mkdir()
    for path in [
        engine_config,
        region_config,
        force_check,
        *([pilot_decision] if pilot_decision is not None else []),
        Path(__file__),
        *implementation_sources,
    ]:
        shutil.copy2(path, source / path.name)
    metadata = {
        "command": command,
        "environment": engine["environment"],
        "sources": manifest,
        "region_audit": audit,
        "pilot_decision": engine["_pilot_decision"],
        "numerical_force_validation_passed": engine["_pilot_decision"] is None,
        "source_box": box.tolist(),
        "box_absolute_tolerance_A_or_degree": 1e-10,
        "analytic_baseline_energy_kcal_per_mol": engine[
            "_baseline_energy_kcal_per_mol"
        ],
        "restraint": {
            "force_constant_kcal_mol_A2": 10.0,
            "all_non_QM_atoms_including_water_H": True,
            "reference": str(restart),
            "atom_count": len(parm.atoms) - len(selected),
        },
        "trajectory_policy": "Installed runmin supports ntwx=1: each force evaluation, not each minimizer iteration.",
        "scope": "At most25cycles; no convergence or ES-state claim; fixed box; no reaction-coordinate restraints.",
    }
    _write(output / "preparation.json", metadata)
    _write(output / "source_manifest.json", manifest)
    _write(
        output / "geometry_acceptance.json", {"accepted": False, "status": "not_run"}
    )
    return metadata


def run_prepared_relaxation(metadata: dict[str, Any]) -> dict[str, Any]:
    """Execute an explicitly requested prepared pilot and retain separate status.

    Args:
        metadata: Fresh record returned by prepare_es_relaxation.

    Returns:
        Process status; scientific endpoint acceptance is recorded separately.

    Raises:
        ValueError: Source mutation or invalid final coordinates.
        RuntimeError: Engine subprocess exits unsuccessfully.
    """
    command = metadata["command"]
    output = Path(command[command.index("-o") + 1]).parent
    if (output / "process.json").exists():
        raise ValueError("Engine attempt must not be rerun in place")
    if any(_sha(Path(row["path"])) != row["sha256"] for row in metadata["sources"]):
        raise ValueError("Preparation source changed before execution")
    status = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "exit_status": None,
    }
    _write(output / "process.json", status)
    env = os.environ.copy()
    env.update(metadata["environment"])
    try:
        with (output / "engine.log").open("w") as log:
            process = subprocess.run(
                command,
                cwd=output,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        status.update(
            exit_status=process.returncode,
            ended_utc=datetime.now(timezone.utc).isoformat(),
        )
        _write(output / "process.json", status)
        if process.returncode:
            raise RuntimeError(f"Bounded ES relaxation exited{process.returncode}")
        topology = Path(command[command.index("-p") + 1])
        final = AmberParm(str(topology), xyz=str(output / "optimized.ncrst"))
        finite = np.isfinite(final.coordinates).all()
        same_box = np.allclose(
            np.asarray(final.box), metadata["source_box"], rtol=0, atol=1e-10
        )
        _write(
            output / "geometry_acceptance.json",
            {
                "accepted": False,
                "coordinates_finite": bool(finite),
                "fixed_box_preserved": bool(same_box),
                "box_absolute_tolerance_A_or_degree": 1e-10,
                "status": "Scientific endpoint review required; exit does not certify convergence or ES state.",
            },
        )
        if not finite or not same_box:
            raise ValueError("Final coordinates nonfinite or periodic box changed")
        return status
    except Exception as exc:
        status.update(error=str(exc), ended_utc=datetime.now(timezone.utc).isoformat())
        _write(output / "process.json", status)
        raise
