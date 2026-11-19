"""Inputs and checks for fully unrestrained classical E-S relaxation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from nyl12_covalent.es_system_preparation import sha256


def minimization_input(max_cycles: int = 10_000, case_label: str = "Model 23") -> str:
    """Return a periodic, fully unrestrained Amber minimization input."""
    if max_cycles < 1:
        raise ValueError("max_cycles must be positive")
    if not case_label.strip():
        raise ValueError("case_label must be nonempty")
    return f"""{case_label} fully unrestrained classical E-S relaxation
&cntrl
 imin=1, maxcyc={max_cycles}, ncyc={max_cycles // 2}, ntmin=1,
 ntb=1, igb=0, cut=10.0, ntc=1, ntf=1,
 ntr=0, nmropt=0, drms=0.01,
 ntpr=25, ntwx=25, ntwr=250, ioutfm=1, ntxo=2,
/
"""


def solvent_relaxation_input(
    max_cycles: int = 10_000, case_label: str = "Model 23"
) -> str:
    """Return an input that relaxes new solvent around fixed solute heavy atoms."""
    if max_cycles < 1:
        raise ValueError("max_cycles must be positive")
    if not case_label.strip():
        raise ValueError("case_label must be nonempty")
    return f"""{case_label} solvent/ion relaxation around the accepted static solute
&cntrl
 imin=1, maxcyc={max_cycles}, ncyc={max_cycles // 2}, ntmin=1,
 ntb=1, igb=0, cut=10.0, ntc=1, ntf=1,
 ntr=1, restraint_wt=10.0, restraintmask=':1-1313 & !@H=', nmropt=0,
 drms=0.01, ntpr=25, ntwx=25, ntwr=250, ioutfm=1, ntxo=2,
/
"""


def prepare_relaxation(
    topology: Path,
    coordinates: Path,
    engine: Path,
    output: Path,
    case_label: str = "Model 23",
    max_cycles: int = 10_000,
) -> dict[str, Any]:
    """Write source-pinned inputs for an unrestrained classical minimization."""
    for path in (topology, coordinates, engine):
        if not path.is_file():
            raise FileNotFoundError(path)
    output.mkdir(parents=True, exist_ok=False)
    mdin = output / "minimize.mdin"
    mdin.write_text(
        minimization_input(max_cycles=max_cycles, case_label=case_label)
    )
    command = [
        str(engine.resolve()),
        "-O",
        "-i",
        str(mdin.resolve()),
        "-p",
        str(topology.resolve()),
        "-c",
        str(coordinates.resolve()),
        "-o",
        str((output / "minimized.mdout").resolve()),
        "-r",
        str((output / "minimized.ncrst").resolve()),
        "-x",
        str((output / "optimization.nc").resolve()),
        "-inf",
        str((output / "mdinfo").resolve()),
    ]
    sources = [topology, coordinates, engine, Path(__file__).resolve(), mdin]
    metadata = {
        "case_label": case_label,
        "max_cycles": max_cycles,
        "command": command,
        "restraints": "none; ntr=0, nmropt=0, all atoms active",
        "trajectory_stride_minimization_cycles": 25,
        "sources": [
            {"path": str(path.resolve()), "sha256": sha256(path)} for path in sources
        ],
    }
    (output / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def prepare_solvent_relaxation(
    topology: Path,
    coordinates: Path,
    engine: Path,
    output: Path,
    case_label: str = "Model 23",
) -> dict[str, Any]:
    """Write a source-pinned solvent-relaxation stage for the static solute."""
    for path in (topology, coordinates, engine):
        if not path.is_file():
            raise FileNotFoundError(path)
    output.mkdir(parents=True, exist_ok=False)
    mdin = output / "minimize.mdin"
    mdin.write_text(solvent_relaxation_input(case_label=case_label))
    command = [
        str(engine.resolve()),
        "-O",
        "-i",
        str(mdin.resolve()),
        "-p",
        str(topology.resolve()),
        "-c",
        str(coordinates.resolve()),
        "-ref",
        str(coordinates.resolve()),
        "-o",
        str((output / "minimized.mdout").resolve()),
        "-r",
        str((output / "solvent_relaxed.ncrst").resolve()),
        "-x",
        str((output / "optimization.nc").resolve()),
        "-inf",
        str((output / "mdinfo").resolve()),
    ]
    sources = [topology, coordinates, engine, Path(__file__).resolve(), mdin]
    metadata = {
        "case_label": case_label,
        "command": command,
        "restraints": {
            "purpose": "Relax LEaP-added solvent/ions without replacing the accepted static solute pose",
            "mask": ":1-1313 & !@H=",
            "weight_kcal_mol_A2": 10.0,
            "released_in_next_DFTB3_stage": True,
        },
        "sources": [
            {"path": str(path.resolve()), "sha256": sha256(path)} for path in sources
        ],
    }
    (output / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata
