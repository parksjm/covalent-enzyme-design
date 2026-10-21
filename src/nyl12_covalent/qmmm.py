"""Reproducible AmberTools26 DFTB3/MM geometry optimization inputs."""

from __future__ import annotations

import json
import logging
import os
import subprocess
from pathlib import Path

from .qm_region import QM_MASK, audit_qm_region

LOGGER = logging.getLogger(__name__)


def write_dftb3_input(output: Path, max_cycles: int = 400) -> Path:
    """Write a restrained acyl-state DFTB3/ff19SB-GAFF minimization input.

    The entire beta Thr1 and retained acyl ligand are quantum atoms, avoiding
    an artificial cut across the new ester. Backbone atoms outside the reactive
    residue are restrained while their side chains and ligand may relax.

    Args:
        output: Acyl model directory with validated topology and coordinates.
        max_cycles: Maximum Amber minimization cycles.

    Returns:
        Amber input path.
    """
    if max_cycles < 1:
        raise ValueError("max_cycles must be positive")
    audit = json.loads((output / "topology_audit.json").read_text())
    if abs(audit["qm_region_partial_charge"] - 1.0) > 0.002:
        raise ValueError("QM region charge audit must pass before optimization")
    mdin = output / "acyl_dftb3.mdin"
    mdin.write_text(
        "Nyl12 acyl enzyme: DFTB3 QM/MM geometry optimization\n"
        " &cntrl\n"
        "  imin=1, ntmin=3, maxcyc=" + str(max_cycles) + ",\n"
        "  drms=0.05, ntb=0, igb=8, cut=999.0,\n"
        "  ntc=1, ntf=1, ntpr=10, ntwr=100,\n"
        "  ntr=1, restraint_wt=5.0,\n"
        "  restraintmask=':1-928,930-1312@N,CA,C,O',\n"
        "  ifqnt=1,\n"
        " /\n"
        " &qmmm\n"
        "  qmmask=':929,1313', qmcharge=1,\n"
        "  qm_theory='DFTB3', scfconv=1.0d-8,\n"
        "  verbosity=1,\n"
        " /\n"
    )
    return mdin


def run_dftb3(output: Path, amber_bin: Path, max_cycles: int = 400) -> Path:
    """Run and record Amber DFTB3/MM acyl-state geometry optimization.

    Args:
        output: Validated acyl model directory.
        amber_bin: AmberTools26 bin directory.
        max_cycles: Maximum minimization cycles.

    Returns:
        Sander output path. A nonzero exit raises with logs retained.
    """
    mdin = write_dftb3_input(output, max_cycles)
    mdout = output / "acyl_dftb3.mdout"
    restart = output / "acyl_dftb3.rst7"
    command = [
        str(amber_bin / "sander"), "-O", "-i", str(mdin), "-o", str(mdout),
        "-p", str(output / "acyl.prmtop"),
        "-c", str(output / "acyl.inpcrd"),
        "-ref", str(output / "acyl.inpcrd"),
        "-r", str(restart),
        "-inf", str(output / "acyl_dftb3.mdinfo"),
    ]
    environment = os.environ.copy()
    environment["AMBERHOME"] = str(amber_bin.parent)
    environment["OMP_NUM_THREADS"] = "4"
    environment["LD_LIBRARY_PATH"] = ":".join(
        part for part in (str(amber_bin.parent / "lib"), environment.get("LD_LIBRARY_PATH"))
        if part
    )
    (output / "qmmm_command.json").write_text(
        json.dumps({"argv": command, "environment": {
            "AMBERHOME": environment["AMBERHOME"],
            "OMP_NUM_THREADS": environment["OMP_NUM_THREADS"],
            "LD_LIBRARY_PATH": environment["LD_LIBRARY_PATH"],
        }}, indent=2) + "\n"
    )
    LOGGER.info("Running DFTB3 optimization in %s", output)
    with (output / "sander.stderr.log").open("w") as stderr:
        result = subprocess.run(
            command, cwd=output, env=environment, stdout=stderr,
            stderr=subprocess.STDOUT, check=False,
        )
    if result.returncode:
        raise RuntimeError(f"Amber DFTB3 failed; inspect {mdout} and sander.stderr.log")
    if not restart.is_file():
        raise ValueError("Amber finished without writing optimized coordinates")
    return mdout


def write_ti1_dftb3_input(output: Path, max_cycles: int = 300) -> Path:
    """Write DFTB3/MM input for the first tetrahedral adduct plus relay water.

    Args:
        output: Directory containing the validated TI1 topology.
        max_cycles: Maximum Amber minimization cycles.

    Returns:
        Amber minimization input path.
    """
    if max_cycles < 1:
        raise ValueError("max_cycles must be positive")
    audit_qm_region(output / "ti1.prmtop")
    mdin = output / "ti1_dftb3.mdin"
    with mdin.open("x") as stream:
        stream.write(
            "Nyl12 first tetrahedral adduct with neutral relay water\n"
            " &cntrl\n"
            f"  imin=1, ntmin=3, maxcyc={max_cycles},\n"
            "  drms=0.05, ntb=0, igb=8, cut=999.0,\n"
            "  ntc=1, ntf=1, ntpr=10, ntwr=100,\n"
            "  ntr=1, restraint_wt=5.0,\n"
            "  restraintmask=':1-928,930-1312@N,CA,C,O',\n"
            "  ifqnt=1,\n"
            " /\n"
            " &qmmm\n"
            f"  qmmask='{QM_MASK}', qmcharge=0,\n"
            "  qm_theory='DFTB3', scfconv=1.0d-8,\n"
            "  verbosity=1,\n"
            " /\n"
        )
    return mdin


def run_ti1_dftb3(
    output: Path,
    amber_bin: Path,
    max_cycles: int = 300,
    *,
    coordinates: Path | None = None,
    reference: Path | None = None,
    attempt: Path | None = None,
    target: str | None = None,
    checkpoint_attempt: Path | None = None,
) -> Path:
    """Run TI1 only with explicit coordinates and a fresh attempt directory.

    Args:
        output: Directory containing the source TI1 topology.
        amber_bin: AmberTools26 bin directory.
        max_cycles: Maximum minimization cycles.
        coordinates: Explicit starting coordinates or source checkpoint.
        reference: Original positional-restraint reference.
        attempt: New attempt directory; existing paths are rejected.
        target: Product-gated substrate identifier.
        checkpoint_attempt: Source attempt for checkpoint continuation.

    Returns:
        Sander output path only after scientific acceptance.

    Raises:
        ValueError: If explicit safe-run arguments are missing.
    """
    from .ti1_attempt import run_attempt

    if any(value is None for value in (coordinates, reference, attempt, target)):
        raise ValueError(
            "TI1 requires explicit coordinates, reference, attempt, and target; "
            "use scripts/run_nyl12_ti1_qmmm.py inside the tmux job launcher"
        )
    return run_attempt(
        topology=output / "ti1.prmtop",
        coordinates=coordinates,
        reference=reference,
        output=attempt,
        amber_bin=amber_bin,
        target=target,
        max_cycles=max_cycles,
        checkpoint_attempt=checkpoint_attempt,
    )
