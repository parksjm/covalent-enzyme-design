"""Classical TIP3P water preparation with the entire nonwater system fixed."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from parmed.amber import AmberParm
from parmed import Structure

from .local_repair_engine import run_checked_engine
from .ti1_checks import check_convergence


def preminimize_water(
    output: Path, executable: Path, environment: dict[str, str], pairs: list[list[int]]
) -> Path:
    """Minimize only TIP3P water; frozen reactive parameters cannot move TI1.

    Args:
        output: Calculation directory with prepared topology and coordinates.
        executable: Explicit Sander executable.
        environment: Recorded Amber runtime environment.
        pairs: Two bridge distance atom pairs in Amber numbering.

    Returns:
        Water-preconditioned restart after engine and frozen-coordinate checks.
    """
    stage = output / "00_tip3p_water_only"
    stage.mkdir()
    before = AmberParm(str(output / "model.prmtop"), xyz=str(output / "prepared.rst7"))
    order = [a.idx for a in before.residues[329].atoms] + [
        a.idx for a in before.atoms if a.residue.idx != 329
    ]
    generic = before.copy(Structure)
    reordered = AmberParm.from_structure(generic[":330"] + generic["!:330"])
    reordered.coordinates = before.coordinates[order]
    for new_atom, old_index in zip(reordered.atoms, order, strict=True):
        old_atom = before.atoms[old_index]
        if new_atom.type != old_atom.type or any(
            abs(getattr(new_atom, field) - getattr(old_atom, field)) > 1e-6
            for field in ("charge", "rmin", "epsilon", "solvent_radius", "screen")
        ):
            raise ValueError(
                "Parameters changed during temporary GB belly atom reorder"
            )
    reordered.write_parm(str(stage / "water_first.prmtop"))
    reordered.save(str(stage / "water_first.rst7"), format="rst7")
    inverse = {old: new + 1 for new, old in enumerate(order)}
    reordered_pairs = [[inverse[a - 1], inverse[b - 1]] for a, b in pairs]
    (stage / "atom_order.json").write_text(
        json.dumps(
            {"original_atom_indices_1based_in_temporary_order": [i + 1 for i in order]},
            indent=2,
        )
        + "\n"
    )
    mdin = stage / "minimize.mdin"
    mdin.write_text(
        "TIP3P water-only coordinate preconditioning; nonwater atoms fixed\n"
        " &cntrl\n imin=1, ntmin=1, ncyc=50, maxcyc=200, drms=0.001,\n"
        " ntb=0, igb=8, cut=999.0, ntc=1, ntf=1, ntpr=10, ntwx=1, ioutfm=0,\n"
        " ibelly=1, bellymask=':1', ifqnt=0, nmropt=1,\n /\n"
        " &wt type='END' /\nDISANG=bridge.RST\n"
    )
    (stage / "bridge.RST").write_text(
        "".join(
            f" &rst iat={a},{b}, r1=2.2, r2=2.6, r3=3.2, r4=3.6, rk2=10.0, rk3=10.0, /\n"
            for a, b in reordered_pairs
        )
    )
    restart, mdout = stage / "optimized.rst7", stage / "minimize.mdout"
    command = [
        str(executable.absolute()),
        "-i",
        str(mdin),
        "-o",
        str(mdout),
        "-p",
        str(stage / "water_first.prmtop"),
        "-c",
        str(stage / "water_first.rst7"),
        "-r",
        str(restart),
        "-inf",
        str(stage / "minimize.mdinfo"),
        "-x",
        str(stage / "optimization.mdcrd"),
    ]
    (stage / "command.json").write_text(
        json.dumps(
            {
                "argv": command,
                "mobile_atoms": "WAT330 O/H1/H2 only",
                "potential": "inherited classical TIP3P/environment parameters; coordinate preparation only; not TI1 energy validation",
            },
            indent=2,
        )
        + "\n"
    )
    process = run_checked_engine(
        command, cwd=stage, env=environment, log_path=stage / "sander.log", mdout=mdout
    )
    (stage / "process.json").write_text(json.dumps(process, indent=2) + "\n")
    if process["exit_status"] or not restart.exists():
        raise RuntimeError("TIP3P water-only pre-minimization failed")
    temporary = AmberParm(str(stage / "water_first.prmtop"), xyz=str(restart))
    mapped = np.empty_like(temporary.coordinates)
    mapped[order] = temporary.coordinates
    after = AmberParm(str(output / "model.prmtop"))
    after.coordinates = mapped
    restored_restart = stage / "original_order.rst7"
    after.save(str(restored_restart), format="rst7")
    fixed = [a.idx for a in before.atoms if a.residue.idx != 329]
    maximum = float(
        np.max(
            np.linalg.norm(before.coordinates[fixed] - after.coordinates[fixed], axis=1)
        )
    )
    report = {
        "nonwater_maximum_displacement_A": maximum,
        "oxygen_displacement_from_AF3_A": float(
            np.linalg.norm(before.coordinates[-3] - after.coordinates[-3])
        ),
        "convergence": check_convergence(mdout, 0.001),
        "used_for": "coordinate preconditioning only; no energetic ranking",
    }
    (stage / "acceptance.json").write_text(json.dumps(report, indent=2) + "\n")
    if maximum > 1e-6:
        raise ValueError(
            "Nonwater atoms moved during fixed-environment preconditioning"
        )
    return restored_restart
