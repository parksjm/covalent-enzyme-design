from __future__ import annotations
import math
import re
import subprocess
from pathlib import Path
from typing import Any
from reference_common import digest, journal, write_json
from reference_embedding import write_direct_gaussian


def run_direct_gaussian(
    source: Path, folder: Path, sander_output: Path, environment: dict[str, str]
) -> dict[str, Any]:
    """Run direct Gaussian Force on the exact QM/link and external-PC export."""
    folder.mkdir(parents=True, exist_ok=False)
    deck = write_direct_gaussian(source, folder)
    scratch = folder / "scratch"
    scratch.mkdir()
    environment = environment | {"GAUSS_SCRDIR": str(scratch.resolve())}
    journal("intent", f"Run direct Gaussian Force at {folder}")
    with (
        deck.open("r") as source_stream,
        (folder / "direct_force.log").open("x") as log,
    ):
        process = subprocess.run(
            ["/usr/local/g16/g16"],
            stdin=source_stream,
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=folder,
            env=environment,
            check=False,
        )
    output = (folder / "direct_force.log").read_text(errors="replace")
    energies = re.findall(r"SCF Done:\s+E\([^)]*\)\s*=\s*([-+\d.EeDd]+)", output)
    if (
        process.returncode != 0
        or not energies
        or "normal termination of gaussian" not in output.lower()
    ):
        raise RuntimeError(f"Direct Gaussian Force failed: {folder}")
    formchk = subprocess.run(
        ["formchk", "direct_force.chk", "direct_force.fchk"],
        cwd=folder,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    fchk = folder / "direct_force.fchk"
    if formchk.returncode != 0 or not fchk.is_file():
        raise RuntimeError(
            f"Gaussian formchk failed in {folder}: {formchk.stdout.decode(errors='replace')}"
        )
    fchk_text = fchk.read_text(errors="replace")
    atoms = re.search(r"Number of atoms\s+I\s+(\d+)", fchk_text)
    grad = re.search(r"Cartesian Gradient\s+R\s+N=\s*(\d+)\s*\n", fchk_text)
    if atoms is None or int(atoms[1]) != 110 or grad is None or int(grad[1]) != 330:
        raise ValueError(
            "Gaussian checkpoint must contain110 atoms and330 gradient values"
        )
    values = fchk_text[grad.end() :].split()[:330]
    if len(values) != 330 or not all(
        math.isfinite(float(v.replace("D", "E"))) for v in values
    ):
        raise ValueError("Gaussian analytic gradient is incomplete or nonfinite")
    self_terms = re.findall(r"Self energy of the charges\s*=\s*([-+\d.EeDd]+)", output)
    if not self_terms:
        raise ValueError("Gaussian did not report external-charge self energy")
    self_energy_hartree = float(self_terms[-1].replace("D", "E"))
    raw_hartree = float(energies[-1].replace("D", "E"))
    corrected_hartree = raw_hartree - self_energy_hartree
    journal(
        "outcome",
        f"Run direct Gaussian Force at {folder}",
        "normal termination;110 atoms,330 gradient values",
    )
    write_json(
        folder / "process.json",
        {
            "status": "finished",
            "exit_status": process.returncode,
            "normal_termination": True,
            "gradient_values": 330,
            "input_sha256": digest(deck),
            "log_sha256": digest(folder / "direct_force.log"),
            "fchk_sha256": digest(fchk),
        },
    )
    return {
        "raw_gaussian_scf_hartree": raw_hartree,
        "point_charge_self_energy_hartree": self_energy_hartree,
        "corrected_gaussian_energy_kcal_mol": corrected_hartree * 627.5094740631,
        "direct_input_sha256": digest(deck),
        "gaussian_log": str(folder / "direct_force.log"),
        "fchk": str(fchk),
        "gradient_artifact": str(fchk),
        "normal_termination": True,
    }
