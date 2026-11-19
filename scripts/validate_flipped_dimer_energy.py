"""Validate the prepared saddle coordinates against the native band energy.

Run this script inside the project tmux launcher. It performs a bounded native
DL-FIND initial energy/gradient evaluation with the original Hamiltonian.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path


def option(text: str, name: str) -> str:
    """Read one required quoted native namelist option."""
    values = re.findall(rf"\b{name}\s*=\s*'([^']*)'", text, re.I)
    if len(values) != 1:
        raise ValueError(f"Expected unique {name}")
    return values[0]


def main() -> None:
    """Evaluate the prepared initial geometry and unlock only a matching job."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--neb-dir", type=Path, required=True)
    args = parser.parse_args()
    prepared, neb = (
        args.prepared.resolve(strict=True),
        args.neb_dir.resolve(strict=True),
    )
    manifest_path = prepared / "preparation.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest["native_cycle"] != 499 or manifest["native_image_1based"] != 8:
        raise ValueError("Expected final cycle499 image8")
    for relative, expected in manifest["prepared_sha256"].items():
        observed = hashlib.sha256((prepared / relative).read_bytes()).hexdigest()
        if observed != expected:
            raise ValueError(f"Prepared input changed: {relative}")
    folder = prepared / "energy_validation"
    folder.mkdir(exist_ok=False)
    original = (prepared / "calculation/dimer.mdin").read_text()
    prefix = original[: original.lower().index("&dlfind")]
    prefix, count = re.subn(r"\bmaxcyc\s*=\s*\d+", "maxcyc=1", prefix, flags=re.I)
    if count != 1:
        raise ValueError("Expected one native maxcyc")
    text = (
        prefix + "&dlfind\n optalg='LBFGS', crdrep='CART',\n"
        f" active='{option(original, 'active')}',\n"
        f" refcrd='{option(original, 'refcrd')}',\n"
        " maxcycle=1, tol=1.d-4, tole=1.d-5, maxstep=0.1, verbosity=2,\n/\n"
    )
    (folder / "initial_energy.mdin").write_text(text)
    old = json.loads((neb / "job/command.json").read_text())
    command = [
        old["command"][0],
        "-O",
        "-i",
        "initial_energy.mdin",
        "-p",
        str(prepared / "inputs/system.prmtop"),
        "-c",
        str(prepared / "inputs/ts_guess.ncrst"),
        "-o",
        "initial_energy.mdout",
        "-r",
        "unused_optimized.ncrst",
        "-inf",
        "mdinfo",
    ]
    environment = dict(os.environ)
    environment.update(old["env_overrides"])
    (folder / "command.json").write_text(
        json.dumps(
            {
                "command": command,
                "cwd": str(folder),
                "env_overrides": old["env_overrides"],
                "script_sha256": hashlib.sha256(
                    Path(__file__).read_bytes()
                ).hexdigest(),
                "input_sha256": hashlib.sha256(text.encode()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    shutil.copy2(Path(__file__), folder / "validation_source.py")
    with (folder / "engine.log").open("w") as handle:
        process = subprocess.run(
            command,
            cwd=folder,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if process.returncode:
        raise RuntimeError(f"Native energy evaluation exited {process.returncode}")
    energies = re.findall(
        r"Energy calculation finished, energy:\s*([-+\d.Ee]+)",
        (folder / "initial_energy.mdout").read_text(),
    )
    if not energies:
        raise ValueError("Native DL-FIND initial Hartree energy missing")
    value = float(energies[0])
    difference = (value - manifest["native_energy_hartree"]) * 627.5094740631
    result = {
        "native_cycle": 499,
        "native_image_1based": 8,
        "expected_energy_hartree": manifest["native_energy_hartree"],
        "recomputed_initial_energy_hartree": value,
        "difference_kcal_mol": difference,
        "tolerance_kcal_mol": 0.01,
        "passed": abs(difference) < 0.01,
        "interpretation": "Native coordinate/energy association; no saddle or method validation",
    }
    (folder / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    if not result["passed"]:
        raise ValueError(f"Native band energy mismatch: {difference} kcal/mol")
    shutil.copy2(manifest_path, folder / "preparation_before_validation.json")
    manifest["native_xyz_cycle_energy_match_verified"] = True
    manifest["native_energy_validation_result"] = str(folder / "result.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    logging.info("Native image energy matched within %.6f kcal/mol", abs(difference))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
