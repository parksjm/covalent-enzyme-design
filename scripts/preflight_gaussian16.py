"""Run a small, recorded Gaussian 16 access and scratch preflight."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shlex
import subprocess
from pathlib import Path

LOGGER = logging.getLogger(__name__)
WATER_INPUT = """%mem=256MB
%nprocshared=1
#p HF/STO-3G SP SCF=Tight

Gaussian installation preflight; not a parameter calculation

0 1
O  0.000000  0.000000  0.000000
H  0.000000  0.757000  0.586000
H  0.000000 -0.757000  0.586000

"""


def run_preflight(install_root: Path, scratch_root: Path, output: Path) -> None:
    """Run a single-point water calculation and require normal termination.

    Args:
        install_root: Parent directory of the Gaussian ``g16`` installation.
        scratch_root: Existing user-writable scratch directory.
        output: New directory for inputs, logs, and provenance.

    Raises:
        PermissionError: If Gaussian or scratch is inaccessible.
        FileExistsError: If the output or its unique scratch directory exists.
        RuntimeError: If Gaussian exits unsuccessfully or lacks normal termination.
    """
    install_root = install_root.resolve(strict=True)
    scratch_root = scratch_root.resolve(strict=True)
    output = output.absolute()
    executable = install_root / "g16" / "g16"
    profile = install_root / "g16" / "bsd" / "g16.profile"
    if not os.access(executable, os.X_OK) or not os.access(profile, os.R_OK):
        raise PermissionError("Gaussian executable/profile is inaccessible")
    if not scratch_root.is_dir() or not os.access(scratch_root, os.W_OK | os.X_OK):
        raise PermissionError(f"Scratch is not writable: {scratch_root}")
    output.mkdir(parents=True, exist_ok=False)
    scratch = scratch_root / f"preflight_{output.name}"
    scratch.mkdir(mode=0o700, exist_ok=False)
    input_path = output / "water.gjf"
    log_path = output / "water.log"
    input_path.write_text(WATER_INPUT)
    shell_script = output / "run_gaussian.sh"
    shell_script.write_text(
        "#!/usr/bin/env bash\nset -eo pipefail\n"
        f"export g16root={shlex.quote(str(install_root))}\n"
        f"source {shlex.quote(str(profile))}\n"
        f"export GAUSS_SCRDIR={shlex.quote(str(scratch))}\n"
        f"cd -- {shlex.quote(str(output))}\n"
        f"{shlex.quote(str(executable))} < water.gjf > water.log 2>&1\n"
    )
    sources = [input_path, shell_script, Path(__file__).resolve(), executable, profile]
    report = {
        "purpose": "installation/access/scratch preflight only; no fitted charges",
        "install_root": str(install_root),
        "scratch": str(scratch),
        "argv": ["bash", str(shell_script)],
        "input_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sources
        },
    }
    provenance = output / "preflight.json"
    provenance.write_text(json.dumps(report, indent=2) + "\n")
    with (output / "launcher.log").open("w") as stream:
        result = subprocess.run(
            report["argv"], cwd=output, stdout=stream,
            stderr=subprocess.STDOUT, check=False,
        )
    normal = log_path.is_file() and "Normal termination of Gaussian 16" in log_path.read_text()
    report.update({"exit_status": result.returncode, "normal_termination": normal})
    provenance.write_text(json.dumps(report, indent=2) + "\n")
    if result.returncode or not normal:
        raise RuntimeError(f"Gaussian preflight failed; inspect {output}")
    LOGGER.info("Gaussian preflight passed; provenance: %s", provenance)


def main() -> None:
    """Parse explicit installation, scratch, and output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--install-root", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    run_preflight(args.install_root, args.scratch_root, args.output)


if __name__ == "__main__":
    main()
