"""Launch explicitly prepared, matched NEB trials through the tmux job runner."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.jobs import launch_job  # noqa: E402

LOGGER = logging.getLogger(__name__)


def render_trial_input(
    base: str, controls: list[Path], images: int, optimizer: str
) -> str:
    """Change only the seed paths and declared numerical NEB test controls.

    Args:
        base: Accepted original native NEB input, including its Hamiltonian.
        controls: Exactly thirteen prepared Amber coordinate files.
        images: Regular band images, excluding the disabled climbing slot.
        optimizer: Supported minimization algorithm for this bounded comparison.

    Returns:
        Native input with unchanged charge, masks, restraints and tolerances.
    """
    if len(controls) != 13 or images not in (24, 48):
        raise ValueError("Expected thirteen controls and 24 or 48 regular images")
    if optimizer not in ("LBFGS", "CG"):
        raise ValueError("Unsupported trial optimizer")
    prefix, delimiter, block = base.partition("&dlfind")
    if not delimiter:
        raise ValueError("Missing native DL-FIND input")
    for key in ("ntr", "nmropt"):
        if re.findall(rf"\b{key}\s*=\s*(\d+)", prefix, re.I) != ["0"]:
            raise ValueError(f"Expected explicit {key}=0")
    replacements = {
        "neb_nimage": str(images + 1),
        "optalg": f"'{optimizer}'",
        "maxstep": "0.02",
        "refcrd": f"'{controls[0]}'",
    }
    for key, value in replacements.items():
        block, count = re.subn(
            rf"\b{key}\s*=\s*(?:'[^']*'|[^,\n]+)",
            f"{key}={value}",
            block,
            flags=re.I,
        )
        if count != 1:
            raise ValueError(f"Expected exactly one {key}")
    for index, path in enumerate(controls[1:], 1):
        block, count = re.subn(
            rf"\bcoords2\({index}\)\s*=\s*'[^']*'",
            f"coords2({index})='{path}'",
            block,
            flags=re.I,
        )
        if count != 1:
            raise ValueError(f"Expected exactly one coords2({index})")
    return prefix + delimiter + block


def launch_trial(config_path: Path) -> Path:
    """Launch one prepared trial with retained sources and explicit environment.

    Args:
        config_path: Trial configuration whose inputs passed preparation gates.

    Returns:
        Recorded job runner script; an existing job directory is an error.
    """
    config = json.loads(config_path.read_text())
    trial = config_path.parent
    preflight = json.loads((trial / "preflight.json").read_text())
    if preflight.get("scientific_preparation_passed") is not True:
        raise ValueError("Trial preparation has not passed its scientific gates")
    if not (trial / "source_snapshot").is_dir():
        raise FileNotFoundError(trial / "source_snapshot")
    calculation = Path(config["calculation"])
    mdin = calculation / "minimize.mdin"
    text = mdin.read_text()
    if config["qm_charge"] != -1 or config["active_atoms"] != 4662:
        raise ValueError("Trial charge or movable region differs from approved case")
    if config["regular_images"] not in (24, 48):
        raise ValueError("Unsupported band size")
    if f"neb_nimage={config['regular_images'] + 1}" not in text:
        raise ValueError("Input band size differs from configuration")
    if config["reactant"] != config["images"][0]:
        raise ValueError("Reactant differs from first seed control")
    sander = Path(config["sander"])
    amber = sander.parent.parent
    command = [
        str(sander),
        "-O",
        "-i",
        "minimize.mdin",
        "-p",
        config["topology"],
        "-c",
        config["reactant"],
        "-o",
        "minimized.mdout",
        "-r",
        "optimized.ncrst",
        "-x",
        "optimization.nc",
        "-inf",
        "mdinfo",
    ]
    inputs = [
        config_path,
        mdin,
        trial / "preflight.json",
        trial / "launch.sh",
        Path(config["topology"]),
        sander,
        amber / "lib/libdlfind.so",
        Path(__file__),
        ROOT / "scripts/run_nyl12_job.py",
        ROOT / "src/nyl12_covalent/jobs.py",
        *[Path(path) for path in config["images"]],
        *sorted((trial / "source_snapshot").glob("*")),
    ]
    return launch_job(
        f"nyl12_neb_{trial.name}_20261119",
        trial / "job",
        calculation,
        command,
        env_overrides={
            "AMBERHOME": str(amber),
            "LD_LIBRARY_PATH": str(amber / "lib"),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "GFORTRAN_UNBUFFERED_ALL": "1",
        },
        input_paths=inputs,
    )


def main() -> None:
    """Launch only the explicitly supplied prepared configurations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("configs", nargs="+", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    for path in args.configs:
        script = launch_trial(path.resolve())
        LOGGER.info("Launched prepared trial %s: %s", path.parent.name, script)


if __name__ == "__main__":
    main()
