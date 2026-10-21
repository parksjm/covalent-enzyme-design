"""Run and cache a Gaussian template Opt/Freq/ESP calculation, never a force-field fit."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.gaussian_checks import assess_gaussian  # noqa: E402
from nyl12_covalent.parameter_cache import ParameterArtifactCache  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, value: object) -> None:
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def run_template(
    source: Path,
    output: Path,
    install_root: Path,
    scratch_root: Path,
    amber_bin: Path,
    cache_root: Path,
) -> None:
    """Execute a prepared template with provenance and scientific acceptance.

    Args:
        source: Prepared directory containing input.gjf, molecule.sdf and atom_map.json.
        output: Fresh attempt directory.
        install_root: Parent of the Gaussian g16 installation directory.
        scratch_root: Existing authorized Gaussian scratch directory.
        amber_bin: Installed AmberTools directory containing espgen.
        cache_root: Validated content-addressed QM-artifact cache.

    Raises:
        RuntimeError: If engine execution or scientific acceptance fails.
        FileExistsError: If output or scratch already exists.
        PermissionError: If Gaussian or the authorized scratch is inaccessible.
    """
    source, install_root, scratch_root, amber_bin = (
        p.resolve(strict=True) for p in (source, install_root, scratch_root, amber_bin)
    )
    output = output.absolute()
    if output.exists():
        raise FileExistsError(output)
    executable = install_root / "g16/g16"
    profile = install_root / "g16/bsd/g16.profile"
    espgen = amber_bin / "espgen"
    if not os.access(executable, os.X_OK) or not os.access(profile, os.R_OK):
        raise PermissionError(
            "Gaussian executable/profile inaccessible; acquire granted g16 group"
        )
    if not os.access(espgen, os.X_OK):
        raise PermissionError(f"espgen inaccessible: {espgen}")
    if not scratch_root.is_dir() or not os.access(scratch_root, os.W_OK | os.X_OK):
        raise PermissionError(f"Authorized scratch inaccessible: {scratch_root}")
    paths = [source / name for name in ("input.gjf", "molecule.sdf", "atom_map.json")]
    identity = {
        "stage": "neutral_L_Thr_NME_QM_optimization_frequency_ESP_not_fitted_parameters_v1",
        "input_sha256": {p.name: _sha(p) for p in paths},
        "engine_sha256": {
            "g16": _sha(executable),
            "profile": _sha(profile),
            "espgen": _sha(espgen),
        },
        "acceptance_source_sha256": _sha(
            ROOT / "src/nyl12_covalent/gaussian_checks.py"
        ),
        "software_versions": {name: version(name) for name in ("rdkit", "numpy")},
    }
    cache = ParameterArtifactCache(cache_root)
    with cache.calculation_lock(identity):
        restored = cache.materialize(identity, output)
        if restored is not None:
            _write(
                output / "cache_reuse.json",
                {"key": cache.key_for(identity), "identity": identity},
            )
            LOGGER.info("Reused validated QM artifacts, not fitted charges: %s", output)
            return
        output.mkdir(parents=True, exist_ok=False)
        for path in paths:
            shutil.copyfile(path, output / path.name)
            if _sha(output / path.name) != identity["input_sha256"][path.name]:
                raise RuntimeError("Prepared input changed during snapshot")
        snapshots = [
            Path(__file__),
            ROOT / "src/nyl12_covalent/gaussian_checks.py",
            ROOT / "src/nyl12_covalent/parameter_cache.py",
            ROOT / "uv.lock",
            ROOT / "src/nyl12_covalent/terminal_template.py",
        ]
        snapshot_dir = output / "source"
        snapshot_dir.mkdir()
        for path in snapshots:
            shutil.copyfile(path, snapshot_dir / path.name)
        scratch = scratch_root / (
            "neutral_thr_" + hashlib.sha256(str(output).encode()).hexdigest()[:16]
        )
        scratch.mkdir(mode=0o700, exist_ok=False)
        shell = output / "run_gaussian.sh"
        shell.write_text(
            "#!/usr/bin/env bash\nset -eo pipefail\n"
            f"export g16root={shlex.quote(str(install_root))}\n"
            f"source {shlex.quote(str(profile))}\n"
            f"export GAUSS_SCRDIR={shlex.quote(str(scratch))}\n"
            f"cd -- {shlex.quote(str(output))}\n"
            f"{shlex.quote(str(executable))} < input.gjf > gaussian.log 2>&1\n"
        )
        _write(
            output / "calculation.json",
            {
                "identity": identity,
                "cache_key": cache.key_for(identity),
                "cache_hit": False,
                "source": str(source),
                "scratch": str(scratch),
                "python": sys.executable,
                "argv": ["bash", str(shell)],
                "source_sha256": {p.name: _sha(p) for p in snapshot_dir.iterdir()},
                "shell_sha256": _sha(shell),
            },
        )
        acceptance = {"accepted": False, "fitted_forcefield_ready": False}
        process = {
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "exit_status": None,
        }
        try:
            with (output / "launcher.log").open("x") as stream:
                result = subprocess.run(
                    ["bash", str(shell)],
                    cwd=output,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            process["exit_status"] = result.returncode
            if result.returncode:
                raise RuntimeError(f"Gaussian exited {result.returncode}; see {output}")
            acceptance.update(
                assess_gaussian(output / "gaussian.log", output / "molecule.sdf")
            )
            if not acceptance["accepted"]:
                raise RuntimeError(
                    "Gaussian completed but geometry/frequency/ESP checks failed"
                )
            command = [
                str(espgen),
                "-i",
                str(output / "gaussian.log"),
                "-o",
                str(output / "esp.dat"),
            ]
            _write(output / "espgen_command.json", {"argv": command})
            with (output / "espgen.log").open("x") as stream:
                subprocess.run(
                    command, stdout=stream, stderr=subprocess.STDOUT, check=True
                )
            header = (output / "esp.dat").read_text().splitlines()[0].split()
            atoms, points = int(header[0]), int(header[1])
            if atoms != len(acceptance["final_coordinates_A"]) or points < 1:
                raise ValueError("ESP atom/grid counts invalid")
            acceptance["esp_atom_count"] = atoms
            acceptance["esp_grid_points"] = points
        except Exception as error:
            acceptance["accepted"] = False
            acceptance["error"] = f"{type(error).__name__}: {error}"
            raise
        finally:
            process["ended_utc"] = datetime.now(timezone.utc).isoformat()
            _write(output / "process.json", process)
            _write(output / "acceptance.json", acceptance)
        artifacts = {
            str(p.relative_to(output)): p for p in output.rglob("*") if p.is_file()
        }
        cache.store(
            identity,
            artifacts,
            label=f"neutral_thr_nme_{source.name}_qm_geometry_frequency_esp",
        )
        LOGGER.info(
            "Validated QM template artifacts cached; RESP/bonded fitting remains: %s",
            output,
        )


def main() -> None:
    """Parse explicit locations and run a prepared Gaussian parameter template."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "source",
        "output",
        "install-root",
        "scratch-root",
        "amber-bin",
        "cache-root",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_template(
        args.source,
        args.output,
        args.install_root,
        args.scratch_root,
        args.amber_bin,
        args.cache_root,
    )


if __name__ == "__main__":
    main()
