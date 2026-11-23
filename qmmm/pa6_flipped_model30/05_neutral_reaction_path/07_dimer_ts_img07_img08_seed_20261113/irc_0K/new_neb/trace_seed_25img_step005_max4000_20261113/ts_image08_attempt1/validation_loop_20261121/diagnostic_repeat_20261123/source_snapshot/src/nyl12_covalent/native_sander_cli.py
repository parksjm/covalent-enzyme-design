"""One-cycle native Sander CLI interface and force-calibration gates."""

from __future__ import annotations

import json
import hashlib
import math
import re
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import Rst7

from scripts.prepare_dimer_continuation import read_xyz_frames, write_restart

FloatArray = np.ndarray
ANGSTROM_PER_BOHR = 0.529177
HARTREE_TO_KCAL_PER_MOL = 627.5094740631
PATH_FORCE_TO_GRADIENT_SCALE = -HARTREE_TO_KCAL_PER_MOL / ANGSTROM_PER_BOHR**2
FORCE_SOURCE_FILES = (
    Path(
        "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/sander/dlfind_module.F90"
    ),
    Path(
        "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/libdlfind/src/dl-find.F90"
    ),
    Path(
        "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/libdlfind/src/dlf_util.F90"
    ),
)
FORCE_SOURCE_REFERENCES = (
    {
        "file": str(FORCE_SOURCE_FILES[0]),
        "lines": "9, 17-19, 1519-1537",
        "detail": "get_active_gradients returns Eh/bohr",
    },
    {
        "file": str(FORCE_SOURCE_FILES[1]),
        "lines": "801",
        "detail": "writes negative xgradient to path_force.xyz",
    },
    {
        "file": str(FORCE_SOURCE_FILES[2]),
        "lines": "806-827",
        "detail": "path_force XYZ output scales by ang_au at line 821",
    },
)


def source_sha256(path: Path) -> str:
    """Return a source file SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_native_input(source: str) -> str:
    """Create a one-cycle 0 K gradient input while preserving the Hamiltonian."""
    prefix, separator, _ = source.partition("&dlfind")
    if not separator:
        raise ValueError("Source input lacks &dlfind")
    prefix, count = re.subn(r"\bmaxcyc\s*=\s*\d+", "maxcyc=1", prefix, flags=re.I)
    if count != 1:
        raise ValueError("Expected one Amber maxcyc in source input")
    prefix, count = re.subn(
        r"\bdftb_telec\s*=\s*[-+\d.EeDd]+", "dftb_telec=0.0", prefix, flags=re.I
    )
    if count != 1:
        raise ValueError("Expected one DFTB electronic temperature setting")
    active = re.findall(r"\bactive\s*=\s*'([^']+)'", source, re.I)
    if len(active) != 1:
        raise ValueError("Expected exactly one quoted active atom mask")
    dlfind = (
        "&dlfind\n optalg='LBFGS', crdrep='CART', dimer=.false.,\n"
        f" active='{active[0]}', maxcycle=1, tol=1.d-4, tole=1.d-5,\n"
        " maxstep=0.1, verbosity=2,\n/\n"
    )
    result = prefix + dlfind
    keys = ("qmmask", "qmcharge", "qm_theory", "qm_ewald", "qm_pme", "qmcut", "active")
    for key in keys:
        left = re.findall(rf"\b{key}\s*=\s*([^,\n]+)", source, re.I)
        right = re.findall(rf"\b{key}\s*=\s*([^,\n]+)", result, re.I)
        if len(left) != 1 or len(right) != 1 or left[0].strip() != right[0].strip():
            raise ValueError(f"Native input changed or lost required setting {key}")
    return result


class NativeSanderEvaluator:
    """Evaluate one exact geometry with an isolated native-Sander CLI call."""

    def __init__(
        self,
        topology: Path,
        source_restart: Rst7,
        active: np.ndarray,
        mdin_text: str,
        executable: Path,
        branch_dir: Path,
        force_to_gradient_scale: float,
    ) -> None:
        self.topology = topology
        self.reference_restart = source_restart
        self.active = active
        self.mdin_text = mdin_text
        self.executable = executable
        self.branch_dir = branch_dir
        self.force_scale = force_to_gradient_scale
        previous = []
        for path in branch_dir.glob("eval_*.mdin"):
            match = re.fullmatch(r"eval_(\d+)\.mdin", path.name)
            if match:
                previous.append(int(match.group(1)))
        self.counter = max(previous, default=0)
        self.evidence: list[dict[str, Any]] = []

    def __call__(self, coordinates: FloatArray) -> tuple[float, FloatArray]:
        """Return native DL-FIND energy and calibrated full-system gradient."""
        self.counter += 1
        tag = f"eval_{self.counter:05d}"
        input_restart = self.branch_dir / f"{tag}.input.ncrst"
        mdin = self.branch_dir / f"{tag}.mdin"
        mdout = self.branch_dir / f"{tag}.mdout"
        output_restart = self.branch_dir / f"{tag}.output.ncrst"
        mdinfo = self.branch_dir / f"{tag}.mdinfo"
        log = self.branch_dir / f"{tag}.log"
        tagged_outputs = [
            input_restart,
            mdin,
            mdout,
            output_restart,
            mdinfo,
            log,
            *(
                self.branch_dir / f"{tag}.{suffix}"
                for suffix in ("path_force.xyz", "path.xyz", "path_active.xyz")
            ),
        ]
        transient_outputs = [
            self.branch_dir / "path_force.xyz",
            self.branch_dir / "path.xyz",
            self.branch_dir / "path_active.xyz",
        ]
        collisions = [
            path for path in (*tagged_outputs, *transient_outputs) if path.exists()
        ]
        if collisions:
            raise FileExistsError(
                "Refusing to overwrite native evaluation artifacts: "
                + ", ".join(str(path) for path in collisions)
            )
        write_restart(input_restart, coordinates, self.reference_restart)
        mdin.write_text(self.mdin_text, encoding="utf-8")
        command = [
            str(self.executable),
            "-O",
            "-i",
            str(mdin),
            "-p",
            str(self.topology),
            "-c",
            str(input_restart),
            "-o",
            str(mdout),
            "-r",
            str(output_restart),
            "-inf",
            str(mdinfo),
        ]
        with log.open("w", encoding="utf-8") as stream:
            subprocess.run(
                command,
                cwd=self.branch_dir,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=True,
            )
        generated: dict[str, Path] = {}
        for suffix in ("path_force.xyz", "path.xyz", "path_active.xyz"):
            source = self.branch_dir / suffix
            if not source.is_file():
                raise FileNotFoundError(f"Native Sander did not write {suffix}")
            destination = self.branch_dir / f"{tag}.{suffix}"
            source.replace(destination)
            generated[suffix] = destination
        force_path = generated["path_force.xyz"]
        force_frames = read_xyz_frames(force_path, int(np.count_nonzero(self.active)))
        gradient = np.zeros_like(coordinates, dtype=float)
        gradient[self.active] = self.force_scale * force_frames[-1]
        energy_text = mdout.read_text(encoding="utf-8", errors="replace")
        matches = re.findall(
            r"Energy calculation finished, energy:\s*([-+\d.EeDd]+)", energy_text
        )
        if not matches:
            raise ValueError(f"No DL-FIND energy found in {mdout}")
        energy = float(matches[-1].replace("D", "E").replace("d", "e"))
        if not math.isfinite(energy) or not np.isfinite(gradient).all():
            raise ValueError("Native Sander returned nonfinite energy or gradient")
        evidence_paths = [input_restart, mdin, mdout, log, *generated.values()]
        self.evidence.append(
            {
                "evaluation": self.counter,
                "energy_hartree": energy,
                "gradient_sha256": hashlib.sha256(gradient.tobytes()).hexdigest(),
                "artifacts_sha256": {
                    path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in evidence_paths
                },
            }
        )
        return energy, gradient


def load_calibration(path: Path) -> tuple[float, dict[str, Any]]:
    """Require a passing finite-difference check of the fixed Amber conversion."""
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("status") != "validated" or record.get("passed") is not True:
        raise ValueError("Force calibration record is not validated and passing")
    observed = float(record["active_gradient"]["observed_force_to_gradient_scale"])
    if not math.isfinite(observed) or not math.isclose(
        observed, PATH_FORCE_TO_GRADIENT_SCALE, rel_tol=0.01, abs_tol=0.0
    ):
        raise ValueError(
            "Finite-difference scale disagrees with Amber path_force source conversion"
        )
    if not record.get("validation_sha256"):
        raise ValueError("Calibration record must identify its validation evidence")
    recorded_sources = {
        item["path"]: item["sha256"]
        for item in record["active_gradient"]["writer_source"]["source_hashes"]
    }
    expected_sources = {
        str(path.resolve()): source_sha256(path) for path in FORCE_SOURCE_FILES
    }
    if recorded_sources != expected_sources:
        raise ValueError(
            "Finite-difference evidence uses different AmberTools source files"
        )
    return PATH_FORCE_TO_GRADIENT_SCALE, record


def read_distance_pairs(path: Path, natom: int) -> list[dict[str, Any]]:
    """Load and validate requested 1-based full-topology distance pairs."""
    pairs = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(pairs, list) or not pairs:
        raise ValueError("Distance-pair JSON must be a nonempty list")
    checked = []
    for pair in pairs:
        i, j = int(pair["atom_i_1based"]) - 1, int(pair["atom_j_1based"]) - 1
        if (
            not pair.get("label")
            or i < 0
            or j < 0
            or i >= natom
            or j >= natom
            or i == j
        ):
            raise ValueError(f"Invalid distance pair: {pair}")
        checked.append({"label": str(pair["label"]), "i": i, "j": j})
    return checked
