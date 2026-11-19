"""Test whether native NFE PMD bias enters the existing DLFind minimizer."""

from __future__ import annotations

import hashlib
import json
import logging
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import parmed as pmd
from parmed.amber.mask import AmberMask

ROOT = Path("/data/jpc/covalent-enzyme-design")
CASE = Path(__file__).resolve().parent / "attempt_2"
BASE = ROOT / "qmmm/pa6_flipped_model30/07_ash274_comparison"
SOURCE = BASE / "06_qmmm_optimization/replica_1_tier_A_seed/water_qm_reoptimization"
AMBER = Path("/data/jpc/code/ambertools26/venv")
NATIVE = Path(
    "/data/jpc/code/ambertools26/build_assets/src/ambertools26_src/AmberTools/src/sander"
)
LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return the SHA-256 identity of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name: str, data: dict) -> None:
    """Save an authoritative JSON record."""
    (CASE / name).write_text(json.dumps(data, indent=2) + "\n")


def main() -> None:
    """Prepare immutable inputs, run two evaluations, and assess bias support."""
    snapshot = CASE / "source_snapshot"
    snapshot.mkdir(parents=True, exist_ok=False)
    top_source = BASE / "00_system/ash274_solvated.prmtop"
    paths = [
        top_source,
        SOURCE / "optimized.ncrst",
        SOURCE / "minimize.mdin",
        SOURCE / "acceptance.json",
        Path(__file__),
        ROOT / "scripts/run_nyl12_job.py",
        ROOT / "src/nyl12_covalent/jobs.py",
    ]
    paths += [
        NATIVE / name
        for name in [
            "sander.F90",
            "force.F90",
            "nfe-pmd-hooks.F90",
            "nfe-cv-N_OF_BONDS.F90",
        ]
    ]
    for path in paths:
        shutil.copy2(path, snapshot / path.name)
    top = snapshot / top_source.name
    coordinates = snapshot / "optimized.ncrst"
    structure = pmd.load_file(str(top), xyz=str(coordinates))
    text = (SOURCE / "minimize.mdin").read_text()
    mask = re.search(r"qmmask='([^']+)'", text).group(1)
    qm = list(AmberMask(structure, mask).Selected())
    hydrogens = [i + 1 for i in qm if structure.atoms[i].atomic_number == 1]
    xyz = np.asarray(structure.coordinates)
    box = np.asarray(structure.box[:3])

    def count(center: int) -> float:
        delta = xyz[np.asarray(hydrogens) - 1] - xyz[center - 1]
        delta -= box * np.rint(delta / box)
        return float(np.sum(1 / (1 + (np.linalg.norm(delta, axis=1) / 1.25) ** 6)))

    cn = count(13905)
    co = count(13916)
    target = cn + 0.20
    pairs = ",".join(f"13905,{h}" for h in hydrogens)
    colvar = (
        "&colvar\n cv_type='N_OF_BONDS',\n"
        f" cv_ni={2 * len(hydrogens)}, cv_i={pairs},\n"
        f" cv_nr=1, cv_r=1.25, anchor_position=0,{target:.12f},{target:.12f},6,\n"
        " anchor_strength=10,10,\n/\n"
    )
    (CASE / "occupancy.colvar").write_text(colvar)
    text = re.sub(r"maxcyc=2000", "maxcyc=1", text)
    text = re.sub(r"maxcycle=2000", "maxcycle=1", text)
    text = text.replace("nmropt=0,", "nmropt=0, infe=INFE,")
    protocol = {
        "hypothesis": "native NFE PMD contributes nonzero bias in DLFind",
        "qm_real_atoms": len(qm),
        "qmcharge": -1,
        "active_atoms": 4662,
        "hydrogens_1based": hydrogens,
        "r0_A": 1.25,
        "CN_N": cn,
        "CN_OG1": co,
        "occupancy_difference_N_minus_OG1": cn - co,
        "CN_N_target": target,
        "expected_bias_kcal_mol": 0.5 * 10 * 0.2**2,
        "native_limitation": "N and O counts can be restrained separately; native PMD cannot combine N-O",
        "future_coordinates": [
            "CN_N(all real QM H)",
            "CN_OG1(all real QM H)",
            "OG1-substrate C distance",
        ],
        "protonation": "ASP272/ASH274(HD2-OD2), neutral Thr233 NH2; unchanged",
        "sampling": "two one-cycle energy/gradient evaluations; identical optimized source; no derivatives",
        "source_hashes": {str(p): digest(p) for p in paths},
        "executable_hashes": {
            str(p): digest(p)
            for p in [
                AMBER / "bin/sander",
                AMBER / "lib/libdlfind.so",
                AMBER / "lib/libsander.so",
            ]
        },
    }
    save("protocol.json", protocol)
    energies = {}
    results = {}
    for label, infe in [("baseline", 0), ("native_pmd", 1)]:
        work = CASE / label
        work.mkdir()
        shutil.copy2(CASE / "occupancy.colvar", work / "occupancy.colvar")
        mdin = text.replace("INFE", str(infe))
        if infe:
            mdin += "&pmd cv_file='occupancy.colvar', output_file='occupancy.dat', output_freq=1 /\n"
        (work / "evaluate.mdin").write_text(mdin)
        command = [
            str(AMBER / "bin/sander"),
            "-O",
            "-i",
            "evaluate.mdin",
            "-o",
            "evaluate.mdout",
            "-p",
            str(top),
            "-c",
            str(coordinates),
            "-r",
            "evaluated.ncrst",
            "-inf",
            "evaluate.mdinfo",
        ]
        (work / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        LOGGER.info("Running %s one-cycle evaluation", label)
        with (work / "stdout.log").open("w") as stream:
            result = subprocess.run(
                command, cwd=work, stdout=stream, stderr=subprocess.STDOUT, check=False
            )
        output = (work / "evaluate.mdout").read_text()
        values = re.findall(r"EPtot\s*=\s*([-\d.E+]+)", output)
        energies[label] = float(values[0]) if values else None
        results[label] = {
            "returncode": result.returncode,
            "EPtot": energies[label],
            "pmd_output_exists": (work / "occupancy.dat").exists(),
            "pmd_initialization_printed": "N_OF_BONDS" in output,
        }
    delta = (
        None
        if None in energies.values()
        else energies["native_pmd"] - energies["baseline"]
    )
    active = (
        delta is not None
        and abs(delta - 0.2) < 0.01
        and results["native_pmd"]["pmd_output_exists"]
    )
    save(
        "process.json",
        {
            "evaluations": results,
            "all_returncodes_zero": all(v["returncode"] == 0 for v in results.values()),
        },
    )
    save(
        "scientific_acceptance.json",
        {
            "accepted": active,
            "energy_delta_kcal_mol": delta,
            "expected_bias_kcal_mol": 0.2,
            "reason": "native_bias_verified"
            if active
            else "unsupported_current_DLFind_hook",
            "limitation": "One-cycle DLFind may internally perform more than one force evaluation; compare first EPtot only.",
        },
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
