"""Prepare an atom-preserving native DFTB3/MM acylenzyme NEB campaign."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.acylenzyme_seed import build_acylenzyme_seed  # noqa: E402

LOGGER = logging.getLogger(__name__)
MASK = "(:113,186,929,930,968,970,1313-1316)<:12.0"


def save_restart(path: Path, xyz: np.ndarray, reference: Rst7) -> None:
    """Save coordinates with the original periodic cell and atom order."""
    restart = Rst7(natom=len(xyz))
    restart.coordinates = xyz
    restart.box = reference.box
    restart.write(str(path), netcdf=True)


def prepare(project: Path) -> Path:
    """Build and provenance-pin inputs while matching the frozen exterior."""
    campaign = project / "qmmm/pa6_resp_gaff2_wt"
    paths = campaign / "11_frame1603_12A_reaction_path_20261109/reaction_paths"
    old = paths / "02_asynchronous_es_to_ti1"
    output = paths / "03_acylenzyme_neb"
    inputs = output / "inputs"
    es = (
        campaign
        / "10_frame1603_nearpass_qmmm_20261106/shell_12A_attempt1/calculation/dftb3_mm"
    )
    topology = Path(
        "/data/jpc/covalent-enzyme-design-pa6-md/models/pa6_modular_md_20261031/campaign_attempt2/systems/wt/md_solvated.prmtop"
    )
    sources = {
        "reactant": es / "optimized.ncrst",
        "proton_transferred": old
        / "01_coupled_coordinate_drive/asynchronous_24/optimized.ncrst",
        "attack": old
        / "03_corrected_attack_from_window24/calculation/windows/attack_03/optimized.ncrst",
        "ti1_like": old
        / "03_corrected_attack_from_window24/calculation/windows/attack_09/optimized.ncrst",
    }
    for path in [topology, es / "minimize.mdin", *sources.values()]:
        if not path.is_file():
            raise FileNotFoundError(path)
    reference = Rst7.open(str(sources["reactant"]))
    xyz = reference.coordinates.reshape(-1, 3)
    parm = AmberParm(str(topology), xyz=str(sources["reactant"]))
    active = np.asarray(AmberMask(parm, MASK).Selection(), dtype=bool)
    seed_reference = Rst7.open(str(sources["ti1_like"]))
    seed, report = build_acylenzyme_seed(
        parm,
        seed_reference.coordinates.reshape(-1, 3),
        ester_distance=1.65,
        leaving_distance=2.3,
        orient_leaving=True,
    )
    moved = (
        np.linalg.norm(seed - seed_reference.coordinates.reshape(-1, 3), axis=1) > 1e-8
    )
    if np.any(moved & ~active):
        raise ValueError("Seed builder moved atoms outside the declared shell")
    if output.exists():
        raise FileExistsError(output)
    inputs.mkdir(parents=True)
    mapped = {}
    adjustments = {}
    for name, source in sources.items():
        rst = Rst7.open(str(source))
        coords = rst.coordinates.reshape(-1, 3).copy()
        if coords.shape != xyz.shape or not np.allclose(rst.box, reference.box):
            raise ValueError("Endpoint atom count or periodic cell mismatch")
        adjustments[name] = float(
            np.max(np.linalg.norm(coords[~active] - xyz[~active], axis=1))
        )
        coords[~active] = xyz[~active]
        target = inputs / f"{name}.ncrst"
        save_restart(target, coords, reference)
        mapped[name] = str(target)
    seed[~active] = xyz[~active]
    seed_path = inputs / "acylenzyme_seed.ncrst"
    save_restart(seed_path, seed, reference)
    parm.coordinates = seed
    parm.save(str(inputs / "acylenzyme_seed.pdb"), overwrite=False)
    shutil.copy2(es / "minimize.mdin", inputs / "base.mdin")
    config = {
        "topology": str(topology),
        "sander": "/data/jpc/code/ambertools26/venv/bin/sander",
        "base_input": str(inputs / "base.mdin"),
        "reactant": mapped["reactant"],
        "seed": str(seed_path),
        "intermediates": [
            mapped[n] for n in ("proton_transferred", "attack", "ti1_like")
        ],
        "method": "DFTB3/3OB periodic QM/MM",
        "qm_charge": -2,
        "active_reference": mapped["reactant"],
        "active_mask": MASK,
        "active_atoms": int(active.sum()),
        "frozen_exterior_reset_max_A": adjustments,
        "endpoint_policy": "FREE in NEB; chemistry-gated unbiased release first",
        "seed_report": report,
        "classical_acylenzyme_topology_valid": False,
        "source_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [topology, es / "minimize.mdin", *sources.values()]
        },
    }
    (output / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    LOGGER.info("Prepared %s with %d movable atoms", output, active.sum())
    return output


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    prepare(ROOT)
