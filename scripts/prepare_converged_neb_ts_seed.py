"""Prepare a full-shell dimer TS seed from the converged direct NEB maximum."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import re
import shutil

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from nyl12_covalent.neb_dimer_seed import local_heavy_atom_clashes

LOGGER = logging.getLogger(__name__)


def _digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _option(text: str, name: str) -> str:
    values = re.findall(rf"\b{name}\s*=\s*'([^']*)'", text, re.I)
    if len(values) != 1:
        raise ValueError(f"Expected one DL-FIND {name} option")
    return values[0]


def _coordinates(path: Path, atom_count: int) -> np.ndarray:
    xyz = np.asarray(Rst7.open(str(path)).coordinates, dtype=float)
    if xyz.shape != (1, atom_count, 3) or not np.isfinite(xyz).all():
        raise ValueError(f"Invalid full-system coordinates: {path}")
    return xyz[0]


def prepare(neb_dir: Path, output: Path) -> dict[str, object]:
    """Write an auditable dimer seed at the accepted NEB maximum.

    Args:
        neb_dir: Converged 25-image direct-route NEB directory.
        output: New attempt directory, which must not exist.

    Returns:
        Input and chemistry manifest for the prepared TS attempt.

    Raises:
        ValueError: NEB acceptance, atom mapping, or input identity fails.
        FileExistsError: The attempt directory already exists.
    """
    neb_dir = neb_dir.resolve(strict=True)
    output = output.resolve()
    if output.exists():
        raise FileExistsError(output)
    bundle = neb_dir / "vmd_final_24img/final_24img"
    acceptance_path = bundle / "scientific_acceptance.json"
    acceptance = json.loads(acceptance_path.read_text())
    if not (
        acceptance["native_optimizer_converged"]
        and acceptance["path_identity_accepted"]
        and acceptance["chemical_endpoint_identity_accepted"]
        and acceptance["peak_image"] == 8
        and acceptance["image_count"] == 24
        and acceptance["qm_atoms"] == 98
        and acceptance["qm_charge"] == -2
        and acceptance["electronic_temperature_K"] == 0.0
    ):
        raise ValueError("NEB chemical/optimizer acceptance differs from plan")
    topology_path = neb_dir / "inputs/system.prmtop"
    source_mdin = neb_dir / "calculation/minimize.mdin"
    reference_path = neb_dir / "inputs/images/image_00.ncrst"
    image_paths = {
        image: bundle / f"source_restarts/optimized.{image:04d}.ncrst"
        for image in (7, 8, 9)
    }
    if _digest(topology_path) != acceptance["topology_sha256"]:
        raise ValueError("NEB topology changed after acceptance")
    mdin = source_mdin.read_text()
    if not re.search(r"\bdftb_telec\s*=\s*0\.0", mdin):
        raise ValueError("NEB Hamiltonian is not the accepted 0 K setup")
    active_mask, rotation_mask = (
        _option(mdin, "active"),
        _option(mdin, "wtmask"),
    )
    topology = AmberParm(str(topology_path))
    frames = {
        image: _coordinates(path, len(topology.atoms))
        for image, path in image_paths.items()
    }
    active = np.flatnonzero(AmberMask(topology, active_mask).Selection())
    rotate = np.flatnonzero(AmberMask(topology, rotation_mask).Selection())
    if len(active) != 4165 or not set(rotate).issubset(set(active)):
        raise ValueError("Active/rotation mask differs from accepted 12 Å shell")
    mapping = acceptance["atom_mapping_zero_based"]
    required = {
        "Thr233_OG1": (929, "OG1"),
        "Thr233_HG1": (929, "HG1"),
        "scissile_C": (1314, "C"),
        "scissile_O": (1314, "O"),
        "leaving_N": (1315, "N"),
    }
    for name, (residue, atom_name) in required.items():
        index = int(mapping[name])
        atom = topology.atoms[index]
        if atom.name != atom_name or atom.residue.idx + 1 != residue:
            raise ValueError(f"Reactive atom identity mismatch: {name}")
    clashes = local_heavy_atom_clashes(
        topology,
        frames[8],
        {int(mapping[name]) for name in required},
    )
    if clashes["severe_contact_count"]:
        raise ValueError("NEB maximum has severe local heavy-atom clashes")
    tangent = frames[9] - frames[7]
    norm = float(np.linalg.norm(tangent[rotate]))
    if not np.isfinite(norm) or norm <= 1e-8:
        raise ValueError("Images 7–9 have no usable reactive tangent")
    mode = frames[8] + tangent * (0.01 / norm)
    output.mkdir(parents=True)
    inputs, calculation = output / "inputs", output / "calculation"
    inputs.mkdir()
    calculation.mkdir()
    sources = {
        "system.prmtop": topology_path,
        "reference.ncrst": reference_path,
        "source_neb.mdin": source_mdin,
        "scientific_acceptance.json": acceptance_path,
        **{f"image_{image:02d}.ncrst": path for image, path in image_paths.items()},
        "prepare_converged_neb_ts_seed.py": Path(__file__).resolve(),
    }
    hashes = {}
    for name, source in sources.items():
        destination = inputs / name
        shutil.copy2(source, destination)
        hashes[name] = _digest(destination)
        if hashes[name] != _digest(source):
            raise RuntimeError(f"Input changed during snapshot: {source}")
    for name, xyz in (("start.ncrst", frames[8]), ("mode.ncrst", mode)):
        restart = Rst7(natom=len(xyz))
        restart.coordinates = xyz
        restart.box = Rst7.open(str(reference_path)).box
        restart.write(str(inputs / name), netcdf=True)
        hashes[name] = _digest(inputs / name)
    prefix = mdin[: mdin.lower().index("&dlfind")]
    dimer = (
        "&dlfind\n"
        " dimer=.true., dimer_mode=1, dimer_tolrot=2., dimer_maxrot=20,\n"
        " optalg='LBFGS', crdrep='CART',\n"
        f" refcrd='{inputs / 'reference.ncrst'}',\n"
        f" coords2(1)='{inputs / 'mode.ncrst'}',\n"
        f" wtmask='{rotation_mask}', active='{active_mask}',\n"
        " tol=1.d-4, tole=1.d-5, maxcycle=1000, maxstep=0.05, verbosity=2,\n/\n"
    )
    (calculation / "dimer.mdin").write_text(prefix + dimer)
    hashes["dimer.mdin"] = _digest(calculation / "dimer.mdin")
    manifest: dict[str, object] = {
        "status": "prepared only; no TS optimization or Hessian run",
        "source_neb": str(neb_dir),
        "peak_image": 8,
        "source_images": [7, 8, 9],
        "mode": "central path tangent (image 9 minus image 7), 0.01 A on wtmask",
        "tangent_norm_A": norm,
        "active_atom_count": len(active),
        "rotation_atom_count": len(rotate),
        "qm_charge": -2,
        "qm_atoms": 98,
        "local_clashes": clashes,
        "input_sha256": hashes,
    }
    (output / "preparation.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Prepared 0 K dimer seed from NEB image 8 in %s", output)
    return manifest


def main() -> None:
    """Expose the audited peak-image seed preparation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neb-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.neb_dir, args.output)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
