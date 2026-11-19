"""Prepare an interpolated local dimer seed from adjacent final NEB images."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

from nyl12_covalent.neb_dimer_seed import (
    build_local_dimer_seed,
    final_neb_rows,
    local_heavy_atom_clashes,
    validate_reactive_atom_mapping,
)
from nyl12_covalent.neb_xyz_snapshot import reconstruct_xyz_band

REQUIRED_ATOMS: dict[int, dict[str, int | str]] = {
    13916: {"name": "OG1", "residue": 929, "atomic_number": 8},
    13917: {"name": "HG1", "residue": 929, "atomic_number": 1},
    19916: {"name": "C", "residue": 1314, "atomic_number": 6},
    19917: {"name": "O", "residue": 1314, "atomic_number": 8},
    19929: {"name": "N", "residue": 1315, "atomic_number": 7},
}


def digest(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def option(text: str, name: str) -> str:
    """Extract exactly one quoted DL-FIND option from the Amber input."""
    values = re.findall(rf"\b{name}\s*=\s*'([^']*)'", text, re.I)
    if len(values) != 1:
        raise ValueError(f"Expected unique {name} option")
    return values[0]


def _geometry(
    coordinates: np.ndarray, atom_indices: dict[int, int]
) -> dict[str, float]:
    """Measure direct acylation distances from the mapped topology atoms."""
    pairs = {
        "OG1_C": (13916, 19916),
        "OG1_HG1": (13916, 13917),
        "HG1_Nsub": (13917, 19929),
        "C_N": (19916, 19929),
        "C_O": (19916, 19917),
    }
    return {
        name: float(
            np.linalg.norm(
                coordinates[atom_indices[first]] - coordinates[atom_indices[second]]
            )
        )
        for name, (first, second) in pairs.items()
    }


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    """Validate final native inputs and write an immutable local dimer seed.

    Args:
        args: Parsed CLI options including NEB directory, map, image pair,
            output location, cycle, and mode displacement.

    Returns:
        The preparation metadata written to the output directory.

    Raises:
        FileExistsError: If the output path already exists.
        ValueError: If source identity, band closure, mapping, or geometry fails.
    """
    neb_dir = args.neb_dir.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    map_path = args.atom_map.resolve(strict=True)
    config_path = neb_dir / "config.json"
    config = json.loads(config_path.read_text())
    calculation = neb_dir / "calculation"
    mdin_path = calculation / "minimize.mdin"
    mdout_path = calculation / "minimized.mdout"
    xyz_path = calculation / "nebpath.xyz"
    top_path = neb_dir / "inputs/system.prmtop"
    mdin, mdout = mdin_path.read_text(), mdout_path.read_text()
    active_mask, rotation_mask = option(mdin, "active"), option(mdin, "wtmask")
    reference_path = Path(option(mdin, "refcrd"))
    if not reference_path.is_absolute():
        reference_path = (calculation / reference_path).resolve(strict=True)
    else:
        reference_path = reference_path.resolve(strict=True)
    if config.get("topology") != str(top_path.resolve()):
        raise ValueError("NEB topology path differs from the recorded config")
    if config.get("active_atoms") != 4165 or config.get("regular_images") != 16:
        raise ValueError("Unexpected native active atom or image count")
    expected_cycle, image_count = args.expected_cycle, config["regular_images"]
    rows = final_neb_rows(mdout, expected_cycle, image_count)
    energies = [row[1] for row in rows]
    peak = int(np.argmax(energies)) + 1
    if peak != args.right_image:
        raise ValueError("Requested right image is not the final NEB energy maximum")
    topology = AmberParm(str(top_path))
    reference = Rst7.open(str(reference_path))
    active_selection = AmberMask(topology, active_mask).Selection()
    active_indices = [
        index
        for index, (selected, atom) in enumerate(zip(active_selection, topology.atoms))
        if selected and atom.mass > 0.1
    ]
    declared_map = json.loads(map_path.read_text()).get("atom_mapping")
    if not isinstance(declared_map, dict):
        raise ValueError("Atom map file has no atom_mapping object")
    atom_indices = validate_reactive_atom_mapping(
        topology, active_indices, declared_map, REQUIRED_ATOMS
    )
    frames, mapped_active = reconstruct_xyz_band(
        xyz_path.read_bytes(),
        topology,
        reference,
        active_mask,
        expected_images=image_count,
        expected_active=len(active_indices),
    )
    if mapped_active != active_indices:
        raise ValueError("Native XYZ active atom mapping differs from topology mask")
    rotation_selection = AmberMask(topology, rotation_mask).Selection()
    rotation_indices = np.array(
        [
            index
            for index, (selected, atom) in enumerate(
                zip(rotation_selection, topology.atoms)
            )
            if selected and atom.mass > 0.1
        ],
        dtype=int,
    )
    seed, mode, tangent_norm = build_local_dimer_seed(
        frames,
        mapped_active,
        rotation_indices,
        args.left_image,
        args.right_image,
        args.displacement_angstrom,
    )
    reactive_indices = {atom_indices[serial] for serial in REQUIRED_ATOMS}
    geometry = _geometry(seed, atom_indices)
    clash_report = local_heavy_atom_clashes(topology, seed, reactive_indices)
    project = Path(__file__).resolve().parents[1]
    sources = [
        config_path,
        mdin_path,
        mdout_path,
        xyz_path,
        top_path,
        reference_path,
        map_path,
        Path(__file__).resolve(),
        project / "src/nyl12_covalent/neb_dimer_seed.py",
        project / "src/nyl12_covalent/neb_xyz_snapshot.py",
    ]
    source_hashes = {str(path): digest(path) for path in sources}
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    output.mkdir(parents=True, exist_ok=False)
    inputs, work = output / "inputs", output / "calculation"
    inputs.mkdir()
    work.mkdir()
    copies = []
    for index, source in enumerate(sources):
        copied = inputs / f"source_{index:02d}_{source.name}"
        shutil.copy2(source, copied)
        if digest(copied) != source_hashes[str(source)]:
            raise RuntimeError(f"Copied source hash mismatch: {source}")
        copies.append(
            {
                "source": str(source),
                "copy": str(copied),
                "sha256": source_hashes[str(source)],
            }
        )
    shutil.copy2(top_path, inputs / "system.prmtop")
    shutil.copy2(reference_path, inputs / "reference.ncrst")
    for name, coordinates in (("ts_guess", seed), ("ts_mode_guess", mode)):
        restart = Rst7(natom=len(coordinates))
        restart.coordinates, restart.box = coordinates, reference.box
        restart.write(str(inputs / f"{name}.ncrst"), netcdf=True)
    dimer_input = (
        "&dlfind\n dimer=.true., dimer_mode=1, dimer_tolrot=2., dimer_maxrot=20,\n"
        " optalg='LBFGS', crdrep='CART',\n"
        f" refcrd='{inputs / 'reference.ncrst'}',\n"
        f" coords2(1)='{inputs / 'ts_mode_guess.ncrst'}',\n"
        f" wtmask='{rotation_mask}', active='{active_mask}',\n"
        " tol=1.d-4, tole=1.d-5, maxcycle=1000, maxstep=0.1, verbosity=2,\n/\n"
    )
    dlfind_start = mdin.lower().find("&dlfind")
    if dlfind_start < 0:
        raise ValueError("Original Amber input has no DL-FIND namelist")
    dimer_path = work / "dimer.mdin"
    dimer_path.write_text(mdin[:dlfind_start] + dimer_input)
    final_gate = re.search(
        rf"Testing convergence\s+in cycle\s+{expected_cycle}\b\s+Energy\s+"
        r"([-+\d.Ee]+)\s+Target:\s*([-+\d.Ee]+)\s+converged\?\s*(\w+)",
        mdout,
    )
    mode_norm = float(np.linalg.norm((mode - seed)[rotation_indices]))
    metadata: dict[str, Any] = {
        "scientific_status": "Prepared dimer seed only; no optimization or job was launched",
        "source_sha256": source_hashes,
        "input_copies": copies,
        "topology_copy_sha256": digest(inputs / "system.prmtop"),
        "reference_copy_sha256": digest(inputs / "reference.ncrst"),
        "native_cycle": expected_cycle,
        "image_count": image_count,
        "left_image_1based": args.left_image,
        "right_image_1based": args.right_image,
        "right_image_is_energy_maximum": peak == args.right_image,
        "image_energies_hartree": [
            {
                "image": image,
                "energy": energy,
                "relative_to_image1_kcal_mol": (energy - energies[0]) * 627.5095,
                "tangential_force": tangent_force,
                "perpendicular_force": perpendicular_force,
            }
            for image, energy, tangent_force, perpendicular_force in rows
        ],
        "native_convergence": {
            "reported_energy_residual": float(final_gate.group(1))
            if final_gate
            else None,
            "reported_energy_target": float(final_gate.group(2))
            if final_gate
            else None,
            "reported_converged": final_gate.group(3).lower() == "yes"
            if final_gate
            else None,
            "climbing_image_enabled": bool(config.get("climbing_enabled")),
        },
        "atom_mapping_1based": {
            str(serial): {
                **REQUIRED_ATOMS[serial],
                "topology_index_0based": atom_indices[serial],
            }
            for serial in REQUIRED_ATOMS
        },
        "seed_interpolation_fraction_from_left": 0.5,
        "seed_reactive_geometry_angstrom": geometry,
        "mode": {
            "vector": "final image 8 minus final image 7",
            "displacement_angstrom_on_wtmask": mode_norm,
            "tangent_norm_angstrom_on_wtmask": tangent_norm,
        },
        "local_heavy_atom_clash_audit": clash_report,
        "active_mask": active_mask,
        "rotation_mask": rotation_mask,
        "outputs": {
            "seed": str(inputs / "ts_guess.ncrst"),
            "mode": str(inputs / "ts_mode_guess.ncrst"),
            "dimer_input": str(dimer_path),
        },
        "ts_validation": "No Hessian or saddle-point validation performed",
    }
    metadata["prepared_sha256"] = {
        str(path.relative_to(output)): digest(path)
        for path in [*inputs.iterdir(), dimer_path]
    }
    (output / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if source_hashes != {str(path): digest(path) for path in sources}:
        raise RuntimeError("A native source changed during seed preparation")
    return metadata


def main() -> None:
    """Parse arguments and prepare one auditable local dimer seed."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neb-dir", type=Path, required=True)
    parser.add_argument("--atom-map", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-cycle", type=int, default=499)
    parser.add_argument("--left-image", type=int, default=7)
    parser.add_argument("--right-image", type=int, default=8)
    parser.add_argument("--displacement-angstrom", type=float, default=0.01)
    args = parser.parse_args()
    metadata = prepare(args)
    logging.info(
        "Prepared local dimer seed at %s; severe local contacts: %d",
        args.output,
        metadata["local_heavy_atom_clash_audit"]["severe_contact_count"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
