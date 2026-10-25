"""Fresh-process, fixed-box SANDER API energy/force consistency checks."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm
from parmed.amber.mask import AmberMask


def _record(path: Path) -> dict[str, Any]:
    before = path.stat()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Source changed while hashing: {path}")
    return {"path": str(path), "sha256": digest, "bytes": after.st_size}


def _write(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _validate(config: dict[str, Any]) -> None:
    for key in ("python_executable", "topology", "coordinates", "parameter_files"):
        if key not in config:
            raise ValueError(f"Missing engine prerequisite: {key}")
    env = config["environment"]
    if not {"AMBERHOME", "LD_LIBRARY_PATH", "OMP_NUM_THREADS"} <= set(env):
        raise ValueError("Explicit Amber environment required")
    if env["OMP_NUM_THREADS"] != "1":
        raise ValueError("Only serial one-thread checks are supported")
    mm, qm = config["mm_options"], config["qm_options"]
    for key, value in {
        "ntb": 1,
        "igb": 0,
        "ntc": 1,
        "ntf": 1,
        "ntr": 0,
        "ifqnt": 1,
    }.items():
        if mm.get(key) != value:
            raise ValueError(f"Explicit periodic unconstrained {key}={value} required")
    if not math.isfinite(mm["cut"]) or mm["cut"] <= 0:
        raise ValueError("Positive explicit MM cutoff required")
    required = {
        "qm_theory",
        "qmcharge",
        "qmmask",
        "adjust_q",
        "qm_ewald",
        "qm_pme",
        "scfconv",
        "qmmm_int",
        "qmcut",
        "lnk_dis",
        "lnk_atomic_no",
        "lnk_method",
    }
    if not required <= set(qm) or qm["qm_theory"] != "DFTB3":
        raise ValueError("Explicit DFTB3 charge, mask and embedding controls required")
    if qm["qmmm_int"] != 1 or qm["qm_ewald"] != 1 or qm["qm_pme"] != 1:
        raise ValueError("Electronic embedding with QM Ewald/PME required")
    if (
        qm["qmcut"] != mm["cut"]
        or qm["lnk_dis"] != 1.09
        or qm["lnk_atomic_no"] != 1
        or qm["lnk_method"] != 1
    ):
        raise ValueError(
            "Explicit matched QM cutoff and standard 1.09A hydrogen links required"
        )
    if not 0 < qm["scfconv"] <= 1e-10 or not math.isfinite(qm["scfconv"]):
        raise ValueError("SCF convergence must be finite and at most 1e-10")
    for key in (
        "force_absolute_tolerance",
        "force_relative_tolerance",
        "energy_reproducibility_tolerance",
    ):
        if not math.isfinite(config[key]) or config[key] < 0:
            raise ValueError(f"Invalid acceptance tolerance: {key}")
    steps = config["step_sizes_angstrom"]
    if (
        len(steps) < 2
        or len(set(steps)) != len(steps)
        or any(not math.isfinite(step) or step <= 0 for step in steps)
    ):
        raise ValueError(
            "At least two distinct positive finite displacement steps required"
        )
    if not config["components"]:
        raise ValueError("Representative force components required")
    if not {"real_atoms", "link_atoms", "electron_count", "boundary_pairs"} <= set(
        config["expected_inventory"]
    ):
        raise ValueError(
            "Explicit expected QM atom, link, electron and boundary inventory required"
        )
    for row in config["components"]:
        if (
            type(row["atom"]) is not int
            or row["atom"] < 1
            or row["axis"] not in (0, 1, 2)
        ):
            raise ValueError("Invalid one-based atom index or Cartesian axis")


def _options(sander: Any, config: dict[str, Any]) -> tuple[Any, Any]:
    mm, qm = sander.pme_input(), sander.QmInputOptions()
    for options, declared in ((mm, config["mm_options"]), (qm, config["qm_options"])):
        for key, value in declared.items():
            if not hasattr(options, key):
                raise ValueError(f"Unsupported installed API option: {key}")
            setattr(options, key, value)
    return mm, qm


def evaluate_request(request_path: Path) -> dict[str, Any]:
    """Evaluate one unchanged or perturbed snapshot with the installed API.

    Args:
        request_path: JSON request with explicit configuration and displacement.

    Returns:
        Full-precision energy, requested forces, fixed box and verified inventory.

    Raises:
        ValueError: Invalid options, inventory, coordinates, or numeric results.
    """
    request = json.loads(request_path.read_text())
    config = request["config"]
    _validate(config)
    structure = AmberParm(config["topology"], xyz=config["coordinates"])
    xyz = np.array(structure.coordinates, dtype=float, copy=True)
    box = np.array(structure.box, dtype=float, copy=True)
    if structure.get_coordinates().shape[0] != 1 or xyz.shape != (
        len(structure.atoms),
        3,
    ):
        raise ValueError("Exactly one complete coordinate frame required")
    if box.shape != (6,) or not np.isfinite(box).all() or not np.isfinite(xyz).all():
        raise ValueError("Finite coordinates and full periodic box required")
    selected = {
        i
        for i, flag in enumerate(
            AmberMask(structure, config["qm_options"]["qmmask"]).Selection()
        )
        if flag
    }
    links = [
        (a.idx, b.idx)
        for bond in structure.bonds
        for a, b in ((bond.atom1, bond.atom2), (bond.atom2, bond.atom1))
        if a.idx in selected and b.idx not in selected
    ]
    inventory = {
        "real_atoms": len(selected),
        "link_atoms": len(links),
        "electron_count": sum(structure.atoms[i].atomic_number for i in selected)
        + len(links)
        - config["qm_options"]["qmcharge"],
        "boundary_pairs": [[a + 1, b + 1] for a, b in links],
        "full_mm_partial_charge": math.fsum(atom.charge for atom in structure.atoms),
        "qm_mm_partial_charge": math.fsum(structure.atoms[i].charge for i in selected),
    }
    for key, expected in config["expected_inventory"].items():
        observed = inventory[key]
        if key == "boundary_pairs":
            observed, expected = sorted(observed), sorted(expected)
        if observed != expected:
            raise ValueError(
                f"QM inventory mismatch: {key}: {inventory[key]} != {expected}"
            )
    if len({a for a, _ in links}) != len(links):
        raise ValueError("Multiple link atoms on one QM parent are unsupported")
    perturbation = request["perturbation"]
    if perturbation is not None:
        xyz[perturbation["atom"] - 1, perturbation["axis"]] += perturbation["delta"]
    sander = importlib.import_module("sander")
    mm, qm = _options(sander, config)
    with sander.setup(config["topology"], xyz, box, mm, qm) as context:
        energy, forces = context.energy_forces()
        forces = np.asarray(forces, dtype=float).reshape(xyz.shape)
        if not np.isfinite(forces).all() or not math.isfinite(float(energy.tot)):
            raise ValueError("Nonfinite engine energy or analytic forces")
        if not np.array_equal(np.asarray(context.box), box):
            raise ValueError("Engine altered the fixed periodic box")
        np.save(request_path.parent / "analytic_forces.npy", forces)
        result = {
            "energy_kcal_per_mol": float(energy.tot),
            "inventory": inventory,
            "box": box.tolist(),
            "perturbation": perturbation,
            "forces": [
                float(forces[row["atom"] - 1, row["axis"]])
                for row in config["components"]
            ],
            "api_sources": [
                _record(Path(sander.__file__)),
                _record(Path(sander._pys.__file__)),
            ],
            "realized_options": {
                name: {
                    key: getattr(options, key)
                    for key in dir(options)
                    if not key.startswith("_")
                }
                for name, options in (("mm", mm), ("qm", qm))
            },
        }
    _write(request_path.parent / "result.json", result)
    return result


def check_es_dftb_forces(config_path: Path, output_directory: Path) -> dict[str, Any]:
    """Check fresh-process forces against central energy differences.

    Args:
        config_path: Approved explicit JSON engine/inventory configuration.
        output_directory: Fresh retained evidence directory.

    Returns:
        Numeric comparisons and separate scientific acceptance.
    """
    from nyl12_covalent.es_dftb_checks import check_es_dftb_forces as run_check

    return run_check(config_path, output_directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", required=True, type=Path)
    evaluate_request(parser.parse_args().worker)
