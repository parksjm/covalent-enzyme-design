"""Check Python-sander support for a converged DFTB3/MM dimer snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
from parmed.amber import AmberParm, Rst7
from parmed.amber.mask import AmberMask

LOGGER = logging.getLogger(__name__)
HARTREE_TO_KCAL_MOL = 627.5095
REFERENCE_HARTREE = -718.134029878265
ENERGY_TOLERANCE_KCAL_MOL = 0.05


def _digest(path: Path) -> dict[str, Any]:
    """Return a path, SHA256 and byte count for an immutable input."""
    data = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
    }


def _value(raw: str) -> Any:
    """Parse the scalar values used by the Amber namelists in this input."""
    text = raw.strip().rstrip(",").strip()
    if text.startswith(("'", '"')) and text.endswith(("'", '"')):
        return text[1:-1]
    lowered = text.lower()
    if lowered in {".true.", "true"}:
        return True
    if lowered in {".false.", "false"}:
        return False
    if re.fullmatch(r"[-+]?\d+", text):
        return int(text)
    try:
        return float(text.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise ValueError(f"Unsupported Amber namelist value: {raw!r}") from exc


def _namelist(path: Path, section: str) -> dict[str, Any]:
    """Parse simple key/value entries from one Amber namelist section."""
    source = path.read_text()
    match = re.search(
        rf"&{re.escape(section)}\b(.*?)/", source, re.IGNORECASE | re.DOTALL
    )
    if match is None:
        raise ValueError(f"Missing &{section} namelist in {path}")
    body = re.sub(r"!.*", "", match.group(1))
    entries: dict[str, Any] = {}
    for item in re.finditer(
        r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)(?=,\s*[A-Za-z][A-Za-z0-9_]*\s*=|$)",
        body,
        re.DOTALL,
    ):
        key, raw = item.group(1).lower(), item.group(2).strip().rstrip(",")
        entries[key] = _value(raw)
    return entries


def _active_mask(path: Path) -> str:
    """Read the quoted DL-FIND active mask without truncating absolute paths."""
    source = path.read_text()
    section = re.search(
        r"&dlfind\b(.*?)^\s*/\s*$",
        source,
        re.IGNORECASE | re.DOTALL | re.MULTILINE,
    )
    match = (
        re.search(r"\bactive\s*=\s*'([^']+)'", section.group(1), re.IGNORECASE)
        if section is not None
        else None
    )
    if match is None:
        raise ValueError("Missing DL-FIND active shell mask")
    return match.group(1)


def _validate_coordinates(coordinates: Any, n_atoms: int) -> np.ndarray:
    """Return finite coordinates with one row per topology atom."""
    xyz = np.asarray(coordinates, dtype=float)
    if xyz.size != n_atoms * 3:
        raise ValueError("Restart coordinate count does not match the topology")
    xyz = xyz.reshape((n_atoms, 3))
    if not np.isfinite(xyz).all():
        raise ValueError("Restart contains nonfinite coordinates")
    return xyz


def _apply_options(
    api: Any, mdin: Path, telec: float
) -> tuple[Any, Any, dict[str, Any]]:
    """Build API option structs from the Hamiltonian settings in the mdin."""
    cntrl, qmmm = _namelist(mdin, "cntrl"), _namelist(mdin, "qmmm")
    qmmm["dftb_telec"] = telec
    mm_options, qm_options = api.pme_input(), api.QmInputOptions()
    mm_keys = {"ntb", "igb", "ntc", "ntf", "cut", "ifqnt", "ntr", "jfastw"}
    missing = [key for key in qmmm if not hasattr(qm_options, key)]
    declared = {key for key in dir(qm_options) if not key.startswith("_")}
    missing_declared = [key for key in qmmm if key not in declared]
    if missing or missing_declared:
        raise AttributeError(
            "Python-sander QmInputOptions lacks declared QM fields: "
            + ", ".join(sorted(set(missing + missing_declared)))
        )
    applied_mm: dict[str, Any] = {}
    for key in mm_keys:
        if key in cntrl:
            if not hasattr(mm_options, key):
                raise AttributeError(f"Python-sander InputOptions lacks {key}")
            setattr(mm_options, key, cntrl[key])
            applied_mm[key] = cntrl[key]
    for key, value in qmmm.items():
        setattr(qm_options, key, value)
    return mm_options, qm_options, {"mm": applied_mm, "qm": qmmm}


def run_probe(
    topology: Path, restart: Path, mdin: Path, output: Path
) -> dict[str, Any]:
    """Check exact-temperature API support, then conditionally evaluate.

    Args:
        topology: Amber topology used by the converged native calculation.
        restart: Converged restart containing the final dimer coordinates.
        mdin: Native converged calculation input.
        output: New JSON result path.

    Returns:
        Compact provenance and feasibility record; no coordinate or force arrays.
    """
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite probe result: {output}")
    api_dir = Path(__import__("sander").__file__).resolve().parent
    import sander

    structure = AmberParm(str(topology))
    n_atoms = len(structure.atoms)
    active = _active_mask(mdin)
    n_active = sum(AmberMask(structure, active).Selection())
    coordinates = Rst7.open(str(restart))
    _validate_coordinates(coordinates.coordinates, n_atoms)
    report: dict[str, Any] = {
        "status": "blocked_before_energy_evaluation",
        "engine": "Python sander.setup / energy_forces(as_numpy=True)",
        "n_atoms": n_atoms,
        "active_mask_atoms": int(n_active),
        "requested_active_shell_angstrom": 12.0,
        "temperature_gate": {
            "reference_hartree": REFERENCE_HARTREE,
            "hartree_to_kcal_per_mol": HARTREE_TO_KCAL_MOL,
            "tolerance_kcal_per_mol": ENERGY_TOLERANCE_KCAL_MOL,
            "300K_energy_kcal_per_mol": None,
            "difference_kcal_per_mol": None,
            "passed": None,
        },
        "analytic_force_summary": None,
        "evaluations": [],
        "source_hashes": [
            _digest(path)
            for path in (
                mdin,
                topology,
                restart,
                Path(sander.__file__),
                Path(sander._pys.__file__),
                Path(__file__),
            )
        ],
        "runtime": {"python": sys.executable, "sander_package": str(api_dir)},
    }
    try:
        mm_options, qm_options, declared_options = _apply_options(sander, mdin, 300.0)
    except AttributeError as exc:
        report["reason"] = str(exc)
        report["api_declared_qm_fields"] = sorted(
            key for key in dir(sander.QmInputOptions()) if not key.startswith("_")
        )
        LOGGER.error("Cannot construct the matched 300 K Hamiltonian: %s", exc)
    else:
        report["declared_options"] = declared_options
        report["evaluations"].append(
            _evaluate(
                sander, topology, restart, mm_options, qm_options, n_atoms, active
            )
        )
        energy = report["evaluations"][-1]["energy_kcal_per_mol"]
        reference = REFERENCE_HARTREE * HARTREE_TO_KCAL_MOL
        difference = energy - reference
        report["temperature_gate"].update(
            {
                "300K_energy_kcal_per_mol": energy,
                "reference_kcal_per_mol": reference,
                "difference_kcal_per_mol": difference,
                "passed": abs(difference) <= ENERGY_TOLERANCE_KCAL_MOL,
            }
        )
        if not report["temperature_gate"]["passed"]:
            report["status"] = "300K_energy_mismatch_requires_diagnosis"
            report["reason"] = (
                "300 K energy differs from the native final energy by more than "
                "0.05 kcal/mol; 0 K evaluation was not run"
            )
        else:
            mm_options, qm_options, declared_options = _apply_options(sander, mdin, 0.0)
            report["declared_0K_options"] = declared_options
            report["evaluations"].append(
                _evaluate(
                    sander, topology, restart, mm_options, qm_options, n_atoms, active
                )
            )
            report["status"] = "matched_300K_and_evaluated_0K"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def _evaluate(
    api: Any,
    topology: Path,
    restart: Path,
    mm_options: Any,
    qm_options: Any,
    n_atoms: int,
    active_mask: str,
) -> dict[str, Any]:
    """Evaluate energy and summarize finite analytic active gradients."""
    structure = AmberParm(str(topology))
    active = np.asarray(AmberMask(structure, active_mask).Selection(), dtype=bool)
    with api.setup(str(topology), str(restart), None, mm_options, qm_options):
        energy, forces = api.energy_forces(as_numpy=True)
        force_array = np.asarray(forces, dtype=float).reshape((n_atoms, 3))
        active_gradient = -force_array[active]
        finite = bool(np.isfinite(float(energy.tot)) and np.isfinite(force_array).all())
        if not finite:
            raise ValueError("SANDER returned a nonfinite energy or analytic force")
        return {
            "dftb_telec_kelvin": float(qm_options.dftb_telec),
            "energy_kcal_per_mol": float(energy.tot),
            "analytic_force_finite": finite,
            "active_gradient_atoms": int(active.sum()),
            "active_gradient_max_abs_kcal_mol_angstrom": float(
                np.max(np.abs(active_gradient))
            ),
            "active_gradient_rms_kcal_mol_angstrom": float(
                np.sqrt(np.mean(active_gradient**2))
            ),
        }


def main() -> None:
    """Run the bounded, single-snapshot API feasibility check."""
    base = Path(
        "qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, default=base / "system.prmtop")
    parser.add_argument(
        "--restart", type=Path, default=base / "calculation/optimized.ncrst"
    )
    parser.add_argument("--mdin", type=Path, default=base / "calculation/dimer.mdin")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    result = run_probe(args.topology, args.restart, args.mdin, args.output)
    print(
        json.dumps(
            {
                key: result[key]
                for key in ("status", "reason", "n_atoms", "active_mask_atoms")
            }
        )
    )


if __name__ == "__main__":
    main()
