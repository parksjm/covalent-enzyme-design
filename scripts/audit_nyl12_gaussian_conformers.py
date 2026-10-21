#!/usr/bin/env python3
"""Audit two neutral Thr-NME Gaussian conformers without refitting charges."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.gaussian_conformer_audit import audit_conformer
from nyl12_covalent.gaussian_conformer_geometry import (
    heavy_atom_kabsch_rmsd,
    sha256,
    write_optimized_sdf,
)

LOGGER = logging.getLogger(__name__)
HARTREE_TO_KCAL_MOL = 627.509474
CONFORMERS = ("conformer_01", "conformer_02")


def _write_structures(
    rows: list[dict[str, Any]], prepared_root: Path, structures_dir: Path
) -> None:
    """Export optimized final geometries as SDF and XYZ with hashes."""
    structures_dir.mkdir(parents=True, exist_ok=False)
    for row in rows:
        name = row["name"]
        xyz = np.asarray(row["optimized_coordinates_A"], dtype=float)
        atom_map = row["atom_map"]
        write_optimized_sdf(
            prepared_root / name / "molecule.sdf",
            xyz,
            structures_dir / f"{name}.sdf",
        )
        lines = [str(len(atom_map)), f"{name}: final Gaussian coordinates; Angstrom"]
        lines.extend(
            f"{atom['element']:2s} {point[0]: .10f} {point[1]: .10f} {point[2]: .10f}"
            for atom, point in zip(atom_map, xyz, strict=True)
        )
        xyz_path = structures_dir / f"{name}.xyz"
        xyz_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        row["optimized_structure_files"] = {
            "sdf": {
                "path": f"{name}.sdf",
                "sha256": sha256(structures_dir / f"{name}.sdf"),
            },
            "xyz": {"path": f"{name}.xyz", "sha256": sha256(xyz_path)},
            "coordinate_units": "angstrom",
        }


def run_audit(
    calculation_root: Path,
    prepared_root: Path,
    output: Path,
    structures_dir: Path,
) -> dict[str, Any]:
    """Audit both conformers and create new report and coordinate exports."""
    if output.exists() or structures_dir.exists():
        raise FileExistsError("Report and structures outputs must not already exist")
    rows = []
    geometries = []
    for name in CONFORMERS:
        row, xyz, elements = audit_conformer(
            name, prepared_root / name, calculation_root / name
        )
        rows.append(row)
        geometries.append((xyz, elements))
    inventories = [
        tuple(
            (atom["index"], atom["name"], atom["residue"], atom["element"])
            for atom in row["atom_map"]
        )
        for row in rows
    ]
    if inventories[0] != inventories[1]:
        raise ValueError("Conformers do not share the same mapped atom inventory")
    reference_energy = min(row["optimized_wb97xd_energy_Eh"] for row in rows)
    for row in rows:
        delta = row["optimized_wb97xd_energy_Eh"] - reference_energy
        row["relative_optimized_wb97xd_energy_Eh"] = delta
        row["relative_optimized_wb97xd_energy_kcal_mol"] = delta * HARTREE_TO_KCAL_MOL
    report = {
        "model": "neutral free-NH2 L-threonine methylamide",
        "geometry_method": "wB97XD/def2-TZVP Opt/Freq",
        "esp_method": "HF/6-31G(d) separate Link1 stage; ESP prerequisite only",
        "energy_reference": "optimized wB97XD electronic SCF energies; independently referenced within this pair",
        "energy_conversion_kcal_mol_per_Eh": HARTREE_TO_KCAL_MOL,
        "heavy_atom_same_order_kabsch_rmsd_A": heavy_atom_kabsch_rmsd(
            geometries[0][0], geometries[1][0], geometries[0][1]
        ),
        "conformers": rows,
        "readiness": {
            "qm_geometry_frequency_esp_prerequisites_passed": True,
            "finite_temperature_populations_estimated": False,
            "charge_model_or_fitted_force_field_ready": False,
            "warning": "Relative electronic energies are not finite-temperature conformer populations; ESP output does not establish charge or force-field readiness.",
        },
    }
    _write_structures(rows, prepared_root, structures_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    LOGGER.info("Wrote conformer audit and optimized structures: %s", output)
    return report


def main() -> None:
    """Run an explicit audit and write the report and optimized structure files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calculation-root", type=Path, required=True)
    parser.add_argument("--prepared-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--structures-dir", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_audit(
        args.calculation_root, args.prepared_root, args.output, args.structures_dir
    )


if __name__ == "__main__":
    main()
