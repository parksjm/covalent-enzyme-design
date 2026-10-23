"""Validate Gaussian geometry/frequency/ESP prerequisites without fitting charges."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import Chem


def final_orientation(text: str) -> tuple[list[int], np.ndarray]:
    """Read the last Gaussian Cartesian orientation in original atom order.

    Args:
        text: Gaussian output text.

    Returns:
        Atomic numbers and finite coordinates in Angstrom.

    Raises:
        ValueError: If the table is absent or malformed.
    """
    starts = list(re.finditer(r"(?:Input|Standard) orientation:", text))
    if not starts:
        raise ValueError("Gaussian output contains no Cartesian orientation")
    table = text[starts[-1].end() :].splitlines()
    separators = 0
    elements, coordinates = [], []
    for line in table:
        if re.fullmatch(r"\s*-{5,}\s*", line):
            separators += 1
            if separators == 3:
                break
            continue
        if separators == 2:
            fields = line.split()
            if len(fields) != 6 or int(fields[0]) != len(elements) + 1:
                raise ValueError("Malformed Gaussian orientation atom order")
            elements.append(int(fields[1]))
            coordinates.append([float(v) for v in fields[3:]])
    xyz = np.asarray(coordinates)
    if (
        separators != 3
        or xyz.shape != (len(elements), 3)
        or not elements
        or not np.isfinite(xyz).all()
    ):
        raise ValueError("Incomplete/nonfinite Gaussian orientation")
    return elements, xyz


def validate_geometry(
    molecule: Chem.Mol, elements: list[int], xyz: np.ndarray
) -> dict[str, Any]:
    """Check atom order, covalent lengths, and original stereocenters after QM.

    Args:
        molecule: Explicit-hydrogen molecular graph with defined stereochemistry.
        elements: Gaussian atomic numbers in input order.
        xyz: Final Gaussian Cartesian coordinates.

    Returns:
        Conservative connectivity and stereochemical screening report.
    """
    if elements != [a.GetAtomicNum() for a in molecule.GetAtoms()]:
        raise ValueError("Gaussian final atom inventory/order differs from input graph")
    if xyz.shape != (len(elements), 3) or not np.isfinite(xyz).all():
        raise ValueError("Invalid final geometry shape/values")
    original = Chem.Mol(molecule)
    Chem.AssignStereochemistry(original, cleanIt=True, force=True)
    expected = {
        a.GetIdx(): a.GetProp("_CIPCode")
        for a in original.GetAtoms()
        if a.HasProp("_CIPCode")
    }
    final = Chem.Mol(original)
    final.RemoveAllConformers()
    conformer = Chem.Conformer(len(elements))
    conformer.Set3D(True)
    for index, point in enumerate(xyz):
        conformer.SetAtomPosition(index, tuple(float(v) for v in point))
    final.AddConformer(conformer)
    Chem.AssignAtomChiralTagsFromStructure(final, replaceExistingTags=True)
    Chem.AssignStereochemistry(final, cleanIt=True, force=True)
    observed = {
        index: final.GetAtomWithIdx(index).GetProp("_CIPCode")
        if final.GetAtomWithIdx(index).HasProp("_CIPCode")
        else None
        for index in expected
    }
    radii = np.asarray([Chem.GetPeriodicTable().GetRcovalent(z) for z in elements])
    bonded = {
        tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
        for b in original.GetBonds()
    }
    failures = []
    for first in range(len(elements)):
        for second in range(first + 1, len(elements)):
            length = float(np.linalg.norm(xyz[first] - xyz[second]))
            radii_sum = float(radii[first] + radii[second])
            lower, upper = 0.65 * radii_sum, 1.30 * radii_sum
            if (first, second) in bonded:
                if not lower <= length <= upper:
                    failures.append(
                        {
                            "pair": [first + 1, second + 1],
                            "distance_A": length,
                            "reason": "input bond outside covalent-length screen",
                        }
                    )
            elif length < 1.10 * radii_sum:
                failures.append(
                    {
                        "pair": [first + 1, second + 1],
                        "distance_A": length,
                        "reason": "new short nonbonded contact/possible proton transfer",
                    }
                )
    return {
        "accepted": expected == observed and not failures,
        "expected_stereocenters": expected,
        "observed_stereocenters": observed,
        "geometry_failures": failures,
    }


def assess_gaussian(log_path: Path, molecule_path: Path) -> dict[str, Any]:
    """Require completed optimization, a vibrational minimum, and ESP output.

    Args:
        log_path: Combined DFT optimization/frequency and HF ESP Gaussian log.
        molecule_path: Original explicit-hydrogen SDF containing one molecule.

    Returns:
        Scientific acceptance of QM prerequisites, never of a fitted force field.
    """
    molecules = list(Chem.SDMolSupplier(str(molecule_path), removeHs=False))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Expected exactly one valid explicit-hydrogen molecule")
    molecule = molecules[0]
    text = log_path.read_text()
    numbers, xyz = final_orientation(text)
    frequencies = [
        float(value)
        for line in re.findall(r"Frequencies --([^\n]+)", text)
        for value in line.split()
    ]
    expected_modes = 3 * molecule.GetNumAtoms() - 6
    frequency_pass = (
        len(frequencies) == expected_modes
        and np.isfinite(frequencies).all()
        and all(value > 0 for value in frequencies)
    )
    geometry = validate_geometry(molecule, numbers, xyz)
    charges = [
        (int(q), int(m))
        for q, m in re.findall(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", text)
    ]
    electronic_state = bool(
        charges
        and all(q == Chem.GetFormalCharge(molecule) and m == 1 for q, m in charges)
    )
    normal_count = text.count("Normal termination of Gaussian 16")
    finished = (
        "Normal termination of Gaussian 16" in text.rstrip().splitlines()[-1]
        and normal_count >= 2
        and "Error termination" not in text
    )
    optimized = "Optimization completed." in text
    esp = "ESP Fit Center" in text and "Electrostatic Properties" in text
    return {
        "accepted": bool(
            finished
            and optimized
            and frequency_pass
            and electronic_state
            and esp
            and geometry["accepted"]
        ),
        "purpose": "QM geometry/frequency/ESP prerequisites only; no validated fitted charges",
        "normal_termination_count": normal_count,
        "completed": finished,
        "optimization_completed": optimized,
        "frequency_minimum_passed": frequency_pass,
        "expected_frequency_count": expected_modes,
        "frequencies_cm_inverse": frequencies,
        "electronic_state_passed": electronic_state,
        "esp_grid_present": esp,
        "geometry": geometry,
        "final_coordinates_A": xyz.tolist(),
    }
