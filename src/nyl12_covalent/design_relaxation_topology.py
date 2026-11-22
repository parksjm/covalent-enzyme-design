"""Build and audit matched four-ligand PA6 relaxation topologies."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import parmed

from .design_charge import protein_formal_charge

LOGGER = logging.getLogger(__name__)
RELEASE_ID = "pa6_l4_modular_resp_gaff2_v1_20261031"
PA6_SEQUENCE = ("P6N", "P6I", "P6I", "P6C")
FIRST_PA6_RESIDUE = 1313
BETA_THR_RESIDUES = (929, 1025, 1121, 1217)
ASP_RESIDUES = (968, 970, 1064, 1066, 1160, 1162, 1256, 1258)
NEUTRAL_THR_CHARGES = {"N": -0.4834, "H1": 0.1224, "H2": 0.1224}


def leap_input(parameters: Path, ions: tuple[int, int] | None) -> str:
    """Return LEaP input for matched ff19SB/TIP3P PA6 L4 preparation.

    Args:
        parameters: Released PA6 parameter directory, or its preserved copy.
        ions: Sodium and chloride counts for ``addIonsRand``, or ``None`` to
            save the solvated, unneutralized system.

    Returns:
        Complete LEaP input text. LEaP is not executed by this function.

    Raises:
        ValueError: If ion counts are negative.
        FileNotFoundError: If the released PA6 LEaP recipe is missing.
    """
    recipe = parameters.resolve() / "leaprc.pa6_shared_p6i"
    if not recipe.is_file():
        raise FileNotFoundError(recipe)
    if ions is not None and (len(ions) != 2 or min(ions) < 0):
        raise ValueError("ions must contain nonnegative sodium and chloride counts")
    edit_lines = [
        f"remove complex complex.{residue}.H3\n"
        f"set complex.{residue}.N charge -0.4834\n"
        f"set complex.{residue}.H1 charge 0.1224\n"
        f"set complex.{residue}.H2 charge 0.1224"
        for residue in BETA_THR_RESIDUES
    ]
    save_name = "noions" if ions is None else "complex"
    ion_line = (
        "" if ions is None else f"addIonsRand complex Na+ {ions[0]} Cl- {ions[1]} 6.0\n"
    )
    return (
        "source leaprc.protein.ff19SB\n"
        "source leaprc.water.tip3p\n"
        f"source {recipe}\n"
        "complex = loadpdb complex_heavy_modular.pdb\n"
        + "\n".join(edit_lines)
        + "\ncheck complex\ncharge complex\n"
        "saveamberparm complex dry.prmtop dry.rst7\n"
        "solvateOct complex TIP3PBOX 12.0\n"
        + ion_line
        + "check complex\ncharge complex\n"
        + f"saveamberparm complex {save_name}.prmtop {save_name}.rst7\nquit\n"
    )


def _mol2_atoms(path: Path) -> dict[str, tuple[str, float]]:
    """Read named atom types and charges from a released MOL2 template."""
    in_atoms = False
    result: dict[str, tuple[str, float]] = {}
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atoms = False
        elif in_atoms and line.strip():
            fields = line.split()
            name, atom_type, charge = fields[1], fields[5], float(fields[8])
            if name in result:
                raise ValueError(f"Duplicate atom {name} in {path}")
            result[name] = (atom_type, charge)
    if not result:
        raise ValueError(f"No atom records in {path}")
    return result


def _residue_charge(residue: Any) -> float:
    """Return a residue's total partial charge."""
    return float(sum(atom.charge for atom in residue.atoms))


def _validate_pa6(residues: list[Any], release: Path) -> tuple[int, float, list[float]]:
    """Check four exact released PA6 copies and their named peptide links."""
    copy_count = 4
    expected_names = [name for _ in range(copy_count) for name in PA6_SEQUENCE]
    ligand_residues = residues[FIRST_PA6_RESIDUE - 1 : FIRST_PA6_RESIDUE - 1 + 16]
    observed = [residue.name.strip() for residue in ligand_residues]
    if observed != expected_names:
        raise ValueError(
            f"Expected four PA6 L4 blocks {expected_names}, found {observed}"
        )
    templates = {
        name: _mol2_atoms(release / f"{name}.mol2") for name in set(PA6_SEQUENCE)
    }
    total_charge = 0.0
    copy_charges: list[float] = []
    for residue, template_name in zip(ligand_residues, expected_names, strict=True):
        observed_atoms = {atom.name.strip(): atom for atom in residue.atoms}
        template = templates[template_name]
        if set(observed_atoms) != set(template):
            raise ValueError(f"Atom names differ for PA6 residue {residue.idx + 1}")
        for atom_name, expected in template.items():
            atom = observed_atoms[atom_name]
            if (
                atom.type.strip() != expected[0]
                or abs(atom.charge - expected[1]) > 1.0e-5
            ):
                raise ValueError(
                    f"Released type/charge mismatch at {residue.idx + 1}:{atom_name}"
                )
        residue_charge = _residue_charge(residue)
        if abs(residue_charge - sum(charge for _, charge in template.values())) > 1e-5:
            raise ValueError(f"PA6 residue charge mismatch at {residue.idx + 1}")
        total_charge += residue_charge
        if (residue.idx + 1 - FIRST_PA6_RESIDUE) % 4 == 3:
            copy_start = residue.idx - (FIRST_PA6_RESIDUE - 1) - 3
            copy_charge = sum(
                _residue_charge(ligand_residues[index])
                for index in range(copy_start, copy_start + 4)
            )
            if abs(copy_charge) > 1.0e-5:
                raise ValueError(
                    f"PA6 L4 copy ending at {residue.idx + 1} has charge {copy_charge}"
                )
            copy_charges.append(copy_charge)
    return (
        sum(len(residue.atoms) for residue in ligand_residues),
        total_charge,
        copy_charges,
    )


def _validate_links(topology: Any) -> int:
    """Require exactly the twelve within-copy carbonyl-C to amide-N links."""
    residue_index = {
        id(residue): index + 1 for index, residue in enumerate(topology.residues)
    }
    pa6_links: list[tuple[int, str, int, str]] = []
    protein_crosslinks: list[tuple[int, int]] = []
    for bond in topology.bonds:
        first, second = bond.atom1, bond.atom2
        first_res = residue_index[id(first.residue)]
        second_res = residue_index[id(second.residue)]
        first_pa6 = FIRST_PA6_RESIDUE <= first_res < FIRST_PA6_RESIDUE + 16
        second_pa6 = FIRST_PA6_RESIDUE <= second_res < FIRST_PA6_RESIDUE + 16
        if first_pa6 != second_pa6:
            protein_crosslinks.append((first_res, second_res))
        if first_pa6 and second_pa6 and first_res != second_res:
            pa6_links.append(
                (first_res, first.name.strip(), second_res, second.name.strip())
            )
    if protein_crosslinks:
        raise ValueError(f"PA6-protein covalent bond found: {protein_crosslinks[0]}")
    expected = {
        (start + offset, "C", start + offset + 1, "N")
        for start in (
            FIRST_PA6_RESIDUE,
            FIRST_PA6_RESIDUE + 4,
            FIRST_PA6_RESIDUE + 8,
            FIRST_PA6_RESIDUE + 12,
        )
        for offset in range(3)
    }
    canonical = {
        (min(r1, r2), n1 if r1 < r2 else n2, max(r1, r2), n2 if r1 < r2 else n1)
        for r1, n1, r2, n2 in pa6_links
    }
    if canonical != expected:
        raise ValueError(
            f"Expected twelve within-copy PA6 C-N links; observed {canonical}"
        )
    return len(canonical)


def _validate_thr_and_asp(residues: list[Any]) -> dict[str, Any]:
    """Validate all four neutral beta Thr termini and eight anionic Asp sites."""
    thr_charges: dict[str, float] = {}
    for number in BETA_THR_RESIDUES:
        residue = residues[number - 1]
        if residue.name.strip() != "THR":
            raise ValueError(f"Expected THR at Amber residue {number}")
        atoms = {atom.name.strip(): atom for atom in residue.atoms}
        if not {"N", "H1", "H2", "HG1"}.issubset(atoms) or "H3" in atoms:
            raise ValueError(f"Incorrect neutral Thr hydrogen inventory at {number}")
        for atom_name, expected in NEUTRAL_THR_CHARGES.items():
            if abs(atoms[atom_name].charge - expected) > 1.0e-5:
                raise ValueError(f"Neutral Thr {number}:{atom_name} charge mismatch")
        charge = _residue_charge(residue)
        if abs(charge) > 1.0e-5:
            raise ValueError(f"Neutral Thr residue {number} has charge {charge}")
        thr_charges[str(number)] = charge
    asp_charges: dict[str, float] = {}
    for number in ASP_RESIDUES:
        residue = residues[number - 1]
        if residue.name.strip() != "ASP":
            raise ValueError(f"Expected ASP at Amber residue {number}")
        charge = _residue_charge(residue)
        if abs(charge + 1.0) > 1.0e-4:
            raise ValueError(f"ASP residue {number} must be -1, found {charge}")
        asp_charges[str(number)] = charge
    return {
        "neutral_thr_residue_charges": thr_charges,
        "asp_residue_charges": asp_charges,
    }


def validate_four_ligand_topology(
    topology: Path, release: Path, solvated: bool
) -> dict[str, Any]:
    """Validate PA6 parameters, links, protonation and system charge.

    Args:
        topology: Amber topology to audit.
        release: Preserved copy of the released modular PA6 parameters.
        solvated: Whether the topology contains water and ions.

    Returns:
        Counts, charge inventory, and validated release identity.

    Raises:
        FileNotFoundError: If topology or required release templates are missing.
        ValueError: If topology composition, parameter identity or charge gates fail.
    """
    if not topology.is_file():
        raise FileNotFoundError(topology)
    if not release.is_dir():
        raise FileNotFoundError(release)
    required = [
        release / "manifest.json",
        release / "leaprc.pa6_shared_p6i",
        *(release / f"{name}.mol2" for name in set(PA6_SEQUENCE)),
        *(release / f"{name}.lib" for name in set(PA6_SEQUENCE)),
        release / "pa6_shared_p6i_gaff2.frcmod",
        release / "pa6_shared_p6i_improper.frcmod",
    ]
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    release_manifest = json.loads((release / "manifest.json").read_text())
    if release_manifest.get("release_id") != RELEASE_ID:
        raise ValueError(
            f"Unexpected PA6 release: {release_manifest.get('release_id')}"
        )
    system = parmed.load_file(str(topology))
    residues = list(system.residues)
    if len(residues) < FIRST_PA6_RESIDUE + 15:
        raise ValueError(f"Expected at least 1328 residues, found {len(residues)}")
    if not solvated and len(residues) != FIRST_PA6_RESIDUE + 15:
        raise ValueError(
            f"Dry topology must contain exactly 1328 residues, found {len(residues)}"
        )
    radii = " ".join(str(value) for value in system.parm_data.get("RADIUS_SET", []))
    if "mbondi3" not in radii.lower():
        raise ValueError(f"Expected mbondi3 radii, found {radii!r}")
    pa6_atom_count, pa6_charge, pa6_copy_charges = _validate_pa6(residues, release)
    link_count = _validate_links(system)
    inventory = _validate_thr_and_asp(residues)
    formal_solute_charge = protein_formal_charge(residues)
    if solvated:
        sodium_count = 0
        chloride_count = 0
        for residue in residues[FIRST_PA6_RESIDUE + 15 :]:
            if residue.name.strip() not in {"Na+", "Cl-"}:
                continue
            expected_ion_charge = 1 if residue.name.strip() == "Na+" else -1
            if abs(_residue_charge(residue) - expected_ion_charge) > 1.0e-4:
                raise ValueError(
                    f"{residue.name.strip()} residue {residue.idx + 1} has "
                    f"incorrect formal ion charge"
                )
            if expected_ion_charge > 0:
                sodium_count += 1
            else:
                chloride_count += 1
        if sodium_count - chloride_count != -formal_solute_charge:
            raise ValueError(
                "Counterion inventory does not neutralize formal solute charge: "
                f"Na={sodium_count}, Cl={chloride_count}, "
                f"solute={formal_solute_charge}"
            )
    total_charge = float(sum(atom.charge for atom in system.atoms))
    expected_charge = 0.0 if solvated else formal_solute_charge
    # Four released L4 copies accumulate about -1.6e-5 e from printed charges.
    # Keep all per-atom/residue chemistry checks tighter; this is system roundoff.
    tolerance = 1.0e-4
    if abs(total_charge - expected_charge) > tolerance:
        raise ValueError(
            f"System charge {total_charge:.8f} differs from expected {expected_charge:.1f}"
        )
    result = {
        "topology": str(topology),
        "release_id": RELEASE_ID,
        "solvated": solvated,
        "residue_count": len(residues),
        "atom_count": len(system.atoms),
        "pa6_residue_count": 16,
        "pa6_atom_count": pa6_atom_count,
        "pa6_charge": pa6_charge,
        "pa6_copy_charges": pa6_copy_charges,
        "pa6_interresidue_link_count": link_count,
        "formal_solute_charge": formal_solute_charge,
        "total_charge": total_charge,
        "radius_set": radii,
        **inventory,
    }
    LOGGER.info("Validated four-ligand topology %s", topology)
    return result
