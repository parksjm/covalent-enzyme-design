from __future__ import annotations
import math
from pathlib import Path
from typing import Any
from reference_common import TOPOLOGY


def parse_external_input(path: Path) -> dict[str, Any]:
    """Parse an Amber Gaussian deck with the fixed 110-atom QM/link model."""
    lines = path.read_text(encoding="utf-8").splitlines()
    charge_index = next(
        (i for i, line in enumerate(lines) if line.strip() == "-1 1"), None
    )
    if charge_index is None:
        raise ValueError(f"Missing charge/multiplicity line in {path}")
    atoms: list[tuple[str, float, float, float]] = []
    index = charge_index + 1
    while index < len(lines) and lines[index].strip():
        fields = lines[index].split()
        if len(fields) != 4:
            raise ValueError(f"Malformed QM/link coordinate: {lines[index]}")
        atoms.append((fields[0], *(float(value) for value in fields[1:])))
        index += 1
    if len(atoms) != 110:
        raise ValueError(f"Expected 102 QM atoms + 8 links, found {len(atoms)}")
    index += 1
    point_charges: list[tuple[float, float, float, float]] = []
    while index < len(lines) and lines[index].strip():
        fields = lines[index].split()
        if len(fields) != 4:
            raise ValueError(f"Malformed external point charge: {lines[index]}")
        point_charges.append(tuple(float(value) for value in fields))
        index += 1
    return {
        "lines": lines,
        "charge": -1,
        "multiplicity": 1,
        "atoms": atoms,
        "point_charges": point_charges,
        "pc_section_end": index,
    }


def write_direct_gaussian(source: Path, output_dir: Path) -> Path:
    """Write a direct Force input, removing only Amber's PC field-query block."""
    embedding = parse_external_input(source)
    lines = embedding["lines"]
    route_indices = [i for i, line in enumerate(lines) if line.lstrip().startswith("#")]
    if len(route_indices) != 1 or "Prop=(Field,Read)" not in lines[route_indices[0]]:
        raise ValueError(
            "Expected one Amber route line with exactly the Field/Read property"
        )
    route = lines[route_indices[0]].replace("Prop=(Field,Read)", "")
    retained = lines[: embedding["pc_section_end"]]
    retained[route_indices[0]] = route
    retained[0] = "%chk=direct_force.chk"
    path = output_dir / "direct_force.gjf"
    path.write_text("\n".join(retained) + "\n\n", encoding="utf-8")
    return path


def write_xtb_inputs(source: Path, output: Path) -> dict[str, Any]:
    """Write the 110-atom GFN2 geometry and exact Sander-adjusted MM charges."""
    embedding = parse_external_input(source)
    xyz_lines = [
        f"{symbol:<2s} {x: .15f} {y: .15f} {z: .15f}"
        for symbol, x, y, z in embedding["atoms"]
    ]
    (output / "qm_link.xyz").write_text(
        "110\nAmber QM atoms plus eight link hydrogens\n" + "\n".join(xyz_lines) + "\n",
        encoding="utf-8",
    )
    pc_lines = [
        f"{q: .16f} {x * 1.8897261254578281: .16f} {y * 1.8897261254578281: .16f} "
        f"{z * 1.8897261254578281: .16f} 999.0"
        for x, y, z, q in embedding["point_charges"]
    ]
    (output / "point_charges.pc").write_text(
        f"{len(pc_lines)}\n" + "\n".join(pc_lines) + "\n", encoding="utf-8"
    )
    (output / "xtb.inp").write_text(
        "$embedding\n  interface=legacy\n  input=point_charges.pc\n$end\n",
        encoding="utf-8",
    )
    return embedding


def validate_pc_membership(
    source: Path,
    membership: dict[str, Any],
    expected_charges: list[float] | None = None,
) -> list[float]:
    """Require exactly the ordered all-MM PC list minus the eight link parents."""
    embedding = parse_external_input(source)
    points = embedding["point_charges"]
    expected_count = int(membership["mm_atom_count"]) - 8
    if len(points) != expected_count:
        raise ValueError(
            f"Expected {expected_count} MM point charges, found {len(points)}"
        )
    if expected_charges is not None:
        charges = [row[3] for row in points]
        if len(charges) != len(expected_charges) or any(
            abs(a - b) > 1e-10 for a, b in zip(charges, expected_charges)
        ):
            raise ValueError(
                "Adjusted point-charge sequence differs from the frozen ES membership"
            )
    return [row[3] for row in points]


def expected_mm_ids(membership: dict[str, Any]) -> list[int]:
    """Derive the fixed MM list by excluding QM atoms and eight bonded link parents."""
    import parmed

    system = parmed.load_file(str(TOPOLOGY))
    qm_ids = {int(value) for value in membership["qm_real_atom_ids"]}
    link_parents: set[int] = set()
    for bond in system.bonds:
        left, right = bond.atom1.idx + 1, bond.atom2.idx + 1
        if left in qm_ids and right not in qm_ids:
            link_parents.add(right)
        elif right in qm_ids and left not in qm_ids:
            link_parents.add(left)
    if len(link_parents) != 8:
        raise ValueError(
            f"Expected exactly eight QM/MM link parents, found {len(link_parents)}"
        )
    ids = [
        int(atom["amber_atom_id"])
        for atom in membership["atoms"]
        if int(atom["amber_atom_id"]) not in link_parents
    ]
    if len(ids) != int(membership["mm_atom_count"]) - 8:
        raise ValueError(
            "Fixed MM membership does not equal non-QM atom IDs minus eight link parents"
        )
    return ids


def verify_baseline_pc_ids(
    source: Path, membership: dict[str, Any], ids: list[int]
) -> None:
    """Map baseline external-PC coordinate order to source Amber atom IDs."""
    embedding = parse_external_input(source)
    raw_by_id = {int(atom["amber_atom_id"]): atom for atom in membership["atoms"]}
    if len(embedding["point_charges"]) != len(ids):
        raise ValueError(
            "Baseline point-charge count does not match derived Amber membership"
        )
    for atom_id, point in zip(ids, embedding["point_charges"]):
        xyz = raw_by_id[atom_id]["xyz_A"]
        distance = math.sqrt(sum((point[i] - xyz[i]) ** 2 for i in range(3)))
        if distance > 1e-5:
            raise ValueError(
                f"PC order/coordinate mismatch at Amber atom {atom_id}: {distance:.6g} A"
            )
