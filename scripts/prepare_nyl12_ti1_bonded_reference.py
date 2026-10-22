"""Prepare a crystal-constrained, capped anionic TI1 Gaussian reference."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.ti1_bonded_model import build_ti1_bonded_model  # noqa: E402
from nyl12_covalent.ti1_reference_coordinates import (  # noqa: E402
    TORSION_LABELS,
    crystal_backbone_reference,
    embed_bonded_reference,
)
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdMolTransforms  # noqa: E402


def prepare_reference(
    cif: Path,
    chain: str,
    residue: int,
    cip: str,
    cap_phi: float | None,
    seed: int,
    cores: int,
    memory_gb: int,
    output: Path,
) -> None:
    """Write immutable molecular, Gaussian and constraint input artifacts.

    Args:
        cif: User-selected protein reference.
        chain: Explicit reference author chain.
        residue: Catalytic Thr author residue number.
        cip: Explicit attacked-carbon configuration R or S.
        cap_phi: Declared cap phi, or None to require a connected source phi.
        seed: ETKDG seed.
        cores: Gaussian CPU allocation.
        memory_gb: Gaussian memory allocation.
        output: New prepared-input directory.
    """
    reference = crystal_backbone_reference(cif, chain, residue)
    phi = reference["phi_degrees"] if cap_phi is None else cap_phi
    if phi is None:
        raise ValueError("Source phi is undefined; supply explicit --cap-phi")
    if cores < 1 or memory_gb < 1:
        raise ValueError("CPU and memory allocations must be positive")
    molecule, metadata = build_ti1_bonded_model(cip)
    indices = metadata["atom_index_by_origin"]
    targets = {
        "phi": phi,
        "psi": reference["psi_degrees"],
        "chi1": reference["chi1_degrees"],
    }
    molecule = embed_bonded_reference(molecule, indices, targets, seed)
    conformer = molecule.GetConformer()
    constraints = []
    for name, labels in TORSION_LABELS.items():
        atoms = [indices[label] + 1 for label in labels]
        value = float(
            rdMolTransforms.GetDihedralDeg(conformer, *[i - 1 for i in atoms])
        )
        constraints.append({"name": name, "atoms_one_based": atoms, "degrees": value})
    route = "#p wB97XD/def2TZVP Int=UltraFine Opt=(Tight,CalcFC,ModRedundant,MaxCycles=200) SCF=(Tight,XQC) NoSymm"
    lines = [
        f"%nprocshared={cores}",
        f"%mem={memory_gb}GB",
        "%chk=reference.chk",
        route,
        "",
        "Capped anionic TI1 bonded surrogate. Constrained reference.",
        "",
        "-1 1",
    ]
    xyz = [str(molecule.GetNumAtoms()), "Capped TI1 initial seed, not optimized"]
    for atom in molecule.GetAtoms():
        point = conformer.GetAtomPosition(atom.GetIdx())
        row = f"{atom.GetSymbol()} {point.x:.10f} {point.y:.10f} {point.z:.10f}"
        lines.append(row)
        xyz.append(row)
    lines.append("")
    for entry in constraints:
        lines.append("D " + " ".join(str(i) for i in entry["atoms_one_based"]) + " F")
    lines.extend(["", ""])
    metadata.update(
        {
            "coordinate_reference": reference,
            "constraints": constraints,
            "cap_phi_override_degrees": cap_phi,
            "seed": seed,
            "geometry_method": "wB97XD/def2TZVP",
            "solvent": "gas phase",
            "cores": cores,
            "memory_gb": memory_gb,
            "full_stationary_minimum_claimed": False,
            "forcefield_validated": False,
        }
    )
    output.mkdir(parents=True, exist_ok=False)
    with Chem.SDWriter(str(output / "molecule.sdf")) as writer:
        writer.write(molecule)
    (output / "coordinates.xyz").write_text("\n".join(xyz) + "\n")
    (output / "atom_map.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (output / "input.gjf").write_text("\n".join(lines))
    for filename, checkpoint, task in (
        ("frequency.gjf", "hessian.chk", "wB97XD/def2TZVP Int=UltraFine Freq"),
        ("esp.gjf", "esp.chk", "HF/6-31G(d) Pop=MK IOp(6/33=2,6/42=6)"),
    ):
        (output / filename).write_text(
            f"%nprocshared={cores}\n%mem={memory_gb}GB\n"
            f"%oldchk=reference.chk\n%chk={checkpoint}\n"
            f"#p {task} Geom=AllCheck Guess=Read SCF=(Tight,XQC) NoSymm\n\n\n"
        )
    logging.info("Prepared constrained %s TI1 bonded reference at %s", cip, output)


def main() -> None:
    """Parse explicit structure, chemistry, resources and output choices."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif", type=Path, required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--residue", type=int, required=True)
    parser.add_argument("--cip", choices=("R", "S"), required=True)
    parser.add_argument("--cap-phi", type=float)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--cores", type=int, required=True)
    parser.add_argument("--memory-gb", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    prepare_reference(
        args.cif,
        args.chain,
        args.residue,
        args.cip,
        args.cap_phi,
        args.seed,
        args.cores,
        args.memory_gb,
        args.output,
    )


if __name__ == "__main__":
    main()
