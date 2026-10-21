"""Build the neutral threonine methylamide charge-template conformers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolTransforms, rdchem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

GeometryMethod = Literal["b3lyp-d3bj-def2tzvp", "wb97xd-def2tzvp"]

_EXPECTED_CIP = {"CA": "S", "CB": "R"}
_ATOM_LABELS = (
    "N",
    "CA",
    "CB",
    "OG1",
    "CG2",
    "C",
    "O",
    "NME_N",
    "NME_C",
)
_HEAVY_BONDS = (
    ("N", "CA", rdchem.BondType.SINGLE),
    ("CA", "CB", rdchem.BondType.SINGLE),
    ("CB", "OG1", rdchem.BondType.SINGLE),
    ("CB", "CG2", rdchem.BondType.SINGLE),
    ("CA", "C", rdchem.BondType.SINGLE),
    ("C", "O", rdchem.BondType.DOUBLE),
    ("C", "NME_N", rdchem.BondType.SINGLE),
    ("NME_N", "NME_C", rdchem.BondType.SINGLE),
)
_HYDROGEN_NAMES = {
    "N": ("H1", "H2"),
    "CA": ("HA",),
    "CB": ("HB",),
    "OG1": ("HG1",),
    "CG2": ("HG21", "HG22", "HG23"),
    "NME_N": ("H",),
    "NME_C": ("HH31", "HH32", "HH33"),
}
_GEOMETRY_ROUTE = (
    "#p B3LYP/def2TZVP EmpiricalDispersion=GD3BJ Int=UltraFine "
    "Opt=(Tight,CalcFC,MaxCycles=200) Freq SCF=(Tight,XQC) NoSymm"
)
_WB97XD_GEOMETRY_ROUTE = (
    "#p wB97XD/def2TZVP Int=UltraFine "
    "Opt=(Tight,CalcFC,MaxCycles=200) Freq SCF=(Tight,XQC) NoSymm"
)
_ESP_ROUTE = (
    "#p HF/6-31G(d) Geom=AllCheck Guess=Read SCF=(Tight,XQC) "
    "NoSymm Pop=MK IOp(6/33=2,6/42=6)"
)


def _build_template() -> tuple[Chem.Mol, dict[str, int]]:
    """Build and validate the explicitly named molecular graph."""
    editable = Chem.RWMol()
    atom_indices: dict[str, int] = {}
    atomic_numbers = {
        "N": 7,
        "CA": 6,
        "CB": 6,
        "OG1": 8,
        "CG2": 6,
        "C": 6,
        "O": 8,
        "NME_N": 7,
        "NME_C": 6,
    }
    for name in _ATOM_LABELS:
        atom = Chem.Atom(atomic_numbers[name])
        atom.SetProp("_Name", name)
        atom.SetProp("atom_name", {"NME_N": "N", "NME_C": "CH3"}.get(name, name))
        atom.SetProp("residue_name", "THR" if not name.startswith("NME_") else "NME")
        atom_indices[name] = editable.AddAtom(atom)
    for first, second, bond_type in _HEAVY_BONDS:
        editable.AddBond(atom_indices[first], atom_indices[second], bond_type)

    molecule = editable.GetMol()
    Chem.SanitizeMol(molecule)
    _set_target_chirality(molecule, atom_indices)
    molecule = Chem.AddHs(molecule)
    _name_hydrogens(molecule, atom_indices)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    _audit_graph_stereochemistry(molecule, atom_indices)
    formula = CalcMolFormula(molecule)
    if formula != "C5H12N2O2":
        raise ValueError(f"Unexpected neutral template formula: {formula}")
    if Chem.GetFormalCharge(molecule) != 0:
        raise ValueError("Neutral template must have formal charge zero")
    return molecule, atom_indices


def _set_target_chirality(molecule: Chem.Mol, atom_indices: dict[str, int]) -> None:
    """Set tetrahedral tags to the requested CIP labels and verify them."""
    for label, target in _EXPECTED_CIP.items():
        atom = molecule.GetAtomWithIdx(atom_indices[label])
        atom.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
        if atom.GetProp("_CIPCode") != target:
            atom.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    _audit_graph_stereochemistry(molecule, atom_indices)


def _name_hydrogens(molecule: Chem.Mol, atom_indices: dict[str, int]) -> None:
    """Assign deterministic names to the explicit hydrogen atoms."""
    next_name = {name: 0 for name in _HYDROGEN_NAMES}
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        parent = atom.GetNeighbors()[0]
        parent_name = parent.GetProp("_Name")
        names = _HYDROGEN_NAMES[parent_name]
        position = next_name[parent_name]
        if position >= len(names):
            raise ValueError(f"Too many hydrogens attached to {parent_name}")
        name = names[position]
        next_name[parent_name] += 1
        atom.SetProp("_Name", name)
        atom.SetProp("atom_name", name)
        atom.SetProp("residue_name", parent.GetProp("residue_name"))
    expected_counts = {name: len(names) for name, names in _HYDROGEN_NAMES.items()}
    if next_name != expected_counts:
        raise ValueError(f"Hydrogen inventory mismatch: {next_name}")


def _audit_graph_stereochemistry(
    molecule: Chem.Mol, atom_indices: dict[str, int]
) -> None:
    """Raise unless the molecular graph carries CA S and CB R stereochemistry."""
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    for label, target in _EXPECTED_CIP.items():
        atom = molecule.GetAtomWithIdx(atom_indices[label])
        if not atom.HasProp("_CIPCode") or atom.GetProp("_CIPCode") != target:
            observed = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "unset"
            raise ValueError(
                f"{label} stereochemistry is {observed}; expected {target}"
            )


def _audit_coordinate_stereochemistry(
    molecule: Chem.Mol, atom_indices: dict[str, int]
) -> None:
    """Raise unless embedded coordinates encode the expected CIP labels."""
    Chem.AssignStereochemistryFrom3D(molecule, confId=0, replaceExistingTags=True)
    _audit_graph_stereochemistry(molecule, atom_indices)


def _write_atom_metadata(molecule: Chem.Mol, path: Path) -> None:
    atoms = []
    for atom in molecule.GetAtoms():
        atoms.append(
            {
                "index": atom.GetIdx(),
                "name": atom.GetProp("atom_name"),
                "residue": atom.GetProp("residue_name"),
                "element": atom.GetSymbol(),
                "formal_charge": atom.GetFormalCharge(),
            }
        )
    path.write_text(json.dumps({"atoms": atoms}, indent=2) + "\n", encoding="utf-8")


def _write_xyz(molecule: Chem.Mol, path: Path) -> None:
    conformer = molecule.GetConformer()
    lines = [str(molecule.GetNumAtoms()), "neutral Thr methylamide conformer"]
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol():2s} {position.x: .8f} {position.y: .8f} {position.z: .8f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _gaussian_input(
    molecule: Chem.Mol, chk_name: str, geometry_method: GeometryMethod
) -> str:
    """Return Gaussian geometry-optimization and MK ESP input text."""
    geometry_routes = {
        "b3lyp-d3bj-def2tzvp": _GEOMETRY_ROUTE,
        "wb97xd-def2tzvp": _WB97XD_GEOMETRY_ROUTE,
    }
    if geometry_method not in geometry_routes:
        raise ValueError(f"Unsupported geometry method: {geometry_method}")
    geometry = geometry_routes[geometry_method]
    lines = [
        "%nprocshared=4",
        "%mem=4GB",
        f"%chk={chk_name}",
        geometry,
        "",
        "neutral beta-N charge-template; not enzyme refinement",
        "",
        "0 1",
    ]
    conformer = molecule.GetConformer()
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol():2s} {position.x: .8f} {position.y: .8f} {position.z: .8f}"
        )
    lines.extend(
        [
            "",
            "--Link1--",
            "%nprocshared=4",
            "%mem=4GB",
            f"%chk={chk_name}",
            _ESP_ROUTE,
            "",
            "",
            "",
        ]
    )
    return "\n".join(lines)


def prepare_templates(
    output: Path,
    geometry_method: GeometryMethod,
    seeds: tuple[int, ...] = (20261020, 20261021),
) -> Path:
    """Create named ETKDGv3 conformers and Gaussian optimization inputs.

    Args:
        output: New directory to create for the generated templates.
        geometry_method: Required geometry method identifier. Select either the
            reviewed wB97XD/def2-TZVP or B3LYP-D3BJ/def2-TZVP route.
        seeds: Deterministic ETKDG random seeds, one per conformer.

    Returns:
        The created output directory.

    Raises:
        FileExistsError: If ``output`` already exists.
        ValueError: If a seed, method, graph, or stereochemical audit is invalid.
    """
    if geometry_method not in ("b3lyp-d3bj-def2tzvp", "wb97xd-def2tzvp"):
        raise ValueError(f"Unsupported geometry method: {geometry_method}")
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("Seeds must be a non-empty tuple of unique integers")
    if output.exists():
        raise FileExistsError(f"Output directory already exists: {output}")

    template, atom_indices = _build_template()
    generated: list[tuple[int, Chem.Mol]] = []
    psi_values = (-40.0, 140.0)
    if len(seeds) > len(psi_values):
        raise ValueError("At most two seeds are supported for the specified psi starts")
    for seed, psi in zip(seeds, psi_values[: len(seeds)], strict=True):
        molecule = Chem.Mol(template)
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = seed
        parameters.numThreads = 1
        if AllChem.EmbedMolecule(molecule, parameters) < 0:
            raise ValueError(f"ETKDGv3 embedding failed for seed {seed}")
        conformer = molecule.GetConformer()
        rdMolTransforms.SetDihedralDeg(
            conformer,
            atom_indices["N"],
            atom_indices["CA"],
            atom_indices["C"],
            atom_indices["NME_N"],
            psi,
        )
        _audit_coordinate_stereochemistry(molecule, atom_indices)
        generated.append((seed, molecule))

    output.mkdir(parents=True)
    conformer_records = []
    for number, (seed, molecule) in enumerate(generated, start=1):
        conformer_dir = output / f"conformer_{number:02d}"
        conformer_dir.mkdir()
        molecule.SetProp("_Name", f"neutral_thr_nme_seed_{seed}")
        molecule.SetProp("SMILES", Chem.MolToSmiles(molecule, isomericSmiles=True))
        molecule.SetProp("formal_charge", str(Chem.GetFormalCharge(molecule)))
        sdf_path = conformer_dir / "molecule.sdf"
        with Chem.SDWriter(str(sdf_path)) as writer:
            writer.SetProps(["SMILES", "formal_charge"])
            writer.write(molecule)
        _write_xyz(molecule, conformer_dir / "coordinates.xyz")
        _write_atom_metadata(molecule, conformer_dir / "atom_map.json")
        (conformer_dir / "input.gjf").write_text(
            _gaussian_input(molecule, "neutral_thr_nme.chk", geometry_method),
            encoding="utf-8",
        )
        conformer_records.append(
            {
                "directory": conformer_dir.name,
                "seed": seed,
                "psi_degrees": psi_values[number - 1],
                "rdkit_version": rdBase.rdkitVersion,
                "geometry_method": geometry_method,
            }
        )

    manifest = {
        "model": "neutral free-NH2 L-threonine methylamide",
        "isomeric_smiles": Chem.MolToSmiles(template, isomericSmiles=True),
        "formula": CalcMolFormula(template),
        "formal_charge": Chem.GetFormalCharge(template),
        "multiplicity": 1,
        "proton_inventory": {
            "N": 2,
            "CA": 1,
            "CB": 1,
            "OG1": 1,
            "CG2": 3,
            "NME_N": 1,
            "NME_C": 3,
        },
        "stereochemistry": _EXPECTED_CIP,
        "rdkit_version": rdBase.rdkitVersion,
        "embedding": "ETKDGv3; no energy minimization",
        "geometry_method": geometry_method,
        "esp_method": "HF/6-31G(d) Merz-Kollman, separate Link1 job step",
        "conformers": conformer_records,
        "parameter_status": "ESP charges are provisional; not a validated force field",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return output
