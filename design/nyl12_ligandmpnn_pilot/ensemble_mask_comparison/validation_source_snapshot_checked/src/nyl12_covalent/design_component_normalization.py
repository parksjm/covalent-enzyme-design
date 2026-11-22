"""Normalize Boltz PA6 labels by an exact heavy-atom graph isomorphism."""

from __future__ import annotations

import hashlib
import json
import math
import pickle
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

_LIGAND_CHAINS = "IJKL"
_EXPECTED_HEAVY_FORMULA = {"C": 24, "N": 4, "O": 5}
_ROLE_BOND_ORDERS = {"C54-O49": "DOUBLE", "C54-N59": "SINGLE"}


def _load_component_molecule(component_pickle: Path) -> Any:
    """Load the trusted Boltz pickle and return its named ``LIG1`` graph."""
    with component_pickle.open("rb") as stream:
        payload = pickle.load(stream)  # noqa: S301 - trusted project Boltz result
    if not isinstance(payload, dict) or "LIG1" not in payload:
        raise ValueError("Boltz molecule pickle must contain the LIG1 component")
    molecule = payload["LIG1"]
    if not hasattr(molecule, "GetAtoms"):
        raise ValueError("Boltz LIG1 component is not an RDKit molecule")
    return molecule


def _load_reference_molecule(reference_mol2: Path) -> Any:
    """Load pinned PA6 MOL2 with its original Tripos atom-name properties."""
    from rdkit import Chem

    molecule = Chem.MolFromMol2File(str(reference_mol2), removeHs=False, sanitize=True)
    if molecule is None:
        raise ValueError(f"Cannot parse reference PA6 MOL2: {reference_mol2}")
    return molecule


def _named_heavy_atoms(molecule: Any, property_name: str, label: str) -> list[Any]:
    """Return heavy atoms after enforcing a complete, unique name property."""
    from rdkit import Chem

    heavy = Chem.RemoveHs(molecule)
    atoms = list(heavy.GetAtoms())
    if len(atoms) != 33 or any(atom.GetAtomicNum() == 1 for atom in atoms):
        raise ValueError(f"{label} must contain exactly 33 heavy atoms")
    if any(not atom.HasProp(property_name) for atom in atoms):
        raise ValueError(f"{label} has an atom without {property_name!r}")
    names = [atom.GetProp(property_name) for atom in atoms]
    if any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError(f"{label} atom names must be non-empty and unique")
    if any(len(name) > 4 for name in names):
        raise ValueError(f"{label} atom names must fit the PDB atom-name field")
    return atoms


def _heavy_formula(atoms: list[Any]) -> dict[str, int]:
    """Count heavy-atom elements in a named graph."""
    return dict(sorted(Counter(atom.GetSymbol() for atom in atoms).items()))


def _exact_graph_match(source: Any, reference: Any, match: tuple[int, ...]) -> bool:
    """Check atom identity, formal charge, bond order, and stereo for one map."""
    if len(match) != source.GetNumAtoms() or len(set(match)) != reference.GetNumAtoms():
        return False
    for source_index, reference_index in enumerate(match):
        source_atom = source.GetAtomWithIdx(source_index)
        reference_atom = reference.GetAtomWithIdx(reference_index)
        if (
            source_atom.GetAtomicNum() != reference_atom.GetAtomicNum()
            or source_atom.GetFormalCharge() != reference_atom.GetFormalCharge()
        ):
            return False
    for left in range(source.GetNumAtoms()):
        for right in range(left + 1, source.GetNumAtoms()):
            source_bond = source.GetBondBetweenAtoms(left, right)
            reference_bond = reference.GetBondBetweenAtoms(match[left], match[right])
            if (source_bond is None) != (reference_bond is None):
                return False
            if source_bond is not None and reference_bond is not None:
                if (
                    source_bond.GetBondType() != reference_bond.GetBondType()
                    or source_bond.GetIsAromatic() != reference_bond.GetIsAromatic()
                    or source_bond.GetStereo() != reference_bond.GetStereo()
                ):
                    return False
    return True


def _unique_graph_match(
    source: Any, reference: Any, matches: tuple[tuple[int, ...], ...]
) -> tuple[int, ...]:
    """Require a single exact graph isomorphism from candidate RDKit matches."""
    exact_matches = [
        match for match in matches if _exact_graph_match(source, reference, match)
    ]
    if len(exact_matches) != 1:
        raise ValueError(
            "Expected one exact PA6 heavy-graph isomorphism; "
            f"found {len(exact_matches)}"
        )
    return exact_matches[0]


def _match_atom_labels(source_mol: Any, reference_mol: Any) -> dict[str, Any]:
    """Return unique source-to-reference atom labels after full graph matching.

    Args:
        source_mol: Boltz ``LIG1`` molecule; atom names use ``name``.
        reference_mol: Pinned whole-PA6 MOL2; atom names use
            ``_TriposAtomName``.

    Returns:
        Exact atom map and graph/charge metadata.

    Raises:
        ValueError: If the heavy graphs differ or have multiple isomorphisms.
    """
    from rdkit import Chem

    source = Chem.RemoveHs(source_mol)
    reference = Chem.RemoveHs(reference_mol)
    source_atoms = _named_heavy_atoms(source, "name", "Boltz component")
    reference_atoms = _named_heavy_atoms(reference, "_TriposAtomName", "PA6 reference")
    if _heavy_formula(source_atoms) != _EXPECTED_HEAVY_FORMULA:
        raise ValueError(
            f"Unexpected Boltz PA6 heavy formula: {_heavy_formula(source_atoms)}"
        )
    if _heavy_formula(reference_atoms) != _EXPECTED_HEAVY_FORMULA:
        raise ValueError(
            f"Unexpected reference PA6 heavy formula: {_heavy_formula(reference_atoms)}"
        )
    source_charge = sum(atom.GetFormalCharge() for atom in source_atoms)
    reference_charge = sum(atom.GetFormalCharge() for atom in reference_atoms)
    if source_charge != 0 or reference_charge != 0:
        raise ValueError(
            f"PA6 heavy-graph formal charges must be zero; found {source_charge} and {reference_charge}"
        )
    if source.GetNumBonds() != reference.GetNumBonds():
        raise ValueError("Boltz and pinned PA6 heavy graphs have different bond counts")
    parameters = Chem.SubstructMatchParameters()
    parameters.uniquify = False
    parameters.useChirality = True
    matches = reference.GetSubstructMatches(source, parameters)
    match = _unique_graph_match(source, reference, matches)
    source_names = [atom.GetProp("name") for atom in source_atoms]
    reference_names = [atom.GetProp("_TriposAtomName") for atom in reference_atoms]
    if not {"C54", "O49", "N59"}.issubset(reference_names):
        raise ValueError("Pinned PA6 reference lacks the C54/O49/N59 amide role names")
    source_to_reference = {
        source_names[source_index]: reference_names[reference_index]
        for source_index, reference_index in enumerate(match)
    }
    reference_index_by_name = {
        name: index for index, name in enumerate(reference_names)
    }
    roles = {name: reference_index_by_name[name] for name in ("C54", "O49", "N59")}
    carbonyl = reference.GetBondBetweenAtoms(roles["C54"], roles["O49"])
    amide = reference.GetBondBetweenAtoms(roles["C54"], roles["N59"])
    if (
        carbonyl is None
        or amide is None
        or str(carbonyl.GetBondType()) != _ROLE_BOND_ORDERS["C54-O49"]
        or str(amide.GetBondType()) != _ROLE_BOND_ORDERS["C54-N59"]
        or reference.GetAtomWithIdx(roles["C54"]).GetAtomicNum() != 6
        or reference.GetAtomWithIdx(roles["O49"]).GetAtomicNum() != 8
        or reference.GetAtomWithIdx(roles["N59"]).GetAtomicNum() != 7
    ):
        raise ValueError("Pinned PA6 C54/O49/N59 do not define the expected amide")
    source_by_name = {atom.GetProp("name"): atom for atom in source_atoms}
    bond_mapping = []
    for bond in source.GetBonds():
        source_first = source_names[bond.GetBeginAtomIdx()]
        source_second = source_names[bond.GetEndAtomIdx()]
        target_first = source_to_reference[source_first]
        target_second = source_to_reference[source_second]
        bond_mapping.append(
            {
                "source_atoms": [source_first, source_second],
                "target_atoms": [target_first, target_second],
                "bond_order": str(bond.GetBondType()),
                "aromatic": bool(bond.GetIsAromatic()),
                "stereo": str(bond.GetStereo()),
            }
        )
    return {
        "source_to_reference": source_to_reference,
        "atom_mapping": [
            {
                "source_atom": source_name,
                "target_atom": target_name,
                "element": source_by_name[source_name].GetSymbol(),
                "formal_charge": source_by_name[source_name].GetFormalCharge(),
            }
            for source_name, target_name in sorted(source_to_reference.items())
        ],
        "graph": {
            "heavy_atom_count": len(source_atoms),
            "heavy_bond_count": source.GetNumBonds(),
            "heavy_formula": _EXPECTED_HEAVY_FORMULA,
            "formal_charge": source_charge,
            "exact_isomorphism_count": 1,
        },
        "bond_mapping": bond_mapping,
        "role_mapping": {
            "C54": source_by_name[
                next(
                    name
                    for name, target in source_to_reference.items()
                    if target == "C54"
                )
            ].GetProp("name"),
            "O49": next(
                name for name, target in source_to_reference.items() if target == "O49"
            ),
            "N59": next(
                name for name, target in source_to_reference.items() if target == "N59"
            ),
            "amide_bonds": _ROLE_BOND_ORDERS,
        },
        "source_elements": {
            name: atom.GetSymbol() for name, atom in source_by_name.items()
        },
        "source_formal_charges": {
            name: atom.GetFormalCharge() for name, atom in source_by_name.items()
        },
    }


def _sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalize_pdb_lines(
    lines: list[str], mapping: dict[str, str], elements: dict[str, str]
) -> tuple[list[str], list[dict[str, Any]]]:
    """Rename exactly the four 33-atom PA6 components while retaining XYZ."""
    atom_lines = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if any(line[16].strip() or line[26].strip() for line in atom_lines):
        raise ValueError("Alternate locations or insertion codes are unsupported")
    if {line[21] for line in atom_lines} != set("ABCDEFGHIJKL"):
        raise ValueError("Boltz PDB must contain exactly chains A-L")
    for line in atom_lines:
        element = line[76:78].strip().upper()
        if element not in {"C", "N", "O", "S", "H"}:
            raise ValueError(f"Unsupported PDB element {element!r}")
        xyz = tuple(float(line[start : start + 8]) for start in (30, 38, 46))
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("Boltz PDB contains non-finite coordinates")
    normalized: list[str] = []
    records: list[dict[str, Any]] = []
    for chain in _LIGAND_CHAINS:
        ligand = [line for line in atom_lines if line[21] == chain]
        if len(ligand) != 33 or any(
            line[76:78].strip().upper() == "H" for line in ligand
        ):
            raise ValueError(
                f"Ligand chain {chain} must contain exactly 33 heavy atoms"
            )
        if len({line[22:27] for line in ligand}) != 1:
            raise ValueError(f"Ligand chain {chain} must occupy one source residue")
        names = [line[12:16].strip() for line in ligand]
        if len(set(names)) != 33 or set(names) != set(mapping):
            raise ValueError(
                f"Ligand chain {chain} source atom labels differ from Boltz LIG1"
            )
        for line in ligand:
            source_name = line[12:16].strip()
            target_name = mapping[source_name]
            element = line[76:78].strip().upper()
            if element != elements[source_name]:
                raise ValueError(
                    f"PDB element mismatch for {chain}:{source_name}: {element}"
                )
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")) or line[21] not in _LIGAND_CHAINS:
            normalized.append(line)
            continue
        source_name = line[12:16].strip()
        target_name = mapping[source_name]
        chain = line[21]
        target_line = (
            f"{line[:12]}{target_name:>4s}{line[16:17]}{line[17:20]} "
            f"{chain}{1:4d}{line[26:]}"
        )
        target_line = f"{target_line[:17]}LIG{target_line[20:]}"
        normalized.append(target_line)
        records.append(
            {
                "source_chain": chain,
                "source_serial": int(line[6:11]),
                "source_residue": line[17:26].strip(),
                "source_atom": source_name,
                "target_chain": chain,
                "target_residue": 1,
                "target_atom": target_name,
                "element": elements[source_name],
            }
        )
    if len(records) != 132:
        raise ValueError(f"Expected 132 normalized PA6 PDB atoms, found {len(records)}")
    return normalized, records


def normalize_model(
    source_pdb: Path,
    component_pickle: Path,
    reference_mol2: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Graph-map Boltz atom labels to pinned PA6 names for chains I-L.

    Args:
        source_pdb: Boltz PDB with canonical protein chains A-H and PA6 I-L.
        component_pickle: Trusted own Boltz results pickle containing ``LIG1``.
        reference_mol2: Pinned whole-L4 MOL2 with Tripos atom names.
        output_dir: New destination directory for normalized PDB and mapping.

    Returns:
        Mapping manifest with source hashes, exact graph evidence, and role map.
    """
    if output_dir.exists():
        raise FileExistsError(output_dir)
    source_molecule = _load_component_molecule(component_pickle)
    reference_molecule = _load_reference_molecule(reference_mol2)
    graph_map = _match_atom_labels(source_molecule, reference_molecule)
    lines = source_pdb.read_text(encoding="utf-8").splitlines()
    normalized_lines, pdb_mapping = _normalize_pdb_lines(
        lines, graph_map["source_to_reference"], graph_map["source_elements"]
    )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_dir.name}.", dir=output_dir.parent
    ) as temp_dir:
        stage = Path(temp_dir)
        normalized_path = stage / "normalized.pdb"
        normalized_path.write_text("\n".join(normalized_lines) + "\n", encoding="utf-8")
        manifest: dict[str, Any] = {
            "status": "accepted_graph_isomorphic_component_normalization",
            "source_pdb": str(source_pdb),
            "component_pickle": str(component_pickle),
            "reference_mol2": str(reference_mol2),
            "input_sha256": {
                str(path): _sha256(path)
                for path in (source_pdb, component_pickle, reference_mol2)
            },
            "graph": graph_map["graph"],
            "role_mapping": graph_map["role_mapping"],
            "atom_mapping": graph_map["atom_mapping"],
            "pdb_atom_mapping": pdb_mapping,
            "normalized_residue": {"name": "LIG", "number": 1},
            "protein_lines_unchanged": True,
            "ligand_coordinates_preserved": True,
            "output": str(output_dir / normalized_path.name),
            "output_sha256": _sha256(normalized_path),
        }
        (stage / "mapping.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        stage.rename(output_dir)
    return manifest
