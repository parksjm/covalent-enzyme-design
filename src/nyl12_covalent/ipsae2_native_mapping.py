"""Strict native Boltz atom/token order qualification before ipSAE2 scoring."""

from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path
from typing import Any

import numpy as np


def artifact(path: Path) -> dict[str, str]:
    """Return the absolute source path and SHA256 of an existing artifact."""
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def compare_native_order(
    expected_atoms: list[tuple[str, int, str, str, int]],
    observed_atoms: list[tuple[str, int, str, str, float]],
    plddt: np.ndarray,
) -> None:
    """Validate every atom identity and its independently serialized token pLDDT.

    Args:
        expected_atoms: Native (chain, residue, residue name, atom name, token).
        observed_atoms: CIF encounter-order identities followed by B factor.
        plddt: Original token confidence array, on the 0–1 scale.

    Raises:
        ValueError: Any identity, ordering, atom count or token confidence differs.
    """
    if len(expected_atoms) != len(observed_atoms):
        raise ValueError("Native/CIF atom count mismatch")
    if not np.isfinite(plddt).all():
        raise ValueError("Nonfinite token pLDDT")
    for index, (expected, observed) in enumerate(zip(expected_atoms, observed_atoms)):
        if expected[:4] != observed[:4]:
            raise ValueError(f"Native/CIF atom order mismatch at atom {index}")
        if not np.isfinite(observed[4]):
            raise ValueError(f"Nonfinite native CIF B factor at atom {index}")
        # Native writer rounds pLDDT*100 to three decimals; float32 arithmetic
        # can shift a final rounding digit by a few millionths.
        if abs(float(plddt[expected[4]]) * 100 - observed[4]) > 0.00052:
            raise ValueError(f"Native/CIF pLDDT token binding mismatch at atom {index}")


def _native_arrays(path: Path, source_root: Path) -> tuple[Any, np.ndarray]:
    from boltz.data.tokenize.boltz2 import tokenize_structure
    from boltz.data.types import StructureV2

    if (
        Path(inspect.getfile(tokenize_structure)).resolve()
        != (source_root / "data/tokenize/boltz2.py").resolve()
    ):
        raise ValueError("Imported tokenizer differs from declared Boltz source")
    names = (
        "atoms",
        "bonds",
        "residues",
        "chains",
        "interfaces",
        "mask",
        "coords",
        "ensemble",
    )
    with np.load(path, allow_pickle=False) as arrays:
        structure = StructureV2(**{name: arrays[name] for name in names})
    if not structure.mask.all():
        raise ValueError("Invalid-chain removal needs separate qualification")
    tokens, _ = tokenize_structure(structure)
    return structure, tokens


def _cif_atoms(path: Path) -> list[tuple[str, int, str, str, float]]:
    import gemmi

    block = gemmi.cif.read_file(str(path)).sole_block()
    fields = [
        "label_asym_id",
        "auth_seq_id",
        "label_comp_id",
        "label_atom_id",
        "B_iso_or_equiv",
        "type_symbol",
        "pdbx_PDB_model_num",
    ]
    rows = block.find("_atom_site.", fields)
    result = []
    for row in rows:
        if row[5].upper() in {"H", "D"} or row[2] in {"HOH", "WAT"}:
            raise ValueError("Unexpected hydrogen/deuterium/water in native CIF")
        if row[6] != "1":
            raise ValueError("Expected a single native model")
        result.append((row[0], int(row[1]), row[2], row[3], float(row[4])))
    if not result:
        raise ValueError("No native CIF atoms")
    return result


def qualify_native_mapping(
    *,
    cif_path: Path,
    pae_path: Path,
    confidence_path: Path,
    plddt_path: Path,
    structure_path: Path,
    boltz_source_root: Path,
) -> dict[str, Any]:
    """Qualify the unmodified prediction triplet against native tokenizer replay.

    Run with the original Boltz Python environment. No coordinates are rewritten.
    This receipt establishes serialization identity, not physical protomer pairing.
    Missing artifacts raise explicitly; scientific/order failures return rejection.

    Args:
        cif_path: Original native prediction CIF.
        pae_path: Matching native PAE NPZ.
        confidence_path: Matching native confidence JSON.
        plddt_path: Matching native per-token pLDDT NPZ.
        structure_path: Original processed StructureV2 NPZ.
        boltz_source_root: Installed Boltz package whose writer produced the files.

    Returns:
        JSON-serializable receipt with acceptance, hashes and per-token mapping.
    """
    paths = dict(
        cif=cif_path,
        pae=pae_path,
        confidence=confidence_path,
        plddt=plddt_path,
        structure=structure_path,
    )
    sources = [
        "data/tokenize/boltz2.py",
        "data/types.py",
        "data/write/writer.py",
        "data/write/mmcif.py",
        "model/loss/confidence.py",
        "model/modules/confidencev2.py",
    ]
    receipt: dict[str, Any] = {
        "schema_version": 1,
        "accepted": False,
        "errors": [],
        "artifacts": {key: artifact(path) for key, path in paths.items()},
        "source_hashes": {
            str(boltz_source_root / p): artifact(boltz_source_root / p)["sha256"]
            for p in sources
        },
        "native_chain_map": [],
        "tokens": [],
        "protein_chains": [],
        "ligand_chains": [],
        "physical_groups_qualified": False,
        "pae_axis_convention": {
            "row": "aligned_frame_token",
            "column": "scored_token",
            "proof": "loss/confidence.py express_coordinate_in_frame: d=b[:,:,None,:,:]-b[:,:,:,None,:]; basis e[:,:,:,None,:]; confidencev2.py aggregates logits without transpose; writer.py saves prediction pae directly",
            "reader_transpose": False,
        },
        "policy": "native writer atom order + native tokenizer replay + per-token pLDDT/B-factor binding; ligand frames inferred by ipSAE2",
    }
    try:
        structure, tokens = _native_arrays(structure_path, boltz_source_root)
        with np.load(pae_path, allow_pickle=False) as data:
            pae = data["pae"]
        with np.load(plddt_path, allow_pickle=False) as data:
            plddt = data["plddt"]
        if pae.shape != (len(tokens), len(tokens)) or plddt.shape != (len(tokens),):
            raise ValueError("PAE/pLDDT/native token dimensions disagree")
        if not np.isfinite(pae).all() or (pae < 0).any():
            raise ValueError("Invalid native PAE values")
        confidence = json.loads(confidence_path.read_text())
        native_ids = [str(int(c["asym_id"])) for c in structure.chains]
        pair = confidence["pair_chains_iptm"]
        if set(pair) != set(native_ids) or any(
            set(row) != set(native_ids) for row in pair.values()
        ):
            raise ValueError("Confidence native chain IDs disagree")
        if set(confidence["chains_ptm"]) != set(native_ids):
            raise ValueError("Confidence chains_ptm IDs disagree")
        atom_tokens: dict[int, int] = {}
        for index, token in enumerate(tokens):
            if int(token["token_idx"]) != index:
                raise ValueError("Nonsequential native token indices")
            for atom_index in range(
                int(token["atom_idx"]), int(token["atom_idx"] + token["atom_num"])
            ):
                if atom_index in atom_tokens:
                    raise ValueError("Overlapping native token atoms")
                atom_tokens[atom_index] = index
        expected = []
        for chain in structure.chains:
            chain_id = str(chain["name"])
            native_id = int(chain["asym_id"])
            selected = np.flatnonzero(tokens["asym_id"] == native_id)
            if not len(selected) or not np.array_equal(
                selected, np.arange(selected[0], selected[-1] + 1)
            ):
                raise ValueError("Noncontiguous native chain tokens")
            receipt["native_chain_map"].append(
                dict(
                    label_asym_id=chain_id,
                    native_asym_id=native_id,
                    mol_type=int(chain["mol_type"]),
                    entity_id=int(chain["entity_id"]),
                    sym_id=int(chain["sym_id"]),
                    token_start=int(selected[0]),
                    token_stop=int(selected[-1] + 1),
                    atom_count=int(chain["atom_num"]),
                    residue_count=int(chain["res_num"]),
                )
            )
            key = {0: "protein_chains", 3: "ligand_chains"}.get(int(chain["mol_type"]))
            if key is None:
                raise ValueError("Unsupported native molecule type")
            receipt[key].append(chain_id)
            residues = structure.residues[
                int(chain["res_idx"]) : int(chain["res_idx"] + chain["res_num"])
            ]
            for residue in residues:
                start, count = int(residue["atom_idx"]), int(residue["atom_num"])
                for atom_index in range(start, start + count):
                    atom = structure.atoms[atom_index]
                    expected.append(
                        (
                            chain_id,
                            int(residue["res_idx"]) + 1,
                            str(residue["name"]),
                            str(atom["name"]),
                            atom_tokens[atom_index],
                        )
                    )
            for index in selected:
                token = tokens[index]
                start, count = int(token["atom_idx"]), int(token["atom_num"])
                receipt["tokens"].append(
                    dict(
                        token_index=int(index),
                        native_asym_id=native_id,
                        label_asym_id=chain_id,
                        residue_index=int(token["res_idx"]) + 1,
                        residue_name=str(token["res_name"]),
                        atom_start=start,
                        atom_count=count,
                        atom_names=structure.atoms["name"][
                            start : start + count
                        ].tolist(),
                    )
                )
        compare_native_order(expected, _cif_atoms(cif_path), plddt)
        receipt["pae_shape"] = list(pae.shape)
        receipt["accepted"] = True
    except (ValueError, KeyError, IndexError) as error:
        receipt["errors"].append(str(error))
    return receipt
