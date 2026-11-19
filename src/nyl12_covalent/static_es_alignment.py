"""Compare a canonical static Nyl12 E-S PDB with the catalytic reference."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np

from nyl12_covalent.catalytic_alignment import (
    BACKBONE,
    FIT_RESIDUES,
    SITE_RESIDUES,
    _cif_atoms,
    _pdb_line,
)
from nyl12_covalent.catalytic_alignment_geometry import fit_coordinates, state_metrics

LOGGER = logging.getLogger(__name__)


def read_canonical_pdb_site(
    path: Path,
    chain: str,
    residues: set[int],
) -> dict[tuple[int, str], tuple[str, np.ndarray]]:
    """Read heavy atoms for selected residues from one canonical PDB chain.

    Args:
        path: Input PDB file.
        chain: Canonical protein chain identifier.
        residues: Canonical residue numbers to retain.

    Returns:
        Heavy atoms keyed by ``(residue number, atom name)``.

    Raises:
        ValueError: If coordinates, alternate locations, or atom identities are
            ambiguous or invalid.
    """
    selected: dict[tuple[int, str], tuple[str, np.ndarray]] = {}
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        if line[:6] not in {"ATOM  ", "HETATM"} or line[21:22] != chain:
            continue
        try:
            residue = int(line[22:26])
        except ValueError as error:
            raise ValueError(
                f"Invalid PDB residue number at line {line_number}"
            ) from error
        if residue not in residues:
            continue
        alternate = line[16:17]
        if alternate not in {" ", "A"}:
            raise ValueError(
                f"Unsupported alternate location {alternate!r} at line {line_number}"
            )
        atom = line[12:16].strip()
        element = line[76:78].strip().upper()
        if not element:
            element = atom.lstrip("0123456789")[:1].upper()
        if element == "H":
            continue
        try:
            xyz = np.array([float(line[start : start + 8]) for start in (30, 38, 46)])
        except ValueError as error:
            raise ValueError(f"Invalid PDB coordinate at line {line_number}") from error
        if not np.isfinite(xyz).all():
            raise ValueError(f"Non-finite PDB coordinate at line {line_number}")
        key = (residue, atom)
        if key in selected:
            raise ValueError(f"Duplicate canonical PDB atom {chain}:{residue}:{atom}")
        selected[key] = (element, xyz)
    if not selected:
        raise ValueError(f"No selected atoms found for chain {chain} in {path}")
    return selected


def align_static_es_model(
    reference_cif: Path,
    model_pdb: Path,
    output_dir: Path,
    model_chain: str = "A",
) -> dict[str, Any]:
    """Align one canonical static E-S model and report Thr233 geometry.

    The fit excludes Thr233 and uses the N, CA, C, and O atoms of Phe113,
    Asn186, Asp272, and Asp274. The same held-out Thr233 and contact metrics
    used for the solvated QM/MM E-S comparison are then evaluated.

    Args:
        reference_cif: Supplied Nyl12 structural reference.
        model_pdb: OpenMM-relaxed, canonical-chain E-S PDB.
        output_dir: New directory for the report and overlay.
        model_chain: Canonical active-site chain to assess.

    Returns:
        JSON-compatible alignment report.
    """
    residue_numbers = {auth for _, auth, _ in SITE_RESIDUES}
    reference = _cif_atoms(reference_cif, "A", residue_numbers)
    model = read_canonical_pdb_site(model_pdb, model_chain, residue_numbers)
    fit_keys = [(auth, atom) for _, auth, _ in FIT_RESIDUES for atom in BACKBONE]
    for key in fit_keys:
        if key not in reference or key not in model:
            raise ValueError(f"Missing required fit atom {key}")
        if reference[key][0] != model[key][0]:
            raise ValueError(f"Element mismatch for fit atom {key}")
    rotation, mobile_center, reference_center, fit_rmsd = fit_coordinates(
        np.array([model[key][1] for key in fit_keys]),
        np.array([reference[key][1] for key in fit_keys]),
    )
    aligned = {
        key: (element, (coordinate - mobile_center) @ rotation + reference_center)
        for key, (element, coordinate) in model.items()
    }
    metrics = state_metrics(reference, aligned)
    output_dir.mkdir(parents=True, exist_ok=True)
    overlay_path = output_dir / "catalytic_site_overlay.pdb"
    report_path = output_dir / "alignment_report.json"
    if overlay_path.exists() or report_path.exists():
        raise FileExistsError(f"Refusing to overwrite outputs in {output_dir}")
    lines = [
        "REMARK 1 Chain A is n12c.cif reference; B is aligned static E-S model.\n",
        "REMARK 2 Thr233 was excluded from the catalytic-scaffold fit.\n",
    ]
    serial = 1
    for output_chain, atoms in (("A", reference), ("B", aligned)):
        for _, residue, residue_name in SITE_RESIDUES:
            for (atom_residue, atom), (element, coordinate) in sorted(atoms.items()):
                if atom_residue != residue:
                    continue
                lines.append(
                    _pdb_line(
                        serial,
                        "ATOM",
                        atom,
                        residue_name,
                        output_chain,
                        residue,
                        coordinate,
                        element,
                    )
                )
                serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    overlay_path.write_text("".join(lines), encoding="utf-8")
    report: dict[str, Any] = {
        "interpretation": (
            "Structural seed screen only; Tier A and native-like Thr233 geometry "
            "do not establish a reaction mechanism or transition state."
        ),
        "inputs": {
            "reference_cif": str(reference_cif),
            "model_pdb": str(model_pdb),
            "model_chain": model_chain,
        },
        "fit": {
            "residues": [auth for _, auth, _ in FIT_RESIDUES],
            "atoms_per_residue": list(BACKBONE),
            "atom_count": len(fit_keys),
            "rmsd_angstrom": fit_rmsd,
        },
        **metrics,
        "outputs": {"overlay_pdb": str(overlay_path), "report_json": str(report_path)},
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    LOGGER.info("Wrote static E-S catalytic alignment to %s", output_dir)
    return report
