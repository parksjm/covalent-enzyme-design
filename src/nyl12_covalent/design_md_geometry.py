"""Fixed-mapping four-site geometry for preliminary PA6 design trajectories."""

from __future__ import annotations

import itertools
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
from parmed import unit
from parmed.geometry import box_lengths_and_angles_to_vectors

ROOT = Path(__file__).resolve().parents[2]
SCREEN = (
    ROOT
    / "design/nyl12_ligandmpnn_pilot/ensemble_mask_comparison/validation_source_snapshot_checked/scripts"
)
sys.path.insert(0, str(SCREEN))
from screen_pa6_binding_orientation import measure  # noqa: E402
from nyl12_covalent.catalytic_alignment import SITE_RESIDUES, _cif_atoms  # noqa: E402


def image_near(points: np.ndarray, origin: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Image points near an origin in a periodic triclinic cell.

    Args:
        points: Cartesian point array of shape N by 3.
        origin: Cartesian reference position.
        box: Six cell lengths and angles.

    Returns:
        Minimum-image Cartesian points near the reference.
    """
    if box.shape != (6,) or not np.isfinite(box).all() or (box <= 0).any():
        raise ValueError("Invalid periodic box")
    cell = np.asarray(
        box_lengths_and_angles_to_vectors(*box).value_in_unit(unit.angstrom)
    )
    frac = (points - origin) @ np.linalg.inv(cell)
    frac -= np.round(frac)
    offsets = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    candidates = (frac[:, None, :] + offsets[None, :, :]) @ cell
    norms = np.sum(candidates * candidates, axis=2)
    return origin + candidates[np.arange(len(points)), np.argmin(norms, axis=1)]


def index_mapping(topology: Any, source_dir: Path) -> list[dict[str, Any]]:
    """Resolve immutable protein and scissile-amide indices for all four sites.

    Args:
        topology: Loaded original Amber topology.
        source_dir: Original relaxation directory containing conversion records.

    Returns:
        Four explicit protein and source-atom index mappings.
    """
    maps = []
    for i, chain in enumerate("IJKL"):
        protein = {}
        for _, auth, _ in SITE_RESIDUES:
            resno = auth + 232 * i if auth <= 232 else 929 + 96 * i + auth - 233
            res = topology.residues[resno - 1]
            for atom in res.atoms:
                if atom.atomic_number != 1:
                    key = (auth, atom.name)
                    if key in protein:
                        raise ValueError("Duplicate protein atom name")
                    protein[key] = (atom.idx, atom.element_name)
        conversion = json.loads(
            (source_dir / "system/components" / chain / "conversion.json").read_text()
        )["atom_mapping"]
        ligand = {}
        for name in ("C54", "O49", "N59", "C62"):
            hits = [r for r in conversion if r["source_atom"] == name]
            if len(hits) != 1:
                raise ValueError("Ambiguous original amide mapping")
            row = hits[0]
            # Each conversion record numbers its single component from 1313;
            # the assembled topology concatenates four four-residue ligands.
            res = topology.residues[row["target_residue"] + 4 * i - 1]
            atoms = [a for a in res.atoms if a.name == row["target_atom"]]
            if len(atoms) != 1:
                raise ValueError("Missing mapped ligand atom")
            ligand[name] = atoms[0].idx
        maps.append({"protein": protein, "ligand": ligand})
    return maps


def frame_rows(
    xyz: np.ndarray, box: np.ndarray, mapping: list[dict[str, Any]], reference: dict
) -> list[dict[str, Any]]:
    """Measure the four fixed assigned sites using the frozen static classifier.

    Args:
        xyz: Full topology-ordered frame coordinates.
        box: Periodic cell lengths and angles.
        mapping: Immutable four-site protein and ligand atom mappings.
        reference: Native CIF atom dictionary.

    Returns:
        Four catalytic and held-out native-Thr geometry records.
    """
    result = []
    for site in mapping:
        p = site["protein"]
        ids = [value[0] for value in p.values()]
        origin = xyz[p[(233, "OG1")][0]]
        positions = image_near(xyz[ids], origin, box)
        protein = {key: (p[key][1], pos) for key, pos in zip(p, positions, strict=True)}
        lig = site["ligand"]
        names = list(lig)
        coords = image_near(xyz[[lig[n] for n in names]], origin, box)
        row = measure(protein, dict(zip(names, coords, strict=True)), reference)
        result.append(row)
    return result


def native_reference(path: Path) -> dict:
    """Read the frozen canonical native alignment and Thr reference.

    Args:
        path: Canonical n12c mmCIF path.

    Returns:
        Canonical residue and atom coordinate mapping.
    """
    return _cif_atoms(path, "A", {r for _, r, _ in SITE_RESIDUES})
