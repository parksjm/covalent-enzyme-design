"""Geometry, ledger, and VMD helpers for the incremental water path exporter."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any

import numpy as np


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_index(structure: Any, residue: int, name: str) -> int:
    """Return the unique zero-based index for a residue atom."""
    matches = [
        atom.idx for atom in structure.residues[residue - 1].atoms if atom.name == name
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one {residue}:{name}; found {len(matches)}")
    return matches[0]


def dist(
    xyz: np.ndarray, structure: Any, a: tuple[int, str], b: tuple[int, str]
) -> float:
    """Calculate an interatomic distance in angstroms."""
    return float(
        np.linalg.norm(xyz[atom_index(structure, *a)] - xyz[atom_index(structure, *b)])
    )


def ang(
    xyz: np.ndarray,
    structure: Any,
    a: tuple[int, str],
    b: tuple[int, str],
    c: tuple[int, str],
) -> float:
    """Calculate the A-B-C angle in degrees."""
    va = xyz[atom_index(structure, *a)] - xyz[atom_index(structure, *b)]
    vc = xyz[atom_index(structure, *c)] - xyz[atom_index(structure, *b)]
    cosine = np.dot(va, vc) / (np.linalg.norm(va) * np.linalg.norm(vc))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def plot_svg(
    path: Path,
    title: str,
    labels: list[str],
    series: dict[str, list[float]],
    ylabel: str,
) -> None:
    """Write a dependency-free SVG line plot."""
    width, height, left, top, pw, ph = 900, 460, 85, 55, 770, 315
    vals = [v for row in series.values() for v in row if np.isfinite(v)]
    low, high = min(vals), max(vals)
    if low == high:
        low, high = low - 1.0, high + 1.0
    pad = (high - low) * 0.08
    low, high = low - pad, high + pad

    def xpos(i: int) -> float:
        return left + pw * i / max(1, len(labels) - 1)

    def ypos(value: float) -> float:
        return top + ph * (high - value) / (high - low)

    colors = ("#1769aa", "#d95f02", "#2b8c4b", "#8b4c9e", "#555555")
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="450" y="28" text-anchor="middle" font-family="sans-serif" font-size="18">{title}</text>',
    ]
    for tick in np.linspace(low, high, 6):
        yy = ypos(float(tick))
        out.extend(
            [
                f'<path d="M{left} {yy:.1f}H{left + pw}" stroke="#ddd"/>',
                f'<text x="{left - 8}" y="{yy + 4:.1f}" text-anchor="end" font-size="11">{tick:.3g}</text>',
            ]
        )
    for n, (name, values) in enumerate(series.items()):
        points = " ".join(
            f"{xpos(i):.1f},{ypos(v):.1f}"
            for i, v in enumerate(values)
            if np.isfinite(v)
        )
        color = colors[n % len(colors)]
        out.append(
            f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>'
        )
        out.append(
            f'<text x="{left + n * 145}" y="425" font-size="12" fill="{color}">{name}</text>'
        )
    for i, label in enumerate(labels):
        if (
            len(labels) < 20
            or i in (0, len(labels) - 1)
            or i % max(1, len(labels) // 10) == 0
        ):
            out.append(
                f'<text x="{xpos(i):.1f}" y="395" text-anchor="middle" font-size="9">{label}</text>'
            )
    out.extend(
        [
            f'<text x="20" y="{top + ph / 2}" transform="rotate(-90 20 {top + ph / 2})" text-anchor="middle" font-size="12">{ylabel}</text>',
            "</svg>",
        ]
    )
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def read_energy(path: Path) -> dict[str, dict[str, float]]:
    """Read unique physical and bias-inclusive energy values for each window."""
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    required = {
        "window",
        "energy_kcal_mol",
        "energy_bias_inclusive_kcal_mol",
    }
    if not rows or not required <= rows[0].keys():
        raise ValueError(f"Energy ledger requires {sorted(required)} columns")
    result = {
        row["window"]: {
            "energy_kcal_mol": float(row["energy_kcal_mol"]),
            "energy_bias_inclusive_kcal_mol": float(
                row["energy_bias_inclusive_kcal_mol"]
            ),
        }
        for row in rows
    }
    if len(result) != len(rows):
        raise ValueError("Duplicate window in energy ledger")
    return result


def vmd_script(
    frame_count: int,
    qm_indices: list[int],
    atoms: dict[str, int],
    *,
    topology: Any,
    manifest: dict[str, Any],
) -> str:
    """Adapt this Asp272 water route to the shared topology-mapped renderer.

    Args:
        frame_count: Positive accepted-frame count.
        qm_indices: QM indices in the explicitly stripped topology.
        atoms: Named stripped indices from the original source topology.
        topology: Original ParmEd topology with bonded substrate connectivity.
        manifest: Explicit stripped-to-original zero-based atom map.

    Returns:
        Portable VMD Tcl with seven route-specific distances and BD angle.

    Raises:
        ValueError: If the frame count or QM102/54H inventory is invalid.
    """
    import sys

    root = next(
        p
        for p in Path(__file__).resolve().parents
        if (p / "src/nyl12_covalent").is_dir()
    )
    sys.path.insert(0, str(root / "src"))
    from nyl12_covalent.relay_progress_vmd import (
        VMD_API_REFERENCES,
        vmd_script as render,
    )

    if frame_count <= 0:
        raise ValueError("At least one accepted frame is required")
    keep = manifest["stripped_to_original_atom_indices_0based"]
    qm_ids = [keep[i] + 1 for i in qm_indices]
    if (
        len(qm_ids) != 102
        or sum(topology.atoms[i - 1].atomic_number == 1 for i in qm_ids) != 54
    ):
        raise ValueError("Expected QM102 including54real hydrogens")
    keys = {
        "thr_og": "og1",
        "thr_n": "thr_n",
        "proton": "hg1",
        "water_o": "water_o",
        "water_h1": "water_h1",
        "water_h2": "water_h2",
        "asp_od1": "asp_od1",
        "carbon": "carbonyl_c",
        "oxygen": "carbonyl_o",
        "leaving_n": "leaving_n",
        "phe_n": "phe113_n",
        "asn_nd2": "asn_nd2",
    }
    case = {
        "qm_atom_ids": qm_ids,
        "atom_map": {key: keep[atoms[value]] + 1 for key, value in keys.items()},
    }
    pairs = (
        ("thr_og", "proton"),
        ("proton", "water_o"),
        ("water_o", "water_h2"),
        ("water_h2", "asp_od1"),
        ("thr_og", "carbon"),
        ("oxygen", "phe_n"),
        ("oxygen", "asn_nd2"),
    )
    renderer = root / "src/nyl12_covalent/relay_progress_vmd.py"
    manifest.setdefault("source_code", {})[str(renderer)] = sha256(renderer)
    manifest["view_rendering"] = {
        "route": "OG1-HG1-waterO-H2-Asp272OD1; seven distances plus BD",
        "distance_atom_ids_1based": [
            [case["atom_map"][a], case["atom_map"][b]] for a, b in pairs
        ],
        "vmd_api_references": VMD_API_REFERENCES,
        "label_contract": "molID/atomID strings; empty add result; integer indices; per-label textoffset",
        "camera_contract": "molinfo attribute/value lists; one nested4x4matrix per attribute",
        "hydrogens": "All54realQMH dynamic1.4A; waterH-H constraint hidden; C-H topology retained",
        "bond_rendering": "BaseElement allvisibleendpoints; orangeMMcarbonoverlay; DynamicBonds cylindersresolution12",
        "substrate": "Connected mapped P6N/P6I/P6C component; explicit original-to-stripped atom indices",
        "validation_limit": "VMD unavailable; strict API-contract Tcl and topology tests",
    }
    return render(manifest, case, topology, distance_pairs=pairs)


def render_manifest(topology: Any, manifest: dict[str, Any]) -> str:
    """Resolve this route's source topology atoms and render its existing bundle.

    Args:
        topology: Original ParmEd topology.
        manifest: Existing or newly generated bundle atom-map metadata.

    Returns:
        Seven-distance Asp272-route viewer without residue-offset selections.
    """
    keep = manifest["stripped_to_original_atom_indices_0based"]
    index_map = {old: new for new, old in enumerate(keep)}
    named = {
        "og1": (929, "OG1"),
        "thr_n": (929, "N"),
        "hg1": (929, "HG1"),
        "water_o": (1576, "O"),
        "water_h1": (1576, "H1"),
        "water_h2": (1576, "H2"),
        "asp_od1": (968, "OD1"),
        "carbonyl_c": (1314, "C"),
        "carbonyl_o": (1314, "O"),
        "leaving_n": (1315, "N"),
        "phe113_n": (113, "N"),
        "asn_nd2": (186, "ND2"),
    }
    atoms = {
        key: index_map[atom_index(topology, *identity)]
        for key, identity in named.items()
    }
    return vmd_script(
        manifest["frames"],
        manifest["qm_indices_stripped_0based"],
        atoms,
        topology=topology,
        manifest=manifest,
    )
