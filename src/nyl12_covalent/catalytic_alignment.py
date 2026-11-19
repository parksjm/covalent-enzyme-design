"""Align an Amber E-S structure to a CIF catalytic scaffold."""

from __future__ import annotations

import json
import logging
import shlex
from pathlib import Path
from typing import Any

import numpy as np
import parmed

from nyl12_covalent.catalytic_alignment_geometry import (
    fit_coordinates,
    state_metrics,
)

LOGGER = logging.getLogger(__name__)
FIT_RESIDUES = (
    (113, 113, "PHE"),
    (186, 186, "ASN"),
    (968, 272, "ASP"),
    (970, 274, "ASP"),
)
SITE_RESIDUES = (
    (113, 113, "PHE"),
    (186, 186, "ASN"),
    (929, 233, "THR"),
    (930, 234, "THR"),
    (968, 272, "ASP"),
    (970, 274, "ASP"),
)
BACKBONE = ("N", "CA", "C", "O")


def _cif_atoms(
    path: Path,
    chain: str,
    residues: set[int],
) -> dict[tuple[int, str], tuple[str, np.ndarray]]:
    """Read selected atom_site rows from a simple one-row-per-atom mmCIF loop."""
    lines = path.read_text(encoding="utf-8").splitlines()
    for start, line in enumerate(lines):
        if (
            line.strip() != "loop_"
            or start + 1 >= len(lines)
            or not lines[start + 1].strip().startswith("_atom_site.")
        ):
            continue
        index = start + 1
        tags: list[str] = []
        while index < len(lines) and lines[index].strip().startswith("_atom_site."):
            tags.append(lines[index].strip())
            index += 1
        fields = {tag: offset for offset, tag in enumerate(tags)}
        required = (
            "_atom_site.auth_asym_id",
            "_atom_site.auth_seq_id",
            "_atom_site.auth_comp_id",
            "_atom_site.auth_atom_id",
            "_atom_site.type_symbol",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.pdbx_PDB_model_num",
        )
        missing = [tag for tag in required if tag not in fields]
        if missing:
            raise ValueError(f"CIF atom_site loop is missing fields: {missing}")
        selected: dict[tuple[int, str], tuple[str, np.ndarray]] = {}
        while (
            index < len(lines)
            and lines[index].strip()
            and not lines[index].startswith(("#", "loop_", "_"))
        ):
            row = shlex.split(lines[index])
            index += 1
            if (
                len(row) != len(tags)
                or row[fields["_atom_site.pdbx_PDB_model_num"]] != "1"
            ):
                continue
            if row[fields["_atom_site.auth_asym_id"]] != chain:
                continue
            auth_id = int(row[fields["_atom_site.auth_seq_id"]])
            if auth_id not in residues:
                continue
            atom_name = row[fields["_atom_site.auth_atom_id"]]
            element = row[fields["_atom_site.type_symbol"]].upper()
            xyz = np.array(
                [
                    float(row[fields[tag]])
                    for tag in (
                        "_atom_site.Cartn_x",
                        "_atom_site.Cartn_y",
                        "_atom_site.Cartn_z",
                    )
                ]
            )
            if not np.isfinite(xyz).all():
                raise ValueError(
                    f"Non-finite CIF coordinate at {chain}:{auth_id}:{atom_name}"
                )
            key = (auth_id, atom_name)
            if key in selected:
                raise ValueError(
                    f"Alternate/duplicate CIF atom requires resolution: {chain}:{auth_id}:{atom_name}"
                )
            selected[key] = (element, xyz)
        return selected
    raise ValueError(f"No atom_site loop found in {path}")


def _topology_atoms(
    structure: Any,
    residues: tuple[tuple[int, int, str], ...],
) -> dict[tuple[int, str], tuple[str, np.ndarray]]:
    """Return requested atom coordinates keyed by reference auth residue/name."""
    xyz = np.asarray(structure.coordinates, dtype=float)
    if xyz.shape != (len(structure.atoms), 3) or not np.isfinite(xyz).all():
        raise ValueError("Restart coordinates must be finite and match the topology")
    mapped: dict[tuple[int, str], tuple[str, np.ndarray]] = {}
    for amber_res, auth_res, expected_name in residues:
        if amber_res < 1 or amber_res > len(structure.residues):
            raise ValueError(f"Amber residue position out of range: {amber_res}")
        residue = structure.residues[amber_res - 1]
        if residue.name.strip().upper() != expected_name:
            raise ValueError(
                f"Amber residue {amber_res} expected {expected_name}, "
                f"found {residue.name}"
            )
        for atom in residue.atoms:
            name = atom.name.strip()
            if atom.atomic_number == 1:
                continue
            key = (auth_res, name)
            if key in mapped:
                raise ValueError(f"Duplicate Amber atom in mapped residue: {key}")
            element = atom.element
            symbol = {6: "C", 7: "N", 8: "O", 16: "S"}.get(
                element, str(element).upper()
            )
            if symbol not in {"C", "N", "O", "S"}:
                raise ValueError(
                    f"Unsupported mapped atom element {element!r} at {amber_res}:{name}"
                )
            mapped[key] = (symbol, xyz[atom.idx].copy())
    return mapped


def _pdb_line(
    serial: int,
    record: str,
    atom: str,
    resname: str,
    chain: str,
    resid: int,
    xyz: np.ndarray,
    element: str,
) -> str:
    """Format one PDB atom record."""
    atom_field = f" {atom:<3}" if len(atom) < 4 else atom[:4]
    return f"{record:<6}{serial:5d} {atom_field:4s} {resname:>3s} {chain:1s}{resid:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"


def align_catalytic_site(
    reference_cif: Path,
    topology_path: Path,
    starting_restart: Path,
    restart_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Align E-S starting and continuation coordinates to CIF chain A residues.

    The Kabsch fit uses backbone N/CA/C/O atoms of Phe113, Asn186, Asp272,
    and Asp274. Thr233 is excluded from fitting and reported as a held-out
    comparison. Coordinates are descriptive and do not classify TI1 retention.

    Args:
        reference_cif: mmCIF reference structure.
        topology_path: ParmEd-readable Amber topology.
        starting_restart: Starting E-S coordinates corresponding to the topology.
        restart_path: Converged continuation E-S coordinates corresponding to the
            topology.
        output_dir: New or existing directory for overlay and JSON report.

    Returns:
        A JSON-compatible report dictionary.
    """
    reference = _cif_atoms(reference_cif, "A", {auth for _, auth, _ in SITE_RESIDUES})
    fit_keys = [(auth, atom) for _, auth, _ in FIT_RESIDUES for atom in BACKBONE]
    for key in fit_keys:
        if key not in reference:
            raise ValueError(
                f"Missing required reference fit atom at auth residue {key[0]} atom {key[1]}"
            )
    ref_fit = np.array([reference[key][1] for key in fit_keys])
    reference_site = {key: value for key, value in reference.items()}
    aligned_states = {}
    for label, path in (
        ("starting", starting_restart),
        ("continuation_endpoint", restart_path),
    ):
        structure = parmed.load_file(str(topology_path), str(path))
        mobile = _topology_atoms(structure, SITE_RESIDUES)
        for key in fit_keys:
            if key not in mobile or reference[key][0] != mobile[key][0]:
                raise ValueError(f"Missing or element-mismatched mapped atom {key}")
        mob_fit = np.array([mobile[key][1] for key in fit_keys])
        rotation, mob_center, ref_center, fit_rmsd = fit_coordinates(
            mob_fit,
            ref_fit,
        )

        def transform(point: np.ndarray) -> np.ndarray:
            return (point - mob_center) @ rotation + ref_center

        model = {
            key: (element, transform(point)) for key, (element, point) in mobile.items()
        }
        aligned_states[label] = {
            "model": model,
            "fit_rmsd_angstrom": fit_rmsd,
            "rotation": rotation,
            "translation": ref_center - mob_center @ rotation,
            **state_metrics(reference_site, model),
        }
    endpoint = aligned_states["continuation_endpoint"]
    model = endpoint["model"]

    output_dir.mkdir(parents=True, exist_ok=True)
    overlay_path = output_dir / "catalytic_site_overlay.pdb"
    report_path = output_dir / "alignment_report.json"
    if overlay_path.exists() or report_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing alignment outputs in {output_dir}"
        )
    lines = [
        "REMARK 1 Chain A is n12c.cif reference; B is starting E-S; "
        "C is continuation endpoint.\n",
        "REMARK 2 Fit: backbone N, CA, C, O of Phe113, Asn186, Asp272, Asp274.\n",
        "REMARK 3 Thr233 excluded from fit; E-S state is descriptive, "
        "not a TI1 claim.\n",
    ]
    serial = 1
    for chain, atoms in (
        ("A", reference_site),
        ("B", aligned_states["starting"]["model"]),
        ("C", model),
    ):
        for _, auth, resname in SITE_RESIDUES:
            for (resid, atom), (element, coord) in sorted(atoms.items()):
                if resid != auth:
                    continue
                lines.append(
                    _pdb_line(
                        serial, "ATOM", atom, resname, chain, auth, coord, element
                    )
                )
                serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    overlay_path.write_text("".join(lines), encoding="utf-8")
    report: dict[str, Any] = {
        "interpretation": (
            "Descriptive alignment of converged E-S coordinates; the separate "
            "unrestrained release reverted from TI1 and this report does not "
            "claim TI1 retention."
        ),
        "inputs": {
            "reference_cif": str(reference_cif),
            "topology": str(topology_path),
            "starting_restart": str(starting_restart),
            "continuation_endpoint_restart": str(restart_path),
        },
        "mapping": {
            "reference_chain": "A",
            "residues": [
                {
                    "reference_auth": auth,
                    "amber_residue_position": amber,
                    "resname": name,
                }
                for amber, auth, name in SITE_RESIDUES
            ],
        },
        "fit": {
            "residues": [auth for _, auth, _ in FIT_RESIDUES],
            "atoms_per_residue": list(BACKBONE),
            "atom_count": len(fit_keys),
        },
        "states": {
            label: {
                "fit_rmsd_angstrom": item["fit_rmsd_angstrom"],
                "thr233_held_out": item["thr233_held_out"],
                "contacts": item["contacts"],
            }
            for label, item in aligned_states.items()
        },
        "thr233_held_out": endpoint["thr233_held_out"],
        "contacts": endpoint["contacts"],
        "outputs": {"overlay_pdb": str(overlay_path), "report_json": str(report_path)},
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    LOGGER.info("Wrote catalytic alignment overlay and report to %s", output_dir)
    return report
