"""Reversible inspection exports around the established four-L4 converter."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import shutil
from typing import Any

from .af3_postprocess import read_atom_sites
from .catalytic_alignment_geometry import wrapped_delta
from .design_pa6_mapping import (
    _pdb_atom_line,
    _protein_chain_residue,
    convert_four_ligands,
)
from .relay_any_amide import evaluate_model, load_graph
from .inspection_metric_audit import serialization_metric_audit


def sha256(path: Path) -> str:
    """Return a file content digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pdb_atoms(path: Path) -> dict[tuple[str, int, str], tuple[str, tuple]]:
    """Read unique heavy-atom identities from a fixed-column inspection PDB."""
    atoms = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[76:78].strip() in ("H", "D"):
            continue
        key = (line[21], int(line[22:26]), line[12:16].strip())
        if key in atoms:
            raise ValueError(f"Duplicate PDB identity: {key}")
        xyz = tuple(float(line[i : i + 8]) for i in (30, 38, 46))
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError("Nonfinite coordinates")
        atoms[key] = (line[76:78].strip(), xyz)
    if not atoms:
        raise ValueError("Empty PDB")
    return atoms


def restore_modular(source: Path, conversion: dict, output: Path) -> None:
    """Restore canonical mature chain-local heavy identities from modular PDB.

    Args:
        source: Sequential Amber-layout PDB, including post-engine exported PDB.
        conversion: Existing convert_four_ligands conversion receipt.
        output: New canonical mature-local heavy PDB destination.

    Raises:
        ValueError: Heavy identities, inventories or residue layout differ.
    """
    if output.exists():
        raise FileExistsError(output)
    inverse = {
        (r["target_residue"], r["target_atom"]): (
            r["source_chain"],
            1,
            r["source_atom"],
        )
        for r in conversion["atom_mapping"]
    }
    if len(inverse) != 132:
        raise ValueError("Expected132unique ligand mappings")
    residue_names = {
        int(line[22:26]): line[17:20].strip()
        for line in source.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    }
    records = []
    observed = set()
    for (_, number, name), (element, xyz) in pdb_atoms(source).items():
        if number <= 1312:
            chain, local = _protein_chain_residue(number)
            key = (chain, local, name)
        else:
            if (number, name) not in inverse:
                raise ValueError(f"Unmapped modular atom: {number}:{name}")
            key = inverse[number, name]
        if key in observed:
            raise ValueError(f"Duplicate restored identity: {key}")
        observed.add(key)
        chain, local, name = key
        if chain in "IJKL":
            residue = "LIG"
        else:
            # Preserve actual residue identities, including engine protonation names.
            residue = residue_names[number]
            if residue in ("HID", "HIE", "HIP"):
                residue = "HIS"
        records.append(
            _pdb_atom_line(
                "ATOM" if chain in "ABCDEFGH" else "HETATM",
                len(records) + 1,
                name,
                residue,
                chain,
                local,
                *xyz,
                element,
            )
        )
    if len([k for k in observed if k[0] in "IJKL"]) != 132:
        raise ValueError("Incomplete ligand restoration")
    output.write_text("\n".join(records + ["END"]) + "\n")


def _metric_error(metric: str, before: float, after: float) -> float:
    """Compare signed phi cyclically; retain linear errors for other metrics."""
    return (
        abs(wrapped_delta(before, after))
        if metric == "herschlag_phi"
        else abs(before - after)
    )


def export_control(
    source: Path, maps: dict, mol2: Path, mapping: Path, release: Path, output: Path
) -> dict[str, Any]:
    """Export native predictions through existing mature-local/modular conventions.

    Args:
        source: Original predictor CIF, never modified.
        maps: Native ligand name to source-MOL2 name per ligand chain.
        mol2: Validated whole-L4 graph.
        mapping: Existing whole-L4 to modular-context mapping.
        release: Existing released P6N/P6I/P6C templates.
        output: New private staging destination.

    Returns:
        Round-trip acceptance and reversible identity receipt.
    """
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(source, output / "raw_prediction.cif")
    sites = read_atom_sites(source)
    local_lines, biological_lines, identities = [], [], []
    for s in sites:
        ligand = s.chain in "IJKL"
        name = maps[s.chain][s.name] if ligand else s.name
        local = 1 if ligand else s.sequence
        precursor = local + 232 if s.chain in "EFGH" else local
        residue = "LIG" if ligand else s.residue
        common = ("HETATM" if ligand else "ATOM", s.serial, name, residue, s.chain)
        local_lines.append(_pdb_atom_line(*common, local, *s.xyz, s.element))
        biological_lines.append(_pdb_atom_line(*common, precursor, *s.xyz, s.element))
        identities.append(
            {
                "native": [s.chain, s.sequence, s.name],
                "mature_local": [s.chain, local, name],
                "biological_precursor": [s.chain, precursor, name],
                "element": s.element,
                "native_xyz_A": s.xyz,
            }
        )
    localpdb = output / "canonical_mature_local.pdb"
    localpdb.write_text("\n".join(local_lines + ["END"]) + "\n")
    (output / "canonical_precursor.pdb").write_text(
        "\n".join(biological_lines + ["END"]) + "\n"
    )
    # Auth identifiers carry biological precursor numbering; label_seq_id stays polymer-local.
    lines = source.read_text().splitlines()
    columns = []
    active = False
    derived = []
    for line in lines:
        if line.startswith("_atom_site."):
            columns.append(line.strip().split(".")[1])
            active = True
            continue
        if active and line.strip() == "#":
            break
        if active and line.strip():
            values = line.split()
            d = dict(zip(columns, values, strict=True))
            chain = d["auth_asym_id"]
            name = d["label_atom_id"]
            if chain in "IJKL":
                d["label_atom_id"] = maps[chain][name]
                if "auth_atom_id" in d:
                    d["auth_atom_id"] = maps[chain][name]
            if chain in "EFGH":
                d["auth_seq_id"] = str(int(d["auth_seq_id"]) + 232)
            derived.append(" ".join(d[c] for c in columns))
    (output / "canonical_precursor.cif").write_text(
        "data_DERIVED_CANONICAL_PRECURSOR\n#\nloop_\n"
        + "\n".join("_atom_site." + c for c in columns)
        + "\n"
        + "\n".join(derived)
        + "\n#\n"
    )
    conversion = convert_four_ligands(
        localpdb, mol2, mapping, release, output / "amber_sequential_input"
    )
    ligand_inverse = {
        (r["source_chain"], r["source_atom"]): (r["target_residue"], r["target_atom"])
        for r in conversion["atom_mapping"]
    }
    for row in identities:
        chain, number, name = row["mature_local"]
        if chain in "IJKL":
            number, name = ligand_inverse[chain, name]
        else:
            number += (
                "ABCD".index(chain) * 232
                if chain in "ABCD"
                else 928 + "EFGH".index(chain) * 96
            )
        row["amber_sequential"] = [number, name]
    restored = output / "roundtrip_mature_local.pdb"
    restore_modular(
        output / "amber_sequential_input/complex_heavy_modular.pdb",
        conversion,
        restored,
    )
    expected = {
        (r["mature_local"][0], r["mature_local"][1], r["mature_local"][2]): (
            r["element"],
            r["native_xyz_A"],
        )
        for r in identities
    }
    actual = pdb_atoms(restored)
    if set(actual) != set(expected):
        raise ValueError("Roundtrip atom identity mismatch")
    delta = max(
        abs(a - b)
        for k in actual
        for a, b in zip(actual[k][1], expected[k][1], strict=True)
    )
    if delta > 0.000501:
        raise ValueError("Roundtrip exceeds PDB rounding tolerance")
    graph = load_graph(mol2)
    native = evaluate_model(
        expected, {c: {v: v for v in maps[c].values()} for c in "IJKL"}, graph
    )
    back = evaluate_model(
        actual, {c: {v: v for v in maps[c].values()} for c in "IJKL"}, graph
    )
    serialized = evaluate_model(
        pdb_atoms(localpdb),
        {c: {v: v for v in maps[c].values()} for c in "IJKL"},
        graph,
    )
    metric_audit = serialization_metric_audit(native, serialized, back)
    cif_atoms = {
        (
            s.chain,
            1 if s.chain in "IJKL" else s.sequence - (232 if s.chain in "EFGH" else 0),
            s.name,
        ): (s.element, s.xyz)
        for s in read_atom_sites(output / "canonical_precursor.cif")
    }
    if cif_atoms != expected:
        raise ValueError("Full-precision canonical CIF identity/coordinate mismatch")
    cif_geometry = evaluate_model(
        cif_atoms, {c: {v: v for v in maps[c].values()} for c in "IJKL"}, graph
    )
    cif_errors = {
        k: max(
            _metric_error(k, a[k], b[k])
            for a, b in zip(
                native["all_amide_site_rows"],
                cif_geometry["all_amide_site_rows"],
                strict=True,
            )
        )
        for k in ("d1", "d2", "d3", "a_bd", "herschlag_phi")
    }
    if any(v != 0.0 for v in cif_errors.values()):
        raise ValueError("Full-precision canonical CIF metric mismatch")
    nearest = [
        min(
            (r for r in native["all_amide_site_rows"] if r["site"] == site),
            key=lambda r: r["d1"],
        )
        for site in ("AE", "BF", "CG", "DH")
    ]
    receipt = {
        "source_sha256": sha256(source),
        "atom_mapping": identities,
        "roundtrip_max_coordinate_component_error_A": delta,
        **metric_audit,
        "full_precision_canonical_CIF_metric_max_errors": cif_errors,
        "raw_atom_count": len(sites),
        "roundtrip_atom_count": len(actual),
        "topology_acceptance": "Not run; heavy sequential input only, no hydrogens/protonation assignment or topology readiness claim",
        "nearest_d1_per_site": nearest,
        "assigned": native["assigned"],
        "all48": native["all_amide_site_rows"],
    }
    (output / "mapping_and_acceptance.json").write_text(
        json.dumps(receipt, indent=2) + "\n"
    )
    return receipt


def publish_subset(stage: Path, destination: Path) -> None:
    """Atomically publish a complete new curated tree with a relative manifest.

    Args:
        stage: Private completed staging tree on the destination filesystem.
        destination: New immutable publication directory.
    """
    if destination.exists():
        raise FileExistsError(destination)
    files = sorted(p for p in stage.rglob("*") if p.is_file())
    if any(p.is_symlink() for p in stage.rglob("*")):
        raise ValueError("Publication contains symlinks")
    (stage / "MANIFEST.sha256").write_text(
        "".join(f"{sha256(f)}  {f.relative_to(stage)}\n" for f in files)
    )
    os.rename(stage, destination)
