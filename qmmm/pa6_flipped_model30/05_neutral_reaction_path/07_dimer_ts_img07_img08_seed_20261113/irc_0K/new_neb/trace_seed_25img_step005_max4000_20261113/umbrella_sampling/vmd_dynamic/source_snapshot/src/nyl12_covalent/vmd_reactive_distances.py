"""Configure native VMD bond-distance labels for each elementary stage."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any


def append_reactive_labels(
    viewer: Path, rows: list[dict[str, Any]], config: dict[str, Any]
) -> None:
    """Add explicit atom selections and frame-stage measurement switching.

    Args:
        viewer: Newly exported VMD loader to annotate.
        rows: Ordered endpoint records containing elementary-stage identities.
        config: Declared stripped atom indices/names and stage bond lists.
    """
    atoms = config["measurement_atoms"]
    definitions = config["stage_measurements"]
    dynamic = config.get("stage_dynamic_bonds", {})
    persistent = config.get("persistent_measurements", {})
    persistent_angles = config.get("persistent_angles", {})
    cutoffs = config.get("dynamic_bond_cutoffs", {})
    proton_atom = config.get("proton_atom")
    stages = [str(row["stage"]) for row in rows]
    if not stages or any(stage not in definitions for stage in stages):
        raise ValueError("Missing measurement definition for reaction stage")
    parts = []
    dynamic_parts = []
    for stage, pairs in definitions.items():
        if not re.fullmatch(r"[A-Za-z0-9_]+", stage) or not pairs:
            raise ValueError("Invalid or empty reaction-stage measurements")
        lines = []
        for pair in pairs:
            a, b = (atoms[key] for key in pair)
            if a["index"] == b["index"]:
                raise ValueError("Reactive measurement needs distinct atoms")
            for atom in (a, b):
                if (
                    not isinstance(atom["index"], int)
                    or atom["index"] < 0
                    or not re.fullmatch(r"[A-Za-z0-9_]+", atom["name"])
                ):
                    raise ValueError("Invalid reactive atom identity")
            lines.append(f"{{{a['index']} {a['name']} {b['index']} {b['name']}}}")
        parts.append(f"{stage} {{{' '.join(lines)}}}")
        dynamic_lines = []
        for pair in dynamic.get(stage, []):
            a, b = (atoms[key] for key in pair)
            _validate_atom_pair(a, b)
            cutoff_key = "-".join(pair)
            cutoff = cutoffs.get(cutoff_key)
            if not isinstance(cutoff, (int, float)) or cutoff <= 0:
                raise ValueError(f"Missing positive dynamic-bond cutoff: {cutoff_key}")
            dynamic_lines.append(
                f"{{{a['index']} {a['name']} {b['index']} {b['name']} {cutoff}}}"
            )
        dynamic_parts.append(f"{stage} {{{' '.join(dynamic_lines)}}}")
    persistent_lines = []
    for label, pair in persistent.items():
        if not re.fullmatch(r"[A-Za-z0-9_]+", label) or len(pair) != 2:
            raise ValueError("Invalid persistent measurement declaration")
        a, b = (atoms[key] for key in pair)
        _validate_atom_pair(a, b)
        persistent_lines.append(
            f"{{{label} {a['index']} {a['name']} {b['index']} {b['name']}}}"
        )
    angle_lines = []
    for label, atom_keys in persistent_angles.items():
        if not re.fullmatch(r"[A-Za-z0-9_]+", label) or len(atom_keys) != 3:
            raise ValueError("Invalid persistent angle declaration")
        angle_atoms = [atoms[key] for key in atom_keys]
        for atom in angle_atoms:
            if (
                not isinstance(atom["index"], int)
                or atom["index"] < 0
                or not re.fullmatch(r"[A-Za-z0-9_]+", atom["name"])
            ):
                raise ValueError("Invalid persistent angle atom identity")
        if len({atom["index"] for atom in angle_atoms}) != 3:
            raise ValueError("Persistent angle needs three distinct atoms")
        angle_lines.append(
            "{"
            + label
            + " "
            + " ".join(f"{atom['index']} {atom['name']}" for atom in angle_atoms)
            + "}"
        )
    if proton_atom is not None:
        proton = atoms[proton_atom]
        if (
            not isinstance(proton["index"], int)
            or proton["index"] < 0
            or not re.fullmatch(r"[A-Za-z0-9_]+", proton["name"])
        ):
            raise ValueError("Invalid proton atom identity")
    viewer.write_text(
        viewer.read_text()
        + "\nsource [file join $::nyl12_dir nyl12_reactive_distances.tcl]\n"
        + "::nyl12_reactive::install $::nyl12_mol "
        + "{"
        + " ".join(stages)
        + "} {"
        + " ".join(parts)
        + "} {"
        + " ".join(dynamic_parts)
        + "} {"
        + " ".join(persistent_lines)
        + "}"
        + (
            f" {proton['index']} {proton['name']}"
            if proton_atom is not None
            else " {} {}"
        )
        + " {"
        + " ".join(angle_lines)
        + "}"
        + "\n"
    )


def _validate_atom_pair(a: dict[str, Any], b: dict[str, Any]) -> None:
    """Validate two distinct zero-based atom identities."""
    for atom in (a, b):
        if (
            not isinstance(atom["index"], int)
            or atom["index"] < 0
            or not re.fullmatch(r"[A-Za-z0-9_]+", atom["name"])
        ):
            raise ValueError("Invalid reactive atom identity")
    if a["index"] == b["index"]:
        raise ValueError("Reactive measurement needs distinct atoms")
