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
    stages = [str(row["stage"]) for row in rows]
    if not stages or any(stage not in definitions for stage in stages):
        raise ValueError("Missing measurement definition for reaction stage")
    parts = []
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
    viewer.write_text(
        viewer.read_text()
        + "\nsource [file join $::nyl12_dir nyl12_reactive_distances.tcl]\n"
        + "::nyl12_reactive::install $::nyl12_mol "
        + "{"
        + " ".join(stages)
        + "} {"
        + " ".join(parts)
        + "}\n"
    )
