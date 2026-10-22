"""Run a real headless PyMOL viewer check without adding or changing atoms."""

from __future__ import annotations

import argparse
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
import logging
import math
from pathlib import Path

import pymol2


def main() -> None:
    """Inspect viewer parsing, all water states, bonds and hydrogen visibility."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--viewer", type=Path, required=True)
    parser.add_argument("--object", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-pass", action="store_true")
    parser.add_argument("--repeat-viewer", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    output = args.output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    viewer = args.viewer.absolute()
    console = io.StringIO()
    with redirect_stdout(console), redirect_stderr(console), pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.cd(str(args.data_dir.absolute()))
        cmd.do("@" + str(viewer), echo=0)
        first_states = cmd.count_states(args.object)
        if args.repeat_viewer:
            cmd.do("@" + str(viewer), echo=0)
        states = cmd.count_states(args.object)
        rows = []
        for state in range(1, states + 1):
            water = cmd.get_model(f"{args.object} and chain J", state=state)
            atoms = {atom.name: atom for atom in water.atom}
            if set(atoms) != {"O", "H1", "H2"}:
                raise ValueError(
                    f"Wrong water inventory at state {state}: {set(atoms)}"
                )
            if [atoms[name].symbol for name in ("O", "H1", "H2")] != ["O", "H", "H"]:
                raise ValueError("Wrong water elements")
            distances = [
                math.dist(atoms["O"].coord, atoms[name].coord) for name in ("H1", "H2")
            ]
            rows.append({"state": state, "OH_A": distances, "bonds": len(water.bond)})
        visible_h = cmd.count_atoms(f"{args.object} and chain J and elem H and visible")
        version = cmd.get_version()[0]
    transcript = console.getvalue()
    (output / "pymol_console.log").write_text(transcript)
    errors = [
        line
        for line in transcript.splitlines()
        if any(word in line for word in ("SyntaxError", "Traceback", "Error:"))
    ]
    report = {
        "pymol_version": version,
        "viewer": str(viewer),
        "viewer_sha256": hashlib.sha256(viewer.read_bytes()).hexdigest(),
        "states": states,
        "water_atoms_per_state": 3 if rows else 0,
        "visible_water_hydrogens": visible_h,
        "water_states": rows,
        "parser_errors": errors,
        "accepted": bool(rows)
        and not errors
        and visible_h == 2
        and all(
            row["bonds"] == 2 and all(0.8 <= d <= 1.2 for d in row["OH_A"])
            for row in rows
        ),
        "coordinates_changed_by_audit": False,
    }
    report["reload_did_not_append_states"] = states == first_states
    report["accepted"] = report["accepted"] and states == first_states
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    logging.info(
        "PyMOL %s: %d states, %d visible water H; accepted=%s",
        version,
        states,
        visible_h,
        report["accepted"],
    )
    if args.require_pass and not report["accepted"]:
        raise ValueError(
            "Viewer/water audit failed; inspect retained report and console"
        )


if __name__ == "__main__":
    main()
