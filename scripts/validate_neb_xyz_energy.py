"""Validate a reconstructed NEB image against its recorded native energy."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import sys


def validate(campaign: Path) -> None:
    """Run one initial energy/gradient and compare to the matched NEB report."""
    sys.path.insert(0, str(campaign / "source_xyz_export"))
    from scripts.run_acylenzyme_neb import input_text, run_engine

    cfg = json.loads((campaign / "config.json").read_text())
    expected = json.loads((campaign / "xyz_validation/expected.json").read_text())
    folder = campaign / "xyz_validation/single_point"
    text = input_text(Path(cfg["base_input"]).read_text(), Path(cfg["reactant"]))
    text = re.sub(r"maxcyc\s*=\s*1000", "maxcyc=1", text)
    text = re.sub(r"maxcycle\s*=\s*1000", "maxcycle=1", text)
    start = campaign / f"xyz_validation/band/optimized.{expected['image']:04d}.ncrst"
    run_engine(cfg, folder, start, text)
    matches = re.findall(
        r"Energy calculation finished, energy:\s*([-+\d.Ee]+)",
        (folder / "minimized.mdout").read_text(),
    )
    if not matches:
        raise ValueError("No initial energy found")
    value = float(matches[0])
    difference = (value - expected["energy_hartree"]) * 627.5094740631
    result = {
        **expected,
        "recomputed_initial_energy_hartree": value,
        "difference_kcal_mol": difference,
        "tolerance_kcal_mol": 0.01,
        "passed": abs(difference) < 0.01,
        "interpretation": "Coordinate/energy consistency only, not Hamiltonian validation",
    }
    (campaign / "xyz_validation/result.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    if not result["passed"]:
        raise ValueError(f"Energy mismatch: {difference} kcal/mol")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    validate(parser.parse_args().campaign.resolve())
