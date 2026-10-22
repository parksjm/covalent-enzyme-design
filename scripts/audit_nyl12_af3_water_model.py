"""Audit and export the completed minimal AF3 PA6 L4 TI1/water model."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.af3_postprocess import read_atom_sites  # noqa: E402
from nyl12_covalent.af3_water_audit import audit_sites  # noqa: E402


def main() -> None:
    """Write an immutable model audit and inspectable unchanged CIF copy."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.output.exists():
        raise FileExistsError(args.output)
    payload = json.loads(args.input.read_text())
    expected = [[["E", 1, "OG1"], ["I", 1, "C54"]]]
    if payload["bondedAtomPairs"] != expected:
        raise ValueError("Input covalent bond differs from the audited PA6 TI1 bond")
    report = audit_sites(read_atom_sites(args.model))
    report["source_inputs"] = {
        label: {
            "path": str(path.absolute()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for label, path in (
            ("model", args.model),
            ("af3_input", args.input),
            ("substrate_ccd", Path(payload["userCCDPath"])),
        )
    }
    report["declared_covalent_bonds"] = payload["bondedAtomPairs"]
    report["input_contact_constraints"] = payload["contactConstraints"]
    report["source_code_sha256"] = {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [
            Path(__file__),
            *sorted((ROOT / "src/nyl12_covalent").glob("*.py")),
        ]
    }
    args.output.mkdir(parents=True)
    (args.output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    shutil.copy2(args.model, args.output / "pa6_l4_alpha_beta_water_test.cif")
    (args.output / "view_af3_water_test.pml").write_text(
        "load pa6_l4_alpha_beta_water_test.cif, af3_water_test\n"
        "hide everything\nshow cartoon, chain A+E\n"
        "show sticks, chain I or (chain E and resi 1+2)\n"
        "show sticks, resi 113+186 and chain A\n"
        "show spheres, chain J\nset sphere_scale, 0.25, chain J\n"
        "set valence, 0\nutil.cnc\n"
        "distance water_to_N, chain J and name O, chain E and resi 1 and name N\n"
        "distance water_to_OG1, chain J and name O, chain E and resi 1 and name OG1\n"
        "zoom chain I+J or (chain E and resi 1), 8\n"
    )
    logging.info(
        "Heavy-atom screen passed: %s; QM/MM acceptance: false",
        report["heavy_atom_screen_passed"],
    )
    logging.info("Bridge distances: %s", report["bridge_distances"])
    logging.info("Active-site metrics: %s", report["active_site"]["metrics"])


if __name__ == "__main__":
    main()
