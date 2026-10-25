"""Audit the energy rise without changing the QM/MM calculation."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import re
import sys

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.path_energies import relay_endpoint_energies  # noqa: E402

RUN = ROOT / "output_nyl12_pa6_dftb_pilot_20261024/asp272_relay_attempt2"


def main() -> None:
    """Separate printed QM, driving and MM restraint changes through window13."""
    endpoints = [{"window": f"stage1_{i:02d}", "accepted": True} for i in range(14)]
    rows = relay_endpoint_energies(RUN, endpoints)
    for row in rows:
        text = (RUN / row["window"] / "minimized.mdout").read_text()
        for label in ("RESTRAINT", "DFTBESCF"):
            row[label] = float(
                re.findall(r"\b" + label + r"\s*=\s*([-+\d.Ee]+)", text)[-1]
            )
        row["MM_positional_restraint_kcal_mol"] = (
            row["RESTRAINT"] - row["driving_bias_kcal_mol"]
        )
        row["other_printed_terms_kcal_mol"] = (
            row["total_kcal_mol"] - row["RESTRAINT"] - row["DFTBESCF"]
        )
    keys = [
        "total_kcal_mol",
        "driving_bias_kcal_mol",
        "MM_positional_restraint_kcal_mol",
        "DFTBESCF",
        "other_printed_terms_kcal_mol",
    ]
    differences = {key: rows[-1][key] - rows[0][key] for key in keys}
    result = {
        "reference": "stage1_00",
        "endpoint": "stage1_13",
        "differences_kcal_mol": differences,
        "rows": rows,
        "caveat": "Printed energy partition, not causal or residue-specific decomposition; same original MM positional reference. No free-energy or pKa inference.",
    }
    (RUN / "relay_energy_component_audit.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    for key, value in differences.items():
        logging.info("%s: %.6f kcal/mol", key, value)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
