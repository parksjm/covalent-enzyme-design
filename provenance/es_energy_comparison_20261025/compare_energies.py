"""Compare printed full-system energies without counting positional bias."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "output_nyl12_pa6_dftb_pilot_20261024"


def last_match(text: str, pattern: str) -> tuple[float, str]:
    """Return a numeric value and its exact last matched source text."""
    values = list(re.finditer(pattern, text, flags=re.M))
    if not values:
        raise ValueError(f"Missing native energy record: {pattern}")
    match = values[-1]
    return float(match[1].replace("D", "E").replace("d", "e")), match[0]


def main() -> None:
    """Retain bias-free energy evidence, separating a live evaluation from a minimum."""
    old = BASE / "analytic_relaxation_continuation_attempt1"
    current = BASE / "asp272_relay_attempt2/phase1_unrestrained"
    old_bytes = (old / "minimized.mdout").read_bytes()
    new_bytes = (current / "minimized.mdout").read_bytes()
    old_text, new_text = old_bytes.decode(), new_bytes.decode()
    original, original_line = last_match(old_text, r"^\s*EAMBER\s*=\s*([-+\d.EeDd]+)")
    bias, bias_line = last_match(old_text, r"RESTRAINT\s*=\s*([-+\d.EeDd]+)")
    latest, latest_line = last_match(new_text, r"^\s*Etot\s*=\s*([-+\d.EeDd]+)\s+EKtot")
    new_bias, new_bias_line = last_match(new_text, r"RESTRAINT\s*=\s*([-+\d.EeDd]+)")
    old_qm, old_qm_line = last_match(old_text, r"DFTBESCF\s*=\s*([-+\d.EeDd]+)")
    new_qm, new_qm_line = last_match(new_text, r"DFTBESCF\s*=\s*([-+\d.EeDd]+)")
    for folder, name, expected in (
        (old, "continuation.mdin", 1),
        (current, "minimize.mdin", 0),
    ):
        values = re.findall(r"\bntr\s*=\s*(\d+)", (folder / name).read_text())
        if values != [str(expected)]:
            raise ValueError("Restraint controls differ from the audited comparison")
    if bias <= 0 or new_bias != 0:
        raise ValueError("Unexpected printed restraint-energy evidence")
    result = {
        "observed_utc": datetime.now(timezone.utc).isoformat(),
        "original_ES_source": str(old),
        "reversion_source": str(current),
        "original_ES_EAMBER_without_position_restraints_kcal_mol": original,
        "original_ES_position_restraint_kcal_mol": bias,
        "original_ES_total_reconstructed_kcal_mol": original + bias,
        "latest_reversion_Etot_kcal_mol": latest,
        "relative_without_position_restraints_kcal_mol": latest - original,
        "printed_DFTBESCF_difference_kcal_mol": new_qm - old_qm,
        "reversion_process": json.loads((current / "process.json").read_text()),
        "scope": "Live last printed evaluation, not necessarily the latest restart or a converged minimum. Same topology/QM/embedding and fixed cell; geometry optimization histories differ.",
        "evidence_lines": [
            original_line,
            bias_line,
            latest_line,
            new_bias_line,
            old_qm_line,
            new_qm_line,
        ],
        "log_snapshot_sha256": {
            str(old / "minimized.mdout"): hashlib.sha256(old_bytes).hexdigest(),
            str(current / "minimized.mdout"): hashlib.sha256(new_bytes).hexdigest(),
        },
        "original_native_energy_excludes_bias": "Native EAMBER is reported separately from the positive RESTRAINT contribution.",
    }
    output = Path(__file__).with_name("comparison.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    logging.info(
        "Original bias-free %.4f; current %.4f; difference %.4f kcal/mol",
        original,
        latest,
        latest - original,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
