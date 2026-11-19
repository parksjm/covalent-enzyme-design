"""Prepare a one-cycle, zero-smearing native SANDER gradient probe."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)


def _record(path: Path) -> dict[str, Any]:
    """Return reproducibility metadata for one input file."""
    payload = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "bytes": len(payload),
    }


def prepare_input(source: str) -> str:
    """Preserve the input Hamiltonian while disabling dimer mode and smearing.

    Args:
        source: Existing converged DL-FIND dimer input.

    Returns:
        An input with only DFTB temperature, dimer mode, and cycle budget changed.
    """
    replacements = (
        (r"\bdftb_telec\s*=\s*300(?:\.0*)?\b", "dftb_telec=0.0", 1),
        (r"\bdimer\s*=\s*\.true\.", "dimer=.false.", 1),
        (r"\bmaxcycle\s*=\s*200\b", "maxcycle=1", 1),
    )
    result = source
    for pattern, replacement, expected in replacements:
        result, count = re.subn(pattern, replacement, result, flags=re.IGNORECASE)
        if count != expected:
            raise ValueError(
                f"Expected {expected} match for {pattern!r}; found {count}"
            )
    for key in ("qmmask", "qmcharge", "qm_theory", "qmcut", "scfconv", "active"):
        if not re.search(rf"\b{key}\s*=", source, re.IGNORECASE):
            raise ValueError(f"Source input lacks required setting {key}")
        before = re.search(rf"\b{key}\s*=\s*([^,\n]+)", source, re.IGNORECASE)
        after = re.search(rf"\b{key}\s*=\s*([^,\n]+)", result, re.IGNORECASE)
        if before is None or after is None or before.group(1) != after.group(1):
            raise ValueError(f"Probe input unexpectedly changed setting {key}")
    return result


def prepare_case(
    source_mdin: Path,
    topology: Path,
    restart: Path,
    sander: Path,
    output_dir: Path,
) -> Path:
    """Create a unique isolated input and source-hash manifest.

    Args:
        source_mdin: Exact converged 300 K dimer input.
        topology: Exact topology used by the converged calculation.
        restart: Converged final candidate coordinates.
        sander: Installed native SANDER executable.
        output_dir: New calculation directory for this bounded attempt.

    Returns:
        Path to the isolated one-cycle input file.
    """
    output_dir.mkdir(parents=True, exist_ok=False)
    mdin = output_dir / "dimer.mdin"
    mdin.write_text(prepare_input(source_mdin.read_text()), encoding="utf-8")
    manifest = {
        "purpose": "one-cycle native SANDER LBFGS gradient feasibility at 0 K",
        "changes": {
            "dftb_telec_kelvin": 0.0,
            "dimer": False,
            "dlfind_maxcycle": 1,
        },
        "preserved": [
            "qmmask",
            "qmcharge",
            "qm_theory",
            "QM/MM embedding",
            "12 A active mask",
        ],
        "inputs": [_record(path) for path in (source_mdin, topology, restart, sander)],
        "prepared_mdin": _record(mdin),
        "preparer": _record(Path(__file__)),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    LOGGER.info("Prepared one-cycle 0 K gradient probe in %s", output_dir)
    return mdin


def main() -> None:
    """Prepare the bounded native SANDER input bundle."""
    base = Path(
        "qmmm/pa6_flipped_model30/05_neutral_reaction_path/07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-mdin", type=Path, default=base / "calculation/dimer.mdin"
    )
    parser.add_argument("--topology", type=Path, default=base / "system.prmtop")
    parser.add_argument(
        "--restart", type=Path, default=base / "calculation/optimized.ncrst"
    )
    parser.add_argument(
        "--sander",
        type=Path,
        default=Path("/data/jpc/code/ambertools26/venv/bin/sander"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    prepare_case(
        args.source_mdin, args.topology, args.restart, args.sander, args.output_dir
    )


if __name__ == "__main__":
    main()
