"""Verify saved Gaussian and xTB embeddings; no new QM calls."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from reference_embedding import parse_external_input  # noqa: E402


def main() -> None:
    """Save explicit PC correspondence and artifact hashes after successful gates."""
    checks = []
    inputs = {}
    for state in ("01_es", "02_activated", "03_attack"):
        folder = HERE.parent / "evaluations" / state
        source = (
            HERE.parent / "es_seed/nonperiodic_full_pc/old.gau_job.inp"
            if state == "01_es"
            else folder / "extern_export/old.gau_job.inp"
        )
        direct = folder / "gaussian_direct/direct_force.gjf"
        pc = folder / "gfn2_xtb/point_charges.pc"
        original = np.asarray(parse_external_input(source)["point_charges"])
        gaussian = np.asarray(parse_external_input(direct)["point_charges"])
        xtb = np.loadtxt(pc, skiprows=1)
        # xTB PC columns q,x_bohr,y_bohr,z_bohr,hardness.
        reordered = np.column_stack((xtb[:, 1:4] / 1.8897261254578281, xtb[:, 0]))
        assert original.shape == gaussian.shape == reordered.shape == (104824, 4)
        assert np.array_equal(original, gaussian)
        assert np.allclose(original, reordered, rtol=0, atol=1e-12)
        checks.append(
            {
                "state": state,
                "PC_count": len(original),
                "DFT_PC_exact_match": True,
                "xTB_PC_max_abs_delta_A_or_e": float(abs(original - reordered).max()),
            }
        )
        for p in (source, direct, pc):
            inputs[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    status = json.loads((HERE / "job/status.json").read_text())
    assert status["exit_status"] == 0
    acceptance = {
        "accepted_utc": datetime.now(timezone.utc).isoformat(),
        "process_complete": True,
        "scientific_mapping_and_embedding_pass": True,
        "embedding_verification": checks,
        "embedding_source_sha256": inputs,
        "coverage": "Focused parser/projection tests passed; line coverage unmeasured.",
        "no_new_QM": True,
        "interpretation": "Local derivatives only; no minima or barriers established.",
        "artifacts_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in HERE.iterdir()
            if p.is_file()
        },
    }
    (HERE / "acceptance.json").write_text(json.dumps(acceptance, indent=2) + "\n")


if __name__ == "__main__":
    main()
