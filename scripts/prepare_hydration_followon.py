#!/usr/bin/env python3
"""Stage the reviewed-source conditional hydration continuation without MD."""

from __future__ import annotations

import difflib
import json
from pathlib import Path
import shutil

from nyl12_covalent.hydration_followon_gates import digest

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "qmmm/pa6_flipped_model30/13_deacylation_pilot/hydration_reconnaissance"


def prepare() -> None:
    """Create exclusive source copies and an immutable conditional plan."""
    snapshot = BASE / "followon_source_snapshot"
    recovery = BASE / "followon_postprocess_source"
    if snapshot.exists() or recovery.exists():
        raise FileExistsError("Followon preparation already exists")
    original = json.loads((BASE / "extension_20ps_config.json").read_text())
    shutil.copytree(
        BASE / "extension_source_snapshot",
        snapshot,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    shutil.copytree(
        BASE / "postprocess_recovery",
        recovery,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    module = snapshot / "nyl12_covalent/deacylation_hydration_extension.py"
    before = module.read_text()
    # Keep only runtime code: preparation is deliberately owned by this plan.
    start = before.index("def prepare_extension(")
    end = before.index("def _load_identity(")
    after = before[:start] + before[end:]
    after = after.replace("range(1, 5)", "range(5, 10)")
    after = after.replace("extension_20ps_config.json", "followon_runtime_config.json")
    after = after.replace("extension_20ps_status.json", "followon_native_status.json")
    after = after.replace("extension_source_snapshot", "followon_source_snapshot")
    after = after.replace(
        "SEEDS = (2_026_112_402, 2_026_112_403, 2_026_112_404, 2_026_112_405)",
        "SEEDS = (2026112406, 2026112407, 2026112408, 2026112409, 2026112410)",
    )
    after = after.replace(
        "command, cwd=segment, check=False, capture_output=True, text=True",
        "command, cwd=segment, check=False, capture_output=True, text=True,\n"
        "                env=dict(os.environ, **config['postprocess_environment']),",
    )
    after = after.replace(
        'if not acceptance["acyl_retained"]:',
        'if not all(acceptance[key] for key in ("acyl_retained", '
        '"scc_verified", "finite_temperature_passed", "native_process_passed")):',
    )
    module.write_text(after)
    (snapshot / "native_runner_changes.diff").write_text(
        "".join(
            difflib.unified_diff(
                before.splitlines(True),
                after.splitlines(True),
                fromfile="original",
                tofile="followon",
            )
        )
    )
    adapter = recovery / "src/nyl12_covalent/hydration_receipt_adapter.py"
    adapter.write_text(
        adapter.read_text()
        .replace("range(1, 5)", "range(1, 10)")
        .replace("extend_04", "extend_09")
    )
    for name in ("hydration_followon_gates.py", "hydration_followon_queue.py"):
        shutil.copyfile(
            ROOT / "src/nyl12_covalent" / name, snapshot / "nyl12_covalent" / name
        )
    shutil.copyfile(
        ROOT / "scripts/run_hydration_followon_queue.py",
        snapshot / "scripts/run_hydration_followon_queue.py",
    )
    recovery_hashes = {str(p): digest(p) for p in recovery.rglob("*.py")}
    for name, expected in recovery_hashes.items():
        dst = snapshot / "postprocess" / Path(name).name
        dst.write_bytes(Path(name).read_bytes())
        assert digest(dst) == expected
    fixed = {
        str(BASE / "extension_20ps_config.json"): digest(
            BASE / "extension_20ps_config.json"
        ),
        str(BASE / "postprocess_recovery_config.json"): digest(
            BASE / "postprocess_recovery_config.json"
        ),
    }
    fixed.update(
        {
            str(BASE / "extension_source_snapshot" / rel): sha
            for rel, sha in original["source_snapshot_sha256"].items()
        }
    )
    fixed.update(
        json.loads((BASE / "postprocess_recovery_config.json").read_text())[
            "source_sha256"
        ]
    )
    fixed.update({str(p): digest(p) for p in snapshot.rglob("*") if p.is_file()})
    fixed.update(recovery_hashes)
    plan = {
        "base": str(BASE),
        "root": str(ROOT),
        "wait_seconds": 172800,
        "poll_seconds": 120,
        "total_limit_ps": 51,
        "initial_required_ps": 26,
        "segment_names": [f"extend_{i:02d}" for i in range(5, 10)],
        "seeds": {f"extend_{i:02d}": 2026112401 + i for i in range(5, 10)},
        "pinned_sha256": fixed,
        "recovery_source_sha256": recovery_hashes,
        "postprocess_argv": [
            str(ROOT / ".venv/bin/python"),
            str(recovery / "scripts/refresh_deacylation_hydration.py"),
            "--base",
            str(BASE),
            "--project-root",
            str(recovery),
        ],
        "postprocess_environment": {
            "PYTHONPATH": f"{recovery}/src:{recovery}",
            "LD_LIBRARY_PATH": "/data/jpc/code/ambertools26/venv/lib",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        },
        "scope": "Conditional five 5 ps blocks; stop at any plausible water or failed gate. All waters remain MM; no protonation/QM/movable/Hamiltonian change. No reheating. Not sampling convergence.",
    }
    (BASE / "followon_25ps_plan.json").write_text(json.dumps(plan, indent=2) + "\n")


if __name__ == "__main__":
    prepare()
