"""Create source-pinned continuation snapshots with explicit patch anchors."""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    """Return a file's SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(source: str, old: str, new: str, label: str) -> str:
    """Apply one explicit source patch, requiring a unique anchor."""
    count = source.count(old)
    if count != 1:
        raise ValueError(f"Patch anchor {label!r} matched {count} times, expected 1")
    return source.replace(old, new, 1)


def patch_native_sources(source: Path, destination: Path) -> dict[str, str]:
    """Clone and specialize the validated queue for extend_10 through extend_19."""
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("__pycache__"))
    package = destination / "nyl12_covalent"
    shutil.copy2(
        ROOT / "src/nyl12_covalent/hydration_continuation_101ps.py",
        package / "hydration_continuation_101ps.py",
    )
    extension = package / "deacylation_hydration_extension.py"
    text = extension.read_text()
    text = replace_once(
        text,
        'SEGMENT_NAMES = tuple(f"extend_{index:02d}" for index in range(5, 10))',
        'SEGMENT_NAMES = tuple(f"extend_{index:02d}" for index in range(10, 20))',
        "10 continuation segment range",
    )
    text = replace_once(
        text,
        "SEEDS = (2026112406, 2026112407, 2026112408, 2026112409, 2026112410)",
        "SEEDS = (2026112610, 2026112611, 2026112612, 2026112613, 2026112614, 2026112615, 2026112616, 2026112617, 2026112618, 2026112619)",
        "unique bounded seeds",
    )
    text = text.replace(
        "followon_runtime_config.json", "continuation_101ps_config.json"
    )
    text = text.replace("followon_source_snapshot", "continuation_101ps_source")
    extension.write_text(text)

    gates = package / "hydration_followon_gates.py"
    text = gates.read_text()
    text = replace_once(
        text,
        '    config_name = (\n        "extension_20ps_config.json" if index <= 4 else "followon_runtime_config.json"\n    )',
        '    config_name = (\n        "extension_20ps_config.json" if index <= 4 else\n        "followon_runtime_config.json" if index <= 9 else\n        "continuation_101ps_config.json"\n    )',
        "config selection by accepted segment index",
    )
    text = replace_once(
        text,
        "from typing import Any\n",
        "from typing import Any\n\nfrom .hydration_continuation_101ps import native_time_matches\n",
        "float32 time helper import",
    )
    text = replace_once(
        text,
        '    if not math.isclose(\n        float(final["time_ps"]) - report["endpoint_time_ps"],\n        6 + count * 5,\n        abs_tol=1e-6,\n    ):',
        '    observed_time = float(final["time_ps"]) - report["endpoint_time_ps"]\n'
        "    if not native_time_matches(observed_time, 6 + count * 5):",
        "precision-aware float32 final-time check",
    )
    gates.write_text(text)

    queue = package / "hydration_followon_queue.py"
    text = queue.read_text()
    old_ready = text[text.index("def initial_ready(") : text.index("\ndef hold(")]
    new_ready = 'def initial_ready(base: Path) -> bool:\n    """Require the complete accepted 51 ps parent and current analysis."""\n    execution = read(base / "followon_queue_resume_job/execution/status.json")\n    if execution.get("ended_utc") is None:\n        return False\n    if execution.get("exit_status") != 0:\n        raise ValueError("51 ps parent continuation process failed")\n    completion = read(base / "followon_completion.json")\n    if completion.get("status") != "completed_51ps_review_required":\n        raise ValueError("The 51 ps parent campaign is not complete")\n    for index in range(1, 10):\n        accepted(base, index)\n    report = read(base / "analysis/report.json")\n    if any(f"extend_{index:02d}" not in report["frame_count_by_stage"]\n           for index in range(1, 10)):\n        raise ValueError("Analysis does not cover all nine accepted extension blocks")\n    return True\n'
    text = replace_once(text, old_ready, new_ready, "51 ps parent readiness gates")
    text = replace_once(
        text,
        "if hold(plan, 4):",
        "if hold(plan, 9):",
        "initial nine-block candidate screen",
    )
    text = replace_once(
        text,
        '        f"followon_candidate_review_after_{count:02d}.json",',
        '        f"continuation_101ps_candidate_review_after_{count:02d}.json",',
        "new candidate receipt name",
    )
    text = replace_once(
        text,
        '    cfg = read(base / "extension_20ps_config.json")',
        '    cfg = read(base / "continuation_101ps_template.json")',
        "new continuation runtime config",
    )
    text = replace_once(
        text,
        '        endpoint_restart_path=str(base / "extend_04/restart.ncrst"),\n        endpoint_restart_sha256=digest(base / "extend_04/restart.ncrst"),',
        '        endpoint_restart_path=str(base / "extend_09/restart.ncrst"),\n        endpoint_restart_sha256=digest(base / "extend_09/restart.ncrst"),',
        "accepted extend09 predecessor",
    )
    text = replace_once(
        text,
        '        if index == 0:\n            (inputs / "restart.ncrst").write_bytes(\n                (base / "extend_04/restart.ncrst").read_bytes()\n            )',
        '        if index == 0:\n            (inputs / "restart.ncrst").write_bytes(\n                (base / "extend_09/restart.ncrst").read_bytes()\n            )',
        "pin first restart to extend09",
    )
    text = replace_once(
        text,
        '                "segment": index + 5,\n                "duration_ps": 5.0,\n                "predecessor": f"extend_{index + 4:02d}",',
        '                "segment": index + 10,\n                "duration_ps": 5.0,\n                "predecessor": f"extend_{index + 9:02d}",',
        "continuation segment lineage",
    )
    text = replace_once(
        text,
        '    snapshot = base / "followon_source_snapshot"',
        '    snapshot = base / "continuation_101ps_source"',
        "new source snapshot identity",
    )
    text = replace_once(
        text,
        '    save(base, "followon_runtime_config.json", cfg)',
        '    save(base, "continuation_101ps_config.json", cfg)',
        "new runtime config receipt",
    )
    text = text.replace(
        "followon_native_status.json", "continuation_101ps_native_status.json"
    )
    text = text.replace(
        "followon_resource_preflight.json", "continuation_101ps_resource_preflight.json"
    )
    text = text.replace(
        "followon_completion.json", "continuation_101ps_completion.json"
    )
    text = replace_once(
        text,
        'completion = read(base / "continuation_101ps_completion.json")',
        'completion = read(base / "followon_completion.json")',
        "read 51 ps parent completion receipt",
    )
    text = text.replace(
        "followon_candidate_review", "continuation_101ps_candidate_review"
    )
    text = text.replace("start=5", "start=10")
    text = text.replace(
        "completed_51ps_review_required", "completed_101ps_review_required"
    )
    text = replace_once(
        text,
        'completion.get("status") != "completed_101ps_review_required"',
        'completion.get("status") != "completed_51ps_review_required"',
        "require parent 51 ps completion status",
    )
    text = text.replace(
        "Five further blocks accepted: 51ps total; stop for review",
        "Ten further blocks accepted: 101ps total; stop for review",
    )
    queue.write_text(text)

    return {
        str(path.relative_to(destination)): sha256(path)
        for path in sorted(destination.rglob("*.py"))
    }


def patch_postprocess_sources(source: Path, destination: Path) -> dict[str, str]:
    """Clone postprocessing, extend receipts through 19, and correct annotation."""
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("__pycache__"))
    package = destination / "src/nyl12_covalent"
    for name in (
        "hydration_vmd_exporter.py",
        "hydration_vmd_viewer.py",
        "hydration_vmd_export_io.py",
        "vmd_movie_policy.py",
        "relay_progress_vmd.py",
        "vmd_chain_groups.py",
    ):
        shutil.copy2(ROOT / "src/nyl12_covalent" / name, package / name)
    adapter = package / "hydration_receipt_adapter.py"
    text = replace_once(
        adapter.read_text(),
        "range(1, 10)",
        "range(1, 20)",
        "accepted receipt range through extend19",
    )
    text = replace_once(
        text,
        "extend_01 through extend_09",
        "extend_01 through extend_19",
        "receipt adapter scope",
    )
    adapter.write_text(text)
    analysis = package / "deacylation_water_analysis.py"
    text = replace_once(
        analysis.read_text(),
        "water-H...Nalpha is descriptive because Nalpha is protonated",
        "water-H...Nalpha distances are descriptive and do not establish protonation; interpret them against each accepted segment's proton inventory",
        "conditional Nalpha criteria annotation",
    )
    analysis.write_text(text)
    reporting = package / "deacylation_water_reporting.py"
    text = replace_once(
        reporting.read_text(),
        '"Water-H···Nalpha values are descriptive only because Nalpha is protonated. "',
        '"Water-H···Nalpha distances are descriptive and do not establish protonation; "\n'
        '            "interpret them against the accepted segment proton inventory. "',
        "conditional Nalpha report text",
    )
    reporting.write_text(text)
    return {
        str(path.relative_to(destination)): sha256(path)
        for path in sorted(destination.rglob("*.py"))
    }
