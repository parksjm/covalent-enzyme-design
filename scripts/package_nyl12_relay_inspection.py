"""Build an inspection bundle from two explicitly accepted relay-repair results."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.prepare_nyl12_inspection_download import (  # noqa: E402
    build_package,
    selected_files,
)

LOGGER = logging.getLogger(__name__)


def package_repairs(
    repo: Path,
    pa6_result: Path,
    maxa_result: Path,
    journal: Path,
    output: Path,
) -> dict:
    """Select accepted repair artifacts and create a fresh integrity-checked bundle.

    Args:
        repo: Repository with retained QM artifacts and versioned documentation.
        pa6_result: Explicit PA6 L4 repair result directory.
        maxa_result: Explicit MAXA|MAXA repair result directory.
        journal: Frozen action-journal snapshot to include.
        output: New bundle directory; existing directories are never overwritten.

    Returns:
        Exact source-file and SHA-256 manifest for the exported bundle.
    """
    resume = repo / "output_nyl12_resume_20261021"
    files = selected_files(
        repo,
        resume / "qm_conformer_audit_v2/report.json",
        resume / "qm_conformer_audit_v2/structures",
        resume / "qm_trajectories",
    )
    files["provenance/agent_actions_at_export.jsonl"] = journal
    for target, result in (("pa6_l4", pa6_result), ("maxa_maxa", maxa_result)):
        acceptance = json.loads((result / "acceptance.json").read_text())
        if acceptance.get("accepted") is not True:
            raise ValueError(f"Cannot publish unaccepted repair for {target}: {result}")
        provenance = json.loads((result / f"{target}_ti1_provenance.json").read_text())
        if provenance.get("target") != target:
            raise ValueError(f"Repair target mismatch: {result}")
        for suffix in (
            "ti1_starting.pdb",
            "residue_map.json",
            "ti1_audit.json",
            "ti1_provenance.json",
        ):
            files[f"protein/{target}_{suffix}"] = result / f"{target}_{suffix}"
        for filename in (
            "search.json",
            "acceptance.json",
            "invariance.json",
            "original_preflight.json",
        ):
            files[f"reports/{target}_relay_repair/{filename}"] = result / filename
        files[f"reports/{target}_original_active_site.json"] = (
            result / "original_active_site.json"
        )
        files[f"reports/{target}_active_site.json"] = result / "active_site.json"
        attempt = result.parent
        for filename in ("command.json", "status.json", "run.sh"):
            files[f"provenance/{target}_repair_job/{filename}"] = (
                attempt / "job" / filename
            )
    for name in (
        "repair_nyl12_relay_water.py",
        "launch_nyl12_relay_repair.py",
        Path(__file__).name,
    ):
        files[f"reproduction/scripts/{name}"] = repo / "scripts" / name
    for path in sorted((repo / "src/nyl12_covalent").glob("relay_repair*.py")):
        files[f"reproduction/src/nyl12_covalent/{path.name}"] = path
    return build_package(files, output)


def main() -> None:
    """Package explicit accepted results without changing source artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--pa6-result", required=True, type=Path)
    parser.add_argument("--maxa-result", required=True, type=Path)
    parser.add_argument("--journal", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    report = package_repairs(
        args.repo, args.pa6_result, args.maxa_result, args.journal, args.output
    )
    LOGGER.info(
        "Exported %s files (%s bytes)", len(report["files"]), report["payload_bytes"]
    )


if __name__ == "__main__":
    main()
