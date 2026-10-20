"""Create an explicit, hashed source inventory for covalent-pipeline extraction."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
import logging
import subprocess
import tarfile
from pathlib import Path
from typing import Sequence

LOGGER = logging.getLogger(__name__)
REQUIRED = {
    **{f"src/nyl12_covalent/{name}.py": "package module" for name in (
        "__init__", "af3", "af3_postprocess", "chemistry", "inspection", "jobs",
        "parameter_cache", "pose", "qm_region", "qmmm", "relay_water", "tetrahedral",
        "topology", "validate_ti1", "validate_topology",
    )},
    **{f"scripts/{name}.py": "pipeline entry point" for name in (
        "build_nyl12_acyl_qmmm", "build_nyl12_ti1_qmmm", "canonicalize_nyl12_af3_ti1",
        "diagnose_nyl12_relay_water_contacts", "prepare_nyl12_af3_ti1", "run_nyl12_job",
        "search_nyl12_relay_water", "audit_nyl12_ti1_qm_region", "preflight_gaussian16",
    )},
    "tests/test_nyl12_inspection.py": "module tests",
    "tests/test_nyl12_jobs.py": "module tests",
    "tests/test_nyl12_parameter_cache.py": "module tests",
    "tests/test_nyl12_qm_region.py": "QM-region module tests",
    "tests/test_diagnose_nyl12_relay_water_contacts.py": "relay-water diagnostic tests",
    "tests/test_inventory_covalent_pipeline.py": "inventory generator tests",
    "tests/test_preflight_gaussian16.py": "Gaussian preflight tests",
    "docs/nyl12_covalent_resume_20261019.md": "workflow handoff",
    "docs/nyl12_pa6_covalent_qmmm_plan_20261019.md": "scientific plan",
    "docs/covalent_pipeline_repository_plan_20261020.md": "extraction plan",
    "tools/skills/scientific-tmux-jobs/SKILL.md": "skill instructions",
    "tools/skills/scientific-tmux-jobs/agents/openai.yaml": "skill metadata",
    "tools/skills/scientific-tmux-jobs/references/pipeline.md": "skill adapter",
    "tools/covalent_repository_template/README.md": "extracted repository guide",
    "tools/covalent_repository_template/AGENTS.md": "extracted repository instructions",
    "tools/covalent_repository_template/pyproject.toml": "extracted repository configuration",
    "tools/covalent_repository_template/.gitignore": "extracted repository ignores",
}


def _git(root: Path, *args: str) -> str:
    """Run a read-only Git query in the source tree."""
    return subprocess.run(["git", *args], cwd=root, check=True, text=True,
                          capture_output=True).stdout.strip()


def build_inventory(root: Path) -> tuple[dict[str, object], dict[str, bytes]]:
    """Build inventory data and exact bytes for an optional archive.

    Args:
        root: Repository root to inventory.

    Returns:
        JSON-ready metadata and allowlisted file contents keyed by safe paths.
    """
    root = root.resolve(strict=True)
    missing = sorted(name for name in REQUIRED if not (root / name).is_file())
    if missing:
        raise FileNotFoundError(f"Required extraction files missing: {missing}")
    allow = dict(REQUIRED)
    self_name = Path(__file__).resolve().relative_to(root).as_posix()
    allow[self_name] = "inventory generator"
    tracked = set(filter(None, _git(root, "ls-files").splitlines()))
    contents: dict[str, bytes] = {}
    entries = []
    for name, role in sorted(allow.items()):
        path = root / name
        resolved = path.resolve(strict=True)
        if path.is_symlink() or not path.is_file() or not resolved.is_relative_to(root):
            raise ValueError(f"Allowlisted path must be a regular file: {name}")
        data = path.read_bytes()
        contents[name] = data
        changed = bool(_git(root, "status", "--porcelain=v1", "--", name))
        state = "untracked" if name not in tracked else "tracked_modified" if changed else "tracked_clean"
        entries.append({"path": name, "role": role, "size_bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest(), "git_status": state})
    return ({"source_root": str(root), "head": _git(root, "rev-parse", "HEAD"),
             "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
             "head_captures_untracked_changes": False, "files": entries}, contents)


def write_snapshot(path: Path, contents: dict[str, bytes]) -> None:
    """Archive exact allowlisted bytes under safe relative names."""
    with tarfile.open(path, "x:gz") as archive:
        for name, data in sorted(contents.items()):
            member = tarfile.TarInfo(name)
            member.size = len(data)
            member.mode = 0o644
            archive.addfile(member, io.BytesIO(data))


def main(argv: Sequence[str] | None = None) -> int:
    """Write inventory JSON and optionally a matching source tarball."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--snapshot", type=Path, help="Optional .tar.gz source snapshot")
    args = parser.parse_args(argv)
    outputs = [args.output, *([args.snapshot] if args.snapshot else [])]
    if len({path.absolute() for path in outputs}) != len(outputs):
        parser.error("output and snapshot paths must differ")
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        parser.error(f"Refusing to overwrite: {existing}")
    inventory, contents = build_inventory(args.root)
    if args.snapshot:
        args.snapshot.parent.mkdir(parents=True, exist_ok=True)
        write_snapshot(args.snapshot, contents)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    LOGGER.info("Wrote %d files to %s", len(contents), args.output)
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
