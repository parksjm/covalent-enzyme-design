"""Compose a NEW export inventory after the root confirms docs/tools are frozen."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE / "tools"))
from student_bundle_support import (  # noqa: E402
    RESERVED,
    SOFTWARE_MAGIC,
    SOFTWARE_SUFFIXES,
    digest,
    regular_source,
    safe_relative,
)

LIMIT = 5 * 1024 * 1024
PREFIXES = ("src", "scripts", "tests", "docs", "manuscript", "configs", "cif")
ROOT_FILES = ("AGENTS.md", "README.md", "pyproject.toml", "uv.lock")
CODE_SUFFIXES = {
    ".py",
    ".sh",
    ".md",
    ".toml",
    ".lock",
    ".json",
    ".yaml",
    ".yml",
    ".tcl",
    ".vmd",
    ".cif",
    ".txt",
    ".csv",
    ".pml",
}
RAW_SUFFIXES = {".nc", ".ncrst", ".rst7", ".prmtop", ".inpcrd", ".dcd", ".chk", ".fchk"}
EXCLUDED_NAMES = {
    "export_sources.json",
    "inventory_supplement_entries.json",
    "historical_mutant_structure_entries.json",
}


def git(*arguments: str) -> str:
    """Read repository metadata without staging or modifying Git state."""
    return subprocess.run(
        ["git", "-C", str(ROOT), *arguments], check=True, capture_output=True, text=True
    ).stdout


def finalize(frozen_path: Path, mutant_path: Path, output: Path) -> dict[str, Any]:
    """Compose and verify an explicit export inventory without copying source data.

    Args:
        frozen_path: Original immutable data_manifest.json.
        mutant_path: Reviewed18-entry historical mutant supplement.
        output: New nonexistent export_sources.json destination.

    Returns:
        Composed manifest with source hashes, omissions and dirty-snapshot metadata.
    """
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    frozen_sha, mutant_sha = digest(frozen_path), digest(mutant_path)
    frozen = json.loads(frozen_path.read_text())
    supplement = json.loads(mutant_path.read_text())
    entries: dict[str, dict[str, Any]] = {}
    omissions: list[dict[str, str]] = []

    def insert(entry: dict[str, Any]) -> None:
        name = str(safe_relative(entry["destination"]))
        if name in RESERVED:
            raise ValueError(f"Reserved bundle-root destination: {name}")
        if name not in entries:
            entries[name] = entry.copy()

    for entry in frozen["entries"]:
        insert(entry)
    for entry in supplement["entries"]:
        insert(entry)

    def add(
        path: Path,
        destination: str,
        profile: str,
        category: str,
        expected: str | None = None,
    ) -> None:
        name = str(safe_relative(destination))
        if name in entries:
            if expected is not None and entries[name]["sha256"] != expected:
                raise ValueError(f"Frozen entry conflicts with VMD checksum: {name}")
            return
        regular_source(path)
        checksum = digest(path)
        if expected is not None and checksum != expected:
            raise ValueError(f"VMD checksum mismatch: {path}")
        insert(
            {
                "source": str(path.absolute()),
                "destination": name,
                "profile": profile,
                "required": True,
                "bytes": path.stat().st_size,
                "sha256": checksum,
                "category": category,
            }
        )

    # Close the existing immutable VMD manifests; do not rewrite their contents.
    done: set[str] = set()
    while True:
        pending = [
            e
            for n, e in entries.items()
            if n not in done
            and e["profile"] != "raw_reference_only"
            and PurePosixPath(n).name == "MANIFEST.sha256"
        ]
        if not pending:
            break
        for entry in pending:
            done.add(entry["destination"])
            path = Path(entry["source"])
            regular_source(path)
            if digest(path) != entry["sha256"]:
                raise ValueError(f"Selected VMD manifest changed: {path}")
            for line in path.read_text().splitlines():
                match = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
                if match is None:
                    raise ValueError(f"Invalid VMD checksum record: {path}")
                relative = safe_relative(match[2])
                destination = str(PurePosixPath(entry["destination"]).parent / relative)
                add(
                    path.parent.joinpath(*relative.parts),
                    destination,
                    entry["profile"],
                    "VMD_checksum_member",
                    match[1],
                )

    def small_source(path: Path, destination: str, category: str) -> None:
        if path.is_symlink() or any(parent.is_symlink() for parent in path.parents):
            omissions.append({"source": str(path), "reason": "symlink omitted"})
            return
        if not path.is_file():
            omissions.append({"source": str(path), "reason": "not a regular file"})
            return
        if path.stat().st_size > LIMIT:
            omissions.append({"source": str(path), "reason": "larger than5MiB"})
            return
        if path.suffix.lower() in SOFTWARE_SUFFIXES:
            omissions.append({"source": str(path), "reason": "software binary suffix"})
            return
        with path.open("rb") as stream:
            magic = stream.read(4)
        if magic in SOFTWARE_MAGIC or magic.startswith(b"MZ"):
            omissions.append(
                {"source": str(path), "reason": "software executable magic"}
            )
            return
        add(path, destination, "core", category)

    for path in sorted(HERE.rglob("*")):
        relative = path.relative_to(HERE)
        if any(
            part in {"__pycache__", ".pytest_cache", ".ruff_cache"}
            for part in relative.parts
        ):
            continue
        if (
            path.name in EXCLUDED_NAMES
            or path.name.startswith("data_manifest")
            or path.suffix == ".pyc"
        ):
            continue
        if path.suffix.lower() in RAW_SUFFIXES:
            continue
        if path.is_file() or path.is_symlink():
            destination = (
                "handoff/README.md"
                if relative.as_posix() == "README.md"
                else relative.as_posix()
            )
            small_source(path, destination, "student_handoff_source_or_document")

    paths = git(
        "ls-files",
        "--cached",
        "--others",
        "--exclude-standard",
        "-z",
        "--",
        *PREFIXES,
        *ROOT_FILES,
    ).split("\0")
    for name in sorted(
        set(filter(None, paths)) | {n for n in ROOT_FILES if (ROOT / n).exists()}
    ):
        relative = safe_relative(name)
        if not (relative.parts[0] in PREFIXES or name in ROOT_FILES):
            continue
        path = ROOT.joinpath(*relative.parts)
        if path.suffix.lower() not in CODE_SUFFIXES:
            continue
        small_source(
            path, f"project_code/{relative}", "selected_working_tree_project_code"
        )

    # Verify every included source against the final explicit inventory.
    for entry in entries.values():
        if entry["profile"] == "raw_reference_only":
            continue  # Never open/hash the excluded75GB raw references.
        path = Path(entry["source"])
        try:
            regular_source(path)
        except FileNotFoundError:
            if entry["required"]:
                raise
            omissions.append(
                {"source": str(path), "reason": "optional frozen source missing"}
            )
            continue
        with path.open("rb") as stream:
            magic = stream.read(4)
        if (
            path.suffix.lower() in SOFTWARE_SUFFIXES
            or magic in SOFTWARE_MAGIC
            or magic.startswith(b"MZ")
        ):
            raise ValueError(f"Software binary in selected inventory: {path}")
        if path.stat().st_size != entry["bytes"] or digest(path) != entry["sha256"]:
            raise ValueError(f"Selected source size/SHA changed: {path}")
    if digest(frozen_path) != frozen_sha or digest(mutant_path) != mutant_sha:
        raise RuntimeError("Input manifest changed during composition")
    status = git("status", "--porcelain", "--untracked-files=normal")
    counts = Counter(e["profile"] for e in entries.values())
    record = {
        **frozen,
        "entries": list(entries.values()),
        "export_composition": {
            "frozen_manifest": str(frozen_path.absolute()),
            "frozen_manifest_sha256": frozen_sha,
            "mutant_supplement_sha256": mutant_sha,
            "git_base_commit": git("rev-parse", "HEAD").strip(),
            "working_tree_dirty": bool(status.strip()),
            "project_snapshot": "Selected current working-tree bytes, including uncommitted/untracked files; NOT a clean-commit snapshot.",
            "profile_entry_counts": dict(counts),
            "all_available_selected_sources_verified": True,
            "VMD_manifests_closed": sorted(done),
            "omissions": omissions,
            "raw_references_included": False,
            "scientific_sources_modified": False,
        },
    }
    with output.open("x") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return record


def main(argv: list[str] | None = None) -> None:
    """Write a new metadata inventory only after docs/tools are explicitly frozen."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=HERE / "data_manifest.json")
    parser.add_argument(
        "--mutants",
        type=Path,
        default=HERE / "historical_mutant_structure_entries.json",
    )
    parser.add_argument("--output", type=Path, default=HERE / "export_sources.json")
    args = parser.parse_args(argv)
    finalize(args.manifest, args.mutants, args.output)


if __name__ == "__main__":
    main()
