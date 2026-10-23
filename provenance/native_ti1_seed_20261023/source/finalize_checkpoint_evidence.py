"""Retain completed test evidence and Git-safe copies without deleting originals."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path("/data/jpc/covalent-enzyme-design")
BASE = ROOT / "output_nyl12_native_reference_20261023"
ARCHIVE = ROOT / "provenance/native_ti1_seed_20261023"


def main() -> None:
    """Preserve successful verification and input bytes for the Git checkpoint."""
    for name in (
        "native_seed_preparation_job",
        "pilot_staging_job",
        "checker_tests_job",
        "repository_tests_job",
        "checkpoint_job",
    ):
        source_dir = BASE / name
        status = json.loads((source_dir / "status.json").read_text())
        if status["exit_status"] != 0 or not status["ended_utc"]:
            raise ValueError(f"Prerequisite did not finish successfully: {name}")
        target = ARCHIVE / "jobs" / name
        target.mkdir(parents=True, exist_ok=False)
        for source in source_dir.iterdir():
            if source.is_file():
                shutil.copy2(
                    source,
                    target
                    / (
                        source.name + ".txt" if source.suffix == ".log" else source.name
                    ),
                )
    shutil.copy2(BASE / "checker_coverage.json", ARCHIVE / "checker_coverage.json")
    shutil.copy2(Path(__file__), ARCHIVE / "source" / Path(__file__).name)
    ignored = {
        ".prmtop",
        ".rst7",
        ".log",
        ".inpcrd",
        ".nc",
        ".dcd",
        ".chk",
        ".fchk",
        ".pyc",
    }
    # Retain originals locally; Git-safe copies preserve identical bytes.
    for source in list(ARCHIVE.rglob("*")):
        if source.is_file() and source.suffix in ignored:
            copy = source.with_name(source.name + ".txt")
            if copy.exists():
                raise FileExistsError(copy)
            shutil.copy2(source, copy)
            if source.read_bytes() != copy.read_bytes():
                raise ValueError(f"Copy differs: {source}")
    manifest = {
        str(p.relative_to(ARCHIVE)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in ARCHIVE.rglob("*")
        if p.is_file()
        and p.suffix not in ignored
        and p != ARCHIVE / "source_manifest.json"
    }
    (ARCHIVE / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
