"""Freeze completed NEB-US data and unchanged validated primary-PMF sources."""
from __future__ import annotations
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.neb_umbrella_analysis import collect, digest


def main() -> None:
    """Copy only accepted blocks, validating identities before and after capture."""
    out = Path(__file__).resolve().parent
    campaign = out.parents[1]
    previous = out.parent / "pmf_refresh_20261121_190812"
    prepared = campaign / "prepared.json"
    cfg, windows = collect(prepared)
    if len(windows) != 63 or any(len(w["blocks"]) != 6 for w in windows):
        raise ValueError("Expected63windows with sixcompletedblocks each")
    shutil.copytree(previous / "source_snapshot", out / "source_snapshot")
    frozen = out / "frozen_inputs"
    frozen.mkdir()
    files = {}
    def capture(source: Path, destination: Path, expected: str) -> None:
        if digest(source) != expected:
            raise ValueError(f"Inputchangedbeforecapture: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        if digest(destination) != expected or digest(source) != expected:
            raise ValueError(f"Inputchangedwhilecapturing: {source}")
        files[str(destination.relative_to(out))] = {"original": str(source), "sha256": expected}
    for name in ("prepared.json", "prepared.sha256"):
        source = campaign / name
        capture(source, frozen / name, digest(source))
    for window in windows:
        for name, expected in window["source_sha256"].items():
            source = Path(name)
            capture(source, frozen / source.relative_to(campaign), expected)
    profile = previous / "frozen_inputs/previous_profile.npz"
    capture(profile, frozen / "previous_profile.npz", digest(profile))
    for source in (out / "source_snapshot").rglob("*.py"):
        files[str(source.relative_to(out))] = {"original": str(previous / source.relative_to(out)), "sha256": digest(source)}
    files[Path(__file__).name] = {"sha256": digest(Path(__file__))}
    manifest = {"frozen_utc": datetime.now(timezone.utc).isoformat(),
        "previous_analysis": str(previous),
        "watcher_generation": json.loads((campaign / "analysis/latest.json").read_text())["generation"],
        "files": files,
        "reference_grid_policy": "Unchanged originalhistogramwidth andfixedESreference; no chemical framecensoring",
        "coverage": {"windows": 63, "blocks": 378, "production_ps_total": 1890,
            "production_ps_min": 30, "production_ps_max": 30,
            "raw_frames": sum(len(b["samples"]) for w in windows for b in w["blocks"])}}
    (out / "input_source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
