"""Checkpoint a verified diagnostic movie and correct status documents."""

from pathlib import Path
import hashlib
import json
import shutil
import sys

ROOT = Path("/data/jpc/covalent-enzyme-design")
sys.path.insert(0, str(ROOT))
from scripts.prepare_nyl12_inspection_download import build_package  # noqa: E402


def main() -> None:
    """Archive small movie artifacts and retain the previous download revision."""
    note = """## Pilot ended: TI1 state loss; trajectory available — 2026-10-23

The63atom donor/water pilot ended17:10:33EDT at20cycle limit, exit1.
No convergence and no TI1 acceptance: OG1-C54=1.90553A, N-H3=1.40361A,
OG1-H3=1.03419A, consistent with reverse proton transfer. WaterO-N=3.30394A
and waterO-OG1=3.53239A; an acceptable pseudo-ring is not established.
No continuation is running. Earlier ACTIVE entries below are historical.

Exported21complete Input orientation tables,63atoms and37explicitH each,
with both water H. Workstation snapshot:
output_nyl12_oxyanion_cluster_20261023/trajectory_snapshot_01/.
The usual rsync pull now includes
qm/ti1_oxyanion_cluster/optimization_progress/view_optimization.pml,
optimization_trajectory.pdb and .xyz, CSV metrics and acceptance report.
Run @view_optimization.pml then mplay from that directory. This is a
snapshot, not a live feed; fixed PDB bonds do not prove connectivity survived.
Source log snapshot/hash, exporter and tmux records preserved. Ruff and
PDB frame/atom/H roundtrip checks pass. No new QM or force-field run launched.

"""
    handoff = ROOT / "docs/standalone_handoff_20261020.md"
    handoff.write_text(note + handoff.read_text())
    guide = ROOT / "docs/laptop_inspection_20261021.md"
    guide.write_text(guide.read_text() + "\n" + note)
    base = ROOT / "output_nyl12_oxyanion_cluster_20261023"
    archive = ROOT / "provenance/ti1_donor_movie_20261023"
    archive.mkdir(exist_ok=False)
    sources = [
        base / "export_progress.py",
        Path(__file__),
        handoff,
        guide,
        base / "svp_pilot_01/scientific_acceptance.json",
        base / "movie_export_job/status.json",
        base / "movie_export_job/command.json",
        base / "movie_export_job/run.sh",
        ROOT / "output_agent_actions/20261021.jsonl",
    ]
    sources += [
        base / "trajectory_snapshot_01" / name
        for name in (
            "trajectory_report.json",
            "trajectory_metrics.csv",
            "view_optimization.pml",
        )
    ]
    for i, source in enumerate(sources):
        shutil.copy2(source, archive / f"{i:02d}_{source.name}")
    (archive / "source_manifest.json").write_text(
        json.dumps(
            {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
            indent=2,
        )
        + "\n"
    )
    current = ROOT / "output_laptop_inspection_20261021"
    prior = json.loads((current / "provenance.json").read_text())
    files = {name: Path(row["source"]) for name, row in prior["files"].items()}
    files["README.md"] = archive / "03_laptop_inspection_20261021.md"
    files["docs/laptop_inspection_20261021.md"] = files["README.md"]
    files["docs/standalone_handoff_20261020.md"] = (
        archive / "02_standalone_handoff_20261020.md"
    )
    files["provenance/agent_actions_at_export.jsonl"] = archive / "09_20261021.jsonl"
    staging = ROOT / "output_laptop_inspection_20261021_donor_movie_docs_staging"
    backup = ROOT / "output_laptop_inspection_20261021_revision_24"
    if backup.exists() or staging.exists():
        raise FileExistsError("Fresh package paths required")
    build_package(files, staging)
    for line in (staging / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        if hashlib.sha256((staging / name).read_bytes()).hexdigest() != digest:
            raise ValueError("Package verification failed")
    current.rename(backup)
    staging.rename(current)


if __name__ == "__main__":
    main()
