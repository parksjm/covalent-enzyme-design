"""Verify ordered post-segment commands without using scientific engines."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts import refresh_deacylation_hydration as hook


@pytest.mark.parametrize("fail_analysis", [False, True])
def test_real_children_stop_before_movie_on_failed_analysis(
    tmp_path: Path, fail_analysis: bool
) -> None:
    """Exercise argument quoting, environment, order and failure propagation."""
    project = tmp_path / "project with spaces"
    scripts = project / "scripts"
    scripts.mkdir(parents=True)
    campaign = tmp_path / "hydration with spaces"
    campaign.mkdir()
    recorder = (
        "import json, os, pathlib, sys\n"
        "with pathlib.Path('called.jsonl').open('a') as stream:\n"
        " stream.write(json.dumps({'argv':sys.argv,'omp':os.environ['OMP_NUM_THREADS'],"
        "'blas':os.environ['OPENBLAS_NUM_THREADS']})+'\\n')\n"
    )
    (scripts / "analyze_deacylation_hydration.py").write_text(
        recorder + f"raise SystemExit({int(fail_analysis)})\n"
    )
    (scripts / "export_hydration_vmd.py").write_text(recorder)
    if fail_analysis:
        with pytest.raises(subprocess.CalledProcessError):
            hook.refresh(campaign, project)
    else:
        hook.refresh(campaign, project)
    rows = [
        json.loads(line) for line in (project / "called.jsonl").read_text().splitlines()
    ]
    assert len(rows) == (1 if fail_analysis else 2)
    assert rows[0]["argv"][1:3] == ["--base", str(campaign)]
    assert all(row["omp"] == row["blas"] == "1" for row in rows)
    if not fail_analysis:
        assert rows[1]["argv"][1:3] == ["--source", str(campaign)]
        assert str(campaign / "analysis/per_water_summary.json") in rows[1]["argv"]


def test_main_uses_explicit_paths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Cover the CLI without launching actual analysis or mail/model services."""
    received = []
    monkeypatch.setattr(hook, "refresh", lambda *args: received.append(args))
    monkeypatch.setattr(
        sys,
        "argv",
        ["refresh", "--base", str(tmp_path), "--project-root", str(tmp_path)],
    )
    hook.main()
    assert received == [(tmp_path, tmp_path)]
