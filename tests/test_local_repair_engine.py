"""Exercise engine exit recording and fail-fast SCF monitoring."""

import os
from pathlib import Path
import sys

from nyl12_covalent.local_repair_engine import run_checked_engine


def test_records_success_without_scf_warning(tmp_path: Path) -> None:
    result = run_checked_engine(
        [sys.executable, "-c", "pass"],
        cwd=tmp_path,
        env=os.environ.copy(),
        log_path=tmp_path / "stdout.log",
        mdout=tmp_path / "mdout",
    )
    assert result == {"exit_status": 0, "stopped_on_scf_failure": False}


def test_stops_process_on_inaccurate_force_warning(tmp_path: Path) -> None:
    program = "from pathlib import Path; import time; Path('mdout').write_text('Convergence could not be achieved in this step'); time.sleep(60)"
    result = run_checked_engine(
        [sys.executable, "-c", program],
        cwd=tmp_path,
        env=os.environ.copy(),
        log_path=tmp_path / "stdout.log",
        mdout=tmp_path / "mdout",
    )
    assert result["exit_status"] != 0
    assert result["stopped_on_scf_failure"]
