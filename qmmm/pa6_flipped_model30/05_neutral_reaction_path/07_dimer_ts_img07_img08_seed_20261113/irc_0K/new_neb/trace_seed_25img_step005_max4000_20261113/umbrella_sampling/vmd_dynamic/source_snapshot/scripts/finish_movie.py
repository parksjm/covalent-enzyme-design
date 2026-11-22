"""Wait for original-host exports, transfer stripped data and verify the movie."""

from __future__ import annotations

import json
import logging
from pathlib import Path
import shlex
import subprocess
import sys
import time

OUT = Path(__file__).resolve().parents[2]
SSH = ["ssh", "-i", "/home/jpc/.ssh/id_ed25519_aster", "-o", "BatchMode=yes"]


def main() -> None:
    """Finalize only after both host jobs complete, failing on explicit errors."""
    for attempt in range(120):
        local = json.loads((OUT / "work/job_violet/status.json").read_text())
        remote = json.loads(
            subprocess.check_output(
                [
                    *SSH,
                    "azalea",
                    "cat " + shlex.quote(str(OUT / "work/job_azalea/status.json")),
                ],
                text=True,
            )
        )
        for host, status in [("violet", local), ("azalea", remote)]:
            if status["exit_status"] not in (None, 0):
                raise RuntimeError(f"{host} export failed: {status['exit_status']}")
        if local["exit_status"] == remote["exit_status"] == 0:
            break
        logging.info(
            "Waiting for source-host exports: violet=%s azalea=%s",
            local["exit_status"],
            remote["exit_status"],
        )
        time.sleep(15)
    else:
        raise TimeoutError("Movie exports did not finish within 30 minutes")
    subprocess.run(
        [
            "rsync",
            "-rt",
            "-e",
            shlex.join(SSH),
            "azalea:" + str(OUT / "work/azalea") + "/",
            str(OUT / "work/azalea") + "/",
        ],
        check=True,
    )
    (OUT / "work/job_azalea").mkdir()
    subprocess.run(
        [
            "rsync",
            "-rt",
            "-e",
            shlex.join(SSH),
            "azalea:" + str(OUT / "work/job_azalea") + "/",
            str(OUT / "work/job_azalea") + "/",
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "coverage",
            "run",
            "--data-file",
            str(OUT / "work/assembly_coverage"),
            "--include",
            str(OUT / "source_snapshot/scripts/assemble_movie.py"),
            str(OUT / "source_snapshot/scripts/assemble_movie.py"),
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "coverage",
            "json",
            "--data-file",
            str(OUT / "work/assembly_coverage"),
            "-o",
            str(OUT / "work/assembly_coverage.json"),
        ],
        check=True,
    )
    subprocess.run(["sha256sum", "-c", "MANIFEST.sha256"], cwd=OUT, check=True)
    logging.info("Dynamic umbrella movie complete: %s", OUT / "view_umbrella.vmd")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    main()
