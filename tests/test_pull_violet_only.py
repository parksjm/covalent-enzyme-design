"""Exercise the entire laptop dry-run without network or nested remote refreshes."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


PULL = Path(__file__).resolve().parents[1] / "scripts/pull_nyl12_inspection.sh"
STUB = """import json
import os
from pathlib import Path
import sys

name = Path(sys.argv[0]).name
args = sys.argv[1:]
with Path(os.environ["PULL_STUB_LOG"]).open("a") as stream:
    stream.write(json.dumps({"command": name, "argv": args}) + "\\n")
if name == "ssh":
    pointer_probe = len(args) == 4 and args[:3] == ["jpc@violet", "readlink", "-f"] and args[3].startswith("/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/") and args[3].endswith(("/latest_verified", "/current"))
    if pointer_probe:
        print(args[3].rsplit("/", 1)[0] + ("/verified_fixture" if args[3].endswith("/latest_verified") else "/checkpoint_fixture"))
    elif len(args) != 4 or args[0] != "jpc@violet" or args[1] != "test" or args[2] not in {"-f", "-d"}:
        sys.exit("Only read-only violet existence probes are allowed")
elif name == "rsync":
    if "--dry-run" not in args:
        sys.exit("Only dry-run transfers are allowed")
    sources = [arg for arg in args if not arg.startswith("-") and ":" in arg]
    if not sources or any(arg.split(":", 1)[0] != "jpc@violet" for arg in sources):
        sys.exit("Every transfer must read existing files from violet")
else:
    sys.exit("Unexpected stub command")
"""


def test_whole_pull_dry_run_uses_only_existing_violet_copies(tmp_path: Path) -> None:
    """No refresh is executed, and final water plus current NEB downloads are reached."""
    commands = tmp_path / "commands"
    commands.mkdir()
    for name in ("ssh", "rsync"):
        executable = commands / name
        executable.write_text(f"#!{sys.executable}\n{STUB}")
        executable.chmod(0o755)
    log = tmp_path / "commands.jsonl"
    destination = tmp_path / "download"
    environment = dict(os.environ)
    environment.update(
        PATH=f"{commands}{os.pathsep}{environment['PATH']}", PULL_STUB_LOG=str(log)
    )
    result = subprocess.run(
        ["bash", str(PULL), "--dry-run", "jpc@violet", str(destination)],
        env=environment,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert result.returncode == 0, result.stderr
    records = [json.loads(line) for line in log.read_text().splitlines()]
    probes = [r["argv"] for r in records if r["command"] == "ssh"]
    transfers = [r["argv"] for r in records if r["command"] == "rsync"]
    assert probes and transfers
    assert all(
        args[:2] in (["jpc@violet", "test"], ["jpc@violet", "readlink"])
        for args in probes
    )
    assert {args[2] for args in probes} == {"-f", "-d"}
    assert sum(args[1:3] == ["readlink", "-f"] for args in probes) == 3
    assert all("--dry-run" in args for args in transfers)
    assert all(
        not any("sync_acylenzyme_progress" in arg for arg in args) for args in probes
    )
    sources = [
        arg
        for args in transfers
        for arg in args
        if not arg.startswith("-") and ":" in arg
    ]
    assert all(source.startswith("jpc@violet:") for source in sources)
    assert any("/06_asynchronous_umbrella_30ps_azalea/" in source for source in sources)
    benchmark = next(
        i for i, source in enumerate(sources) if "/final_analysis/" in source
    )
    final_water = next(
        i
        for i, source in enumerate(sources)
        if "/07_water_to_asp272_ti1_like/vmd_progress/" in source
    )
    assert final_water > benchmark
    relay = next(
        args
        for args in transfers
        if any(arg.endswith("/09_water_relay_acylation/") for arg in args)
    )
    assert "--include=/04_neb/calculation/minimize.mdin" in relay
    assert "--include=/04_neb/progress/latest/***" in relay
    assert any("/new_neb/plots/" in source for source in sources)
    assert not destination.exists()
