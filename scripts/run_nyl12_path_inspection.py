"""Update endpoint movie and distance diagnostics after every completed window."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.path_inspection import publish_inspection  # noqa: E402

LOGGER = logging.getLogger(__name__)


def _terminal_ready(
    state: str,
    child_present: bool,
    child_entered_active_state: bool,
    child_finished: bool,
) -> bool:
    return state in {"completed", "stopped"} and (
        not child_present or child_entered_active_state or child_finished
    )


def main() -> None:
    """Inspect an existing path or supervise a new controller with inspection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--poll-seconds", type=float, default=15)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    if args.poll_seconds <= 0:
        raise ValueError("Polling interval must be positive")
    command = args.command
    if command[:1] == ["--"]:
        command = command[1:]
    if command and not args.watch:
        raise ValueError("A supervised path command requires --watch")
    # Require a prepared run and all inspection prerequisites before calculations.
    config = json.loads(args.config.read_text())
    for name in ("cpptraj", "plot_python", "viewer_assets"):
        if not Path(config[name]).exists():
            raise FileNotFoundError(config[name])
    status_path = args.run / "status.json"
    initial_state = json.loads(status_path.read_text())
    child = subprocess.Popen(command, cwd=ROOT) if command else None
    terminal_states = {"completed", "stopped"}
    child_entered_active_state = (
        child is None or initial_state["state"] not in terminal_states
    )
    signature = None
    try:
        while True:
            text = status_path.read_text()
            state = json.loads(text)
            if state["state"] not in terminal_states:
                child_entered_active_state = True
            if text != signature:
                if state["completed_windows"]:
                    output = publish_inspection(
                        args.run.resolve(),
                        args.config.resolve(),
                        ROOT / "scripts/plot_reaction_distances.py",
                    )
                    LOGGER.info("Inspection updated: %s", output)
                signature = text
            child_finished = child is not None and child.poll() is not None
            terminal_ready = _terminal_ready(
                state["state"],
                child is not None,
                child_entered_active_state,
                child_finished,
            )
            if not args.watch or terminal_ready:
                break
            if child_finished:
                raise RuntimeError("Path controller exited without a terminal status")
            time.sleep(args.poll_seconds)
        if child and child.wait() != 0:
            raise RuntimeError("Path controller failed; diagnostics retained")
    finally:
        if child and child.poll() is None:
            child.terminate()
            child.wait()


if __name__ == "__main__":
    main()
