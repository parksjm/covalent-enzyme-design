"""Run one prepared native N-terminal relay release."""

from __future__ import annotations
import argparse
from pathlib import Path
from nyl12_covalent.relay_release_run import run

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path, required=True)
    run(parser.parse_args().seed.resolve())
