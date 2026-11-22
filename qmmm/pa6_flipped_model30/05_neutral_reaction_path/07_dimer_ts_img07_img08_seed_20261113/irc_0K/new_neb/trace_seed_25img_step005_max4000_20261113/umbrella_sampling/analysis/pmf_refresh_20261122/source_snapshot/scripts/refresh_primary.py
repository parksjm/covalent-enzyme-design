"""Refresh frozen accepted-block MBAR and the previously validated q_acyl marginal."""
from pathlib import Path
import sys
import requested_path_pmf as base
import primary_q_acyl


def main() -> None:
    """Analyze immutable inputs without touching the running sampling campaign."""
    out = Path(__file__).resolve().parents[2]
    base.refresh(out / "frozen_inputs/prepared.json", out / "mbar", plot_python=Path(sys.executable))
    primary_q_acyl.main()


if __name__ == "__main__":
    main()
