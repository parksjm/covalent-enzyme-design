"""Run the bounded native-Sander low-mode probe."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(PROJECT_ROOT / "src"), str(PROJECT_ROOT)]
from nyl12_covalent.native_hessian_probe import run_probe  # noqa: E402


def main() -> None:
    """Parse arguments and execute the finite-difference spectrum probe."""
    base = Path(
        "qmmm/pa6_flipped_model30/05_neutral_reaction_path/"
        "07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113"
    )
    calibration = base / (
        "native_gradient_0K_onecycle_20261113_attempt2/scientific_acceptance.json"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, default=base / "system.prmtop")
    parser.add_argument(
        "--restart", type=Path, default=base / "calculation/optimized.ncrst"
    )
    parser.add_argument(
        "--mode", type=Path, default=base / "calculation/path_tsmode.xyz"
    )
    parser.add_argument("--mdin", type=Path, default=base / "calculation/dimer.mdin")
    parser.add_argument("--calibration", type=Path, default=calibration)
    parser.add_argument(
        "--sander",
        type=Path,
        default=Path("/data/jpc/code/ambertools26/venv/bin/sander"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-iterations", type=int, default=80)
    parser.add_argument("--count", type=int, default=2)
    parser.add_argument("--max-hvp-per-step", type=int, default=200)
    parser.add_argument("--subspace-size", type=int, default=12)
    parser.add_argument(
        "--solver", choices=("arpack", "block_lobpcg"), default="arpack"
    )
    parser.add_argument("--tolerance", type=float, default=1.0e-3)
    parser.add_argument("--step", type=float, nargs="+", default=[0.01])
    parser.add_argument(
        "--initial-spectrum",
        type=Path,
        help="Prior spectrum JSON with two Ritz modes for block_lobpcg warm start",
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_probe(
        topology_path=args.topology,
        restart_path=args.restart,
        mdin_path=args.mdin,
        mode_path=args.mode,
        calibration_path=args.calibration,
        sander_path=args.sander,
        output_dir=args.output_dir,
        step_sizes=tuple(args.step),
        max_iterations=args.max_iterations,
        tolerance=args.tolerance,
        max_hvp_per_step=args.max_hvp_per_step,
        subspace_size=args.subspace_size,
        solver=args.solver,
        initial_spectrum_path=args.initial_spectrum,
        count=args.count,
    )


if __name__ == "__main__":
    main()
