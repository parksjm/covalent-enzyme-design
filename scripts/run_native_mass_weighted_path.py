#!/usr/bin/env python3
"""Run or extend both native-Sander numerical descent branches from the PA6 TS."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from nyl12_covalent.native_irc_driver import run_two_branches

BASE = Path(
    "qmmm/pa6_flipped_model30/05_neutral_reaction_path/"
    "07_dimer_ts_img07_img08_seed_20261113/continuation_cycle200_20261113"
)
TS_CALC = BASE / "calculation"


def main() -> None:
    """Parse immutable scientific inputs and execute the requested bounded path."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, default=BASE / "system.prmtop")
    parser.add_argument("--restart", type=Path, default=TS_CALC / "optimized.ncrst")
    parser.add_argument("--mdin", type=Path, default=TS_CALC / "dimer.mdin")
    parser.add_argument("--path-xyz", type=Path, default=TS_CALC / "path.xyz")
    parser.add_argument("--mode-xyz", type=Path, default=TS_CALC / "path_tsmode.xyz")
    parser.add_argument(
        "--sander",
        type=Path,
        default=Path("/data/jpc/code/ambertools26/venv/bin/sander"),
    )
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--distance-pairs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--accepted-steps",
        type=int,
        default=20,
        help="Total accepted steps desired per branch; an existing run can be extended",
    )
    parser.add_argument("--initial-step-size", type=float, default=0.02)
    parser.add_argument("--minimum-step-size", type=float, default=1e-5)
    parser.add_argument(
        "--gradient-rms-stop",
        type=float,
        help="Optional active mass-weighted gradient RMS threshold in kcal/mol/angstrom/sqrt(amu)",
    )
    parser.add_argument(
        "--energy-plateau-window",
        type=int,
        default=0,
        help="Stop after this many recent attempts span <= plateau tolerance; 0 disables",
    )
    parser.add_argument("--energy-plateau-tolerance-hartree", type=float, default=1e-8)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    result = run_two_branches(
        topology_path=args.topology,
        restart_path=args.restart,
        mdin_path=args.mdin,
        path_xyz=args.path_xyz,
        mode_path=args.mode_xyz,
        sander=args.sander,
        output_dir=args.output_dir,
        calibration_path=args.calibration,
        distance_pairs_path=args.distance_pairs,
        accepted_steps=args.accepted_steps,
        initial_step_size=args.initial_step_size,
        minimum_step_size=args.minimum_step_size,
        gradient_rms_stop=args.gradient_rms_stop,
        energy_plateau_window=args.energy_plateau_window,
        energy_plateau_tolerance_hartree=args.energy_plateau_tolerance_hartree,
    )
    logging.info("Path branches written to %s", result)


if __name__ == "__main__":
    main()
