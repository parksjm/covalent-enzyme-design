"""Run guarded relay chemical dimer seeds or conditional two-sided downhill tests."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import time

import numpy as np

from nyl12_covalent.nterminal_relay_ts import _dimer_input
from nyl12_covalent.relay_release import minimization_input
from nyl12_covalent.relay_ts_guarded import run_guarded
from nyl12_covalent.stationary_native import NativeState, digest, write_json
from scripts.prepare_dimer_continuation import write_restart

LOGGER = logging.getLogger(__name__)


def chemical_direction(state: NativeState, kind: str) -> np.ndarray:
    """Build a declared bond-distance-gradient mode without MM tangent mixing."""
    direction = np.zeros_like(state.xyz)
    pairs = (
        [(19916, 13916, 1.0)]
        if kind == "attack"
        else [(13917, 13905, 1.0), (13917, 13916, -1.0)]
    )
    for first, second, sign in pairs:
        i, j = first - 1, second - 1
        if not state.active[i] or not state.active[j]:
            raise ValueError("Reactive seed atoms outside original active region")
        vector = state.xyz[i] - state.xyz[j]
        vector /= np.linalg.norm(vector)
        direction[i] += sign * vector
        direction[j] -= sign * vector
    direction /= np.linalg.norm(direction)
    return direction[state.active].reshape(-1)


def main() -> None:
    """Execute one immutable native stage from explicit CLI inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--restart", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--kind", choices=("attack", "proton_transfer", "downhill"), required=True
    )
    parser.add_argument("--screen", type=Path)
    parser.add_argument("--sign", type=int, choices=(-1, 1))
    parser.add_argument("--chemical-sign-gate", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    model, restart, output = (
        p.resolve() for p in (args.model, args.restart, args.output)
    )
    output.mkdir(parents=True, exist_ok=False)
    state = NativeState(model, restart)
    energy, gradient = state.gradient(state.xyz, output / "seed_gradient")
    if args.kind == "downhill":
        if args.screen is None or args.sign is None:
            raise ValueError("Downhill stage requires explicit screen and sign")
        screen = args.screen.resolve()
        started = time.monotonic()
        receipt = screen / "scientific_acceptance.json"
        while not (
            receipt.exists() and (screen / "lifted_cartesian_direction.npy").exists()
        ):
            if time.monotonic() - started > 3 * 3600:
                raise TimeoutError(
                    "Partial curvature stage not ready within three hours"
                )
            time.sleep(15)
        # Receipt is written before the final half-step checks; wait for those.
        while True:
            result = json.loads(receipt.read_text())
            if "relative_half_step_difference" in result:
                break
            if time.monotonic() - started > 3 * 3600:
                raise TimeoutError("Directional half-step gate unfinished")
            time.sleep(15)
        if args.chemical_sign_gate is not None:
            gate = json.loads(args.chemical_sign_gate.read_text())
            if (
                not gate["passed"]
                or gate["restart_sha256"] != digest(restart)
                or gate["mode_sha256"]
                != digest(screen / "lifted_cartesian_direction.npy")
                or gate["source_precision_receipt_sha256"] != digest(receipt)
            ):
                raise ValueError(
                    "Explicit chemical-sign gate identity/acceptance failed"
                )
        elif not result["partial_negative_direction_verified"]:
            raise RuntimeError("Conditional chemical negative-curvature gate failed")
        direction = np.load(screen / "lifted_cartesian_direction.npy")
        if direction.shape != (len(state.indices) * 3,):
            raise ValueError("Partial mode/full-active inventory mismatch")
        accepted = None
        probes = []
        for displacement in (0.05, 0.1, 0.2, 0.4):
            xyz = state.xyz.copy()
            xyz[state.active] += args.sign * displacement * direction.reshape((-1, 3))
            trial, grad = state.gradient(xyz, output / f"displacement_{displacement}")
            outward = float(args.sign * np.dot(grad.reshape(-1), direction))
            probes.append(
                dict(
                    displacement_A=displacement,
                    energy_hartree=trial,
                    outward_slope=outward,
                )
            )
            if trial < energy - 1e-7 and outward < 0:
                accepted = xyz
                break
        write_json(
            output / "downhill_start.json",
            {
                "probes": probes,
                "passed": accepted is not None,
                "full_active_index_verified": False,
                "mode_origin": "Partial QM block; conditional connectivity test",
            },
        )
        if accepted is None:
            raise RuntimeError(
                "No demonstrably downhill both-energy-and-slope displacement"
            )
        start = output / "downhill_start.ncrst"
        write_restart(start, accepted, state.restart)
        mdin = minimization_input(state.mdin_source, maxcycles=12000)
        mdin = mdin.replace("maxstep=0.05", "maxstep=0.02")
        run_guarded(model, start, output / "native", mdin, allowance=5.0)
        return
    direction = chemical_direction(state, args.kind)
    np.save(output / "chemical_direction.npy", direction, allow_pickle=False)
    probes = []
    for index, step in enumerate((0.0025, 0.00125)):
        slopes = []
        for sign, label in ((-1, "minus"), (1, "plus")):
            xyz = state.xyz.copy()
            xyz[state.active] += sign * step * direction.reshape((-1, 3))
            value, grad = state.gradient(xyz, output / f"preflight_{index}" / label)
            slopes.append(float(np.dot(grad.reshape(-1), direction)))
            probes.append(dict(step_A=sign * step, energy_hartree=value))
        LOGGER.info(
            "Chemical seed curvature h=%s: %s Eh/A2",
            step,
            (slopes[1] - slopes[0]) / (2 * step),
        )
    write_json(
        output / "seed_preflight.json",
        {
            "kind": args.kind,
            "probes": probes,
            "initial_energy_hartree": energy,
            "maximum_gradient_Eh_per_bohr": float(np.max(np.abs(gradient)) * 0.529177),
            "mode_definition": "C-OG distance gradient"
            if args.kind == "attack"
            else "gradient(d(HG,N)-d(HG,OG))",
            "full_active_index_verified": False,
        },
    )
    xyz = state.xyz.copy()
    xyz[state.active] += 0.01 * direction.reshape((-1, 3))
    mode = output / "mode.ncrst"
    write_restart(mode, xyz, state.restart)
    mdin = _dimer_input(state.mdin_source, restart, mode)
    mdin = mdin.replace("maxstep=0.05", "maxstep=0.01")
    run_guarded(model, restart, output / "native", mdin)


if __name__ == "__main__":
    main()
