"""Prepare, but never launch, the bounded N-terminal relay umbrella pilot."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from nyl12_covalent.neb_umbrella_inputs import STAGE_STEPS, build_mdin
from nyl12_covalent.relay_umbrella_cv import RelayBiasWindow

PILOT_FORCE_CONSTANT = 25.0
PILOT_OUTER_EXTENT_A = 10.0
PILOT_SEED_BASE = 2_026_112_500
STAGES = (
    ("01_heating_1ps", "heat", 2_000),
    ("02_equilibration_1ps", "equil", 2_000),
    ("03_production_2ps", "production", 4_000),
)
MOBILE_MASK = ":27-32,34,43-56,61,70-81,83,108-124,141,144-170,181-192,199-201,211,215,222-230,232,277,279,736-738,750-758,760-765,783,787,790,792-805,929-934,946,948-974,1006-1019,1037,1039-1040,1042-1051,1053-1054,1080-1097,1099-1100,1313-1316,1580,1639,1962,2010,2779,2942,3012,3152,3327,3528,4023,4079,4165,4529,4883,4917,4938,5084,5706,5709,5795,6363,6578,6687,6850,7060,7157,7321,7526,7766,8129,8194,8313,8357,8519,8530,8602,8607,8626,8632,8643,8695,8801,8826,8842,9270,9368,9381,9575,9845,9914,9980,10148,10260,10591,10594,10630,10664-10665,10669,10834,10881,11029,11088,11091,11108,11121,11137,11141,11523,11543,11609,11803,11818,12165,12507,12521,12889,13192,13712,13717,14449,14521,14601,14619,14694,14706,14860,15096,15728,15741,15909,15945,16110,16175,16273,16299,16403,16467,16499,16512,16585,16640,16651,16662,16670,16680,16783,16919,17370,17435,17448,17454-17455,17462,17486,17491,17497,17515,17564,17588-17589,17789,17827,18240,18323,18789,19041,19620,19919,20400,20505,21023,21288,21305,21674,21682,21697,21702,21713,21744,21792,21847,21873,21879,22201,22246,22254,22259,22352,22355,22373,22377,22382,22392,22396,22412,22414,22422,22473-22474,22504,22513,22535,22591,22622,22646,22658,22665,22676,22696,22708,22763,23143,23173,23182,23204,23217,23241,23266,23286,23335,23424,25439,25441,25664,25896,25954,26056,26160,26298,26361,26543,26714,26814,26993,26996,27090,27118,27154,27186,27219,27222,27231,27662,28155,28979,29147,29459,29695"
QM_MASK = "@1621-1655,2727-2734,13905-13931,14555-14560,14574-14579,19914-19917,19925-19931,19937-19941"
ATOM_IDS = {"carbonyl_C": 19916, "Nleave": 19929, "OG1": 13916, "Nalpha": 13905, "HG1": 13917}


@dataclass(frozen=True)
class PilotWindow:
    """One explicit CV target, selected seed, and Amber bias pair."""

    label: str
    chart: str
    bias: RelayBiasWindow
    seed_source: Path
    seed_evidence: Path
    seed_sha256: str
    seed_id: str
    seed_window: str
    seed_distance_A: float
    seed_offset_A: tuple[float, float]


def digest(path: Path) -> str:
    """Return streaming SHA-256 for a file."""
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def pilot_stage_input(
    old_mdin: str,
    stage: str,
    seed: int,
    restraint_path: Path,
    active_mask: str = MOBILE_MASK,
) -> str:
    """Reuse the validated native input builder with bounded pilot durations."""
    steps_by_stage = {kind: steps for _, kind, steps in STAGES}
    if stage not in steps_by_stage:
        raise ValueError("Unsupported relay pilot stage")
    return native_stage_input(
        old_mdin,
        stage,
        seed,
        restraint_path,
        active_mask=active_mask,
        steps=steps_by_stage[stage],
    )


def native_stage_input(
    old_mdin: str,
    stage: str,
    seed: int,
    restraint_path: Path,
    active_mask: str = MOBILE_MASK,
    steps: int | None = None,
) -> str:
    """Render one existing native stage with an explicit step count."""
    if STAGE_STEPS.get(stage) is None:
        raise ValueError("Native builder does not define the requested stage")
    text = build_mdin(old_mdin, stage, seed, active_mask, restraint_path)
    old_steps = STAGE_STEPS[stage]
    replacement_steps = old_steps if steps is None else steps
    if type(replacement_steps) is not int or replacement_steps < 1:
        raise ValueError("Explicit stage steps must be a positive integer")
    marker = f"nstlim={old_steps},"
    if text.count(marker) != 1:
        raise ValueError("Expected one native nstlim card for bounded replacement")
    return text.replace(marker, f"nstlim={replacement_steps},", 1)


def geometry_features(acceptance: Mapping[str, Any]) -> dict[str, float]:
    """Map an accepted relay geometry receipt to its three distance features."""
    if acceptance.get("native_converged") is not True:
        raise ValueError("Seed geometry must have native convergence acceptance")
    if acceptance.get("electronic_history", {}).get("passed") is not True:
        raise ValueError("Seed geometry must pass electronic-history acceptance")
    distances = acceptance.get("geometry", {}).get("distances_A", {})
    required = ("C_NL", "C_OG", "HG_OG", "HG_N", "HG_NL")
    if not all(name in distances for name in required):
        raise ValueError("Seed acceptance lacks required mapped distances")
    cn, co, oh, nh, nl_h = (float(distances[name]) for name in required)
    values = {"q_acyl": cn - co, "q_act": oh - nh, "q_del": nh - nl_h}
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("Seed CV values must be finite")
    return values


def seed_candidates(branch_directory: Path) -> list[dict[str, Any]]:
    """Read accepted optimized restart candidates from one retained branch."""
    candidates = []
    for evidence in sorted((branch_directory / "01_coordinate_drive").glob("window_*/acceptance.json")):
        acceptance = json.loads(evidence.read_text())
        values = geometry_features(acceptance)
        restart = evidence.parent / "optimized.ncrst"
        if not restart.is_file() or Path(acceptance["restart"]).resolve() != restart.resolve():
            raise FileNotFoundError(f"Missing or mismapped accepted restart: {restart}")
        candidates.append(
            {
                "seed_id": f"{branch_directory.name}_{evidence.parent.name}",
                "window": evidence.parent.name,
                "restart": restart,
                "evidence": evidence,
                "sha256": digest(restart),
                "features": values,
            }
        )
    if not candidates:
        raise ValueError(f"No accepted restart candidates under {branch_directory}")
    return candidates


def select_seed(
    candidates: Sequence[Mapping[str, Any]],
    feature_names: tuple[str, str],
    centers_A: tuple[float, float],
) -> tuple[Mapping[str, Any], float]:
    """Select the nearest accepted restart in the exact two biased features."""
    if len(feature_names) != 2 or len(centers_A) != 2 or not candidates:
        raise ValueError("Seed selection requires two CVs and nonempty candidates")
    target = np.asarray(centers_A, dtype=float)
    if not np.isfinite(target).all():
        raise ValueError("Seed target centers must be finite")
    best: tuple[float, Mapping[str, Any]] | None = None
    for candidate in candidates:
        values = candidate["features"]
        point = np.asarray([values[name] for name in feature_names], dtype=float)
        if not np.isfinite(point).all():
            raise ValueError("Seed feature values must be finite")
        distance = float(np.linalg.norm(point - target))
        if best is None or distance < best[0]:
            best = (distance, candidate)
    if best is None:
        raise ValueError("No finite seed candidate was available")
    return best[1], best[0]


def pilot_window_specs() -> list[tuple[str, tuple[str, str], tuple[float, float]]]:
    """Return the immutable 18 target definitions approved for input prep."""
    q_acyl = (-0.45, -0.25, -0.05)
    q_act = (0.90, 1.10, 1.30)
    q_del = (-1.05, -0.90, -0.75)
    specs = []
    for chart, names, x_centers, y_centers in (
        ("activation", ("q_acyl", "q_act"), q_acyl, q_act),
        ("delivery", ("q_acyl", "q_del"), q_acyl, q_del),
    ):
        index = 0
        for y_center in y_centers:
            for x_center in x_centers:
                index += 1
                specs.append(
                    (f"{chart}_{index:02d}", names, (x_center, y_center))
                )
    return specs


def pilot_window_records(
    base: Path,
) -> list[PilotWindow]:
    """Build both charts and resolve every target to an accepted branch seed."""
    branches = {
        "activation": seed_candidates(base / "01_activation_reverse"),
        "delivery": seed_candidates(base / "02_acylation_forward"),
    }
    records = []
    for label, names, centers in pilot_window_specs():
        chart = label.split("_", maxsplit=1)[0]
        seed, distance = select_seed(branches[chart], names, centers)
        records.append(
            PilotWindow(
                label=label,
                chart=chart,
                bias=RelayBiasWindow(
                    feature_names=names,
                    centers_A=centers,
                    force_constants=(PILOT_FORCE_CONSTANT, PILOT_FORCE_CONSTANT),
                    outer_extent_A=PILOT_OUTER_EXTENT_A,
                ),
                seed_source=seed["restart"],
                seed_evidence=seed["evidence"],
                seed_sha256=seed["sha256"],
                seed_id=seed["seed_id"],
                seed_window=seed["window"],
                seed_distance_A=distance,
                seed_offset_A=tuple(
                    float(seed["features"][name] - center)
                    for name, center in zip(names, centers, strict=True)
                ),
            )
        )
    return records


def write_json(path: Path, value: Any) -> None:
    """Write deterministic finite JSON with one trailing newline."""
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
