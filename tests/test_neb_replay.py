"""Same-cycle replay comparison separates precision agreement from science gates."""

from pathlib import Path

import pytest

from nyl12_covalent.neb_replay import (
    compare_replay,
    parse_replay_profiles,
    printed_quantum,
)
from scripts.verify_neb_cg_replay import watch_replay


def block(
    cycle: int,
    *,
    energy_shift: float = 0.0,
    maximum: str = "1.2345E-03",
    rms: str = "1.1111E-05",
) -> bytes:
    """Return one complete24image native report and matching convergence block."""
    lines = ["NEB Report\n"]
    lines += [
        f"Img {image:4d}{-717.0 + image * 0.001 + energy_shift:15.7f}\n"
        for image in range(1, 25)
    ]
    lines += [
        f"Total path length 1.0\nTesting convergence in cycle {cycle}\n",
        f" Max grad {maximum} Target: 1.0E-4 converged? no\n",
        f" RMS grad {rms} Target: 6.0E-5 converged? no\n",
    ]
    return "".join(lines).encode()


def test_parser_associates_samecycle_gradients_and_ignores_partial_tail() -> None:
    payload = block(1) + block(2) + block(3).split(b" RMS grad")[0]
    parsed = parse_replay_profiles(payload)
    assert list(parsed) == [1, 2]
    assert len(parsed[1].energies) == 24
    assert parsed[1].max_gradient == pytest.approx(0.0012345)
    assert parsed[1].end_byte < len(payload)
    with pytest.raises(ValueError, match="Duplicate"):
        parse_replay_profiles(block(1) + block(1))


def test_printedprecision_allows_one_energy_quantum_and_matching_gradient_rounding() -> (
    None
):
    pilot = parse_replay_profiles(block(1))
    new = parse_replay_profiles(
        block(1, energy_shift=1e-7, maximum="1.2346E-03", rms="1.1112E-05")
    )
    report = compare_replay(pilot, new)
    assert report["initial_24image_profile_verified"]
    assert report["replay_matched_pilot"] and not report["replay_passed_pilot"]
    assert report["scientific_success_claim"] is False
    assert printed_quantum("1.2345D-03") == pytest.approx(1e-7)
    assert printed_quantum("1.234") == pytest.approx(0.001)


@pytest.mark.parametrize("defect", ["energy", "max", "rms"])
def test_replay_divergence_is_explicit(defect: str) -> None:
    options = (
        {"energy_shift": 1e-5}
        if defect == "energy"
        else ({"maximum": "1.3345E-03"} if defect == "max" else {"rms": "1.2111E-05"})
    )
    result = compare_replay(
        parse_replay_profiles(block(1)), parse_replay_profiles(block(1, **options))
    )
    assert result["status"] == "diverged" and result["failed_cycles"] == [1]


def test_missing_intermediate_pilot_cycle_fails() -> None:
    with pytest.raises(ValueError, match="consecutive"):
        compare_replay(parse_replay_profiles(block(1) + block(3)), {})


def test_watcher_retains_pilot_and_checkpoint79_and_flags_only_replay_progress(
    tmp_path: Path,
) -> None:
    pilot = tmp_path / "pilot.mdout"
    replay = tmp_path / "replay.mdout"
    pilot_bytes = b"".join(block(i) for i in range(1, 80))
    pilot.write_bytes(pilot_bytes)
    replay.write_bytes(pilot_bytes + block(80))
    output = tmp_path / "verification"
    report = watch_replay(pilot, replay, output)
    assert report["replay_passed_pilot"] and report["compared_cycles"] == 79
    assert report["replay_last_complete_cycle"] == 80
    assert (output / "pilot_reference.mdout").read_bytes() == pilot_bytes
    for cycle in (1, 50, 79):
        assert (output / f"checkpoint_cycle_{cycle:03d}/comparison.json").is_file()
    assert report["stationarity_or_TS_verified"] is False


def test_missing_replay_is_pending_not_a_pass(tmp_path: Path) -> None:
    pilot = tmp_path / "pilot.mdout"
    pilot.write_bytes(block(1))
    result = watch_replay(pilot, tmp_path / "not_started.mdout", tmp_path / "out")
    assert result["status"] == "pending_initial_profile"
    assert not result["initial_24image_profile_verified"]
    with pytest.raises(ValueError):
        watch_replay(pilot, tmp_path / "none", tmp_path / "bad", interval=0)
