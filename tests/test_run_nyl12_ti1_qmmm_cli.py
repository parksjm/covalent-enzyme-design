from __future__ import annotations

from pathlib import Path

import pytest

from scripts import run_nyl12_ti1_qmmm as cli


def test_cli_forwards_explicit_inputs_and_defaults(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    received: dict[str, object] = {}

    def fake_run_attempt(**kwargs: object) -> Path:
        received.update(kwargs)
        return Path("attempt")

    monkeypatch.setattr(cli, "run_attempt", fake_run_attempt)

    result = cli.main(
        [
            "--topology",
            "system.prmtop",
            "--coordinates",
            "system.rst7",
            "--reference",
            "reference.rst7",
            "--output",
            "attempt",
            "--amber-bin",
            "/opt/amber/bin",
            "--target",
            "pa6_l2",
        ]
    )

    assert result == 0
    assert received == {
        "topology": Path("system.prmtop"),
        "coordinates": Path("system.rst7"),
        "reference": Path("reference.rst7"),
        "output": Path("attempt"),
        "amber_bin": Path("/opt/amber/bin"),
        "target": "pa6_l2",
        "max_cycles": 300,
        "checkpoint_attempt": None,
    }


def test_cli_forwards_optional_cycle_limit_and_checkpoint(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    received: dict[str, object] = {}
    monkeypatch.setattr(
        cli, "run_attempt", lambda **kwargs: received.update(kwargs) or Path("attempt")
    )

    cli.main(
        [
            "--topology",
            "system.prmtop",
            "--coordinates",
            "system.rst7",
            "--reference",
            "reference.rst7",
            "--output",
            "attempt",
            "--amber-bin",
            "/opt/amber/bin",
            "--target",
            "poly_xama_maxa",
            "--max-cycles",
            "125",
            "--checkpoint-attempt",
            "previous-attempt",
        ]
    )

    assert received["target"] == "poly_xama_maxa"
    assert received["max_cycles"] == 125
    assert received["checkpoint_attempt"] == Path("previous-attempt")


def test_cli_requires_all_input_and_output_paths() -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main([])

    assert exc_info.value.code == 2


def test_cli_rejects_unknown_target() -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(
            [
                "--topology",
                "system.prmtop",
                "--coordinates",
                "system.rst7",
                "--reference",
                "reference.rst7",
                "--output",
                "attempt",
                "--amber-bin",
                "/opt/amber/bin",
                "--target",
                "unknown",
            ]
        )

    assert exc_info.value.code == 2
