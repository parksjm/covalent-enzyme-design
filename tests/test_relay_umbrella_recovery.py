"""Reviewed heat recovery preserves physical inputs and original lineage."""

from __future__ import annotations
import json
from pathlib import Path
import pytest
from nyl12_covalent.amber_md_records import instantaneous_md_text
from nyl12_covalent.neb_umbrella_campaign import digest
from nyl12_covalent.relay_umbrella_recovery import heat_adoption


@pytest.mark.parametrize(
    "delimiter",
    [
        "A V E R A G E S OVER 2000 STEPS",
        "AVERAGES OVER 2000 STEPS",
        "R M S F L U C T U A T I O N S",
    ],
)
def test_instantaneous_parser_excludes_statistical_temperatures(delimiter: str) -> None:
    text = "NSTEP=2000 TEMP(K)=300\n" + delimiter + "\nNSTEP=2000 TEMP(K)=50"
    assert "TEMP(K)=50" not in instantaneous_md_text(text)
    assert "TEMP(K)=300" in instantaneous_md_text(text)


def test_reviewed_heat_adoption_requires_exact_native_and_old_plan(
    tmp_path: Path,
) -> None:
    old = tmp_path / "oldplan.json"
    old.write_text("original")
    stage = tmp_path / "heat"
    stage.mkdir()
    native = stage / "restart.ncrst"
    native.write_text("originalheatcheckpoint")
    receipt = tmp_path / "heat_recovery.json"
    record = dict(
        window="window",
        stage=str(stage),
        temperature_acceptance_role="completed_nonstationary_heating_ramp",
        native_source_sha256={str(native): digest(native)},
        original_input_plan=str(old),
        original_input_plan_sha256=digest(old),
    )
    receipt.write_text(json.dumps(record))
    plan = dict(
        adopted_heat=dict(window=dict(receipt=str(receipt), sha256=digest(receipt)))
    )
    assert heat_adoption(plan, "window", stage) == (receipt, old)
    assert heat_adoption(plan, "other", stage) is None
    native.write_text("changed")
    with pytest.raises(ValueError, match="native source changed"):
        heat_adoption(plan, "window", stage)
