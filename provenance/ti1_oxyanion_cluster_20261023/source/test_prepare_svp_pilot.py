"""Check pilot inputs, physical constraints and fail-closed preparation."""

from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent


@pytest.mark.parametrize("bad_seed", [False, True])
def test_stages_only_accepted_seed(tmp_path: Path, bad_seed: bool) -> None:
    """Reject invalid seeds and keep all reactive bonds and hydrogen coordinates free."""
    shutil.copytree(
        BASE / "seed_01", tmp_path / "seed_01", ignore=shutil.ignore_patterns("source")
    )
    for name in (
        "assess_pilot.py",
        "prepare_oxyanion_cluster.py",
        "donor_fragments.py",
    ):
        shutil.copy2(BASE / name, tmp_path / name)
    spec = importlib.util.spec_from_file_location(
        "svp_preparer", BASE / "prepare_svp_pilot.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.BASE = tmp_path
    report_path = tmp_path / "seed_01/preflight.json"
    if bad_seed:
        report = json.loads(report_path.read_text())
        report["geometry_screen_passed"] = False
        report_path.write_text(json.dumps(report))
        with pytest.raises(ValueError, match="Unacceptable donor seed"):
            module.main()
        assert not (tmp_path / "svp_pilot_01").exists()
        return
    module.main()
    out = tmp_path / "svp_pilot_01"
    policy = json.loads((out / "pilot_policy.json").read_text())
    atoms = json.loads((out / "atom_map.json").read_text())["atoms"]
    frozen = policy["frozen_atom_indices_0based"]
    assert len(frozen) == 21
    assert all(row["index"] not in frozen for row in atoms if row["element"] == "H")
    mobile_heavy = {
        row["origin_id"]
        for row in atoms
        if row["element"] != "H" and row["index"] not in frozen
    }
    assert mobile_heavy == {"THR:OG1", "PA6:C54", "PA6:O49", "PA6:N59", "WAT:O"}
    text = (out / "input.gjf").read_text()
    assert "wB97XD/def2SVP" in text
    assert "EmpiricalDispersion" not in text
    assert text.count("(Frozen)") == 21 * 3 + 1
    assert "MaxCycles=20" in text
    assert policy["reactive_bond_restraints"] is False
    assert policy["all_hydrogens_mobile"] is True
    assert (
        policy["prospective_checks"]["ammonium_proton_limits"]["N_H_maximum_A"] == 1.25
    )
    assert (out / "native_seed.xyz").read_bytes() == (
        BASE / "seed_01/native_seed.xyz"
    ).read_bytes()
    assert len(list((out / "source").iterdir())) == 10
    with pytest.raises(FileExistsError):
        module.main()
