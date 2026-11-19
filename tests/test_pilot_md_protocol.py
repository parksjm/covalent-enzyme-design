"""Scientific controls and process-failure checks for matched MD pilots."""

from pathlib import Path

import pytest

from nyl12_covalent.pilot_md_protocol import mdin_files, protocol_manifest
from scripts.run_model30_md_replica import stage_acceptance


def test_thermal_ramp_enabled_and_production_unrestrained() -> None:
    """Require TEMP0 weights to be active and remove all production restraints."""
    files = mdin_files(1234)
    heating = files["01_heat_nvt"]
    assert "nmropt=1" in heating
    assert "type='TEMP0'" in heating
    assert "value1=0.0, value2=300.0" in heating
    production = files["05_production_10ns"]
    assert "nmropt=0" in production and "ntr=0" in production
    assert "restraintmask" not in production
    assert "nstlim=5000000" in production
    assert "ntc=2, ntf=2" in production
    assert "barostat=2" in production


def test_matched_duration_and_independent_seed() -> None:
    """Equal stage lengths preserve comparisons while seeds change independently."""
    one, two = protocol_manifest(1234), protocol_manifest(4321)
    assert [s["duration_ps"] for s in one["stages"]] == [100, 100, 100, 200, 10000]
    assert [s["heavy_restraint_kcal_mol_A2"] for s in one["stages"]] == [
        5,
        2,
        0.5,
        None,
        None,
    ]
    assert set(s["seed"] for s in one["stages"]).isdisjoint(
        s["seed"] for s in two["stages"]
    )
    with pytest.raises(ValueError, match="positive"):
        mdin_files(0)


def test_stage_completion_and_temperature(tmp_path: Path) -> None:
    """A returned engine must have reached the intended step and temperature."""
    output = tmp_path / "stage.out"
    output.write_text(
        "NSTEP = 50000 TIME(PS) = 100.000 TEMP(K) = 299.34\n"
        "A V E R A G E S O V E R 50 S T E P S\n"
        "NSTEP = 50000 TIME(PS) = 100.000 TEMP(K) = 149.75\n"
        "R M S F L U C T U A T I O N S\n"
        "NSTEP = 50000 TIME(PS) = 100.000 TEMP(K) = 85.76\n"
    )
    assert stage_acceptance(output, 50000)["final_temperature_K"] == 299.34
    with pytest.raises(ValueError, match="Missing final step"):
        stage_acceptance(output, 100000)
    output.write_text("NSTEP = 50000 TIME(PS) = 100.000 TEMP(K) = 650.0\n")
    with pytest.raises(ValueError, match="temperature"):
        stage_acceptance(output, 50000)
    output.write_text("NSTEP = 50000 TIME(PS) = 100.000 TEMP(K) = 300.0\n Etot=NaN\n")
    with pytest.raises(ValueError, match="Numeric/engine"):
        stage_acceptance(output, 50000)


def test_90ns_extension_preserves_restart_and_unrestrained_controls() -> None:
    """Extend sampling from a checkpoint while retaining the matched Hamiltonian."""
    from scripts.continue_model30_md import continuation_input

    source = mdin_files(731110)["05_production_10ns"]
    extension = continuation_input(source, 781110)
    assert "nstlim=45000000" in extension
    assert "ig=781110" in extension
    assert "irest=1" in extension and "ntx=5" in extension
    assert "ntr=0" in extension and "nmropt=0" in extension
    assert "dt=0.002" in extension
    with pytest.raises(ValueError, match="ntr"):
        continuation_input(source.replace("ntr=0", "ntr=1"), 781110)
    with pytest.raises(ValueError, match="irest"):
        continuation_input(source.replace("irest=1", "irest=0"), 781110)
