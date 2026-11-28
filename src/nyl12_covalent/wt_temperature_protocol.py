"""Bounded full-system WT E-S thermal equilibration after accepted 300 K MD."""

from __future__ import annotations

import math
import re
from typing import Any

STAGES = (
    ("01_ramp_nvt", 100_000),
    ("02_density_npt", 500_000),
    ("03_check_npt", 1_000_000),
)
TARGET_K = 348.15


def mdin_files(seed: int) -> dict[str, str]:
    """Return 200 ps NVT ramp, 1 ns NPT and 2 ns unrestrained NPT inputs.

    Args:
        seed: First of three distinct positive signed-32-bit Langevin seeds.

    Returns:
        Native input texts preserving the accepted WT PME/SHAKE protocol.
    """
    if type(seed) is not int or not 0 < seed <= 2147483645:
        raise ValueError("Three valid consecutive explicit seeds required")
    inputs = {}
    for index, (name, steps) in enumerate(STAGES):
        ramp = index == 0
        text = f"""WT E-S full-system thermal equilibration: {name}
&cntrl
 imin=0, ntx=5, irest=1,
 nstlim={steps}, dt=0.002, ntc=2, ntf=2,
 ntt=3, gamma_ln=1.0, temp0=348.15, ig={seed + index},
 ntb={1 if ramp else 2}, ntp={0 if ramp else 1}, cut=10.0,
 ntpr=1000, ntwx=1000, ntwr=50000,
 ioutfm=1, ntxo=2, iwrap=0,
 ntr=0, ibelly=0, ifqnt=0, nmropt={int(ramp)},
"""
        if not ramp:
            text += " barostat=2, pres0=1.0, taup=2.0,\n"
        text += "/\n"
        if ramp:
            text += "&wt type='TEMP0', istep1=0, istep2=100000, value1=300.0, value2=348.15 /\n&wt type='END' /\n"
        validate_input(text, name, seed + index)
        inputs[name] = text
    return inputs


def validate_input(text: str, stage: str, seed: int) -> None:
    """Reject unintended restraints, frozen regions, QM and wrong stage controls."""
    controls = {
        key.lower(): value.strip(" '\"")
        for key, value in re.findall(r"(\w+)\s*=\s*([^,/\n]+)", text)
    }
    steps = dict(STAGES)[stage]
    expected = {
        "irest": "1",
        "ntx": "5",
        "ntr": "0",
        "ibelly": "0",
        "ifqnt": "0",
        "dt": "0.002",
        "ntc": "2",
        "ntf": "2",
        "nstlim": str(steps),
        "ntt": "3",
        "gamma_ln": "1.0",
        "cut": "10.0",
        "temp0": "348.15",
        "ig": str(seed),
        "ntb": "1" if stage == STAGES[0][0] else "2",
        "ntp": "0" if stage == STAGES[0][0] else "1",
    }
    if any(controls.get(key) != value for key, value in expected.items()):
        raise ValueError("Thermal input differs from reviewed full-system protocol")
    if any(
        word in text.lower()
        for word in (
            "restraint",
            "disang",
            "&rst",
            "bellymask",
            "qmmask",
            "&qmmm",
            "&ewald",
        )
    ):
        raise ValueError("Unexpected restraint/selection/extra Hamiltonian block")
    if stage == STAGES[0][0]:
        if (
            controls.get("nmropt") != "1"
            or controls.get("value1") != "300.0"
            or controls.get("value2") != "348.15"
            or controls.get("istep2") != "100000"
        ):
            raise ValueError("Incorrect gradual thermal ramp")
    elif (
        controls.get("nmropt") != "0"
        or controls.get("barostat") != "2"
        or controls.get("pres0") != "1.0"
        or controls.get("taup") != "2.0"
    ):
        raise ValueError("Incorrect unchanged NPT controls")


def audit_output(text: str, steps: int) -> dict[str, Any]:
    """Check native completion and thermodynamic finiteness, without claiming equilibration."""
    if re.search(r"\b(?:NaN|Infinity)\b|FATAL|ERROR:|SHAKE cannot", text, re.I):
        raise ValueError("Native numeric or engine failure")
    rows = re.findall(
        r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([-+\d.EeDd]+)\s+TEMP\(K\)\s*=\s*([-+\d.EeDd]+)",
        text.split("A V E R A G E S", 1)[0],
    )
    energies = re.findall(r"EPtot\s*=\s*([-+\d.EeDd]+)", text)

    def number(value: str) -> float:
        return float(value.replace("D", "E").replace("d", "e"))

    if not rows or int(rows[-1][0]) != steps or not energies:
        raise ValueError("Missing exact native completion or potential energies")
    if not all(math.isfinite(number(x)) for row in rows for x in row[1:]) or not all(
        math.isfinite(number(e)) for e in energies
    ):
        raise ValueError("Nonfinite native thermodynamics")
    final_temperature = number(rows[-1][2])
    if not TARGET_K - 25 <= final_temperature <= TARGET_K + 25:
        raise ValueError("Final temperature outside target +/-25 K process gate")
    return {
        "completed_steps": steps,
        "final_time_ps": number(rows[-1][1]),
        "final_temperature_K": final_temperature,
        "finite_native_process_gate": True,
        "thermal_equilibration_accepted": False,
        "sampling_converged": False,
        "review_required": "Temperature/density/volume/block stability and protein/substrate geometry require separate review",
    }
