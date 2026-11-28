"""Validate accepted native endpoints and launch the authorized100ps prefix queue."""
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
from parmed.amber import Rst7

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
PARENT = BASE / "continuation_60ps_20261126/plan.json"
TARGET = BASE / "continuation_100ps_20261127"
sys.path.insert(0, str(PARENT.parent / "source_snapshot/src"))
from nyl12_covalent.relay_umbrella_runtime import verify_plan, native_acceptance
from relay_umbrella_continuation import prepare_continuation


def main() -> None:
    """Check frozen inputs/endpoints and launch unchanged accepted-prefix workers."""
    logging.basicConfig(level=logging.INFO)
    processes = subprocess.run(["pgrep", "-x", "sander"], capture_output=True, text=True)
    if processes.returncode != 1:
        raise RuntimeError("Sander occupancy changed; no launch")
    plan, cfg, base = verify_plan(PARENT)
    rows = []
    for row in cfg["windows"]:
        stage = base / row["label"] / "14_production_5ps"
        accepted = native_acceptance(stage, 10000, plan["temperature_gate_K"])
        rst = Rst7.open(str(stage / "restart.ncrst"))
        for name in ["coordinates", "velocities", "box"]:
            value = np.asarray(getattr(rst, name))
            if not value.size or not np.isfinite(value).all():
                raise ValueError("Invalid restart field: " + name)
        if rst.natom != 104922 or np.asarray(rst.velocities).shape != (104922, 3):
            raise ValueError("Unexpected atom/velocity inventory")
        rows.append(dict(window=row["label"], native=accepted, natom=rst.natom,
                         restart_time_ps=float(rst.time), velocities_retained=True,
                         finite_coordinates_velocities_box=True))
        logging.info("Accepted endpoint %s", row["label"])
    receipt = dict(status="accepted60ps_endpoints_verified", checked_utc=datetime.now(timezone.utc).isoformat(),
                   windows=rows, no_new_equilibration=True, unchanged_topology=True)
    (HERE / "native_prefix_acceptance.json").write_text(json.dumps(receipt, indent=2) + "\n")
    prepare_continuation(PARENT, TARGET, target_blocks=20, workers=57)
    script = TARGET / "source_snapshot/scripts/run_relay_umbrella.py"
    subprocess.run([plan["python"], str(script), "launch", "--plan", str(TARGET / "plan.json")], check=True)
    logging.info("Launched57 persistent native workers toward100ps")


if __name__ == "__main__":
    main()
