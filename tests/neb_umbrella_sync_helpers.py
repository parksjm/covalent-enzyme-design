"""Shared temporary scientific payload fixtures for compact synchronization tests."""

from pathlib import Path
import numpy as np
import pytest
from nyl12_covalent import neb_umbrella_sync as s
from nyl12_covalent.neb_umbrella_dispatch import partition

STAGE = "03_production_block_01_5ps"
WINDOW = "window_017"


def hashes(folder: Path, names: tuple = s.CRITICAL) -> dict:
    """Hash a compact test payload."""
    return {name: s.digest(folder / name) for name in names}


def block(folder: Path, cfg: dict) -> dict:
    """Write one valid small production result without native MD."""
    (folder / "observables").mkdir(parents=True)
    csv = folder / "observables/observables.csv"
    csv.write_text("q_acyl,q_h,cn_progress,cn_nalpha_hpool\n1,2,3,4\n2,3,4,5\n")
    s.write(
        folder / "observables/observables.json",
        {"frames": 2, "frames_censored": 0, "csv_sha256": s.digest(csv)},
    )
    s.write(
        folder / "acceptance.json",
        {"finite_energy_temperature": True, "completed_steps": 10000},
    )
    s.write(folder / "process.json", {"exit_status": 0, "ended_utc": "now"})
    shape = (len(cfg["windows"]), 2)
    np.savez_compressed(
        folder / "samples_and_cross_bias.npz",
        cv_values=[[1.0, 2.0], [2.0, 3.0]],
        bias_kcal_mol=np.zeros(shape),
        u_kn_bias_only=np.zeros(shape),
    )
    s.write(
        folder / "analysis_request.json",
        {
            "samples": 2,
            "bias_matrix_shape": list(shape),
            "all_window_bias_definitions": [w["bias"] for w in cfg["windows"]],
        },
    )
    return hashes(folder)


def view(folder: Path) -> dict:
    """Write a manifest-complete small endpoint bundle."""
    folder.mkdir(parents=True)
    for name in ("view.pdb", "view_umbrella.vmd", "manifest.json"):
        (folder / name).write_text(name)
    payload = hashes(folder, ("view.pdb", "view_umbrella.vmd", "manifest.json"))
    (folder / "MANIFEST.sha256").write_text(
        "".join(f"{v}  {k}\n" for k, v in payload.items())
    )
    return hashes(folder, (*payload, "MANIFEST.sha256"))


@pytest.fixture
def campaign(tmp_path: Path) -> tuple:
    """Build a two-host plan with all 63 windows declared once."""
    root = tmp_path / "campaign"
    root.mkdir()
    labels = [f"window_{i:03d}" for i in range(1, 64)]
    cfg = {
        "windows": [
            {
                "label": label,
                "bias": {"q": i},
                "stages": [{"label": STAGE, "kind": "production"}],
            }
            for i, label in enumerate(labels)
        ]
    }
    prepared = root / "prepared.json"
    s.write(prepared, cfg)
    (root / "prepared.sha256").write_text(s.digest(prepared))
    groups = partition(labels, labels[:8])
    plan = {
        "local_host": "violet",
        "remote_host": "azalea",
        "windows": labels,
        "assignments": {"violet": groups["local"], "azalea": groups["remote"]},
        "hosts": {
            h: {"prepared": str(prepared), "sha256": s.digest(prepared)}
            for h in ("violet", "azalea")
        },
    }
    folder = root / "dispatch"
    folder.mkdir()
    path = folder / "plan.json"
    s.write(path, plan)
    (folder / "plan.sha256").write_text(s.digest(path))
    incoming = tmp_path / "remote_payload"
    values = block(incoming, cfg)
    expected = s.identity(path, "azalea", str(root / WINDOW / STAGE), values)
    return path, cfg, incoming, values, expected
