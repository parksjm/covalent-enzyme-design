"""Immutable preparation of dense two-CV umbrellas from an audited old NEB band."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import shutil
from typing import Any

import numpy as np

from .neb_umbrella_bias import HarmonicWindow, amber_restraints
from .neb_umbrella_inputs import build_mdin

ROOT = Path(__file__).resolve().parents[2]
STAGES = [
    ("01_heating_1ps", "heat"),
    ("02_equilibration_5ps", "equil"),
    *[(f"03_production_block_{i:02d}_5ps", "production") for i in range(1, 7)],
]


def digest(path: Path) -> str:
    """Return streaming SHA-256 of an immutable input/source."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            result.update(block)
    return result.hexdigest()


def write_json(path: Path, value: Any) -> None:
    """Atomically write stage/status JSON with finite numeric values."""
    pending = path.with_name(path.name + ".pending")
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def dense_centers(
    rows: list[dict[str, Any]], maximum_step: float
) -> list[dict[str, Any]]:
    """Retain every original center and linearly subdivide only the two CV targets.

    Args:
        rows: Audited image records with q_acyl,q_h,image,full_restart and sha256.
        maximum_step: Positive maximum Euclidean target spacing in Angstrom.

    Returns:
        Ordered targets with nearest ORIGINAL seed, ancestry and seed-CV distance.
    """
    if len(rows) < 2 or not math.isfinite(maximum_step) or maximum_step <= 0:
        raise ValueError("Need multiple auditedcenters andpositivefiniteCVstep")
    original = np.array([[r["q_acyl"], r["q_h"]] for r in rows], dtype=float)
    if not np.isfinite(original).all():
        raise ValueError("CVcentersmustbefinite")
    values = [(original[0], rows[0]["image"])]
    for index, (a, b) in enumerate(zip(original, original[1:])):
        parts = max(1, int(np.ceil(np.linalg.norm(b - a) / maximum_step)))
        for j in range(1, parts + 1):
            values.append(
                (
                    b.copy() if j == parts else a + (b - a) * j / parts,
                    rows[index + 1]["image"] if j == parts else None,
                )
            )
    result = []
    for index, (center, origin) in enumerate(values):
        lengths = np.linalg.norm(original - center, axis=1)
        nearest = int(np.argmin(lengths))
        result.append(
            {
                "label": f"window_{index + 1:03d}",
                "q_acyl_center": float(center[0]),
                "q_h_center": float(center[1]),
                "original_neb_image": origin,
                "seed_image": rows[nearest]["image"],
                "source_restart": rows[nearest]["full_restart"],
                "source_restart_sha256": rows[nearest]["sha256"],
                "seed_cv_distance_A": float(lengths[nearest]),
            }
        )
    return result


def priority_order(windows: list[dict[str, Any]], peak_image: int) -> list[str]:
    """Prioritize endpoints, the source peak and its neighbors, then spread in CV space."""
    points = np.array([[r["q_acyl_center"], r["q_h_center"]] for r in windows])
    peak = next(
        i for i, row in enumerate(windows) if row["original_neb_image"] == peak_image
    )
    chosen = list(
        dict.fromkeys(
            [
                0,
                len(windows) - 1,
                peak,
                max(0, peak - 1),
                min(len(windows) - 1, peak + 1),
            ]
        )
    )
    while len(chosen) < len(windows):
        distance = np.min(
            np.linalg.norm(points[:, None] - points[chosen], axis=2), axis=1
        )
        distance[chosen] = -1
        chosen.append(int(np.argmax(distance)))
    return [windows[i]["label"] for i in chosen]


def source_paths() -> list[Path]:
    """Return the complete small runtime source tree retained with this campaign."""
    names = [
        "scripts/run_neb_umbrella.py",
        "scripts/run_qmmm_umbrella_campaign.py",
        "src/nyl12_covalent/__init__.py",
        "src/nyl12_covalent/jobs.py",
        "src/nyl12_covalent/qmmm_umbrella.py",
        "src/nyl12_covalent/restrained_qmmm_md.py",
        "src/nyl12_covalent/relay_progress_vmd.py",
        *[
            f"src/nyl12_covalent/neb_umbrella_{name}.py"
            for name in (
                "campaign",
                "runtime",
                "diagnostics",
                "observables",
                "bias",
                "inputs",
            )
        ],
    ]
    return [ROOT / name for name in names]


def prepare(config_path: Path) -> Path:
    """Freeze audited full seeds and native stages; never launch MD.

    Args:
        config_path: Explicit audit/output/Sander/Python/seed/force/temperature config.

    Returns:
        Prepared campaign config path after all native inputs have been rendered.
    """
    cfg = json.loads(config_path.read_text())
    audit_path = Path(cfg["seed_audit"]).resolve(strict=True)
    audit = json.loads(audit_path.read_text())
    if (
        audit["passed"] is not True
        or len(audit["windows"]) != 24
        or audit["active_atoms"] != 4165
        or audit["qm_atoms"] != 98
        or audit["qm_charge"] != -2
    ):
        raise ValueError("Requirepassedold24image/4165active/98QM/charge−2seed audit")
    if type(cfg["workers"]) is not int or not 1 <= cfg["workers"] <= 8:
        raise ValueError("One toeightCPUworkersrequired")
    if type(cfg["seed_base"]) is not int or cfg["seed_base"] <= 0:
        raise ValueError("Positiveexplicitseedbase required")
    if cfg["temperature_K"] != 300 or cfg["temperature_gate_K"] != [250, 350]:
        raise ValueError("Explicitreviewed300K/250–350gate required")
    if cfg["heating_temperature_gate_K"] != [25, 350]:
        raise ValueError("Explicit reviewed25–350K heating gate required")
    output = Path(cfg["output"]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "prepared.json").exists() or (output / "inputs").exists():
        raise FileExistsError("Preparedcampaignalreadyexists")
    windows = dense_centers(audit["windows"], cfg["maximum_cv_step_A"])
    roles = dict(audit["roles_1based"])
    for old, new in (("og1", "thr_og"), ("hg1", "proton"), ("nalpha", "thr_n")):
        if old in roles:
            roles[new] = roles[old]
    sources = [
        config_path.resolve(),
        audit_path,
        Path(audit["topology"]),
        Path(audit["source_input"]),
        Path(cfg["sander"]),
        Path(cfg["python"]),
        *[Path(r["full_restart"]) for r in audit["windows"]],
        *source_paths(),
    ]
    hashes = {str(p): digest(p) for p in sources}
    for row in audit["windows"]:
        if digest(Path(row["full_restart"])) != row["sha256"]:
            raise ValueError("AuditedNEBseedchanged")
    for name, expected in audit["source_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Audited source changed: {name}")
    inputs = output / "inputs"
    inputs.mkdir()
    seeds = inputs / "seeds"
    seeds.mkdir()
    topology = inputs / "system.prmtop"
    shutil.copy2(audit["topology"], topology)
    old = Path(audit["source_input"]).read_text()
    (inputs / "source_neb.mdin").write_text(old)
    seed_paths = {}
    for row in audit["windows"]:
        target = seeds / f"image_{row['image']:03d}.ncrst"
        shutil.copy2(row["full_restart"], target)
        seed_paths[row["image"]] = target
    frozen = output / "source_snapshot"
    for source in source_paths():
        destination = frozen / source.relative_to(ROOT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    cfg["frozen_source_sha256"] = {str(p): digest(p) for p in frozen.rglob("*.py")}
    used = set()
    for i, row in enumerate(windows):
        window = HarmonicWindow(
            row["q_acyl_center"], row["q_h_center"], cfg["k_acyl"], cfg["k_h"]
        )
        row["bias"] = asdict(window)
        row["seed"] = str(seed_paths[row["seed_image"]])
        row["stages"] = []
        folder = output / row["label"]
        folder.mkdir()
        for j, (label, kind) in enumerate(STAGES):
            seed = cfg["seed_base"] + i * 100 + j
            if seed in used or seed > 2147483647:
                raise ValueError("Duplicate/outofrangeLangevinseed")
            used.add(seed)
            stage = folder / label
            stage.mkdir()
            restraint = stage / "umbrella.RST"
            restraint.write_text(amber_restraints(window, roles))
            (stage / "run.mdin").write_text(
                build_mdin(old, kind, seed, audit["active_mask"], Path("umbrella.RST"))
            )
            row["stages"].append({"label": label, "kind": kind, "seed": seed})
    cfg.update(
        topology=str(topology),
        roles_1based=roles,
        active_mask=audit["active_mask"],
        windows=windows,
        source_sha256=hashes,
        runtime_script=str(frozen / "scripts/run_neb_umbrella.py"),
        priority_schedule=priority_order(windows, cfg["source_peak_image"]),
        execution_sha256={
            str(Path(cfg[k])): digest(Path(cfg[k])) for k in ("sander", "python")
        },
        active_atoms=4165,
        qm_atoms=98,
        qm_charge=-2,
        coordinate_interpolation="CVtargets only;fullseedsalwaysnearestvalidoriginalNEBimage",
        bias_convention="Amberrk:centralE=k*delta²;nativeouterlineartails;nohalf",
        analysis_policy="2DMBARpendingcorrelation/overlapreview;allvalidsamplesretainedincludingprotonbranches",
    )
    cfg["prepared_input_sha256"] = {
        str(p): digest(p)
        for p in [
            topology,
            *seeds.iterdir(),
            *output.glob("window_*/*/run.mdin"),
            *output.glob("window_*/*/umbrella.RST"),
        ]
    }
    prepared = output / "prepared.json"
    write_json(prepared, cfg)
    (output / "prepared.sha256").write_text(digest(prepared) + "\n")
    write_json(
        output / "preflight.json",
        {
            "prepared": True,
            "window_count": len(windows),
            "original_neb_centers_retained": 24,
            "workers": cfg["workers"],
            "unique_stage_seeds": len(used),
            "minimum_production_ps_per_window": 30,
            "launch_requires_root_review": True,
            "density_review_required_over_100_windows": len(windows) > 100,
        },
    )
    return prepared


def frozen(path: Path) -> dict[str, Any]:
    """Reject changed inputs/runtime sources before any native stage is started."""
    if digest(path) != (path.parent / "prepared.sha256").read_text().strip():
        raise ValueError("Preparedcampaignconfigchanged")
    cfg = json.loads(path.read_text())
    for name, expected in {
        **cfg["frozen_source_sha256"],
        **cfg["execution_sha256"],
        **cfg["prepared_input_sha256"],
    }.items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Frozenumbrella source/inputchanged:{name}")
    return cfg
