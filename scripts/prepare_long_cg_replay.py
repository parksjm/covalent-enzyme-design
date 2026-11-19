#!/usr/bin/env python3
"""Prepare an identical-start CG replay with only two native budget changes."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import re
import shlex
import shutil

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = (
    ROOT
    / "qmmm/pa6_flipped_model30/07_ash274_comparison/09_water_relay_acylation/04_neb/spacing_tests/04_uniform_24_cg"
)
LOGGER = logging.getLogger(__name__)


def digest(path: Path) -> str:
    """Return an input/source SHA-256 without loading large binaries into memory."""
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def replace_budget(text: str, namelist: str, key: str, old: int, new: int) -> str:
    """Replace exactly one declared integer budget while preserving other bytes.

    Args:
        text: Complete native Amber input.
        namelist: Required named block, excluding ampersand.
        key: Exactly one integer key within that block.
        old: Required original value.
        new: Approved replacement value.

    Returns:
        Original text with only the matched integer token replaced.

    Raises:
        ValueError: If blocks, key counts or the original value differ.
    """
    blocks = list(re.finditer(rf"(?ims)^[ \t]*&{namelist}\b.*?^[ \t]*/[ \t]*$", text))
    if len(blocks) != 1:
        raise ValueError(f"Expected exactly one &{namelist} block")
    block = blocks[0]
    entries = list(re.finditer(rf"(?i)\b{key}\s*=\s*([^,\s]+)", block.group()))
    if len(entries) != 1 or entries[0][1] != str(old):
        raise ValueError(f"Expected exactly one &{namelist}.{key}={old}")
    entry = entries[0]
    start, end = block.start() + entry.start(1), block.start() + entry.end(1)
    return text[:start] + str(new) + text[end:]


def replay_input(text: str) -> str:
    """Apply and reversibly verify the two authorized budget-only changes.

    Args:
        text: Original pilot input with maxcycle250 and maxcyc100000.

    Returns:
        Long replay input with maxcycle32000 and maxcyc300000.

    Raises:
        ValueError: If expected budgets, CG/step or control inventory differ.
    """
    for pattern, expected in (
        (r"\boptalg\s*=\s*'([^']+)'", ["CG"]),
        (r"\bmaxstep\s*=\s*([^,\s]+)", ["0.02"]),
        (r"\bneb_nimage\s*=\s*(\d+)", ["25"]),
    ):
        if re.findall(pattern, text, re.I) != expected:
            raise ValueError("Expected unchanged CG/maxstep0.02/neb_nimage25")
    controls = re.findall(r"\bcoords2\((\d+)\)\s*=", text, re.I)
    if (
        controls != [str(i) for i in range(1, 13)]
        or len(re.findall(r"\brefcrd\s*=", text, re.I)) != 1
    ):
        raise ValueError("Expected exactly thirteen ordered original seed controls")
    changed = replace_budget(text, "dlfind", "maxcycle", 250, 32000)
    changed = replace_budget(changed, "cntrl", "maxcyc", 100000, 300000)
    restored = replace_budget(changed, "dlfind", "maxcycle", 32000, 250)
    restored = replace_budget(restored, "cntrl", "maxcyc", 300000, 100000)
    if restored != text:
        raise ValueError("Changes extend beyond authorized budget tokens")
    return changed


def prepare(pilot: Path, destination: Path) -> Path:
    """Freeze the pilot lineage and prepare a new replay child without launching.

    Args:
        pilot: Completed short CG trial with accepted seed-preparation records.
        destination: New long-run directory; an existing directory is an error.

    Returns:
        New configuration path for root review and explicit launch.

    Raises:
        ValueError: If inherited scientific gates/configuration or control paths differ.
    """
    config_path = pilot / "config.json"
    mdin = pilot / "calculation/minimize.mdin"
    cfg = json.loads(config_path.read_text())
    preflight_path = pilot / "preflight.json"
    preflight = json.loads(preflight_path.read_text())
    if (
        preflight["scientific_preparation_passed"] is not True
        or cfg["regular_images"] != 24
        or cfg["active_atoms"] != 4662
        or cfg["qm_charge"] != -1
        or cfg["optimizer"] != "CG"
        or len(cfg["qm_atom_ids"]) != 102
        or len(cfg["images"]) != 13
    ):
        raise ValueError(
            "Pilot must retain acceptedCG/QM102/-1/24images/4662active/13controls"
        )
    original = mdin.read_text()
    changed = replay_input(original)
    controls = re.findall(
        r"\b(?:refcrd|coords2\(\d+\))\s*=\s*'([^']+)'", original, re.I
    )
    if controls != cfg["images"] or cfg["reactant"] != controls[0]:
        raise ValueError("Original native control references differ from pilotconfig")
    seed = Path(cfg["source_seed_manifest"])
    if digest(seed) != cfg["trial_source_sha256"][str(seed)]:
        raise ValueError("Original seed manifest identity changed")
    sander = Path(cfg["sander"])
    amber = sander.parent.parent
    sources = [
        config_path,
        mdin,
        preflight_path,
        seed,
        Path(cfg["case_path"]),
        Path(cfg["topology"]),
        sander,
        amber / "lib/libdlfind.so",
        ROOT / "scripts/run_nyl12_job.py",
        ROOT / "src/nyl12_covalent/jobs.py",
        Path(__file__).resolve(),
        *map(Path, controls),
    ]
    for name in ("startup_precision_review.json", "startup_verification.json"):
        sources.append(pilot / name)
    hashes = {str(p): digest(p) for p in sources}
    destination.mkdir(parents=True, exist_ok=False)
    calculation = destination / "calculation"
    calculation.mkdir()
    snapshot = destination / "source_snapshot"
    snapshot.mkdir()
    copies = {}
    for index, path in enumerate(sources):
        if path in (sander, amber / "lib/libdlfind.so"):
            continue  # Engine binaries are identity-pinned and rehashed by the job runner.
        copy = snapshot / f"{index:02d}_{path.name}"
        shutil.copy2(path, copy)
        copies[str(path)] = str(copy)
    (calculation / "minimize.mdin").write_text(changed)
    cfg.update(
        calculation=str(calculation),
        base_input=str(calculation / "minimize.mdin"),
        parent_trial_config=str(config_path),
        trial_name="04_uniform_24_cg_long_replay",
        session_name="nyl12_neb_cg_long_20261119",
        dlfind_maxcycle=32000,
        cntrl_maxcyc=300000,
        native_energy_evaluation_budget=224000,
        replay_kind="Identical original seed replay; not an in-memory optimizer restart",
        long_replay_source_sha256=hashes,
    )
    output = destination / "config.json"
    output.write_text(json.dumps(cfg, indent=2) + "\n")
    record = {
        "scientific_preparation_passed": True,
        "budget_only_delta_verified": True,
        "pilot_config": str(config_path),
        "pilot_preflight_sha256": digest(preflight_path),
        "source_seed_manifest": str(seed),
        "source_seed_manifest_sha256": digest(seed),
        "source_sha256": hashes,
        "snapshots": copies,
        "controls_unchanged": controls,
        "edits": {"dlfind.maxcycle": [250, 32000], "cntrl.maxcyc": [100000, 300000]},
        "evaluation_budget": 224000,
        "pilot_budget": 1750,
        "budget_multiplier": 128,
        "force_converged": None,
        "chemical_acceptance": None,
        "replay": cfg["replay_kind"],
        "new_mdin_sha256": digest(calculation / "minimize.mdin"),
    }
    (destination / "preflight.json").write_text(json.dumps(record, indent=2) + "\n")
    command = [
        str(amber / "bin/python"),
        str(ROOT / "scripts/run_nyl12_job.py"),
        "launch",
        "--session",
        cfg["session_name"],
        "--output-dir",
        str(destination / "job"),
        "--cwd",
        str(calculation),
    ]
    for key, value in {
        "AMBERHOME": str(amber),
        "LD_LIBRARY_PATH": str(amber / "lib"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "GFORTRAN_UNBUFFERED_ALL": "1",
    }.items():
        command.extend(["--env", f"{key}={value}"])
    for path in [
        *sources,
        output,
        destination / "preflight.json",
        calculation / "minimize.mdin",
    ]:
        command.extend(["--input", str(path)])
    command.extend(
        [
            "--",
            str(sander),
            "-O",
            "-i",
            "minimize.mdin",
            "-p",
            cfg["topology"],
            "-c",
            cfg["reactant"],
            "-o",
            "minimized.mdout",
            "-r",
            "optimized.ncrst",
            "-x",
            "optimization.nc",
            "-inf",
            "mdinfo",
        ]
    )
    launcher = destination / "launch.sh"
    launcher.write_text(
        "#!/bin/bash\nset -euo pipefail\ncd "
        + shlex.quote(str(ROOT))
        + "\nexec "
        + shlex.join(command)
        + "\n"
    )
    launcher.chmod(0o750)
    for source in snapshot.iterdir():
        source.chmod(0o444)
    if digest(mdin) != hashes[str(mdin)]:
        raise ValueError("Pilot input changed during preparation")
    LOGGER.info("Prepared %s; budget224000; pilot unchanged; not launched", output)
    return output


def main() -> None:
    """Prepare the explicit pilot child; never execute its launch script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, default=DEFAULT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    pilot = args.pilot.resolve()
    prepare(pilot, args.output.resolve() if args.output else pilot / "long_run")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
