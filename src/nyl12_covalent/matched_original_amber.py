"""Run source-bound matched original-backbone predictions through Violet Amber."""

from __future__ import annotations

from collections import defaultdict
import importlib.util
from itertools import product
import json
import math
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
from typing import Any

from .af3_postprocess import read_atom_sites
from .control_inspection import export_control, pdb_atoms, sha256
from .design_case_postprocessing import annotate_amide_registers, process_design_case
from .design_model_relaxation import relax_model
from .negative_amide_gate import evaluate_signed_amide_gate
from .relay_broad_prediction_io import audit_prediction

ENGINE = Path("/data/jpc/code/pmemd26/bin/pmemd.cuda_SPFP")
ENGINE_SHA256 = "7ed3f439d09d6beff94210bb3f356585e6918cd3f500dcc58046067ac8fd9602"
ROOT = Path(os.environ.get("NYL12_PROJECT_ROOT", Path(__file__).resolve().parents[2]))


def _write(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def _journal(phase: str, action: str, result: str) -> None:
    if not (ROOT / "scripts/record_agent_action.py").is_file():
        raise FileNotFoundError("Set NYL12_PROJECT_ROOT to the active journal project")
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "matched_original_amber",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            result,
            "--rollback",
            "Preserve native predictions and failed outputs; new attempt uses new output directory",
        ],
        check=True,
    )


def validate_plan(plan: dict[str, Any]) -> None:
    """Validate a complete matched two-design plus WT, two-seed frozen plan.

    Args:
        plan: Explicit cases, scientific inputs, source pins and Amber configuration.

    Raises:
        ValueError: Case matching, positive sign, executable or a frozen input fails.
    """
    if (
        plan["schema_version"] != 1
        or type(plan["target_phi_sign"]) is not int
        or plan["target_phi_sign"] != 1
    ):
        raise ValueError("Require schema 1 and explicit target_phi_sign=+1")
    cases = plan["cases"]
    candidates = {c["candidate"] for c in cases}
    expected = {(c, s) for c in candidates for s in (20261126, 20261127)}
    if (
        len(candidates) != 3
        or "WT" not in candidates
        or len(cases) != 6
        or {(c["candidate"], c["seed"]) for c in cases} != expected
    ):
        raise ValueError("Require two designs plus WT, each with seeds 20261126/27")
    names = [c["case"] for c in cases]
    if len(set(names)) != 6 or any(
        not n or Path(n).name != n or n in {".", ".."} for n in names
    ):
        raise ValueError("Require six unique safe case identities")
    cfg = plan["relaxation_config"]
    if Path(cfg["engine"]) != ENGINE or sha256(ENGINE) != ENGINE_SHA256:
        raise ValueError("Violet validated Amber executable path/hash changed")
    required = [
        Path(plan[k])
        for k in (
            "approved_reference_manifest",
            "canonical_adapter",
            "legacy_processor",
        )
    ]
    required += [Path(cfg[k]) for k in ("mol2", "mapping", "engine", "tleap")]
    required += sorted(p for p in Path(cfg["release"]).rglob("*") if p.is_file())
    if not Path(cfg["release"]).is_dir() or not required:
        raise ValueError("Require released modular RESP/GAFF2 parameters")
    required += [Path(c["case_input"]) for c in cases]
    pins = plan["input_sha256"]
    for path in required:
        if not path.is_absolute() or str(path) not in pins:
            raise ValueError(f"Missing absolute input/source hash pin: {path}")
    for name, digest in pins.items():
        if not Path(name).is_absolute() or sha256(Path(name)) != digest:
            raise ValueError(f"Frozen input changed: {name}")
    reference = json.loads(Path(plan["approved_reference_manifest"]).read_text())
    if reference["target_phi_sign"] != 1:
        raise ValueError("Approved reference must bind positive phi")


def severe_clashes(source: Path) -> dict[str, Any]:
    """Count interchain heavy pairs strictly below 1 Å with a spatial grid.

    Args:
        source: Canonical mature-local heavy-atom PDB.

    Returns:
        Pair count and bounded examples; any such pair fails assembly preflight.
    """
    cells: dict[tuple[int, ...], list[tuple[Any, Any]]] = defaultdict(list)
    count, examples = 0, []
    atoms = pdb_atoms(source)
    if not atoms:
        raise ValueError("Empty canonical source")
    for key, (element, xyz) in atoms.items():
        if element.upper() in {"H", "D"}:
            continue
        if not all(math.isfinite(v) for v in xyz):
            raise ValueError("Non-finite canonical coordinates")
        cell = tuple(math.floor(v) for v in xyz)
        for offset in product((-1, 0, 1), repeat=3):
            neighbor = tuple(a + b for a, b in zip(cell, offset, strict=True))
            for other, position in cells[neighbor]:
                distance2 = sum(
                    (a - b) ** 2 for a, b in zip(xyz, position, strict=True)
                )
                if key[0] != other[0] and distance2 < 1.0:
                    count += 1
                    if len(examples) < 20:
                        examples.append(
                            {"atoms": [key, other], "distance_A": math.sqrt(distance2)}
                        )
        cells[cell].append((key, xyz))
    return {
        "threshold_A": 1.0,
        "interchain_heavy_pairs_below_threshold": count,
        "examples": examples,
        "accepted": count == 0,
    }


def _canonicalize(case: dict[str, Any], plan: dict[str, Any], out: Path) -> Path:
    receipt = json.loads(Path(case["native_receipt"]).read_text())
    for key in ("case", "candidate", "seed"):
        if receipt[key] != case[key]:
            raise ValueError(f"Native receipt identity mismatch: {key}")
    if (
        receipt["status"] != "native_prediction_completed"
        or receipt["exit_status"] != 0
        or receipt["target_phi_sign"] != 1
    ):
        raise ValueError("Require completed native prediction with positive target")
    raw = Path(receipt["prediction_cif"])
    if sha256(raw) != receipt["prediction_cif_sha256"]:
        raise ValueError("Native prediction hash changed")
    if Path(receipt["approved_reference_manifest"]) != Path(
        plan["approved_reference_manifest"]
    ):
        raise ValueError("Native receipt approved reference mismatch")
    for name, digest in receipt["input_sha256"].items():
        if sha256(Path(name)) != digest:
            raise ValueError(f"Native prediction input changed: {name}")
    shutil.copy2(raw, out / "original_native_prediction.cif")
    shutil.copy2(case["native_receipt"], out / "native_prediction.json")
    inputs = json.loads(Path(case["case_input"]).read_text())
    if any(inputs[k] != case[k] for k in ("case", "candidate", "seed")):
        raise ValueError("Case input identity mismatch")
    seq = inputs["sequences"]
    graph = Path(receipt["ligand_graph_manifest"])
    if str(graph) not in plan["input_sha256"]:
        raise ValueError("Ligand graph manifest needs frozen hash pin")
    maps = json.loads(graph.read_text())["ligand_atom_maps"]["boltz"]
    _write(
        out / "native_identity_acceptance.json",
        audit_prediction(raw, seq["alpha"], seq["beta"], maps["I"]),
    )
    spec = importlib.util.spec_from_file_location(
        "matched_canonical_adapter", plan["canonical_adapter"]
    )
    if spec is None or spec.loader is None:
        raise ImportError(plan["canonical_adapter"])
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    atoms = {
        (s.chain, s.sequence, s.name): (s.element, s.xyz) for s in read_atom_sites(raw)
    }
    pairs, evidence = adapter.pairing(atoms, adapter.reference_atoms())
    _write(out / "physical_pairing.json", evidence)
    derived = out / "canonical_pair_input.cif"
    _write(
        out / "canonicalization_receipt.json",
        adapter.canonical_cif(raw, derived, pairs),
    )
    cfg = plan["relaxation_config"]
    report = export_control(
        derived,
        maps,
        Path(cfg["mol2"]),
        Path(cfg["mapping"]),
        Path(cfg["release"]),
        out / "canonical",
    )
    rows = annotate_amide_registers(
        report["all48"],
        whole_context_mapping=Path(cfg["mapping"]),
        source_mol2=Path(cfg["mol2"]),
    )
    _write(
        out / "native_all48_positive_gate.json",
        {
            "all_amide_site_rows": rows,
            **evaluate_signed_amide_gate(rows, target_phi_sign=1),
        },
    )
    return out / "canonical/canonical_mature_local.pdb"


def run_case(
    case: dict[str, Any], plan: dict[str, Any], output: Path
) -> dict[str, Any]:
    """Preserve, screen, minimize and publish one case without a method fallback.

    Args:
        case: Exact candidate/seed native receipt and case-input binding.
        plan: Validated frozen campaign plan.
        output: New isolated case artifact directory.

    Returns:
        Separate native assembly, force and positive-pose acceptance records.
    """
    output.mkdir(parents=True, exist_ok=False)
    _write(output / "case_input.json", {"case": case, "plan": plan})
    source = _canonicalize(case, plan, output)
    clash = severe_clashes(source)
    _write(output / "assembly_precheck.json", clash)
    if not clash["accepted"]:
        return {
            "case": case["case"],
            "status": "failed_assembly_precheck",
            "amber_launched": False,
            "assembly": clash,
        }
    native = relax_model(source, output / "relaxation", plan["relaxation_config"])
    base = output / "relaxation"
    report = process_design_case(
        case=case["case"],
        target_phi_sign=1,
        topology=base / "minimization/endpoint_complex.prmtop",
        before_restart=base / "system/dry.rst7",
        after_restart=base / "minimization/endpoint_complex.rst7",
        conversion=base / "system/conversion.json",
        output_dir=output / "postprocess",
        legacy_processor=Path(plan["legacy_processor"]),
        whole_context_mapping=Path(plan["relaxation_config"]["mapping"]),
        source_mol2=Path(plan["relaxation_config"]["mol2"]),
    )
    return {
        "case": case["case"],
        "status": "completed",
        "amber_launched": True,
        "final_rms_gradient": native["final_rms_gradient"],
        "force_converged": native["force_converged"],
        "positive_ab_before": report["snapshots"]["before"]["assigned_ab_count"],
        "positive_ab_after": report["snapshots"]["after"]["assigned_ab_count"],
        "native_thr_policy": "diagnostic only",
        "postprocess": str(output / "postprocess"),
    }


def run_campaign(
    plan: dict[str, Any], output: Path, *, case_name: str | None = None
) -> dict[str, Any]:
    """Run selected cases once; retain failures and publish a failure census.

    Args:
        plan: Source-bound matched six-case configuration.
        output: Artifact root; existing case directories are never overwritten.
        case_name: Optional single case for an externally reserved GPU worker.

    Returns:
        Per-case process/scientific outcomes, including preserved failures.
    """
    validate_plan(plan)
    if socket.gethostname().split(".")[0].lower() != "violet":
        raise ValueError("This exact executable is qualified only on Violet")
    cases = [c for c in plan["cases"] if case_name is None or c["case"] == case_name]
    if not cases:
        raise ValueError("Unknown case")
    _journal("intent", "Initialize matched Amber artifacts", str(output))
    output.mkdir(parents=True, exist_ok=True)
    _journal("outcome", "Initialize matched Amber artifacts", str(output))
    summary: dict[str, Any] = {"target_phi_sign": 1, "cases": {}, "failure_count": 0}
    census = output / (
        f"failure_census_{case_name}.json" if case_name else "failure_census.json"
    )
    if census.exists() or any((output / c["case"]).exists() for c in cases):
        raise FileExistsError(
            "Existing campaign outputs require a new attempt directory"
        )
    for case in cases:
        validate_plan(plan)
        action = "Matched original Amber followthrough " + case["case"]
        _journal(
            "intent",
            action,
            "Reserved GPU supplied by explicit frozen plan; no process signaling or fallback",
        )
        try:
            result = run_case(case, plan, output / case["case"])
            validate_plan(plan)
        except Exception as exc:
            result = {
                "case": case["case"],
                "status": "failed",
                "error": repr(exc),
                "evidence": str(output / case["case"]),
            }
        summary["cases"][case["case"]] = result
        summary["failure_count"] += result["status"] != "completed"
        _write(census, summary)
        _journal("outcome", action, json.dumps(result, allow_nan=False))
    return summary
