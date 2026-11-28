"""Resume only the frozen benchmark's Mistral model on an isolated endpoint."""

from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

BENCHMARK_ROOT = Path(
    "/data/jpc/mechanism-agent/benchmarks/nemotron_mistral_nano_20261125"
)
ISOLATED_PORT = 18112
LOGGER = logging.getLogger(__name__)


def load_module(name: str, path: Path) -> ModuleType:
    """Load a frozen benchmark module by source path."""
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load benchmark module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def run(root: Path, output: Path, port: int = ISOLATED_PORT) -> dict[str, Any]:
    """Run Mistral on an isolated port and retain a separate completion receipt."""
    if port != ISOLATED_PORT:
        raise ValueError(f"Recovery port must be {ISOLATED_PORT}")
    root = root.resolve()
    output = output.resolve()
    if output == root:
        raise ValueError("Completion output must be a new directory")
    if output.is_relative_to(root):
        relative_parts = output.relative_to(root).parts
        if relative_parts[0] in {"results", "source_snapshot", "inputs"}:
            raise ValueError("Completion output overlaps protected campaign data")
    output.mkdir(parents=True, exist_ok=False)

    source = root / "source_snapshot"
    sys.path.insert(0, str(source))
    try:
        protocol = load_module("protocol", source / "protocol.py")
        runner = load_module("mistral_recovery_runner", source / "runner.py")
        cases, campaign, cases_hash = runner.load_inputs(root)
    finally:
        sys.path.remove(str(source))

    # Preserve every frozen setting; alter only the transport port in memory.
    campaign = json.loads(json.dumps(campaign))
    original_port = campaign["models"]["mistral"]["port"]
    campaign["models"]["mistral"]["port"] = port
    # The frozen runner imports these transport guards from protocol.py.
    protocol.PORTS["mistral"] = port
    protocol.EXPECTED_MODELS["mistral"]["port"] = port
    if (
        campaign["models"]["nemotron"]["port"]
        != protocol.EXPECTED_MODELS["nemotron"]["port"]
    ):
        raise ValueError("Nemotron configuration differs from the frozen protocol")

    # Exercise the first frozen case as a canary before permitting the full resume.
    if not cases:
        raise ValueError("Frozen campaign has no cases")
    canary_ok = runner.run_model("mistral", root, cases[:1], campaign, cases_hash, True)
    if not canary_ok:
        raise RuntimeError("Mistral first-case canary did not complete")
    progress_path = root / "results" / "mistral" / "progress.json"
    progress = json.loads(progress_path.read_bytes())
    first_ref = progress.get("cases", {}).get(cases[0]["case_id"])
    if first_ref is None:
        raise RuntimeError("Mistral first-case canary did not publish progress")
    first_receipt = json.loads(
        (root / "results" / "mistral" / first_ref["receipt"]).read_bytes()
    )
    if (
        first_receipt.get("process_status") != "complete"
        or first_receipt.get("finish_reason") != "stop"
    ):
        raise RuntimeError(
            "Mistral first-case canary was not a successful stop response"
        )

    ok = runner.run_model("mistral", root, cases, campaign, cases_hash, True)
    progress_path = root / "results" / "mistral" / "progress.json"
    progress = json.loads(progress_path.read_bytes())
    case_entries = progress.get("cases", {})
    successful = 0
    for reference in case_entries.values():
        receipt = json.loads(
            (root / "results" / "mistral" / reference["receipt"]).read_bytes()
        )
        successful += (
            receipt.get("finish_reason") == "stop"
            and receipt.get("process_status") == "complete"
        )
    failed = len(cases) - successful
    request_hashes = {
        case_id: json.loads(
            (root / "results" / "mistral" / ref["receipt"]).read_bytes()
        ).get("request_sha256")
        for case_id, ref in case_entries.items()
    }
    result = {
        "model": "mistral",
        "process": "complete" if ok and successful == len(cases) else "failed",
        "cases_sha256": cases_hash,
        "original_port": original_port,
        "recovery_port": port,
        "successful_stop_responses": successful,
        "failed_or_incomplete": failed,
        "request_sha256_by_case": request_hashes,
        "finished_utc": protocol.utc_now(),
    }
    (output / "mistral_recovery_completion.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    """Parse command-line options and run the isolated Mistral resume."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=BENCHMARK_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--port", type=int, default=ISOLATED_PORT)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    result = run(args.root, args.output, args.port)
    LOGGER.info("Mistral recovery finished: %s", result["process"])
    return 0 if result["process"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
