"""Check immutable-source and no-overwrite contracts of the uncertainty CLI."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SPEC = importlib.util.spec_from_file_location(
    "uncertainty_cli",
    Path(__file__).resolve().parents[1] / "scripts/analyze_relay_pmf_uncertainty.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_cli_keeps_source_identity_and_never_overwrites(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request, report, receipt = (
        tmp_path / name for name in ("request.json", "report.json", "receipt.json")
    )
    receipt.write_text("{}")
    request.write_text(
        json.dumps(
            dict(
                temperature_K=300,
                physical_hamiltonian_sha256="a" * 64,
                diagnostic_columns=[],
                accepted_block_receipts=[str(receipt)],
            )
        )
    )
    sha = hashlib.sha256(request.read_bytes()).hexdigest()
    report.write_text(json.dumps(dict(request_sha256=sha)))
    options = SimpleNamespace(
        request=request,
        report=report,
        output=tmp_path / "output.json",
        bootstrap_replicates=200,
        bootstrap_workers=2,
    )
    monkeypatch.setattr(MODULE.argparse.ArgumentParser, "parse_args", lambda _: options)
    monkeypatch.setattr(MODULE, "load_block", lambda *a, **kw: dict(validated=True))
    monkeypatch.setattr(
        MODULE,
        "report_uncertainty",
        lambda *a, **kw: dict(
            identity_checked=True, bootstrap_replicates=kw["bootstrap_replicates"]
        ),
    )
    MODULE.main()
    result = json.loads(options.output.read_text())
    assert result["request_sha256"] == sha
    assert result["identity_checked"]
    assert result["bootstrap_replicates"] == 200
    with pytest.raises(FileExistsError):
        MODULE.main()
    options.output = tmp_path / "other.json"
    request.write_text(request.read_text() + " ")
    with pytest.raises(ValueError, match="hash"):
        MODULE.main()
