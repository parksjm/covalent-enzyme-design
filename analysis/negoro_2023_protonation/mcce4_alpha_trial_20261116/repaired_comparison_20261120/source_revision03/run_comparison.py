"""Run a profile-driven, upstream-gated native MCCE diagnostic comparison."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace
from typing import Any

LOGGER = logging.getLogger(__name__)


def sha(path: Path) -> str:
    """Return a retained artifact digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path: Path, value: dict) -> None:
    """Atomically publish status or evidence."""
    pending = path.with_suffix(".pending")
    pending.write_text(json.dumps(value, indent=2) + "\n")
    pending.replace(path)


class Controller:
    """Keep one structure's native process and scientific records separate."""

    def __init__(self, profile_path: Path) -> None:
        self.profile = json.loads(profile_path.read_text())
        self.base = Path(self.profile["output"])
        self.root = Path(self.profile["root"])
        self.base.mkdir(exist_ok=False)
        self.status("starting")
        shutil.copy2(profile_path, self.base / "profile.json")

    def journal(self, phase: str, action: str) -> None:
        """Append each controller action before and after execution."""
        subprocess.run(
            [
                sys.executable,
                str(self.root / "scripts/record_agent_action.py"),
                "--log",
                str(self.root / "output_agent_actions/20261021.jsonl"),
                "--agent",
                "nylc_mcce_queue",
                "--phase",
                phase,
                "--action",
                action,
                "--permission-mode",
                "escalated",
                "--reason",
                "Authorized native MCCE comparison",
                "--rollback",
                "Retain all immutable inputs, outputs and failures",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )

    def status(self, state: str, **values: Any) -> None:
        """Publish controller state without claiming scientific pKa validity."""
        write(
            self.base / "controller_status.json",
            {
                "state": state,
                "utc": datetime.now(timezone.utc).isoformat(),
                **values,
            },
        )

    def wait(self, path: Path) -> dict:
        """Wait at most twelve hours for exactly the declared upstream job."""
        self.journal("intent", f"Wait for declared upstream {path}")
        deadline = time.monotonic() + 43200
        while time.monotonic() < deadline:
            if not path.exists():
                self.status(
                    "waiting_upstream", upstream=str(path), status_file_exists=False
                )
                time.sleep(30)
                continue
            status = json.loads(path.read_text())
            if status[self.profile["upstream_exit_key"]] is not None:
                write(self.base / f"{path.parent.name}_upstream_status.json", status)
                if status[self.profile["upstream_exit_key"]] != 0:
                    raise ValueError(f"Upstream failed: {path}")
                self.journal("outcome", f"Upstream completed successfully {path}")
                return status
            self.status("waiting_upstream", upstream=str(path))
            time.sleep(30)
        raise TimeoutError(f"Upstream wait expired: {path}")

    def stage(self, label: str, argv: list[str], seconds: int) -> None:
        """Run an owned GNU timeout group and retain its actual exit status."""
        command = [
            "/usr/bin/timeout",
            "--signal=TERM",
            "--kill-after=30s",
            str(seconds),
            *argv,
        ]
        record = {
            "command": command,
            "started_utc": datetime.now(timezone.utc).isoformat(),
        }
        self.journal("intent", f"{label}: {json.dumps(command)}")
        write(self.base / f"{label}_process_intent.json", record)
        self.status("running", stage=label)
        with Path(f"{label}.log").open("x") as stream:
            result = subprocess.run(
                command, stdout=stream, stderr=subprocess.STDOUT, check=False
            )
        Path(f"{label}.exit").write_text(str(result.returncode) + "\n")
        record.update(
            exit_status=result.returncode,
            ended_utc=datetime.now(timezone.utc).isoformat(),
        )
        write(self.base / f"{label}_process_outcome.json", record)
        self.journal("outcome", f"{label}: exit {result.returncode}")
        text = Path(f"{label}.log").read_text(errors="replace")
        if result.returncode or re.search(
            r"Error!|\bFATAL\b|\[ERROR\]|^ERROR[: ]|Traceback", text, re.I | re.M
        ):
            raise RuntimeError(f"{label} failed; inspect retained native log")

    def run(self) -> None:
        """Adopt accepted heavy atoms, then run all four native MCCE stages."""
        p = self.profile
        for name, expected in p["source_sha256"].items():
            if sha(Path(name)) != expected:
                raise ValueError(f"Source changed: {name}")
        upstream = self.wait(Path(p["upstream_status"]))
        if "refinement" in p:
            r = p["refinement"]
            os.chdir(self.base)
            self.stage(
                "refinement",
                [
                    r["python"],
                    r["script"],
                    "--candidate",
                    r["candidate"],
                    "--original",
                    r["original"],
                    "--outputdir",
                    r["output"],
                    "--termini",
                    p["termini_path"],
                ],
                43200,
            )
        receipt = json.loads(Path(p["receipt"]).read_text())
        if receipt["state"] != "accepted":
            raise ValueError("Refinement scientific acceptance is not accepted")
        validate_adoption(p, receipt, upstream, self.base)
        source = Path(p["heavy"])
        if sha(source) != receipt["output_sha256"]:
            raise ValueError("Accepted heavy input digest mismatch")
        self.journal(
            "intent", f"Adopt accepted heavy input {source}, canonicalize atom names"
        )
        write(self.base / "upstream_scientific_acceptance.json", receipt)
        write(
            self.base / "adoption.json",
            {"heavy_sha256": sha(source), "upstream_status": upstream},
        )
        from nyl12_covalent.mcce_preparation import prepare_input
        from nyl12_covalent.mcce_validation import audit_preparation

        folder = self.base / "prepared"
        prepare_input(source, folder, receipt["output_sha256"], p["expected_termini"])
        self.journal("outcome", f"Canonical input ready {folder}")
        software = Path(p["software"])
        os.environ["PATH"] = (
            str(software / "bin")
            + ":"
            + str(Path(p["python"]).parent)
            + ":/usr/bin:/bin"
        )
        if shutil.which("delphi") != str(software / "bin/delphi"):
            raise ValueError("Declared installed DelPhi not resolved")
        sys.path.insert(0, str(software / "bin"))
        import step1
        import step2
        import step4

        os.chdir(folder)
        self.journal("intent", "Write native step1 preprocessing and run.prm")
        step1.write_runprm(
            SimpleNamespace(noter=False, d=4.0, dry=True, u="", load_runprm="")
        )
        lines = step1.identify_ligands(step1.fix_format("protein.pdb", noter=False))
        Path("step0_out.pdb").write_text("".join(lines))
        shutil.copy2("run.prm", self.base / "step1_run.prm")
        self.journal("outcome", "Native step1 preprocessing written")
        self.stage("step1", [str(software / "bin/mcce")], 600)
        acceptance = audit_preparation(
            folder, folder / "protein.pdb", p["expected_termini"]
        )
        write(self.base / "preparation_acceptance.json", acceptance)
        if not acceptance["accepted"]:
            raise ValueError(
                f"MCCE step1 scientific gate failed: {acceptance['errors']}"
            )
        self.journal("intent", "Write native step2 l1 d4 parameters")
        step2.write_runprm(SimpleNamespace(d=4.0, l=1, u="", load_runprm="", ftpl=""))
        shutil.copy2("run.prm", self.base / "step2_run.prm")
        self.journal("outcome", "Native step2 parameters written")
        self.stage("step2", [str(software / "bin/mcce")], 7200)
        if not Path("step2_out.pdb").is_file():
            raise FileNotFoundError("step2_out.pdb")
        self.stage(
            "step3",
            [
                p["python"],
                str(software / "bin/step3.py"),
                "-s",
                "delphi",
                "-d",
                "4",
                "-do",
                "80",
                "-salt",
                "0.15",
                "-p",
                "2",
                "-t",
                "pbe_tmp",
            ],
            43200,
        )
        if (
            Path("failed_conformers.txt").exists()
            and Path("failed_conformers.txt").stat().st_size
        ):
            raise ValueError("Native PBE conformers failed")
        rows = [
            line.split()
            for line in Path("head3.lst").read_text().splitlines()[1:]
            if line.strip()
        ]
        if not rows or any(len(r) < 17 or r[16] not in {"t", "d"} for r in rows):
            raise ValueError("Incomplete electrostatics conformer table")
        if any(
            not (Path("energies") / f"{r[1]}.opp").is_file()
            for r in rows
            if r[16] == "t"
        ):
            raise ValueError("Missing native pair-energy artifact")
        self.journal(
            "intent",
            "Write native six-run 2000-iteration Monte Carlo, MS_OUT false storage only",
        )
        step4.write_runprm(
            SimpleNamespace(
                t="ph",
                i="0.0",
                d="1.0",
                n="15",
                xts=True,
                ms=False,
                u=f"MONTE_SEED={p['seed']},MONTE_T=298.15",
                load_runprm="",
            )
        )
        shutil.copy2("run.prm", self.base / "step4_run.prm")
        self.journal("outcome", "Native step4 parameters retained")
        self.stage("step4", [str(software / "bin/mcce")], 7200)
        validate_occupancies(Path("fort.38"), [r[1] for r in rows])
        names = ("head3.lst", "pK.out", "fort.38", "run.prm.record")
        write(
            self.base / "native_completion.json",
            {
                "state": "native_completed_pka_unqualified",
                "parameter_model": p["parameter_model"],
                "pka_qualification": "pending_known_proline_FTPL_legacy_inconsistency",
                "assigned_states": False,
                "sampling_convergence": False,
                "monte_runs": 6,
                "MONTE_NITER": 2000,
                "temperature_K": 298.15,
                "seed": p["seed"],
                "MS_OUT": False,
                "MS_OUT_effect": "storage_only_sampling_unchanged",
                "artifact_sha256": {n: sha(Path(n)) for n in names},
            },
        )
        self.status("native_completed_pka_unqualified")


def validate_adoption(p: dict, receipt: dict, upstream: dict, base: Path) -> None:
    """Verify exact heavy residue sequences and explicit legacy source identity."""
    if receipt["geometry"]["failures"]:
        raise ValueError("Refinement geometry failures remain")
    if "declared_termini" in receipt:
        if receipt["declared_termini"] != p["expected_termini"]:
            raise ValueError("Refinement termini declaration differs")
    elif p["adoption_contract"] != "legacy_refinement_without_termini_field":
        raise ValueError("Missing explicit refinement termini declaration")
    observed: dict[str, dict[int, str]] = defaultdict(dict)
    oxt: dict[str, list[int]] = defaultdict(list)
    for line in Path(p["heavy"]).read_text().splitlines():
        if line.startswith("ATOM  "):
            chain, number = line[21], int(line[22:26])
            if line[26] != " ":
                raise ValueError("Unexpected insertion code")
            observed[chain][number] = line[17:20]
            if line[12:16].strip() == "OXT":
                oxt[chain].append(number)
    sequences: dict[str, list[str]] = defaultdict(list)
    for line in Path(p["original"]).read_text().splitlines():
        if line.startswith("SEQRES"):
            sequences[line[11]].extend(line[19:70].split())
    expected = p["expected_termini"]
    if set(observed) != set(expected) or set(sequences) != set(expected):
        raise ValueError("Adopted chain inventory differs")
    for chain, ends in expected.items():
        numbers = list(range(ends["n"], ends["c"] + 1))
        if sorted(observed[chain]) != numbers or oxt[chain] != [ends["c"]]:
            raise ValueError(f"Adopted residue/OXT inventory differs: {chain}")
        if [observed[chain][n] for n in numbers] != sequences[chain]:
            raise ValueError(f"Adopted SEQRES sequence differs: {chain}")
    if p["adoption_contract"] == "legacy_refinement_without_termini_field":
        snapshot = Path(p["receipt"]).parent / "refine_rebuilt_protein_atoms.py"
        expected_hash = upstream["input_sha256"][
            str(Path(p["root"]) / "scripts/refine_rebuilt_protein_atoms.py")
        ]
        if sha(snapshot) != expected_hash:
            raise ValueError("Legacy refinement source snapshot identity differs")
    write(
        base / "adoption_validation.json",
        {
            "accepted": True,
            "contract": p["adoption_contract"],
            "expected_termini": expected,
            "observed_residues": observed,
            "observed_OXT": oxt,
            "sequences_match_original_SEQRES": True,
            "native_upstream_exit_status": 0,
            "refinement_source_hashes": receipt["source_hashes"],
        },
    )


def validate_occupancies(path: Path, conformers: list[str]) -> None:
    """Require exact conformers, finite probabilities and residue normalization."""
    lines = path.read_text().splitlines()
    if [float(v) for v in lines[0].split()[1:]] != list(range(15)):
        raise ValueError("Unexpected titration pH grid")
    seen = []
    sums: dict[str, list[float]] = defaultdict(lambda: [0.0] * 15)
    counts: dict[str, int] = defaultdict(int)
    for line in lines[1:]:
        fields = line.split()
        if not fields:
            continue
        name, values = fields[0], [float(v) for v in fields[1:]]
        if len(values) != 15 or not all(
            math.isfinite(v) and 0 <= v <= 1 for v in values
        ):
            raise ValueError("Invalid native occupancy")
        seen.append(name)
        residue = name[:3] + name[5:11]
        counts[residue] += 1
        sums[residue] = [a + b for a, b in zip(sums[residue], values)]
    if len(set(seen)) != len(seen) or set(seen) != set(conformers):
        raise ValueError("Native occupancy conformer inventory mismatch")
    if any(
        abs(v - 1) > 0.0005 * counts[r] + 1e-9 for r, row in sums.items() for v in row
    ):
        raise ValueError(
            "Native per-residue occupancy does not sum to one at output precision"
        )


def main() -> None:
    """Run the declared profile, preserving explicit failure state."""
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, required=True)
    args = parser.parse_args()
    controller = Controller(args.profile.resolve())
    try:
        controller.run()
    except BaseException as error:
        controller.status("failed", error=f"{type(error).__name__}: {error}")
        controller.journal(
            "outcome", f"Controller failed: {type(error).__name__}: {error}"
        )
        raise


if __name__ == "__main__":
    main()
