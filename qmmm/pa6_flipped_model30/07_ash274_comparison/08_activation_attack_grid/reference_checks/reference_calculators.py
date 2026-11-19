from __future__ import annotations
import re
import shutil
import subprocess
from pathlib import Path
from reference_common import (
    BASELINE,
    SANDER,
    XTB,
    TOPOLOGY,
    digest,
    journal,
    write_json,
)
from reference_embedding import parse_external_input, write_xtb_inputs


def make_mdin(
    template: str, theory: str, executable_override: str | None = None
) -> str:
    """Build DFTB3 or EXTERN single-point input from the frozen ES template."""
    if theory == "EXTERN":
        return (
            template
            if executable_override is None
            else template.replace(
                "executable='/usr/local/g16/g16'", f"executable='{executable_override}'"
            )
        )
    if theory != "DFTB3":
        raise ValueError(f"Unsupported Sander theory: {theory}")
    text = template.replace("qm_theory='EXTERN'", "qm_theory='DFTB3'").replace(
        "verbosity=1", "verbosity=2"
    )
    start = text.index("&gau\n")
    end = text.index("&dlfind", start)
    return text[:start] + text[end:]


def run_sander(
    folder: Path, restart: Path, theory: str, template: str, environment: dict[str, str]
) -> float:
    """Run one frozen-geometry Sander force evaluation and return its energy."""
    folder.mkdir(parents=True, exist_ok=False)
    shutil.copy2(restart, folder / "input.ncrst")
    (folder / "minimize.mdin").write_text(make_mdin(template, theory), encoding="utf-8")
    command = [
        str(SANDER),
        "-O",
        "-i",
        "minimize.mdin",
        "-p",
        str(TOPOLOGY),
        "-c",
        "input.ncrst",
        "-o",
        "single_point.mdout",
        "-r",
        "single_point.ncrst",
        "-inf",
        "mdinfo",
    ]
    record = {
        "status": "running",
        "method": theory,
        "command": command,
        "input_sha256": {
            name: digest(folder / name) for name in ("minimize.mdin", "input.ncrst")
        },
        "topology_sha256": digest(TOPOLOGY),
    }
    write_json(folder / "process.json", record)
    journal("intent", f"Run {theory} reference at {folder}")
    with (folder / "launcher.log").open("x") as log:
        process = subprocess.run(
            command,
            cwd=folder,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    output = (folder / "single_point.mdout").read_text(errors="replace")
    field = "DFTBESCF" if theory == "DFTB3" else "EXTERNESCF"
    values = re.findall(field + r"\s*=\s*([-+\d.EeDd]+)", output)
    complete = (
        process.returncode == 0
        and bool(values)
        and "Energy calculation finished" in output
    )
    record.update(
        {
            "status": "finished" if complete else "failed",
            "exit_status": process.returncode,
            "energy_kcal_mol": float(values[-1].replace("D", "E")) if values else None,
            "gradient_artifact": "Sander/Gaussian raw output retained; projection pending",
            "outputs_sha256": {
                p.name: digest(p)
                for p in folder.iterdir()
                if p.is_file() and p.name != "process.json"
            },
        }
    )
    write_json(folder / "process.json", record)
    journal("outcome", f"Run {theory} reference at {folder}", record["status"])
    if not complete:
        raise RuntimeError(
            f"{theory} failed or lacked its energy/gradient output: {folder}"
        )
    return float(values[-1].replace("D", "E"))


def run_xtb(
    folder: Path,
    gaussian_input: Path,
    environment: dict[str, str],
    baseline_charges: list[float],
) -> float:
    """Run GFN2-xTB on Sander's exact QM/link and adjusted-PC export."""
    embedding = write_xtb_inputs(gaussian_input, folder)
    charges = [charge[3] for charge in embedding["point_charges"]]
    if len(charges) != len(baseline_charges) or any(
        abs(a - b) > 1e-10 for a, b in zip(charges, baseline_charges)
    ):
        raise ValueError(
            "External point-charge count/order/charges differ from the ES baseline"
        )
    command = [
        str(XTB),
        "qm_link.xyz",
        "--gfn",
        "2",
        "--chrg",
        "-1",
        "--uhf",
        "0",
        "--input",
        "xtb.inp",
        "--namespace",
        "gfn2",
        "--grad",
    ]
    with (folder / "gfn2.out").open("x") as log:
        process = subprocess.run(
            command,
            cwd=folder,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    output = (folder / "gfn2.out").read_text(errors="replace")
    values = re.findall(r"\| TOTAL ENERGY\s+([-+\d.Ee]+)\s+Eh", output)
    if (
        process.returncode != 0
        or not values
        or "normal termination of xtb" not in output
    ):
        raise RuntimeError(f"GFN2-xTB failed: {folder}")
    return float(values[-1]) * 627.5094740631


def capture_gaussian_export(
    folder: Path, restart: Path, template: str, environment: dict[str, str]
) -> Path:
    """Ask Sander to write its EXTERN deck, stopping before a field-query Gaussian run."""
    folder.mkdir(parents=True, exist_ok=False)
    shutil.copy2(restart, folder / "input.ncrst")
    wrapper = folder / "capture.sh"
    shutil.copy2(BASELINE / "gau_job.tpl", folder / "gau_job.tpl")
    wrapper.write_text(
        "#!/bin/bash\nset -euo pipefail\n/bin/cp gau_job.inp old.gau_job.inp\nexit 99\n",
        encoding="utf-8",
    )
    wrapper.chmod(0o755)
    (folder / "minimize.mdin").write_text(
        make_mdin(template, "EXTERN", "./capture.sh"), encoding="utf-8"
    )
    command = [
        str(SANDER),
        "-O",
        "-i",
        "minimize.mdin",
        "-p",
        str(TOPOLOGY),
        "-c",
        "input.ncrst",
        "-o",
        "single_point.mdout",
        "-r",
        "single_point.ncrst",
        "-inf",
        "mdinfo",
    ]
    journal(
        "intent",
        f"Export Sander-adjusted EXTERN point charges at {folder}; stop before Gaussian field evaluation",
    )
    with (folder / "launcher.log").open("x") as log:
        subprocess.run(
            command,
            cwd=folder,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    deck = folder / "old.gau_job.inp"
    if not deck.is_file():
        raise FileNotFoundError(deck)
    parse_external_input(deck)
    write_json(
        folder / "process.json",
        {
            "status": "external_input_captured",
            "expected_stopped_before_g16": True,
            "command": command,
            "input_sha256": digest(folder / "minimize.mdin"),
            "external_input": deck.name,
            "external_input_sha256": digest(deck),
            "sander_log": "single_point.mdout",
        },
    )
    journal(
        "outcome",
        f"Export Sander-adjusted EXTERN point charges at {folder}",
        "captured before G16 field query",
    )
    return deck
