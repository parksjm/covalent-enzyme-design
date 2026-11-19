#!/usr/bin/env python3
"""Prepare and evaluate matched coupled-path Hamiltonian single points."""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.coupled_ti1_drive import coupled_geometry  # noqa: E402
from nyl12_covalent.enzyme_path_method_benchmark import (  # noqa: E402
    _mdin,
    sha256,
    write_xtb_embedding,
)

LOGGER = logging.getLogger(__name__)
METHODS = {"dftb3_mm": "DFTB3", "pm6_dhplus_mm": "PM6-DH+", "wb97xd_mm": "EXTERN"}


def prepare(project: Path, output: Path, topology: Path) -> None:
    """Copy six provenance-pinned geometries and prepare unbiased evaluations."""
    campaign = project / "qmmm/pa6_resp_gaff2_wt"
    paths = campaign / "11_frame1603_12A_reaction_path_20261109/reaction_paths"
    asynchronous = paths / "02_asynchronous_es_to_ti1"
    attack = asynchronous / "03_corrected_attack_from_window24/calculation/windows"
    sources = {
        "01_es": campaign
        / "10_frame1603_nearpass_qmmm_20261106/shell_12A_attempt1/calculation/dftb3_mm/optimized.ncrst",
        "02_direct_pt": paths
        / "01_proton_transfer/01_coordinate_drive/pt_19/optimized.ncrst",
        "03_coupled_pt": asynchronous
        / "01_coupled_coordinate_drive/asynchronous_24/optimized.ncrst",
        "04_attack": attack / "attack_03/optimized.ncrst",
        "05_late_attack": attack / "attack_06/optimized.ncrst",
        "06_ti1_like": attack / "attack_09/optimized.ncrst",
    }
    records = []
    output.mkdir(parents=True, exist_ok=False)
    for name, source in sources.items():
        geometry = coupled_geometry(topology, source)
        expected_n_hg1 = name != "01_es"
        if any(geometry[f"N_{h}_A"] > 1.25 for h in ("H1", "H2")):
            raise ValueError(f"Invalid terminal proton identity: {name}")
        if expected_n_hg1 and geometry["N_HG1_A"] > 1.25:
            raise ValueError(f"HG1 not on N: {name}")
        if not expected_n_hg1 and geometry["OG1_HG1_A"] > 1.25:
            raise ValueError("ES does not retain OG1-HG1")
        folder = output / "windows" / name
        folder.mkdir(parents=True)
        shutil.copy2(source, folder / "input.ncrst")
        for method, theory in METHODS.items():
            directory = folder / method
            directory.mkdir()
            (directory / "minimize.mdin").write_text(_mdin(theory))
        (folder / "gfn2_xtb").mkdir()
        (folder / "wb97xd_mm/gau_job.tpl").write_text(
            "#P wB97XD/def2SVP Int=UltraFine SCF=(XQC,Tight,MaxCycle=512)\n\n"
        )
        records.append(
            {
                "window": name,
                "source": str(source),
                "sha256": sha256(source),
                "geometry_A": geometry,
            }
        )
    manifest = {
        "records": records,
        "topology": str(topology),
        "topology_sha256": sha256(topology),
        "methods": list(METHODS) + ["gfn2_xtb"],
        "interpretation": "Fixed geometry diagnostic; no stationary-point claim. Compare method double differences, not raw finite-cutoff profiles.",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def run(output: Path, scratch: Path) -> None:
    """Evaluate each state sequentially and retain engine logs and energy table."""
    manifest = json.loads((output / "manifest.json").read_text())
    topology = Path(manifest["topology"])
    if sha256(topology) != manifest["topology_sha256"]:
        raise ValueError("Topology changed")
    amber = Path("/data/jpc/code/ambertools26/venv")
    xtb = "/data/jpc/code/enzyme-substrate-modeling/.pixi/envs/analysis/bin/xtb"
    environment = dict(
        os.environ,
        AMBERHOME=str(amber),
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        LD_LIBRARY_PATH=str(amber / "lib"),
    )
    rows = []
    for record in manifest["records"]:
        folder = output / "windows" / record["window"]
        restart = folder / "input.ncrst"
        if sha256(restart) != record["sha256"]:
            raise ValueError("Restart changed")
        work_scratch = scratch / record["window"]
        work_scratch.mkdir(parents=True, exist_ok=False)
        environment["GAUSS_SCRDIR"] = str(work_scratch)
        row = {"state": record["window"], "energies_kcal_mol": {}}
        for method in METHODS:
            directory = folder / method
            with (directory / "launcher.log").open("x") as stream:
                subprocess.run(
                    [
                        str(amber / "bin/sander"),
                        "-O",
                        "-i",
                        "minimize.mdin",
                        "-p",
                        str(topology),
                        "-c",
                        str(restart),
                        "-o",
                        "single_point.mdout",
                        "-r",
                        "single_point.ncrst",
                        "-inf",
                        "mdinfo",
                    ],
                    cwd=directory,
                    env=environment,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            text = (directory / "single_point.mdout").read_text()
            field = {
                "dftb3_mm": "DFTBESCF",
                "pm6_dhplus_mm": "PM6ESCF",
                "wb97xd_mm": "EXTERNESCF",
            }[method]
            matches = re.findall(field + r"\s*=\s*([-+\d.EeDd]+)", text)
            if not matches or "Energy calculation finished" not in text:
                raise ValueError(f"Incomplete {method} for {record['window']}")
            row["energies_kcal_mol"][method] = float(matches[-1].replace("D", "E"))
        directory = folder / "gfn2_xtb"
        write_xtb_embedding(folder / "wb97xd_mm/old.gau_job.inp", directory)
        with (directory / "gfn2.out").open("x") as stream:
            subprocess.run(
                [
                    xtb,
                    "qm_link.xyz",
                    "--gfn",
                    "2",
                    "--chrg",
                    "-2",
                    "--uhf",
                    "0",
                    "--input",
                    "xtb.inp",
                    "--namespace",
                    "gfn2",
                    "--grad",
                ],
                cwd=directory,
                env=environment,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=True,
            )
        text = (directory / "gfn2.out").read_text()
        matches = re.findall(r"\| TOTAL ENERGY\s+([-+\d.Ee]+)\s+Eh", text)
        if not matches or "normal termination of xtb" not in text:
            raise ValueError("Incomplete xTB calculation")
        row["energies_kcal_mol"]["gfn2_xtb"] = float(matches[-1]) * 627.5094740631
        rows.append(row)
        (output / "energies.json").write_text(json.dumps(rows, indent=2) + "\n")
        LOGGER.info("Completed all four methods for %s", record["window"])


def main() -> None:
    """Dispatch preparation or persistent execution."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["prepare", "run"])
    parser.add_argument("--project", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--topology", type=Path)
    parser.add_argument("--scratch", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.action == "prepare":
        prepare(args.project, args.output, args.topology)
    else:
        run(args.output, args.scratch)


if __name__ == "__main__":
    main()
