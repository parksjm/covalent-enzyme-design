"""Prepare three matched WT states without minimization or dynamics."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
from parmed.amber import AmberParm

from nyl12_covalent.wt_matched_support import (
    BETA,
    digest,
    distant_sodium,
    restore_and_check,
    write_full_pdb,
)

STATES = ("nh2_asp272_asp274", "nh2_asp272_ash274", "nh3_asp272_asp274")
RELEASE = Path(
    "/data/jpc/polyamide-amber-ff/parameters/releases/pa6_l4_modular_resp_gaff2_v1_20261031"
)
AMBER = Path("/data/jpc/code/ambertools26/venv")


def prepare(root: Path, output: Path) -> None:
    """Build and audit matched states from immutable Model30 parent coordinates."""
    topology = root / "qmmm/pa6_flipped_model30/00_system/md_solvated.prmtop"
    restart = root / "qmmm/pa6_flipped_model30/01_mm_minimization/min.rst7"
    if output.exists() and set(p.name for p in output.iterdir()) != {
        "preparation_history"
    }:
        raise FileExistsError(output)
    output.mkdir(parents=True, exist_ok=True)
    frozen = output / "source_inputs"
    frozen.mkdir()
    sources = [
        topology,
        restart,
        Path(__file__),
        Path(__file__).with_name("wt_matched_support.py"),
        root / "scripts/prepare_wt_protonation_md.py",
    ]
    hashes = {}
    for source in sources:
        shutil.copy2(source, frozen / source.name)
        hashes[str(source)] = digest(source)
    params = output / "parameters"
    shutil.copytree(RELEASE, params)
    hashes.update({str(p): digest(p) for p in RELEASE.iterdir() if p.is_file()})
    parent = AmberParm(str(topology), xyz=str(restart))
    removed = distant_sodium(parent)
    manifest = {
        "sources_sha256": hashes,
        "removed_sodium_parent_residues": removed,
        "neutral_recipe": {
            "N": -0.4834,
            "H1": 0.1224,
            "H2": 0.1224,
            "declaration": "Custom whole-neutral NTHR-derived recipe; N3 and other terms retained; N charge +0.001e versus historical recipe",
        },
        "states": {},
    }
    (output / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for state in STATES:
        charged = state.startswith("nh3")
        ash = state.endswith("ash274")
        folder = output / state
        folder.mkdir()
        mapping = write_full_pdb(
            parent, folder / "input_full.pdb", removed if charged or ash else [], ash
        )
        edits = (
            ""
            if charged
            else "\n".join(
                f"remove complex complex.{mapping.index(r) + 1}.H3\n"
                f"set complex.{mapping.index(r) + 1}.N charge -0.4834\n"
                f"set complex.{mapping.index(r) + 1}.H1 charge 0.1224\n"
                f"set complex.{mapping.index(r) + 1}.H2 charge 0.1224"
                for r in BETA
            )
        )
        text = (
            "source leaprc.protein.ff19SB\nsource leaprc.water.tip3p\n"
            f"source {params}/leaprc.pa6_shared_p6i\n"
            "complex = loadpdb input_full.pdb\n" + edits + "\n"
            "check complex\ncharge complex\n"
            "saveamberparm complex leap.prmtop leap.rst7\nquit\n"
        )
        (folder / "prepare.in").write_text(text)
        env = os.environ.copy()
        env["AMBERHOME"] = str(AMBER)
        env["LD_LIBRARY_PATH"] = (
            str(AMBER / "lib") + ":" + env.get("LD_LIBRARY_PATH", "")
        )
        with (folder / "tleap.log").open("w") as log:
            subprocess.run(
                [str(AMBER / "bin/tleap"), "-I", str(params), "-f", "prepare.in"],
                cwd=folder,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        log = (folder / "tleap.log").read_text()
        if "FATAL" in log:
            raise ValueError(f"LEaP failed: {folder}")
        product = AmberParm(str(folder / "leap.prmtop"), xyz=str(folder / "leap.rst7"))
        evidence = restore_and_check(parent, product, mapping, charged, ash)
        product.box = np.array(parent.box, copy=True)
        product.parm_data["POINTERS"][27] = parent.parm_data["POINTERS"][27]
        product.pointers["IFBOX"] = parent.ptr("ifbox")
        product.save(str(folder / "topology.prmtop"), overwrite=False)
        product.save(str(folder / "start.rst7"), format="ncrst", overwrite=False)
        check = AmberParm(
            str(folder / "topology.prmtop"), xyz=str(folder / "start.rst7")
        )
        if check.ptr("ifbox") != parent.ptr("ifbox") or not np.array_equal(
            check.box, parent.box
        ):
            raise ValueError("Periodic box bookkeeping changed")
        # ASCII restart precision must preserve all parent values at native precision.
        checked_xyz = check.coordinates
        parent_xyz = parent.coordinates
        for residue, original in zip(check.residues, mapping, strict=True):
            old = {a.name: a for a in parent.residues[original - 1].atoms}
            for atom in residue.atoms:
                if atom.name in old and not np.array_equal(
                    checked_xyz[atom.idx], parent_xyz[old[atom.name].idx]
                ):
                    raise ValueError(
                        "Shared coordinate changed after restart serialization"
                    )
        evidence.update(
            process_complete=True,
            preparation_accepted=True,
            minimization_claimed=False,
            state=state,
            topology_sha256=digest(folder / "topology.prmtop"),
            coordinates_sha256=digest(folder / "start.rst7"),
            residue_parent_mapping=mapping,
            parent_ifbox=parent.ptr("ifbox"),
            box=parent.box.tolist(),
        )
        (folder / "preparation_acceptance.json").write_text(
            json.dumps(evidence, indent=2) + "\n"
        )
        manifest["states"][state] = evidence
        logging.info("Prepared and accepted %s", state)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
