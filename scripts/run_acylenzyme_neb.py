"""Relax an atom-matched acylenzyme endpoint and run native DL-FIND NEB."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
from parmed.amber import Rst7

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.unrestrained_relay import native_converged  # noqa: E402

LOGGER = logging.getLogger(__name__)
PAIRS = {
    "OG1_C": (13916, 19916),
    "C_O": (19916, 19917),
    "C_Nleave": (19916, 19929),
    "Nleave_HG1": (19929, 13917),
    "N_H1": (13905, 13906),
    "N_H2": (13905, 13907),
    "N_HG1": (13905, 13917),
    "OG1_HG1": (13916, 13917),
}


def write_json(path: Path, data: dict) -> None:
    """Write a finite JSON process or scientific record."""
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def geometry(path: Path) -> dict[str, float]:
    """Measure chemically diagnostic distances in the original atom mapping."""
    xyz = Rst7.open(str(path)).coordinates.reshape(-1, 3)
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite restart coordinates")
    return {
        k: float(np.linalg.norm(xyz[a - 1] - xyz[b - 1])) for k, (a, b) in PAIRS.items()
    }


def acyl_retained(g: dict[str, float]) -> bool:
    """Return whether ester and neutral leaving-amine geometry is retained."""
    return (
        1.2 < g["OG1_C"] < 1.75
        and 1.1 < g["C_O"] < 1.35
        and g["C_Nleave"] > 2.0
        and 0.8 < g["Nleave_HG1"] < 1.3
        and all(0.8 < g[k] < 1.3 for k in ("N_H1", "N_H2"))
        and g["N_HG1"] > 1.4
        and g["OG1_HG1"] > 1.4
    )


def input_text(
    base: str,
    reference: Path,
    *,
    restrained: bool = False,
    images: list[Path] | None = None,
) -> str:
    """Create fixed-reference endpoint or NEB input without positional bias."""
    if "nmropt=0" not in base.replace(" ", "") or "DISANG" in base:
        raise ValueError("Expected unbiased base input")
    text = re.sub(r"maxcyc\s*=\s*\d+", "maxcyc=1000", base)
    text = re.sub(r"maxcycle\s*=\s*\d+", "maxcycle=1000", text)
    text = re.sub(r"ntwx\s*=\s*\d+", "ntwx=5", text)
    text = re.sub(r"maxstep\s*=\s*[\d.]+", "maxstep=0.2", text)
    text = text.replace("&dlfind", f"&dlfind\n refcrd='{reference}',")
    if restrained:
        text = text.replace("nmropt=0", "nmropt=1")
        text = text.replace(
            "&dlfind", "&wt type='END' /\nDISANG=restraints.RST\n&dlfind"
        )
    if images is not None:
        if not 1 <= len(images) <= 12:
            raise ValueError("Native DL-FIND supports at most 12 seed frames")
        controls = "\n neb='FREE', neb_nimage=17, nebk=0.01,\n"
        controls += " neb_climb_test=3.0, neb_freeze_test=0.0,\n"
        controls += " wtmask='@1621-1655,2727-2734,13905-13931,14555-14560,14574-14579,19914-19917,19925-19931,19937-19941',\n"
        controls += "".join(f" coords2({i})='{p}',\n" for i, p in enumerate(images, 1))
        text = text.replace("&dlfind", "&dlfind" + controls)
        text = text.replace("verbosity=1", "verbosity=2")
    return text


def restraints() -> str:
    """Return temporary endpoint-building restraints, never NEB restraints."""
    targets = {
        "OG1_C": 1.43,
        "C_O": 1.23,
        "C_Nleave": 2.8,
        "Nleave_HG1": 1.03,
        "N_H1": 1.03,
        "N_H2": 1.03,
    }
    return "".join(
        f"&rst iat={PAIRS[k][0]},{PAIRS[k][1]}, r1={v - 0.5:.3f}, "
        f"r2={v - 0.05:.3f}, r3={v + 0.05:.3f}, r4={v + 0.5:.3f}, "
        "rk2=100, rk3=100, /\n"
        for k, v in targets.items()
    )


def run_engine(config: dict, folder: Path, start: Path, text: str) -> Path:
    """Run Sander with retained commands and separate process status."""
    folder.mkdir(exist_ok=False)
    (folder / "minimize.mdin").write_text(text)
    if "DISANG" in text:
        (folder / "restraints.RST").write_text(restraints())
    result = folder / "optimized.ncrst"
    cmd = [
        config["sander"],
        "-O",
        "-i",
        "minimize.mdin",
        "-p",
        config["topology"],
        "-c",
        str(start),
        "-o",
        "minimized.mdout",
        "-r",
        str(result),
        "-x",
        "optimization.nc",
        "-inf",
        "mdinfo",
    ]
    record = {
        "command": cmd,
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "start_sha256": hashlib.sha256(start.read_bytes()).hexdigest(),
    }
    write_json(folder / "process.json", record)
    amber = Path(config["sander"]).parent.parent
    env = {
        **os.environ,
        "AMBERHOME": str(amber),
        "LD_LIBRARY_PATH": str(amber / "lib"),
        "GFORTRAN_UNBUFFERED_ALL": "1",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    with (folder / "engine.log").open("w") as log:
        run = subprocess.run(
            cmd, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, check=False
        )
    record.update(
        exit_status=run.returncode, ended_utc=datetime.now(timezone.utc).isoformat()
    )
    write_json(folder / "process.json", record)
    if run.returncode or "SANDER BOMB" in (folder / "minimized.mdout").read_text():
        raise RuntimeError(f"Sander failed: {folder}")
    return result


def run(config_path: Path) -> None:
    """Prepare, release and chemically gate the product before free-end NEB."""
    cfg = json.loads(config_path.read_text())
    root = config_path.parent
    base = Path(cfg["base_input"]).read_text()
    reference = Path(cfg["reactant"])
    state = {"phase": "endpoint_preparation", "scientific_acceptance": False}
    write_json(root / "status.json", state)
    try:
        endpoint = Path(cfg["seed"])
        for name, restrained in (
            ("01_endpoint_preparation", True),
            ("02_endpoint_release", False),
        ):
            state["phase"] = name
            write_json(root / "status.json", state)
            folder = root / name
            endpoint = run_engine(
                cfg,
                folder,
                endpoint,
                input_text(base, reference, restrained=restrained),
            )
            g = geometry(endpoint)
            accepted = acyl_retained(g)
            write_json(
                folder / "acceptance.json",
                {
                    "geometry_A": g,
                    "acylenzyme_retained": accepted,
                    "native_converged": native_converged(folder / "minimized.mdout"),
                    "reaction_restraints": restrained,
                    "stationary_point_claim": False,
                },
            )
            if not accepted:
                raise RuntimeError(f"Acylenzyme identity lost in {name}: {g}")
        state["phase"] = "03_neb"
        state["endpoint_policy"] = (
            "FREE: endpoints also optimize; release convergence recorded separately"
        )
        write_json(root / "status.json", state)
        images = [Path(p) for p in cfg["intermediates"]] + [endpoint]
        run_engine(
            cfg, root / "03_neb", reference, input_text(base, reference, images=images)
        )
        state["phase"] = "finished"
        state["native_converged"] = native_converged(root / "03_neb/minimized.mdout")
        write_json(root / "status.json", state)
    except Exception as exc:
        state.update(phase="stopped", reason=str(exc))
        write_json(root / "status.json", state)
        raise


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    run(parser.parse_args().config.resolve())
