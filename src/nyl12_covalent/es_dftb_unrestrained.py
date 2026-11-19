"""Prepare a source-audited, fully unrestrained DFTB3/MM E-S optimization."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from nyl12_covalent.es_qm_region import audit_es_qm_region
from nyl12_covalent.es_system_preparation import sha256


def compact_mask(indices: list[int]) -> str:
    """Return a deterministic Amber atom mask from sorted one-based indices."""
    values = sorted(set(indices))
    if not values or len(values) != len(indices):
        raise ValueError("QM indices must be nonempty and unique")
    spans: list[str] = []
    start = end = values[0]
    for value in values[1:]:
        if value == end + 1:
            end = value
            continue
        spans.append(str(start) if start == end else f"{start}-{end}")
        start = end = value
    spans.append(str(start) if start == end else f"{start}-{end}")
    return "@" + ",".join(spans)


def dftb_input(
    mask: str,
    max_cycles: int = 10_000,
    active_mask: str = "@*",
    case_label: str = "Model 23",
) -> str:
    """Return the unrestrained periodic DFTB3/MM DL-FIND input."""
    if (
        not mask.startswith("@")
        or max_cycles < 1
        or not active_mask.strip()
        or not case_label.strip()
    ):
        raise ValueError("An Amber atom mask and positive max_cycles are required")
    return f"""{case_label} noncovalent E-S DFTB3/MM optimization
&cntrl
 nmropt=0,
 imin=1, maxcyc={max_cycles}, ntmin=5, ncyc=10, dx0=0.02, drms=0.05,
 ntb=1, igb=0, ntc=1, ntf=1, cut=10.0, ifqnt=1, jfastw=4,
 ntr=0,
 ntpr=25, ntwx=25, ntwr=100, ioutfm=1, ntxo=2,
/
&qmmm
 qmmask='{mask}', qmcharge=-2, qm_theory='DFTB3',
 qmmm_int=1, adjust_q=2, qm_ewald=1, qm_pme=1, qmcut=10.0,
 scfconv=1.0d-10, lnk_dis=1.09, lnk_atomic_no=1, lnk_method=1,
/
&dlfind
 optalg='LBFGS', crdrep='CART', active='{active_mask}', tol=1.d-4, tole=1.d-5,
 maxcycle={max_cycles}, maxstep=0.5, verbosity=1,
/
"""


def prepare_dftb_optimization(
    region_config: Path,
    sander: Path,
    output: Path,
    active_mask: str = "@*",
    case_label: str = "Model 23",
) -> dict[str, Any]:
    """Audit the Model 23 region and write source-pinned DFTB3/MM inputs."""
    if not region_config.is_file() or not sander.is_file():
        raise FileNotFoundError("Region configuration or Sander is missing")
    payload = json.loads(region_config.read_text())
    variants = [
        row for row in payload["variants"] if row["name"] == "primary_106_nuclei"
    ]
    if len(variants) != 1:
        raise ValueError("Exactly one approved primary QM region is required")
    topology = Path(payload["topology"])
    coordinates = Path(payload["coordinates"])
    audit = audit_es_qm_region(topology, coordinates, variants[0])
    if (
        audit["real_atom_count"],
        audit["link_H_count"],
        audit["total_QM_nuclei"],
        audit["qmcharge"],
    ) != (98, 8, 106, -2):
        raise ValueError("Expected the reviewed 98-real/8-link, charge -2 QM region")
    output.mkdir(parents=True, exist_ok=False)
    mask = compact_mask([row["index"] for row in audit["qm_atoms"]])
    mdin = output / "minimize.mdin"
    mdin.write_text(dftb_input(mask, active_mask=active_mask, case_label=case_label))
    command = [
        str(sander.resolve()),
        "-O",
        "-i",
        str(mdin.resolve()),
        "-p",
        str(topology.resolve()),
        "-c",
        str(coordinates.resolve()),
        "-o",
        str((output / "minimized.mdout").resolve()),
        "-r",
        str((output / "optimized.ncrst").resolve()),
        "-x",
        str((output / "optimization.nc").resolve()),
        "-inf",
        str((output / "mdinfo").resolve()),
    ]
    sources = [
        region_config,
        topology,
        coordinates,
        sander,
        Path(__file__).resolve(),
        mdin,
    ]
    metadata = {
        "case_label": case_label,
        "active_mask": active_mask,
        "command": command,
        "environment": {
            "AMBERHOME": str(sander.resolve().parent.parent),
            "LD_LIBRARY_PATH": str(sander.resolve().parent.parent / "lib"),
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        },
        "restraints": "none; ntr=0, nmropt=0; distant inactive atoms are fixed by DL-FIND active mask",
        "trajectory_stride_optimizer_cycles": 25,
        "region_audit": audit,
        "sources": [
            {"path": str(path.resolve()), "sha256": sha256(path)} for path in sources
        ],
    }
    (output / "preparation.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    )
    return metadata
