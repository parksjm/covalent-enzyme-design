"""Saved-gradient projections; no electronic-structure calls or optimization."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from pathlib import Path
import re
import sys

import numpy as np
import parmed
from parmed.amber.mask import AmberMask

HERE = Path(__file__).resolve().parent
REFERENCE = HERE.parent
ROOT = HERE.parents[5]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.proton_coordination import coordination_transfer  # noqa: E402

BOHR = 0.529177210903
HARTREE = 627.5094740631
NATIVE_SCALE = -HARTREE / 0.529177**2
LOGGER = logging.getLogger(__name__)


def fchk_array(path: Path, label: str) -> np.ndarray:
    """Read a declared numeric array from a Gaussian formatted checkpoint."""
    with path.open() as stream:
        for line in stream:
            if line.startswith(label) and "N=" in line:
                count = int(line.split("N=")[1])
                values = []
                while len(values) < count:
                    values.extend(
                        float(x.replace("D", "E")) for x in next(stream).split()
                    )
                if len(values) != count:
                    raise ValueError("fchk array length mismatch")
                return np.asarray(values)
    raise ValueError(f"Missing fchk array {label}")


def xyz(path: Path) -> tuple[list[str], np.ndarray]:
    """Read one XYZ frame, requiring exact declared row count."""
    lines = path.read_text().splitlines()
    n = int(lines[0])
    rows = [line.split() for line in lines[2:] if line.strip()]
    if len(rows) != n:
        raise ValueError(f"Expected exactly one frame in {path}")
    return [r[0] for r in rows], np.array([r[1:4] for r in rows], float)


def deck(path: Path) -> tuple[list[str], np.ndarray]:
    """Read the fixed 110 QM/link atoms in original Sander order."""
    lines = path.read_text().splitlines()
    i = lines.index("-1 1") + 1
    rows = []
    while lines[i].strip():
        rows.append(lines[i].split())
        i += 1
    if len(rows) != 110:
        raise ValueError("Expected 110 QM/link atoms")
    return [r[0] for r in rows], np.array([r[1:4] for r in rows], float)


def engrad(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read xTB energy-gradient format (Eh/bohr; coordinates bohr)."""
    rows = [
        x.strip()
        for x in path.read_text().splitlines()
        if x.strip() and not x.startswith("#")
    ]
    n = int(rows[0])
    gradient = np.array(rows[2 : 2 + 3 * n], float).reshape(n, 3)
    atoms = np.array([r.split() for r in rows[2 + 3 * n :]], float)
    if atoms.shape != (n, 4):
        raise ValueError("Malformed engrad coordinate section")
    return atoms[:, 0].astype(int), atoms[:, 1:] * BOHR, gradient


def distance(x: np.ndarray, a: int, b: int) -> tuple[float, np.ndarray]:
    """Return a distance in Å and its Cartesian derivative."""
    delta = x[a] - x[b]
    r = float(np.linalg.norm(delta))
    g = np.zeros_like(x)
    g[a], g[b] = delta / r, -delta / r
    return r, g


def dual_directions(rows: np.ndarray) -> np.ndarray:
    """Return minimum-norm dx/dCV directions with J D.T = identity."""
    j = rows.reshape(len(rows), -1)
    if np.linalg.matrix_rank(j) != len(rows):
        raise ValueError("Dependent CV gradients")
    return np.linalg.solve(j @ j.T, j).reshape(rows.shape)


def sha(path: Path) -> str:
    """Return SHA-256 for retained immutable input evidence."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    """Verify shared geometries and write projections, arrays and provenance."""
    topology = REFERENCE.parent.parent / "00_system/ash274_solvated.prmtop"
    # reference_checks -> grid -> comparison; system belongs to comparison.
    topology = REFERENCE.parent.parent / "00_system/ash274_solvated.prmtop"
    system = parmed.load_file(str(topology))
    membership = REFERENCE.parent / "raw_mm_membership.json"
    ids = json.loads(membership.read_text())["qm_real_atom_ids"]
    index = {v: i for i, v in enumerate(ids)}
    cut_qm = sorted(
        {
            a.idx + 1
            for a in system.atoms
            if a.idx + 1 in index
            and any(b.idx + 1 not in index for b in a.bond_partners)
        }
    )
    all_rows, verification, arrays = [], [], {}
    hashes = {
        str(p): sha(p)
        for p in (
            Path(__file__),
            topology,
            membership,
            ROOT / "src/nyl12_covalent/proton_coordination.py",
        )
    }
    for state in ("01_es", "02_activated", "03_attack"):
        folder = REFERENCE / "evaluations" / state
        source = (
            REFERENCE / "es_seed/nonperiodic_full_pc/old.gau_job.inp"
            if state == "01_es"
            else folder / "extern_export/old.gau_job.inp"
        )
        fchk = folder / "gaussian_direct/direct_force.fchk"
        xtb = folder / "gfn2_xtb/gfn2.qm_link.engrad"
        symbols, x = deck(source)
        dft_x = fchk_array(fchk, "Current cartesian coordinates").reshape(-1, 3) * BOHR
        dft_z = fchk_array(fchk, "Atomic numbers").astype(int)
        xtb_z, xtb_x, xtb_g = engrad(xtb)
        dft_g = fchk_array(fchk, "Cartesian Gradient").reshape(-1, 3)
        expected_z = np.array(
            [system.atoms[i - 1].atomic_number for i in ids] + [1] * 8
        )
        if not np.array_equal(dft_z, expected_z) or not np.array_equal(
            xtb_z, expected_z
        ):
            raise ValueError("QM/link atomic order differs")
        errors = {
            "DFT_export_max_A": float(abs(dft_x - x).max()),
            "xTB_export_max_A": float(abs(xtb_x - x).max()),
        }
        if max(errors.values()) > 1e-6:
            raise ValueError(f"Geometry mismatch: {errors}")
        hydrogens = np.flatnonzero(expected_z[:102] == 1)
        if len(hydrogens) != 54:
            raise ValueError("Expected 54 real QM hydrogen atoms")
        real_x = x[:102]
        n, o, h, c = [index[i] for i in (13905, 13916, 13917, 19920)]
        occupancy = coordination_transfer(real_x, o, n, hydrogens)
        attack, attack_j = distance(real_x, o, c)
        oh, oh_j = distance(real_x, o, h)
        nh, nh_j = distance(real_x, n, h)
        cvs = {
            "q_occ": (occupancy.value, occupancy.gradient),
            "r_OH": (oh, oh_j),
            "r_NH": (nh, nh_j),
        }
        mdin = folder / "dftb3/minimize.mdin"
        mask = re.search(r"active='([^']+)'", mdin.read_text())[1]
        active_ids = np.flatnonzero(AmberMask(system, mask).Selection()) + 1
        native_path = folder / "dftb3/path_force.xyz"
        native_xyz = folder / "dftb3/path.xyz"
        labels, native_raw = xyz(native_path)
        geometry_labels, native_x = xyz(native_xyz)
        active_index = {v: i for i, v in enumerate(active_ids)}
        real_active = [active_index[v] for v in ids]
        if len(active_ids) != 4662 or labels != geometry_labels:
            raise ValueError("Native active atom order/count mismatch")
        if labels != [system.atoms[i - 1].element_name for i in active_ids]:
            raise ValueError("Native labels differ from topology")
        errors["native_real_export_max_A"] = float(
            abs(native_x[real_active] - real_x).max()
        )
        if errors["native_real_export_max_A"] > 1e-5:
            raise ValueError("Native geometry differs")
        gradients = {
            "DFT_fixed_PC_QM": dft_g[:102] * HARTREE / BOHR,
            "GFN2_fixed_PC_QM": xtb_g[:102] * HARTREE / BOHR,
            "DFTB3_native_TOTAL": native_raw[real_active] * NATIVE_SCALE,
        }
        for cv, (value, cv_j) in cvs.items():
            directions = dual_directions(np.array([cv_j, attack_j]))
            # Boundary QM atoms and all external charges stay fixed: links do too.
            if np.any(directions[:, [index[i] for i in cut_qm]] != 0):
                raise ValueError("Moving cut parent requires explicit link chain rule")
            arrays[f"{state}_{cv}_directions_A_per_CV"] = directions
            arrays[f"{state}_{cv}_jacobian"] = np.array([cv_j, attack_j])
            for axis, name in enumerate((cv, f"attack_hold_{cv}")):
                for method, g in gradients.items():
                    all_rows.append(
                        {
                            "state": state,
                            "method": method,
                            "direction": name,
                            "derivative": float(np.sum(g * directions[axis])),
                            "units": "kcal/mol per dimensionless q"
                            if name == "q_occ"
                            else "kcal/mol/A",
                            "q_occ": occupancy.value,
                            "attack_A": attack,
                            "r_OH_A": oh,
                            "r_NH_A": nh,
                        }
                    )
        arrays[f"{state}_coordinates_A"] = x
        arrays[f"{state}_DFT_gradient_Eh_bohr"] = dft_g
        arrays[f"{state}_GFN2_gradient_Eh_bohr"] = xtb_g
        arrays[f"{state}_native_raw_path_force"] = native_raw
        arrays[f"{state}_native_total_gradient_kcal_A"] = native_raw * NATIVE_SCALE
        verification.append({"state": state, **errors, "real_H": 54})
        for p in (source, fchk, xtb, mdin, native_path, native_xyz):
            hashes[str(p)] = sha(p)
    arrays["real_QM_amber_ids"] = np.array(ids)
    arrays["native_active_amber_ids"] = active_ids
    np.savez_compressed(HERE / "arrays.npz", **arrays)
    result = {
        "verification": verification,
        "cut_QM_amber_ids": cut_qm,
        "sources_sha256": hashes,
        "projections": all_rows,
        "conventions": {
            "q_occ": "CN(N,H)-CN(OG1,H), all 54 real H; r0=1.25 A",
            "directions": "D=(J J^T)^-1 J; J D^T=I; Euclidean Cartesian minimum norm",
            "positive": "positive derivative means energy rises with increasing CV",
            "DFT_xTB": "embedded QM/link energy gradients; external PCs and cut/link atoms fixed",
            "native": "total QM/MM gradient; finite nonperiodic all-PC reference, not original periodic path",
            "native_conversion": NATIVE_SCALE,
            "Eh_bohr_to_kcal_A": HARTREE / BOHR,
        },
        "limitations": [
            "DFTB3 total QM/MM derivatives cannot be ranked against isolated embedded DFT/xTB derivatives.",
            "No stationary minima, transition states, barriers or PMFs follow from these single-point derivatives.",
            "q_occ spans all real H and can mix distinct proton exchange mechanisms.",
        ],
    }
    (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    with (HERE / "projections.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    LOGGER.info(
        "Saved %d derivative rows; all coordinate/mapping gates passed", len(all_rows)
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
