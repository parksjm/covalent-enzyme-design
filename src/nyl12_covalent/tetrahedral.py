"""Build the first covalent Nyl12 tetrahedral adduct and relay-water model."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .chemistry import Substrate, product_report
from .pose import pdb_coordinates
from .topology import QM_BRIDGE, _leap_base, _normalize_library, _run, _write_ligand_pdb

# ff19SB NTHR has +1 charge. Removing H3 leaves +0.8066, so placing the
# compensating charge on N yields a neutral NH2 terminus. This provisional
# QM/MM-only charge assignment requires a RESP fit before classical MD.
NEUTRAL_BETA_N_CHARGE = -0.6254


def _unit(vector: np.ndarray) -> np.ndarray:
    """Return a normalized nonzero Cartesian vector."""
    norm = float(np.linalg.norm(vector))
    if norm < 1.0e-8:
        raise ValueError("Cannot normalize a zero-length geometry vector")
    return vector / norm


def _protein_points(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Read protein heavy atoms and the beta Thr1 catalytic site."""
    points = []
    catalytic: dict[str, np.ndarray] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[12:16].strip().startswith("H"):
            continue
        name = line[12:16].strip()
        point = np.array([float(line[i : i + 8]) for i in (30, 38, 46)])
        points.append(point)
        if line[21] == "E" and int(line[22:26]) == 1 and name in {"N", "OG1"}:
            catalytic[name] = point
    if set(catalytic) != {"N", "OG1"}:
        raise ValueError("Source protein lacks beta E1 N or OG1")
    return np.asarray(points), catalytic


def place_relay_water(
    protein_pdb: Path, ligand_pdb: Path, oxyanion_name: str, output: Path
) -> dict[str, object]:
    """Place one neutral relay water between Thr OG1 and beta N-terminus.

    This is a geometric seed. The subsequent QM/MM geometry optimization must
    decide its hydrogen bonds; the water is not evidence for a reaction path.

    Args:
        protein_pdb: Protein coordinates with chain E residue 1.
        ligand_pdb: Full substrate PDB with stable atom names.
        oxyanion_name: Former carbonyl oxygen in the TI1 state.
        output: Water PDB file.

    Returns:
        Water coordinates and nearest-neighbor diagnostics.
    """
    protein, catalytic = _protein_points(protein_pdb)
    ligand = pdb_coordinates(ligand_pdb)
    nitrogen, og1 = catalytic["N"], catalytic["OG1"]
    axis = _unit(og1 - nitrogen)
    midpoint = (nitrogen + og1) / 2.0
    half_distance = float(np.linalg.norm(og1 - nitrogen) / 2.0)
    if half_distance >= 2.8:
        raise ValueError("Catalytic Thr N/OG1 separation cannot support a relay water")
    offset = math.sqrt(2.8**2 - half_distance**2)
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, reference))) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    first = _unit(np.cross(axis, reference))
    second = _unit(np.cross(axis, first))
    other_protein = protein[
        (np.linalg.norm(protein - nitrogen, axis=1) > 0.01)
        & (np.linalg.norm(protein - og1, axis=1) > 0.01)
    ]
    other_ligand = np.asarray(list(ligand.values()))
    other = np.vstack([other_protein, other_ligand])
    oxyanion = ligand[oxyanion_name]
    best: tuple[float, np.ndarray] | None = None
    for index in range(72):
        angle = 2.0 * math.pi * index / 72
        point = midpoint + offset * (math.cos(angle) * first + math.sin(angle) * second)
        distances = np.linalg.norm(other - point, axis=1)
        clash = float(np.square(np.maximum(2.4 - distances, 0)).sum())
        hb_offset = (float(np.linalg.norm(oxyanion - point)) - 2.8) ** 2
        score = clash + 0.25 * hb_offset
        if best is None or score < best[0]:
            best = (score, point)
    if best is None:
        raise ValueError("No relay-water location was generated")
    oxygen = best[1]
    direction = _unit(oxyanion - oxygen)
    orthogonal = _unit(np.cross(direction, axis))
    bond = 0.9572
    theta = math.radians(104.52)
    h1 = oxygen + bond * direction
    h2 = oxygen + bond * (math.cos(theta) * direction + math.sin(theta) * orthogonal)
    entries = [("O", oxygen, "O"), ("H1", h1, "H"), ("H2", h2, "H")]
    lines = [
        f"HETATM{index:5d} {name:4s} WAT A   1    "
        f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}"
        f"  1.00  0.00          {element:>2s}"
        for index, (name, point, element) in enumerate(entries, start=1)
    ]
    output.write_text("\n".join([*lines, "END"]) + "\n")
    return {
        "oxygen_xyz": oxygen.tolist(),
        "hydrogen_1_xyz": h1.tolist(),
        "hydrogen_2_xyz": h2.tolist(),
        "n_to_water_o": float(np.linalg.norm(nitrogen - oxygen)),
        "og1_to_water_o": float(np.linalg.norm(og1 - oxygen)),
        "oxyanion_to_water_o": float(np.linalg.norm(oxyanion - oxygen)),
        "placement_penalty": best[0],
        "method": "72-point circle search, N/OG1 2.8 A, nonbonded clash penalty",
    }


def _poly_library(mol2: Path, frcmod: Path, output: Path, amber_bin: Path) -> Path:
    """Create a stable LIG library from the intact poly(X,AMA) MOL2."""
    library = output / "intact_poly_xama.lib"
    script = output / "intact_library.leap"
    script.write_text(
        "\n".join([
            *_leap_base(frcmod), f"lig = loadmol2 {mol2}",
            f"saveoff lig {library}", "quit",
        ]) + "\n"
    )
    _run([str(amber_bin / "tleap"), "-f", str(script)], output, "intact_library.leap.log")
    _normalize_library(library)
    return library


def _assemble(
    substrate: Substrate,
    protein_pdb: Path,
    ligand_library: Path,
    frcmod: Path,
    output: Path,
    amber_bin: Path,
    oxyanion_charge: float | None,
) -> None:
    """Build a four-coordinate carbonyl carbon with intact amide C–N."""
    bridge = output / "ti1_qm_bridge.frcmod"
    bridge.write_text(QM_BRIDGE)
    commands = [
        *_leap_base(frcmod),
        "source leaprc.water.tip3p",
        f"loadamberparams {bridge}",
        f"loadoff {ligand_library}",
        f"prot = loadpdb {protein_pdb}",
        f"lig = loadpdb {output / 'intact.pdb'}",
        f"water = loadpdb {output / 'relay_water.pdb'}",
        "remove prot prot.929.HG1",
        "set prot.929.OG1 type os",
        "set prot.929.OG1 charge -0.4650",
    ]
    for residue in (1025, 1121, 1217):
        commands.extend([
            f"remove prot prot.{residue}.H3",
            f"set prot.{residue}.N charge {NEUTRAL_BETA_N_CHARGE:.4f}",
        ])
    if oxyanion_charge is not None:
        commands.append(
            f"set lig.1.{substrate.carbonyl_o} charge {oxyanion_charge:.8f}"
        )
    commands.extend([
        "complex = combine {prot lig water}",
        f"bond complex.929.OG1 complex.1313.{substrate.carbonyl_c}",
        "check complex", "charge complex",
        f"saveamberparm complex {output / 'ti1.raw.prmtop'} {output / 'ti1.inpcrd'}",
        f"savepdb complex {output / 'ti1.pdb'}", "quit",
    ])
    script = output / "ti1_assembly.leap"
    script.write_text("\n".join(commands) + "\n")
    _run([str(amber_bin / "tleap"), "-f", str(script)], output, "ti1_assembly.leap.log")
    log = (output / "ti1_assembly.leap.log").read_text()
    if "Errors = 0" not in log or not (output / "ti1.raw.prmtop").stat().st_size:
        raise ValueError("LEaP did not produce a fully parameterized TI1 topology")


def build_ti1(
    substrate: Substrate,
    protein_pdb: Path,
    ligand_pose: Path,
    ligand_mol2: Path,
    frcmod: Path,
    ligand_library: Path | None,
    output: Path,
    amber_bin: Path,
) -> dict[str, object]:
    """Build and audit the product-gated first tetrahedral covalent state.

    Args:
        substrate: Full intact L4 or poly(X,AMA) oligomer graph.
        protein_pdb: Mature Nyl12 alpha/beta protein starting coordinates.
        ligand_pose: Named full-substrate bound pose.
        ligand_mol2: Full-substrate Amber MOL2 with atom types and charges.
        frcmod: Ligand GAFF parameters.
        ligand_library: PA6 LIG library, or None for a new polyamide MOL2.
        output: Dedicated TI1 directory.
        amber_bin: AmberTools26 executable directory.

    Returns:
        Product, water, and topology audit.
    """
    output.mkdir(parents=True, exist_ok=True)
    product = product_report(substrate)
    _write_ligand_pdb(ligand_pose, output / "intact.pdb")
    water = place_relay_water(
        protein_pdb, output / "intact.pdb", substrate.carbonyl_o,
        output / "relay_water.pdb",
    )
    library = ligand_library or _poly_library(ligand_mol2, frcmod, output, amber_bin)
    _assemble(substrate, protein_pdb, library, frcmod, output, amber_bin, None)
    validator = Path(__file__).with_name("validate_ti1.py")
    _run(
        [str(amber_bin / "python"), str(validator), str(output),
         substrate.carbonyl_c, substrate.carbonyl_o, substrate.leaving_n],
        output, "ti1_audit.log",
    )
    initial = json.loads((output / "ti1_audit.json").read_text())
    correction = -initial["qm_region_partial_charge"]
    _assemble(
        substrate, protein_pdb, library, frcmod, output, amber_bin,
        initial["oxyanion_partial_charge"] + correction,
    )
    _run(
        [str(amber_bin / "python"), str(validator), str(output),
         substrate.carbonyl_c, substrate.carbonyl_o, substrate.leaving_n],
        output, "ti1_audit.log",
    )
    topology = json.loads((output / "ti1_audit.json").read_text())
    if abs(topology["qm_region_partial_charge"]) > 0.002:
        raise ValueError("TI1 QM region charge failed to normalize to zero")
    report = {
        "state": "tetrahedral TI1, intact scissile amide C-N, oxyanion",
        "product_gate": product,
        "relay_water": water,
        "topology": topology,
        "charge_correction_on_qm_oxyanion": correction,
        "force_field_scope": "QM/MM only; cross-link and reactive-center parameters are not MD-ready",
        "other_beta_n_termini": (
            "F/G/H neutral NH2; ff19SB NTHR H3 removed and compensating "
            "charge assigned to N in LEaP; RESP refinement needed for MD"
        ),
    }
    (output / "ti1_provenance.json").write_text(json.dumps(report, indent=2) + "\n")
    return report
