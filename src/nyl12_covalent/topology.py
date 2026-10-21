"""Assemble assayed Nyl12 acyl-enzyme Amber topologies with provenance."""

from __future__ import annotations

import json
import hashlib
import logging
import subprocess
from pathlib import Path

import numpy as np
import rdkit
from rdkit import Chem

from .chemistry import Substrate, retained_names
from .parameter_cache import ParameterArtifactCache

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[2]

QM_BRIDGE = """QM-region-only ester junction analogs; not validated for classical MD
MASS

BOND
3C-os  308.6  1.4316  GAFF c3-os analog; both atoms in QM region

ANGLE
3C-os-c  63.3  115.98  GAFF c3-os-c analog; all QM
os-3C-CX  68.0  107.97  GAFF os-c3-c3 analog; all QM
os-3C-H1  50.8  109.78  GAFF os-c3-h1 analog; all QM
os-3C-CT  68.0  107.97  GAFF os-c3-c3 analog; all QM

DIHE
X -3C-os-X   3  1.150  0.0  3.0  GAFF X-c3-os-X analog; all QM
os-3C-CT-HC  1  0.000  0.0  3.0  Placeholder; all QM

IMPROPER

NONBON
"""


def _renumber_ligand_line(line: str, serial: int) -> str:
    """Put a selected ligand atom in residue LIG A 1 for LEaP."""
    return f"HETATM{serial:5d}{line[11:21]}A   1{line[26:]}"


def _run(command: list[str], workdir: Path, log_name: str) -> None:
    """Run an AmberTools command and preserve its complete output."""
    LOGGER.info("Running %s", " ".join(command))
    result = subprocess.run(command, cwd=workdir, text=True, capture_output=True, check=False)
    (workdir / log_name).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{' '.join(command)} failed; see {workdir / log_name}")


def _write_retained_pdb(source: Path, retained: frozenset[str], output: Path) -> None:
    """Write only the retained heavy atoms from a named bound ligand pose."""
    lines = []
    found = set()
    for line in source.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[17:20] != "LIG":
            continue
        name = line[12:16].strip()
        if (
            name in retained
            and not name.startswith("H")
            and line[76:78].strip().upper() not in {"H", "D"}
        ):
            lines.append(_renumber_ligand_line(line, len(lines) + 1))
            found.add(name)
    required = {name for name in retained if not name.startswith("H")}
    if found != required:
        raise ValueError(f"Bound pose lacks retained atoms: {sorted(required - found)}")
    output.write_text("\n".join([*lines, "END"]) + "\n")


def rename_mol2_atoms(mol: Chem.Mol, input_path: Path, output_path: Path) -> None:
    """Map Amber atom types and charges to stable graph names by coordinates.

    Args:
        mol: The exact RDKit molecule used to make the input SDF.
        input_path: Antechamber MOL2 with its assigned atom types and charges.
        output_path: MOL2 with graph-stable atom names.
    """
    lines = input_path.read_text().splitlines()
    marker = lines.index("@<TRIPOS>ATOM")
    stop = lines.index("@<TRIPOS>BOND")
    conformer = mol.GetConformer()
    source_points = np.asarray(conformer.GetPositions())
    unused = set(range(mol.GetNumAtoms()))
    for line_number in range(marker + 1, stop):
        fields = lines[line_number].split()
        if len(fields) < 9:
            continue
        point = np.array([float(value) for value in fields[2:5]])
        candidates = sorted(
            (
                (float(np.linalg.norm(source_points[index] - point)), index)
                for index in unused
                if mol.GetAtomWithIdx(index).GetSymbol().upper()
                == fields[5][0].upper()
            )
        )
        if not candidates or candidates[0][0] > 0.06:
            raise ValueError(f"Cannot map Amber atom at {point} to RDKit graph")
        index = candidates[0][1]
        unused.remove(index)
        fields[1] = mol.GetAtomWithIdx(index).GetProp("name")
        lines[line_number] = " ".join(fields)
    if unused:
        raise ValueError(f"Antechamber omitted {len(unused)} ligand atoms")
    output_path.write_text("\n".join(lines) + "\n")


def parameterize_poly_xama(
    mol: Chem.Mol, output: Path, amber_bin: Path, charge_method: str = "bcc"
) -> tuple[Path, Path]:
    """Assign GAFF types and requested provisional charges to the oligomer.

    Args:
        mol: Bound-pose oligomer with explicit hydrogens.
        output: Target directory.
        amber_bin: AmberTools26 bin directory.
        charge_method: ``bcc`` for AM1-BCC or ``gas`` for QM-only Gasteiger.

    Returns:
        Named MOL2 and GAFF frcmod paths.
    """
    sdf = output / "intact_poly_xama.sdf"
    if charge_method not in {"bcc", "gas"}:
        raise ValueError(f"Unsupported charge method {charge_method}")
    raw_mol2 = output / "antechamber_raw.mol2"
    named_mol2 = output / "intact_poly_xama.mol2"
    frcmod = output / "poly_xama.frcmod"
    tools = {name: amber_bin / name for name in ("antechamber", "parmchk2", "sqm")}
    gaff_data = amber_bin.parent / "dat" / "antechamber" / "gaff.dat"

    def sha256(path: Path) -> str:
        checksum = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                checksum.update(block)
        return checksum.hexdigest()

    identity = {
        "graph": Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True, isomericSmiles=True),
        "atom_names": [mol.GetAtomWithIdx(i).GetProp("name") for i in range(mol.GetNumAtoms())],
        "coordinates": np.asarray(mol.GetConformer().GetPositions()).round(8).tolist(),
        "state": "intact_neutral_poly_xama",
        "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        "charge_method": charge_method,
        "rdkit_version": rdkit.__version__,
        "tools": {name: sha256(path) for name, path in tools.items()},
        "gaff_data": sha256(gaff_data),
    }
    cache = ParameterArtifactCache(ROOT / ".cache" / "nyl12_covalent" / "parameters")
    artifacts = {
        "intact_poly_xama.sdf": sdf,
        "antechamber_raw.mol2": raw_mol2,
        "intact_poly_xama.mol2": named_mol2,
        "poly_xama.frcmod": frcmod,
        "antechamber.log": output / "antechamber.log",
        "parmchk2.log": output / "parmchk2.log",
    }
    with cache.calculation_lock(identity):
        if cache.materialize(identity, output) is None:
            writer = Chem.SDWriter(str(sdf))
            writer.write(mol)
            writer.close()
            _run(
                [str(tools["antechamber"]), "-i", str(sdf), "-fi", "sdf",
                 "-o", str(raw_mol2), "-fo", "mol2", "-at", "gaff",
                 "-c", charge_method, "-nc", "0", "-pf", "y"],
                output, "antechamber.log",
            )
            rename_mol2_atoms(mol, raw_mol2, named_mol2)
            _run(
                [str(tools["parmchk2"]), "-i", str(named_mol2), "-f", "mol2",
                 "-o", str(frcmod), "-s", "gaff"], output, "parmchk2.log",
            )
            cache.store(
                identity,
                artifacts,
                label=f"poly_xama_intact_neutral_{charge_method}_gaff",
            )
    return named_mol2, frcmod


def _leap_base(frcmod: Path) -> list[str]:
    """Return pinned protein, ligand, and additive parameter commands."""
    return [
        "source leaprc.protein.ff19SB",
        "source leaprc.gaff",
        "set default PBradii bondi",
        f"loadamberparams {frcmod}",
    ]


def _write_ligand_pdb(source: Path, output: Path) -> None:
    """Extract the intact ligand from a bound complex without its protein."""
    lines = [
        _renumber_ligand_line(line, index)
        for index, line in enumerate(
            (
                line for line in source.read_text().splitlines()
                if line.startswith(("ATOM  ", "HETATM")) and line[17:20] == "LIG"
            ),
            start=1,
        )
    ]
    if not lines:
        raise ValueError(f"No LIG atoms in {source}")
    output.write_text("\n".join([*lines, "END"]) + "\n")


def _normalize_library(path: Path) -> None:
    """Name a saved fragment unit LIG so LEaP recognizes its PDB residue."""
    contents = path.read_text()
    if ' "lig"' not in contents or "!entry.lig." not in contents:
        raise ValueError("Unexpected LEaP fragment library format")
    contents = contents.replace(' "lig"', ' "LIG"', 1).replace("!entry.lig.", "!entry.LIG.")
    contents = contents.replace(' "MOL" 1 ', ' "LIG" 1 ', 1)
    path.write_text(contents)


def _fragment_library(
    substrate: Substrate,
    source_mol2: Path,
    source_library: Path | None,
    frcmod: Path,
    output: Path,
    amber_bin: Path,
) -> Path:
    """Remove the departing oligomer by graph-defined atom names."""
    retained = retained_names(substrate)
    departing = [
        atom.GetProp("name")
        for atom in substrate.mol.GetAtoms()
        if atom.GetProp("name") not in retained
    ]
    commands = _leap_base(frcmod)
    if source_library is None:
        commands.append(f"lig = loadmol2 {source_mol2}")
    else:
        commands.extend([f"loadoff {source_library}", f"lig = loadpdb {output / 'intact.pdb'}"])
    commands.extend(f"remove lig lig.1.{name}" for name in departing)
    library = output / "retained.lib"
    commands.extend([f"saveoff lig {library}", "quit"])
    leap_input = output / "fragment.leap"
    leap_input.write_text("\n".join(commands) + "\n")
    _run([str(amber_bin / "tleap"), "-f", str(leap_input)], output, "fragment.leap.log")
    _normalize_library(library)
    return library


def _assemble(
    substrate: Substrate,
    protein_pdb: Path,
    frcmod: Path,
    library: Path,
    output: Path,
    amber_bin: Path,
    carbonyl_charge: float | None,
) -> None:
    """Build the Thr1 ester and neutral beta N-terminus in one active site."""
    commands = _leap_base(frcmod)
    bridge = output / "acyl_qm_bridge.frcmod"
    bridge.write_text(QM_BRIDGE)
    commands.append(f"loadamberparams {bridge}")
    commands.extend(
        [
            f"loadoff {library}",
            f"prot = loadpdb {protein_pdb}",
            f"lig = loadpdb {output / 'retained.pdb'}",
            "remove prot prot.929.H3",
            "remove prot prot.929.HG1",
            "set prot.929.N charge -0.4844",
            "set prot.929.H1 charge 0.1224",
            "set prot.929.H2 charge 0.1224",
            "set prot.929.OG1 type os",
            "set prot.929.OG1 charge -0.4650",
        ]
    )
    if carbonyl_charge is not None:
        commands.append(f"set lig.1.{substrate.carbonyl_c} charge {carbonyl_charge:.8f}")
    commands.extend(
        [
            "complex = combine {prot lig}",
            f"bond complex.929.OG1 complex.1313.{substrate.carbonyl_c}",
            "check complex",
            "charge complex",
            f"saveamberparm complex {output / 'acyl.raw.prmtop'} {output / 'acyl.inpcrd'}",
            f"savepdb complex {output / 'acyl.pdb'}",
            "quit",
        ]
    )
    script = output / "assembly.leap"
    script.write_text("\n".join(commands) + "\n")
    _run([str(amber_bin / "tleap"), "-f", str(script)], output, "assembly.leap.log")
    log = (output / "assembly.leap.log").read_text()
    if "Errors = 0" not in log or not (output / "acyl.raw.prmtop").stat().st_size:
        raise ValueError("LEaP reported missing parameters or an empty topology")


def assemble_topology(
    substrate: Substrate,
    source_dir: Path,
    output: Path,
    amber_bin: Path,
    ligand_mol2: Path,
    frcmod: Path,
    ligand_library: Path | None,
    intact_pose: Path,
) -> dict[str, object]:
    """Produce a charge-consistent acyl topology and its machine-readable audit.

    Args:
        substrate: Product-gated substrate graph.
        source_dir: Selected mature Nyl12 WT model directory.
        output: New model directory.
        amber_bin: AmberTools26 bin directory.
        ligand_mol2: Full ligand MOL2 with GAFF types and charges.
        frcmod: Parameters for the full ligand.
        ligand_library: Original PA6 library, or None for poly(X,AMA).
        intact_pose: Bound-pose PDB used for atom coordinates.

    Returns:
        Topology audit including charge correction and exclusion repair.
    """
    output.mkdir(parents=True, exist_ok=True)
    _write_ligand_pdb(intact_pose, output / "intact.pdb")
    _write_retained_pdb(intact_pose, retained_names(substrate), output / "retained.pdb")
    library = _fragment_library(
        substrate, ligand_mol2, ligand_library, frcmod, output, amber_bin
    )
    _assemble(substrate, source_dir / "protein.pdb", frcmod, library, output, amber_bin, None)
    validator = Path(__file__).with_name("validate_topology.py")
    _run(
        [str(amber_bin / "python"), str(validator), str(output), substrate.carbonyl_c],
        output, "topology_audit.log",
    )
    initial = json.loads((output / "topology_audit.json").read_text())
    correction = round(initial["total_charge"]) - initial["total_charge"]
    baseline = initial["carbonyl_partial_charge"]
    _assemble(
        substrate, source_dir / "protein.pdb", frcmod, library,
        output, amber_bin, baseline + correction,
    )
    _run(
        [str(amber_bin / "python"), str(validator), str(output), substrate.carbonyl_c],
        output, "topology_audit.log",
    )
    final = json.loads((output / "topology_audit.json").read_text())
    if abs(final["total_charge"] - round(final["total_charge"])) > 0.002:
        raise ValueError("Acyl topology has noninteger total charge")
    if abs(final["qm_region_partial_charge"] - 1.0) > 0.002:
        raise ValueError("Acyl QM region does not carry the expected +1 charge")
    final["charge_correction_on_qm_carbonyl_c"] = correction
    final["charge_model"] = "provisional QM/MM only; do not run classical MD"
    (output / "topology_audit.json").write_text(json.dumps(final, indent=2) + "\n")
    return final
