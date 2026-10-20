"""Prepare AF3 covalent tetrahedral-intermediate predictions."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from rdkit import Chem

from .chemistry import Substrate, atom_by_name, product_report
from .pose import pdb_coordinates


def _sha256(path: Path) -> str:
    """Hash a scientific input for model provenance."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_ti1_ccd(substrate: Substrate, pose: Path, output: Path) -> str:
    """Write a custom CCD with the tetrahedral C–O and intact amide bonds.

    The external Thr OG1–carbonyl bond is specified separately in AF3 JSON.
    This CCD includes the full substrate, so the L2 or MAXA product register
    remains unambiguous. The carbonyl O is a formal oxyanion in this state.

    Args:
        substrate: Product-gated intact substrate graph.
        pose: Named bound ligand pose used for ideal coordinates.
        output: New user CCD path.

    Returns:
        Unique chemical component identifier.
    """
    product_report(substrate)
    component = {"pa6_l2": "PA6TI1", "poly_xama_maxa": "XAMATI1"}[substrate.key]
    coordinates = pdb_coordinates(pose)
    heavy = [atom for atom in substrate.mol.GetAtoms() if atom.GetAtomicNum() > 1]
    names = {atom.GetProp("name") for atom in heavy}
    if set(coordinates) != names and not names.issubset(coordinates):
        raise ValueError(f"CCD pose lacks atoms {sorted(names - set(coordinates))}")
    atoms = [
        f"data_{component}", "#", f"_chem_comp.id {component}",
        "_chem_comp.name 'Nyl12 tetrahedral covalent adduct substrate'",
        "_chem_comp.type non-polymer", "_chem_comp.formula ?",
        "_chem_comp.mon_nstd_parent_comp_id ?", "_chem_comp.pdbx_synonyms ?",
        "_chem_comp.formula_weight ?", "#", "loop_",
        "_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol", "_chem_comp_atom.charge",
        "_chem_comp_atom.pdbx_leaving_atom_flag",
        "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
        "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
    ]
    for atom in heavy:
        name = atom.GetProp("name")
        point = coordinates[name]
        charge = -1 if name == substrate.carbonyl_o else atom.GetFormalCharge()
        atoms.append(
            f"{component} {name} {atom.GetSymbol()} {charge} N "
            f"{point[0]:.3f} {point[1]:.3f} {point[2]:.3f}"
        )
    atoms.extend(
        ["#", "loop_", "_chem_comp_bond.atom_id_1",
         "_chem_comp_bond.atom_id_2", "_chem_comp_bond.value_order",
         "_chem_comp_bond.pdbx_aromatic_flag"]
    )
    carbon = atom_by_name(substrate.mol, substrate.carbonyl_c)
    oxygen = atom_by_name(substrate.mol, substrate.carbonyl_o)
    scissile = frozenset((carbon.GetIdx(), oxygen.GetIdx()))
    bond_orders = {
        Chem.BondType.SINGLE: "SING",
        Chem.BondType.DOUBLE: "DOUB",
        Chem.BondType.TRIPLE: "TRIP",
        Chem.BondType.AROMATIC: "AROM",
    }
    for bond in substrate.mol.GetBonds():
        first, second = bond.GetBeginAtom(), bond.GetEndAtom()
        if first.GetAtomicNum() == 1 or second.GetAtomicNum() == 1:
            continue
        pair = frozenset((first.GetIdx(), second.GetIdx()))
        order = "SING" if pair == scissile else bond_orders[bond.GetBondType()]
        aromatic = "Y" if bond.GetIsAromatic() else "N"
        atoms.append(
            f"{first.GetProp('name')} {second.GetProp('name')} {order} {aromatic}"
        )
    output.write_text("\n".join([*atoms, "#"]) + "\n")
    return component


def prepare_af3_ti1(
    substrate: Substrate,
    pose: Path,
    msa_source: Path,
    output: Path,
    seed: int = 20261019,
    gpu: int = 0,
) -> dict[str, object]:
    """Create a covalent AF3 input and local inference runner.

    Args:
        substrate: Product-gated intact substrate graph.
        pose: Atom-named PA6 or poly(X,AMA) bound-pose PDB.
        msa_source: Existing AF3 data JSON for the same Nyl12 octamer.
        output: Dedicated predictor job directory.
        seed: AF3 diffusion seed.
        gpu: CUDA device number.

    Returns:
        Input, software, and chemistry provenance.
    """
    output.mkdir(parents=True, exist_ok=True)
    ccd = output / "tetrahedral_substrate.cif"
    component = write_ti1_ccd(substrate, pose, ccd)
    prior = json.loads(msa_source.read_text())
    proteins = [item for item in prior["sequences"] if "protein" in item]
    if len(proteins) != 2:
        raise ValueError("Expected one alpha and one beta Nyl12 protein entity")
    if proteins[0]["protein"]["id"] != ["A", "B", "C", "D"]:
        raise ValueError("Alpha-chain MSA source IDs changed")
    if proteins[1]["protein"]["id"] != ["E", "F", "G", "H"]:
        raise ValueError("Beta-chain MSA source IDs changed")
    job_name = f"nyl12_{substrate.key}_ti1_af3_seed{seed}"
    payload = {
        "name": job_name,
        "dialect": "alphafold3",
        "version": 3,
        "modelSeeds": [seed],
        "sequences": [*proteins, {"ligand": {"id": "I", "ccdCodes": [component]}}],
        "bondedAtomPairs": [
            [["E", 1, "OG1"], ["I", 1, substrate.carbonyl_c]]
        ],
        "userCCDPath": str(ccd.resolve()),
    }
    input_path = output / "af3_ti1_input.json"
    input_path.write_text(json.dumps(payload, indent=2) + "\n")
    af3_root = Path("/data/jpc/code/alphafold3")
    environment = af3_root / ".conda/envs/af3"
    runner = output / "run_af3_ti1.sh"
    runner.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        f"export CUDA_VISIBLE_DEVICES={gpu}\n"
        f"export CONDA_PREFIX={environment}\n"
        'export PATH="$CONDA_PREFIX/bin:$PATH"\n'
        'export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"\n'
        'export XLA_PYTHON_CLIENT_PREALLOCATE=false\n'
        'export XLA_FLAGS="${XLA_FLAGS:-} --xla_gpu_enable_triton_gemm=false"\n'
        f"cd {af3_root}\n"
        f"{environment / 'bin/python'} {af3_root / 'run/run_af3.py'} "
        f"--json_path {input_path.resolve()} "
        f"--output_dir {(output / 'prediction').resolve()} "
        "--force_output_dir 1 --gpu_device 0 "
        "--num_diffusion_samples 1 --num_recycles 5 "
        f"--jax_compilation_cache_dir {(output / 'jax_cache').resolve()} "
        f"> {(output / 'af3_run.log').resolve()} 2>&1\n"
    )
    runner.chmod(0o755)
    manifest: dict[str, object] = {
        "state": "first tetrahedral covalent intermediate, amide C-N intact",
        "carbonyl_bond_order": "single",
        "oxyanion_atom": substrate.carbonyl_o,
        "protein_ligand_bond": payload["bondedAtomPairs"][0],
        "assayed_product_gate": product_report(substrate),
        "ccd": str(ccd.resolve()),
        "input": str(input_path.resolve()),
        "runner": str(runner.resolve()),
        "msa_source": str(msa_source.resolve()),
        "msa_source_sha256": _sha256(msa_source),
        "bound_pose": str(pose.resolve()),
        "bound_pose_sha256": _sha256(pose),
        "af3_root": str(af3_root),
        "seed": seed,
        "gpu": gpu,
        "limitation": "AF3 does not specify the beta N-terminal protonation or relay-water mechanism; inspect and refine the predicted covalent geometry by QM/MM.",
    }
    (output / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
