"""Build product-gated Nyl12 acyl-enzyme models and DFTB3 inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.chemistry import (  # noqa: E402
    build_poly_xama_maxa,
    load_pa6,
    product_report,
)
from nyl12_covalent.pose import embed_on_pa6_register, write_named_pdb  # noqa: E402
from nyl12_covalent.qmmm import run_dftb3, write_dftb3_input  # noqa: E402
from nyl12_covalent.topology import (  # noqa: E402
    assemble_topology,
    parameterize_poly_xama,
)

SOURCE = ROOT / (
    "output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/1_relaxation/"
    "pa6_l4z_l2_nyl12_wt_model_41"
)
AMBER_BIN = Path("/data/jpc/code/ambertools26/venv/bin")
LOGGER = logging.getLogger(__name__)


def _digest(path: Path) -> str:
    """Return a SHA-256 digest for one provenance file."""
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def build(
    target: str,
    output: Path,
    source_dir: Path = SOURCE,
    amber_bin: Path = AMBER_BIN,
    conformers: int = 24,
    poly_charge_method: str = "bcc",
    run_qmmm: bool = False,
) -> dict[str, object]:
    """Build one assayed Nyl12 acyl state with complete provenance.

    Args:
        target: ``pa6_l2`` or ``poly_xama_maxa``.
        output: Dedicated new model directory.
        source_dir: Matched source PA6 complex and Amber parameters.
        amber_bin: AmberTools26 bin directory.
        conformers: Oligomer conformers for poly(X,AMA) pose selection.
        poly_charge_method: Amber ``bcc`` or QM-only ``gas`` charge assignment.
        run_qmmm: Execute DFTB3 after topology validation.

    Returns:
        Product, source, topology, and QM/MM model audit.
    """
    if target not in {"pa6_l2", "poly_xama_maxa"}:
        raise ValueError(f"Unsupported product target {target}")
    output = output.resolve()
    source_dir = source_dir.resolve()
    amber_bin = amber_bin.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Model output already exists and is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    source_mol2 = source_dir / "lig_sub/lig_h.mol2"
    pose_audit: dict[str, object] | None = None
    if target == "pa6_l2":
        substrate = load_pa6(source_mol2)
        ligand_mol2 = source_mol2
        frcmod = source_dir / "lig_sub.frcmod"
        library: Path | None = source_dir / "lig_sub.lib"
        pose = source_dir / "complex.pdb"
    else:
        substrate = build_poly_xama_maxa()
        embedded, pose_audit = embed_on_pa6_register(
            substrate, source_dir / "complex.pdb", conformers=conformers
        )
        pose = output / "intact_poly_xama_pose.pdb"
        write_named_pdb(embedded, pose)
        ligand_mol2, frcmod = parameterize_poly_xama(
            embedded, output, amber_bin, poly_charge_method
        )
        library = None
    product = product_report(substrate)
    audit = assemble_topology(
        substrate, source_dir, output, amber_bin,
        ligand_mol2, frcmod, library, pose,
    )
    mdin = write_dftb3_input(output)
    sources = [
        source_dir / "complex.pdb", source_dir / "protein.pdb",
        ligand_mol2, frcmod,
    ]
    if library is not None:
        sources.append(library)
    manifest: dict[str, object] = {
        "target": target,
        "experimental_readout": {
            "pa6_l2": "PA6 L2, [M+H]+ m/z 245.19, IDOT/OPSI-MS MRM 245.19>114.09",
            "poly_xama_maxa": "MAXA from poly(X,AMA), IDOT/OPSI-MS in Dishner 2026 JACS",
        }[target],
        "product_gate": product,
        "pose_audit": pose_audit,
        "poly_charge_method": poly_charge_method if pose_audit is not None else None,
        "source_files_sha256": {str(path): _digest(path) for path in sources},
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "amber_bin": str(amber_bin),
        "topology_audit": audit,
        "qmmm_input": str(mdin),
        "limitations": [
            "The source PA6 pose is a Boltz affinity-mode model; atom names are taken from its saved MOL2/PDB pair.",
            "Cross-residue ester partial charges are provisional and QM-region-only; classical MD requires new validated full parameters.",
            "No reaction barrier or catalytic rate follows from a single acyl-state minimum.",
        ],
    }
    if pose_audit is not None:
        manifest["limitations"].append(
            "The poly(X,AMA) conformer is a register-aligned candidate with an explicit clash score, not a validated bound pose."
        )
        if poly_charge_method == "gas":
            manifest["limitations"].append(
                "Gasteiger charges are used only to construct an all-ligand-QM model; AM1-BCC or RESP is required for classical MD."
            )
    (output / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if run_qmmm:
        if pose_audit is not None and pose_audit["protein_clash_penalty"] > 5.0:
            raise ValueError("Poly(X,AMA) pose has unresolved steric clashes; QM/MM not started")
        run_dftb3(output, amber_bin)
    LOGGER.info("Built %s in %s", target, output)
    return manifest


def main() -> None:
    """Parse CLI arguments and run the requested reproducible build."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("pa6_l2", "poly_xama_maxa"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--amber-bin", type=Path, default=AMBER_BIN)
    parser.add_argument("--conformers", type=int, default=24)
    parser.add_argument("--poly-charge-method", choices=("bcc", "gas"), default="bcc")
    parser.add_argument("--run-qmmm", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    output = args.output or ROOT / "output_nyl12_covalent_qmmm_20261019" / args.target
    build(
        args.target, output, args.source_dir, args.amber_bin,
        args.conformers, args.poly_charge_method, args.run_qmmm,
    )


if __name__ == "__main__":
    main()
