"""Build the first Nyl12 tetrahedral covalent intermediate with relay water."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nyl12_covalent.chemistry import (  # noqa: E402
    build_poly_xama_maxa,
    load_pa6,
)
from nyl12_covalent.qmmm import write_ti1_dftb3_input  # noqa: E402
from nyl12_covalent.tetrahedral import build_ti1  # noqa: E402

SOURCE = ROOT / (
    "output_nyl12_pa6_pa66_ssl_sites_1_4_ds1_aff1/1_relaxation/"
    "pa6_l4z_l2_nyl12_wt_model_41"
)
POLY_MODEL = ROOT / "output_nyl12_covalent_qmmm_20261019/poly_xama_maxa_gas2"
AMBER_BIN = Path("/data/jpc/code/ambertools26/venv/bin")
LOGGER = logging.getLogger(__name__)


def _sha(path: Path) -> str:
    """Hash an exact input file used to generate a model."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(
    target: str,
    output: Path,
    source_dir: Path = SOURCE,
    poly_model: Path = POLY_MODEL,
    amber_bin: Path = AMBER_BIN,
    run_qmmm: bool = False,
    protein_pdb: Path | None = None,
    ligand_pdb: Path | None = None,
) -> dict[str, object]:
    """Build a product-gated TI1 topology, water, and DFTB3 specification.

    Args:
        target: ``pa6_l2`` or ``poly_xama_maxa``.
        output: New TI1 directory.
        source_dir: Original matched Nyl12 WT bound complex.
        poly_model: Parameterized intact poly(X,AMA) candidate directory.
        amber_bin: Installed AmberTools26 binary directory.
        run_qmmm: Legacy flag; rejected before writes. Use the explicit attempt CLI.
        protein_pdb: Optional canonical AF3 protein coordinates.
        ligand_pdb: Optional canonical AF3 full-substrate coordinates.

    Returns:
        Machine-readable chemistry and topology provenance.
    """
    if run_qmmm:
        raise ValueError(
            "Build and screen TI1 first; run scripts/run_nyl12_ti1_qmmm.py with "
            "explicit coordinates, reference, and a fresh output via tmux"
        )
    if target not in {"pa6_l2", "poly_xama_maxa"}:
        raise ValueError(f"Unsupported target {target}")
    output = output.resolve()
    source_dir = source_dir.resolve()
    poly_model = poly_model.resolve()
    amber_bin = amber_bin.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"TI1 output already exists and is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    protein = (protein_pdb or source_dir / "protein.pdb").resolve()
    if target == "pa6_l2":
        substrate = load_pa6(source_dir / "lig_sub/lig_h.mol2")
        ligand_pose = (ligand_pdb or source_dir / "complex.pdb").resolve()
        ligand_mol2 = source_dir / "lig_sub/lig_h.mol2"
        frcmod = source_dir / "lig_sub.frcmod"
        library: Path | None = source_dir / "lig_sub.lib"
    else:
        substrate = build_poly_xama_maxa()
        ligand_pose = (ligand_pdb or poly_model / "intact_poly_xama_pose.pdb").resolve()
        ligand_mol2 = poly_model / "intact_poly_xama.mol2"
        frcmod = poly_model / "poly_xama.frcmod"
        library = None
    report = build_ti1(
        substrate, protein, ligand_pose, ligand_mol2, frcmod,
        library, output, amber_bin,
    )
    mdin = write_ti1_dftb3_input(output)
    sources = [protein, ligand_pose, ligand_mol2, frcmod]
    if library is not None:
        sources.append(library)
    manifest: dict[str, object] = {
        **report,
        "input_sha256": {str(path): _sha(path) for path in sources},
        "qmmm_input": str(mdin),
        "amber_bin": str(amber_bin),
        "proton_inventory": {
            "initial_ES_beta_N": "neutral NH2 (separate state, not this topology)",
            "TI1_beta_N": "cationic NH3+ with H1,H2,H3",
            "TI1_OG1": "covalent O-C bond; HG1 absent",
            "relay_water": "neutral H2O before and after transfer",
            "atom_tracking": "H3 on beta N represents the OG1 proton after water-mediated transfer",
        },
    }
    (output / "ti1_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Built %s first tetrahedral intermediate at %s", target, output)
    return manifest


def main() -> None:
    """Parse CLI options and build the requested TI1 model."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", choices=("pa6_l2", "poly_xama_maxa"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--poly-model", type=Path, default=POLY_MODEL)
    parser.add_argument("--amber-bin", type=Path, default=AMBER_BIN)
    parser.add_argument("--run-qmmm", action="store_true")
    parser.add_argument("--protein-pdb", type=Path)
    parser.add_argument("--ligand-pdb", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    output = args.output or ROOT / "output_nyl12_covalent_qmmm_20261019" / f"{args.target}_ti1_graph"
    build(
        args.target, output, args.source_dir, args.poly_model,
        args.amber_bin, args.run_qmmm, args.protein_pdb, args.ligand_pdb,
    )


if __name__ == "__main__":
    main()
