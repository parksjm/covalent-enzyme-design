"""Snapshot and launch the explicit small PA6 L4 TI1/water AF3 contact test."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path
import shlex
import shutil
import sys

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from nyl12_covalent.jobs import launch_job  # noqa: E402


def prepare_input(original: Path, output: Path) -> Path:
    """Retain only alpha A, beta E, covalent substrate I, and water J."""
    payload = json.loads(original.read_text())
    proteins = [item for item in payload["sequences"] if "protein" in item]
    ligands = [item for item in payload["sequences"] if "ligand" in item]
    if len(proteins) != 2 or len(ligands) != 1:
        raise ValueError("Expected two protein entities and one covalent substrate")
    expected = [[["E", 1, "OG1"], ["I", 1, "C54"]]]
    if payload["bondedAtomPairs"] != expected:
        raise ValueError("Input must be the audited PA6 L4 TI1 covalent graph")
    for item, ids, keep, length in zip(
        proteins,
        (["A", "B", "C", "D"], ["E", "F", "G", "H"]),
        ("A", "E"),
        (232, 96),
        strict=True,
    ):
        protein = item["protein"]
        if protein["id"] != ids or len(protein["sequence"]) != length:
            raise ValueError("Unexpected Nyl12 mature chain identity or length")
        for field in ("unpairedMsa", "pairedMsa", "templates"):
            if field not in protein:
                raise ValueError(f"Missing precomputed {field}")
        protein["id"] = keep
    ccd = Path(payload.pop("userCCDPath")).resolve(strict=True)
    shutil.copy2(ccd, output / "tetrahedral_substrate.cif")
    payload.update(
        name="nyl12_pa6_l4_ti1_alpha_beta_water_contact", modelSeeds=[20261021]
    )
    payload["userCCDPath"] = str(output / "tetrahedral_substrate.cif")
    payload["sequences"] = [
        *proteins,
        *ligands,
        {"ligand": {"id": "J", "ccdCodes": ["HOH"]}},
    ]
    payload["contactConstraints"] = [
        {
            "token1": ["E", 1, atom],
            "token2": ["J", 1, "O"],
            "minDistance": 2.6,
            "maxDistance": 3.2,
            "weight": 1.0,
            "mobileEndpoint": "token2",
        }
        for atom in ("N", "OG1")
    ]
    path = output / "input.json"
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def main() -> None:
    """Preserve sources and inputs, then launch preflight followed by inference."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-input", required=True, type=Path)
    parser.add_argument("--fork", required=True, type=Path)
    parser.add_argument("--python", required=True, type=Path)
    parser.add_argument("--purelib", required=True, type=Path)
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--session", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    output = args.output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(args.original_input, output / "original_input.json")
    input_path = prepare_input(args.original_input, output)
    source = output / "source"
    source.mkdir()
    shutil.copytree(
        args.fork / "src",
        source / "src",
        symlinks=True,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    for origin in (
        args.fork / "run_alphafold.py",
        ROOT / "scripts/check_af3_water_features.py",
        Path(__file__),
        ROOT / "src/nyl12_covalent/jobs.py",
    ):
        shutil.copy2(origin, source / origin.name)
    manifest = {
        "purpose": "One-seed one-sample three-recycle AF3 water/contact smoke test; not TI1 validation",
        "source_fork": str(args.fork),
        "original_input": str(args.original_input),
        "python": str(args.python),
        "model_dir": str(args.model_dir),
        "atom_mapping": "Full-sequence Thr233 is E1; alpha A is residues1-232; beta E residues233-328",
        "limitations": "AF3 predicts heavy atoms; NH3 protonation and water H directions require downstream rebuilding and validation",
        "files_sha256": {
            str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(output.rglob("*"))
            if p.is_file()
        },
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    check = [
        str(args.python),
        "-S",
        str(source / "check_af3_water_features.py"),
        "--input",
        str(input_path),
        "--output",
        str(output / "feature_acceptance.json"),
    ]
    inference = [
        str(args.python),
        "-S",
        str(source / "run_alphafold.py"),
        "--json_path",
        str(input_path),
        "--output_dir",
        str(output / "prediction"),
        "--model_dir",
        str(args.model_dir),
        "--run_data_pipeline=false",
        "--num_diffusion_samples=1",
        "--num_recycles=3",
        "--buckets=512",
        "--gpu_device=0",
        "--contact_guidance_weight=0.5",
        "--contact_guidance_final_weight_scale=1.0",
        "--contact_guidance_interval=1",
        "--contact_guidance_steps=20",
        "--jax_compilation_cache_dir",
        str(output / "jax_cache"),
    ]
    workflow = output / "workflow.sh"
    workflow.write_text(
        "#!/bin/bash\nset -euo pipefail\n"
        + shlex.join(check)
        + "\n"
        + shlex.join(inference)
        + "\n"
    )
    env = {
        "CUDA_VISIBLE_DEVICES": "0",
        "PYTHONPATH": f"{source / 'src'}:{args.purelib}",
        "PATH": f"{args.python.parent}:/usr/local/bin:/usr/bin:/bin",
        "XLA_PYTHON_CLIENT_PREALLOCATE": "false",
        "OMP_NUM_THREADS": "4",
        "OPENBLAS_NUM_THREADS": "1",
    }
    launch_job(
        args.session,
        output / "job",
        output,
        ["/bin/bash", str(workflow)],
        env,
        [input_path, manifest_path, workflow],
    )
    logging.info("Launched %s; feature acceptance must precede inference", args.session)


if __name__ == "__main__":
    main()
