"""Prepare a provenance-pinned unbiased DFT/MM TI1 endpoint release."""

import hashlib
import json
from pathlib import Path
import shutil


def prepare(project: Path) -> None:
    """Reuse the audited DFT interface with 28 physical CPU threads."""
    parent = (
        project
        / "qmmm/pa6_resp_gaff2_wt/11_frame1603_12A_reaction_path_20261109/reaction_paths/02_asynchronous_es_to_ti1"
    )
    source = (
        parent
        / "03_corrected_attack_from_window24/calculation/windows/attack_09/optimized.ncrst"
    )
    old = parent / "03_window24_wb97xd_def2svp_mm/optimization"
    output = parent / "04_ti1_unbiased_dft_mm"
    output.mkdir(exist_ok=False)
    shutil.copy2(source, output / "input.ncrst")
    text = (old / "minimize.mdin").read_text()
    text = text.replace("Window 24", "Corrected attack_09 TI1-like")
    text = text.replace("num_threads=16", "num_threads=28")
    text = text.replace("mem='32GB'", "mem='64GB'")
    text = text.replace("ntr=0,", "ntr=0, nmropt=0,")
    if any(
        token in text.lower() for token in ("disang", "restraintmask", "restraint_wt")
    ):
        raise ValueError("Unexpected restraint in DFT release input")
    (output / "minimize.mdin").write_text(text)
    shutil.copy2(old / "gau_job.tpl", output / "gau_job.tpl")
    files = [
        source,
        old / "minimize.mdin",
        output / "input.ncrst",
        output / "minimize.mdin",
        output / "gau_job.tpl",
    ]
    record = {
        "starting_state": "restrained attack_09 TI1-like endpoint",
        "reaction_restraints": False,
        "positional_restraints": False,
        "movable_region": "existing 12 A shell; outside atoms fixed",
        "stationary_point_claim": False,
        "sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    }
    (output / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    prepare(Path(__file__).resolve().parents[1])
