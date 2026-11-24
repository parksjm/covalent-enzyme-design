"""Launch the frozen, balanced MD queues on their explicitly assigned host."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy
import parmed

PROFILES = {
    "violet": {
        "engine": "/data/jpc/code/pmemd26/bin/pmemd.cuda_SPFP",
        "sha256": "7ed3f439d09d6beff94210bb3f356585e6918cd3f500dcc58046067ac8fd9602",
        "config_prefix": "",
        "queue_suffix": "",
        "receipt": "",
    },
    "azalea": {
        "engine": "/data/jpc/code/pmemd26_blackwell_build_20261005/pmemd26_src/build_sysgcc_cuda128_v2/src/pmemd/src/pmemd.cuda_SPFP",
        "sha256": "f4b0b32eda07d1e6b4fbe597e8c3bd4d39cc492a2b01c19c13db7ebe92ff52b1",
        "config_prefix": "azalea_blackwell_",
        "queue_suffix": "_blackwell",
        "receipt": "blackwell_retry",
    },
}


def main() -> None:
    """Prepare missing host-local plans and launch two serial GPU workers."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", required=True, choices=("violet", "azalea"))
    parser.add_argument("--prepare-plans", action="store_true")
    args = parser.parse_args()
    if platform.node().split(".")[0] != args.host:
        raise ValueError("Requested host differs from the execution workstation")
    campaign = Path(__file__).resolve().parent
    repo = campaign / "source_snapshot"
    profile = PROFILES[args.host]
    engine = Path(profile["engine"])
    if hashlib.sha256(engine.read_bytes()).hexdigest() != profile["sha256"]:
        raise ValueError(
            "Host-specific Amber executable hash differs; read runtime guide"
        )
    configs = [
        campaign / "configs" / f"{profile['config_prefix']}{label}.json"
        for label in ("WT", "broader_mask07", "focused_mask03", "V111C")
    ]
    for config in configs:
        selected = Path(json.loads(config.read_text())["cuda_engine"]).resolve(
            strict=True
        )
        if selected != engine.resolve(strict=True):
            raise ValueError("Frozen plan selects the wrong workstation CUDA engine")
    sys.path.insert(0, str(repo / "src"))
    from nyl12_covalent.jobs import launch_job
    from nyl12_covalent.prepared_md_queue_plan import runtime

    runner = runtime(repo)
    if args.prepare_plans:
        for config in configs:
            runner.plan(config, repo)
    information = {
        "host": platform.node(),
        "python": sys.executable,
        "numpy": numpy.__version__,
        "parmed": parmed.__version__,
        "cuda_engine": str(engine),
        "cuda_engine_sha256": profile["sha256"],
        "gpu": subprocess.run(
            ["nvidia-smi", "--query-gpu=index,name,uuid", "--format=csv,noheader"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout,
    }
    (campaign / "queues").mkdir(exist_ok=True)
    with (campaign / "queues" / f"{args.host}_environment.json").open("x") as stream:
        json.dump(information, stream, indent=2)
        stream.write("\n")
    for gpu in ("0", "1"):
        config = (
            campaign
            / "configs"
            / f"queue_{args.host}_gpu{gpu}{profile['queue_suffix']}.json"
        )
        job = campaign / "queues" / f"{args.host}_gpu{gpu}" / profile["receipt"] / "job"
        environment = {
            "AMBERHOME": "/data/jpc/code/ambertools26/venv",
            "LD_LIBRARY_PATH": "/data/jpc/code/cuda128_findcuda_compat/lib64:/data/jpc/code/pmemd26/lib:/data/jpc/code/ambertools26/venv/lib"
            if args.host == "azalea"
            else "/data/jpc/code/pmemd26/lib:/data/jpc/code/ambertools26/venv/lib",
            "CUDA_VISIBLE_DEVICES": gpu,
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
        launch_job(
            f"nyl12_design_followup_{args.host}{profile['queue_suffix']}_gpu{gpu}_20261123",
            job,
            repo,
            [
                sys.executable,
                str(repo / "scripts/run_prepared_md_sequence.py"),
                "--config",
                str(config),
                "--repo",
                str(repo),
                "--gpu",
                gpu,
                "--status",
                str(job.parent / "status.json"),
            ],
            env_overrides=environment,
            input_paths=[
                config,
                campaign / "selection.json",
                Path(__file__),
                repo / "scripts/run_prepared_md_sequence.py",
            ],
        )


if __name__ == "__main__":
    main()
