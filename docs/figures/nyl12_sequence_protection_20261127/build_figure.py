"""Build the exact canonical protection inventory and publication SI figure.

Run with the validated nyl12-pymol interpreter, OMP_NUM_THREADS=1, nice -n 15.
No coordinates are optimized; the display uses a proper protein-defined rotation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from prepare_relay_broad_pilot import CENSUS_PARENTS  # noqa: E402
from nyl12_covalent.relay_broad_pilot import protection_census  # noqa: E402

OUT = Path(__file__).resolve().parent
MAGENTA = "#D400B5"
GRAY = "#E6E6E6"
TEAL = "#009E73"
MECHANISTIC_CORE = frozenset({113, 156, 186, 233, 272, 274})


def digest(path: Path) -> str:
    """Return an exact-byte SHA-256 digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    """Write a new JSON artifact, refusing replacement."""
    mode = (
        "w"
        if "--recolor" in sys.argv
        and path.name in {"render_manifest.json", "manifest.json"}
        else "x"
    )
    with path.open(mode) as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def inventory(receipt: dict) -> list[dict]:
    """Expand the policy union into all 328 canonical sequence positions."""
    policy = receipt["protection_policy"]
    sequence = (
        receipt["source_sequences"]["alpha"] + receipt["source_sequences"]["beta"]
    )
    rows = []
    for position, aa in enumerate(sequence, 1):
        alpha = position <= 232
        family = "alpha" if alpha else "beta_local"
        local = position if alpha else position - 232
        reasons = []
        for key in ("hard_core", "recurring_inter_repeat_interface", "disulfides"):
            if local in policy[key][family]:
                reasons.append(key)
        fixed = position in receipt["fixed"]["precursor_positions"]
        if fixed != bool(reasons):
            raise ValueError(f"Policy-union mismatch at {position}")
        rows.append(
            dict(
                canonical_position=position,
                wt=aa,
                wt_label=f"{aa}{position}",
                family=family,
                local_position=local,
                protected=fixed,
                copies="ABCD" if alpha else "EFGH",
                reasons=";".join(reasons) if reasons else "mutable",
            )
        )
    assert len(rows) == 328 and sum(r["protected"] for r in rows) == 101
    return rows


def audit_candidates(receipt: dict) -> dict:
    """Audit preserved generated sequences without changing running jobs."""
    pilot = (
        ROOT
        / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/protection_corrected_pilot_20261127"
    )
    pool = pilot / "generation_run01/deduplicated_pool.json"
    wt = receipt["source_sequences"]["alpha"] + receipt["source_sequences"]["beta"]
    results = []
    for candidate in json.loads(pool.read_text()):
        sequence = candidate["alpha_sequence"] + candidate["beta_sequence"]
        chains = candidate["sequence"].split(":")
        if (
            len(chains) != 8
            or chains
            != [candidate["alpha_sequence"]] * 4 + [candidate["beta_sequence"]] * 4
        ):
            raise ValueError("Candidate sequence copy symmetry failed")
        violations = [
            f"{wt[p - 1]}{p}{sequence[p - 1]}"
            for p in receipt["fixed"]["precursor_positions"]
            if wt[p - 1] != sequence[p - 1]
        ]
        identifier = candidate.get("candidate_id", candidate.get("id"))
        if identifier is None:
            raise ValueError(f"Missing candidate identifier: {list(candidate)}")
        results.append(
            dict(
                candidate_id=identifier,
                residue156=sequence[155],
                violations=violations,
                eligible=not violations,
            )
        )
    processes = {}
    for name in ("generate_01", "generate_02", "generate_03", "generate_04"):
        path = pilot / f"generation_run01/{name}.process.json"
        processes[name] = json.loads(path.read_text())
    score = pilot / "scoring_recovery_attempt01/scores_run01"
    return dict(
        pool=str(pool),
        pool_sha256=digest(pool),
        candidates=results,
        eligible_ids=[r["candidate_id"] for r in results if r["eligible"]],
        generation_processes=processes,
        completed_score_receipts=len(list(score.glob("completed_relay_corr01_*.json"))),
        audit_is_snapshot=True,
        jobs_modified=False,
    )


def render(receipt: dict) -> dict:
    """Render complete WT TS-like tetramer with the approved protein camera."""
    import pymol2

    source = CENSUS_PARENTS[0]
    data = {
        (line[21], int(line[22:26]), line[12:16].strip()): np.array(
            [float(line[s : s + 8]) for s in (30, 38, 46)]
        )
        for line in source.read_text().splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    }
    phe, asn, thr = [
        data[k] for k in (("A", 113, "CA"), ("A", 186, "CA"), ("E", 233, "CA"))
    ]
    x = asn - phe
    x /= np.linalg.norm(x)
    y = (phe + asn) / 2 - thr
    y -= x * np.dot(y, x)
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    rotation = np.column_stack((x, y, z))
    assert np.isclose(np.linalg.det(rotation), 1)
    fixed = receipt["fixed"]
    selection = (
        "(chain A+B+C+D and resi "
        + "+".join(map(str, fixed["alpha"]))
        + ") or (chain E+F+G+H and resi "
        + "+".join(str(p + 232) for p in fixed["beta_local"])
        + ")"
    )
    commands = [
        f"load {source}, model",
        "bg_color white",
        "set ray_opaque_background, 1",
        "set max_threads, 1",
        "set valence, 0",
        "set stick_radius, 0.2",
        "hide everything",
        "show cartoon, chain A+B+C+D+E+F+G+H",
        "set_color protected_magenta, [0.83137,0,0.70980]",
        "set_color mechanistic_teal, [0,0.61961,0.45098]",
        "set_color variable_gray, [0.90196,0.90196,0.90196]",
        "set ambient, 0.65",
        "set direct, 0.35",
        "set reflect, 0.15",
        "set specular, 0",
        "set ray_shadows, off",
        "set depth_cue, off",
        "color variable_gray, all",
        f"color protected_magenta, {selection}",
        "color mechanistic_teal, (chain A+B+C+D and resi 113+156+186) or (chain E+F+G+H and resi 233+272+274)",
        "show sticks, chain I+J+K+L",
        "color orange, chain I+J+K+L and elem C",
        "color blue, chain I+J+K+L and elem N",
        "color red, chain I+J+K+L and elem O",
        "set orthoscopic, on",
        "set antialias, 2",
        "orient all",
        "zoom all, 3",
    ]
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        for command in commands:
            cmd.do(command)
        if cmd.count_atoms(f"({selection}) and name CA") != 404:
            raise ValueError(
                "Rendered protection selection does not contain 404 CA atoms"
            )
        if (
            cmd.count_atoms(
                "((chain A+B+C+D and resi 113+156+186) or (chain E+F+G+H and resi 233+272+274)) and name CA"
            )
            != 24
        ):
            raise ValueError("Mechanistic subset must contain 24 physical residues")
        view = list(cmd.get_view())
        view[:9] = rotation.ravel().tolist()
        cmd.set_view(view)
        cmd.png(str(OUT / "structure.png"), width=1800, height=1400, dpi=300, ray=1)
        cmd.save(str(OUT / "protected_positions.pse"))
        version = cmd.get_version()[0]
    commands += ["set_view (" + ",".join(map(str, view)) + ")"]
    (OUT / "protected_positions.pml").write_text("\n".join(commands) + "\n")
    return dict(
        source=str(source),
        source_sha256=digest(source),
        pymol_version=version,
        camera=view,
        rotation_determinant=float(np.linalg.det(rotation)),
        camera_landmarks="A:Phe113 CA, A:Asn186 CA, E:Thr233 CA; approved protein frame",
        protein_color_exception="Protection status replaces chain coloring",
        max_threads=1,
        palette={
            "other_protected": MAGENTA,
            "mechanistic_core": TEAL,
            "variable": GRAY,
        },
        mechanistic_core_canonical=sorted(MECHANISTIC_CORE),
        mechanistic_core_interpretation="Six-position mechanistic subset; includes homolog-informed K156; not all experimentally established in Nyl12",
    )


def compose(rows: list[dict]) -> None:
    """Compose a molecular overview and exact sequence-level mask."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig = plt.figure(figsize=(10, 9), facecolor="white")
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1.2], hspace=0.03)
    ax = fig.add_subplot(gs[0])
    ax.imshow(plt.imread(OUT / "structure.png"))
    ax.axis("off")
    ax.set_title(
        "A   Protected sequence identities in the WT TS-like tetramer",
        loc="left",
        fontsize=14,
    )
    ax = fig.add_subplot(gs[1])
    for row in rows:
        p = row["canonical_position"]
        line = (p - 1) // 82
        col = (p - 1) % 82
        ax.add_patch(
            plt.Rectangle(
                (col, line),
                0.9,
                0.7,
                color=TEAL
                if p in MECHANISTIC_CORE
                else MAGENTA
                if row["protected"]
                else GRAY,
                lw=0,
            )
        )
    ax.set_xlim(-1, 83)
    ax.set_ylim(4, -0.6)
    ax.set_xticks([0, 19, 39, 59, 81], ["1", "20", "40", "60", "82"])
    ax.set_yticks([0.35, 1.35, 2.35, 3.35], ["1–82", "83–164", "165–246", "247–328"])
    ax.set_xlabel("Position within each 82-residue row", fontsize=10)
    ax.set_title(
        "B   Canonical precursor mask (328 positions; four symmetry-tied copies)",
        loc="left",
        fontsize=12,
    )
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    ax.tick_params(length=0)
    fig.legend(
        handles=[
            Patch(color=MAGENTA, label="Other protected: 95"),
            Patch(color=TEAL, label="Mechanistic core: 6"),
            Patch(color=GRAY, label="Variable: 227"),
            Patch(color="orange", label="PA6 carbon"),
        ],
        loc="lower center",
        ncol=2,
        frameon=False,
    )
    fig.savefig(
        OUT / "protected_vs_variable.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        OUT / "protected_vs_variable.pdf", bbox_inches="tight", facecolor="white"
    )
    plt.close(fig)


def main() -> None:
    """Create fresh policy, full inventory, candidate audit and SI artwork."""
    import subprocess

    logging.basicConfig(level=logging.INFO)
    if "--render-only" in sys.argv:
        receipt = json.loads(
            (OUT / "canonical_protection_restored_K156_20261127.json").read_text()
        )
        write_json(OUT / "render_manifest.json", render(receipt))
        return
    if "--recolor" in sys.argv:
        receipt = json.loads(
            (OUT / "canonical_protection_restored_K156_20261127.json").read_text()
        )
        rows = inventory(receipt)
        subprocess.run(
            [
                "/data/jpc/code/miniforge3/envs/nyl12-pymol/bin/python",
                str(Path(__file__)),
                "--render-only",
                "--recolor",
            ],
            check=True,
        )
        compose(rows)
        manifest = json.loads((OUT / "render_manifest.json").read_text())
        manifest.update(
            script_sha256=digest(Path(__file__)),
            prior_version="prior_v1",
            coordinate_optimization=False,
            mask_sha256=digest(
                OUT / "canonical_protection_restored_K156_20261127.json"
            ),
            artifacts={
                p.name: digest(p)
                for p in OUT.iterdir()
                if p.is_file()
                and p.name not in {"manifest.json", "acceptance.json", "build.log"}
            },
        )
        write_json(OUT / "manifest.json", manifest)
        return
    receipt = protection_census(CENSUS_PARENTS)
    if receipt["counts"]["fixed_total"] != 101:
        raise ValueError("Expected restored K156 policy with 101 protected positions")
    rows = inventory(receipt)
    write_json(OUT / "canonical_protection_restored_K156_20261127.json", receipt)
    write_json(OUT / "canonical_position_inventory.json", rows)
    with (OUT / "canonical_position_inventory.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json(OUT / "candidate_eligibility_audit.json", audit_candidates(receipt))
    subprocess.run(
        [
            "/data/jpc/code/miniforge3/envs/nyl12-pymol/bin/python",
            str(Path(__file__)),
            "--render-only",
        ],
        check=True,
    )
    manifest = json.loads((OUT / "render_manifest.json").read_text())
    compose(rows)
    manifest.update(
        python=sys.executable,
        numpy_version=np.__version__,
        coordinate_optimization=False,
        script_sha256=digest(Path(__file__)),
        artifacts={
            p.name: digest(p)
            for p in OUT.iterdir()
            if p.is_file() and p.name != "manifest.json"
        },
    )
    write_json(OUT / "manifest.json", manifest)
    logging.info("Built SI figure and exact %s", receipt["counts"])


if __name__ == "__main__":
    main()
