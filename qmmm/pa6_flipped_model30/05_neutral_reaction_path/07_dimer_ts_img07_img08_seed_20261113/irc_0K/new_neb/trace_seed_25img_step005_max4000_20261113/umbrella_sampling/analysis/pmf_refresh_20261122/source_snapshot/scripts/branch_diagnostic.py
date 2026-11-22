"""Extract bounded uncensored branch geometry from existing frozen MBAR observables."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np


def main() -> None:
    """Write biased raw counts, exact-thinning branch masses and frame locators."""
    p = Path(__file__).resolve().parents[2]
    old = p.parent / "requested_mbar_20261121"
    analysis = json.loads((old / "mbar/latest.json").read_text())
    weights = np.load(old / "path_pmf_bootstrap.npz")["weights"]
    counts = []
    locations = []
    selected_owner = []
    selected_ti = []
    for summary in analysis["window_summaries"]:
        rows = []
        for file in sorted((p / "frozen_inputs" / summary["label"]).glob("*/observables/observables.csv")):
            for row in csv.DictReader(file.open()):
                row["locator"] = str(file.relative_to(p / "frozen_inputs"))
                rows.append(row)
        owner = [r["hg1_owner"] == "13905" for r in rows]
        ti = [o and float(r["r_c_og1"]) < 1.7 and float(r["r_c_n"]) < 1.7
              for o, r in zip(owner, rows, strict=True)]
        counts.append(dict(window=summary["label"], raw_frames=len(rows),
                           hg1_on_thr_n_biased_frames=sum(owner), ti1_like_biased_frames=sum(ti)))
        stride = summary["decorrelation_stride"]
        selected_owner.extend(owner[::stride])
        selected_ti.extend(ti[::stride])
        chosen = [r for r, match in zip(rows, ti, strict=True) if match]
        for row in chosen[::max(1, len(chosen) // 3)][:3]:
            locations.append({k: row[k] for k in ["locator", "frame", "time_ps", "r_c_n", "r_c_og1", "r_og1_hg1", "r_nleave_hg1", "hg1_owner", "nalpha_h1_owner", "nalpha_h2_owner"]})
    assert len(weights) == len(selected_owner) == len(selected_ti) == sum(analysis["mbar_N_k"])
    result = dict(
        criterion="HG1 nearest owner=Thr-N atom13905 AND C-OG<1.7A AND C-N<1.7A",
        carbonyl_CO_length="Unavailable in frozen scalar observables; not tested",
        interpretation="Geometry-screened TI1-like branch candidates only; neither a verified protonation/connectivity assignment nor a stable minimum",
        raw_frame_counts_are_biased=True,
        total_biased_frames=sum(r["raw_frames"] for r in counts),
        hg1_on_thr_n_biased_frames=sum(r["hg1_on_thr_n_biased_frames"] for r in counts),
        ti1_like_biased_frames=sum(r["ti1_like_biased_frames"] for r in counts),
        hg1_on_thr_n_MBAR_mass=float(weights[np.array(selected_owner)].sum()),
        ti1_like_MBAR_mass=float(weights[np.array(selected_ti)].sum()),
        weighted_mass_limit="Exact deterministic sample selection and common-Hamiltonian MBAR weights; exploratory only with targetESS28.28 and unestablished slow-branch sampling",
        window_counts=counts, frame_locators=locations)
    (p / "ti1_like_branch_diagnostic.json").write_text(json.dumps(result, indent=2) + "\n")
    with (p / "ti1_like_branch_per_window.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(counts[0]))
        writer.writeheader()
        writer.writerows(counts)
    print(json.dumps({k:v for k,v in result.items() if k not in ["window_counts", "frame_locators"]},indent=2))


if __name__ == "__main__":
    main()
