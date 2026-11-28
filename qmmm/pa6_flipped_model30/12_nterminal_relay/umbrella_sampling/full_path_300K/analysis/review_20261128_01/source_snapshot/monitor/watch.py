"""Bounded point-only current MBAR publisher; all scientific sampling unchanged."""

from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from collect import inventory, collect
from nyl12_covalent.relay_pooled_estimator import point
from nyl12_covalent.relay_pooled_analysis import write
from nyl12_covalent.relay_umbrella_analysis_io import digest

WORK = Path('/data/jpc/covalent-enzyme-design/qmmm/pa6_flipped_model30/12_nterminal_relay/umbrella_sampling/full_path_300K/analysis/main100_continuation_20261127/current_ext100')
ROOT = Path("/data/jpc/covalent-enzyme-design")


def journal(phase: str, action: str) -> None:
    """Append action intent/outcome, stopping dependent work if logging fails."""
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/record_agent_action.py"),
            "--log",
            str(ROOT / "output_agent_actions/20261021.jsonl"),
            "--agent",
            "relay_current_eligible_monitor",
            "--phase",
            phase,
            "--action",
            action,
            "--permission-mode",
            "escalated",
            "--reason",
            "Authorized dynamic current eligible MBAR only; no native changes",
            "--rollback",
            "Immutable generations retained; stop only derived monitor if needed",
        ],
        check=True,
    )


def publish(selected: dict, number: int) -> None:
    """Fit existing global MBAR and atomically expose a complete immutable generation."""
    eligible = [x for x in selected["supplemental_lineages"] if x["eligible"]]
    folder = (
        WORK
        / f"checkpoint_{number:03d}_main{selected['main_horizon_ps']}ps_supplement{len(eligible)}"
    )
    folder.mkdir()
    data = collect(selected, folder)
    report, definitions, reduced, logweights = point(data)
    write(folder / "point.json", report)
    write(
        folder / "sources.json",
        {k: v for k, v in data.items() if k not in ("q", "counts")},
    )
    np.savez_compressed(
        folder / "weights.npz",
        q=data["q"],
        counts=data["counts"],
        logweights=logweights,
        f_k=report["f_k"],
    )
    primary = json.loads(Path(selected["primary_report"]).read_text())
    p = primary["primary"]
    x = np.asarray(p["centers"])
    f = np.asarray(p["pmf_kcal_mol"], float)
    pooled = np.asarray(report["projections"]["hard"]["pmf"], float)
    ci = primary["bootstrap"]["5ps"]
    lower = np.asarray(ci["lower_95_kcal_mol"], float)
    upper = np.asarray(ci["upper_95_kcal_mol"], float)
    supported = np.asarray(ci["support_fraction"]) == 1
    lower[~supported] = np.nan
    upper[~supported] = np.nan
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11, 8),
        sharex=True,
        constrained_layout=True,
        height_ratios=[3, 1],
    )
    axes[0].fill_between(
        x,
        lower,
        upper,
        color="gray",
        alpha=0.25,
        label="Original57 conditional pointwise95%; 5psblocks",
    )
    axes[0].plot(
        x, f, color="black", label=f"Original57 × {selected['main_horizon_ps']}ps"
    )
    axes[0].plot(
        x,
        pooled,
        color="tab:orange",
        linestyle="--",
        label=f"Current pool: {report['states']} independent lineages; POINT ONLY",
    )
    scope = "; ".join(
        f"{g}: "
        + ",".join(
            str(h)
            for h in sorted(
                {lineage["horizon_ps"] for lineage in eligible if lineage["group"] == g}
            )
        )
        + "ps"
        for g in sorted({lineage["group"] for lineage in eligible})
    )
    axes[0].set(
        title=f"Current eligible relay MBAR · main57@{selected['main_horizon_ps']}ps\n{scope}",
        ylabel="Projected ΔF [kcal/mol]",
    )
    axes[0].legend(fontsize=9)
    axes[0].text(
        0.98,
        0.03,
        "Fixed 5ps discard ONCE per lineage. Conditional CIs do not establish mixing.\nPooled point only; final100 conditional bootstrap monitor remains separate.",
        transform=axes[0].transAxes,
        ha="right",
        fontsize=9,
    )
    axes[1].semilogy(x, p["weight_ess_per_bin"], color="black")
    axes[1].semilogy(
        x, report["projections"]["hard"]["ESS"], color="tab:orange", linestyle="--"
    )
    axes[1].set(
        xlabel="Declared hard-projected progress s",
        ylabel="Weight ESS\n(not temporal ESS)",
        xlim=(-0.075, 1.075),
    )
    fig.savefig(folder / "latest_pmf_comparison.png", dpi=150)
    plt.close(fig)
    peak = int(np.nanargmax(pooled))
    write(
        folder / "latest.json",
        dict(
            main_horizon_ps=selected["main_horizon_ps"],
            independent_lineages=report["states"],
            supplemental_lineages=len(eligible),
            retained_frames=report["frames"],
            unique_bias_potentials=len(
                {json.dumps(s["bias"], sort_keys=True) for s in data["states"]}
            ),
            supplemental_inventory=selected["supplemental_lineages"],
            discarded_ps_once_per_lineage=5,
            pooled_maximum_kcal_mol=float(pooled[peak]),
            pooled_maximum_s=float(x[peak]),
            CI_status="POINT ONLY; pooled uncertainty not computed; original57 band explicitly conditional",
            primary_source=selected["primary_report"],
            sampling_converged=False,
            primary_replacement=False,
        ),
    )
    write(
        folder / "acceptance.json",
        dict(
            process_complete=True,
            native_parity_validated=True,
            connected=report["connected"],
            point_sha256=digest(folder / "point.json"),
            sources_sha256=digest(folder / "sources.json"),
            sampling_converged=False,
        ),
    )
    payloads = [p for p in folder.iterdir() if p.is_file()]
    (folder / "SHA256SUMS.txt").write_text(
        "".join(f"{digest(p)}  {p.name}\n" for p in sorted(payloads))
    )
    temporary = WORK / "current.pending"
    temporary.symlink_to(folder.name, target_is_directory=True)
    temporary.replace(WORK / "current")


def main() -> None:
    """Check native receipts every180s for72h and fit only changed input signatures."""
    pinned = json.loads((WORK / "source_manifest.json").read_text())
    deadline = time.monotonic() + 72 * 3600
    previous = None
    number = 0
    while time.monotonic() < deadline:
        journal(
            "intent",
            "Verify pinned source/prepared identities and inspect accepted lineage prefixes",
        )
        for path, expected in pinned.items():
            if digest(Path(path)) != expected:
                raise ValueError("Pinned updater source/prepared changed:" + path)
        selected = inventory()
        if selected["main_horizon_ps"] < 65:
            write(WORK / "monitor_state.json", dict(status="waiting_for_qualified_extended_main_fit", minimum_horizon_ps=65, old60_outputs_unchanged=True))
            journal("outcome", "Extended pool waits for qualified65ps fixed-protocol fit; no newPMF")
            time.sleep(180)
            continue
        signature = hashlib.sha256(
            json.dumps(selected, sort_keys=True).encode()
        ).hexdigest()
        if signature != previous:
            number += 1
            journal(
                "intent",
                "Fit new current eligible global MBAR generation; source-bound native parity, one math CPU",
            )
            publish(selected, number)
            previous = signature
            journal(
                "outcome",
                "Current eligible point generation complete and atomically published; no native changes",
            )
        write(
            WORK / "monitor_state.json",
            dict(
                status="watching",
                poll_seconds=180,
                signature=signature,
                generation=number,
                main_horizon_ps=selected["main_horizon_ps"],
                included_supplements=sum(
                    lineage["eligible"] for lineage in selected["supplemental_lineages"]
                ),
                native_changes=False,
            ),
        )
        journal(
            "outcome",
            "Receipt inspection complete; unchanged signatures skip fitting; next check in180seconds",
        )
        time.sleep(180)
    raise TimeoutError("Authorized72h monitor bound reached")


if __name__ == "__main__":
    main()
