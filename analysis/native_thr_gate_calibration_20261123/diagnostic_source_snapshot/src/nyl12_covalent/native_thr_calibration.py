"""Reference calibration of the unchanged held-out native Thr233 gate."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

RESIDUES = {113: "PHE", 186: "ASN", 233: "THR", 234: "THR", 272: "ASP", 274: "ASP"}
FIT = (113, 186, 272, 274)
BACKBONE = ("N", "CA", "C", "O")
LIMITS = {"og1_shift": 0.75, "nterm_shift": 1.25, "psi_delta": 35.0, "chi1_delta": 20.0}
AtomMap = dict[tuple[int, str], tuple[str, np.ndarray]]


def frozen_geometry(path: Path) -> Any:
    """Load an explicitly supplied immutable metric implementation."""
    spec = importlib.util.spec_from_file_location("frozen_geometry", path)
    if spec is None or spec.loader is None:
        raise ValueError(f"Cannot load metric implementation: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_sites(path: Path) -> tuple[dict[str, AtomMap], dict[str, Any]]:
    """Read author-numbered sites and retain author/label identity metadata.

    Alternate conformers are resolved by highest occupancy, then blank/A/name.
    The policy is explicit and selections are retained in returned metadata.
    """
    data = MMCIF2Dict(str(path))
    cols = {
        k.removeprefix("_atom_site."): v
        for k, v in data.items()
        if k.startswith("_atom_site.")
    }
    sites: dict[str, AtomMap] = {}
    choices: dict[tuple[str, int, str], tuple[float, str]] = {}
    metadata: dict[str, Any] = {}
    for i, atom in enumerate(cols["auth_atom_id"]):
        if cols["pdbx_PDB_model_num"][i] != "1" or cols["group_PDB"][i] != "ATOM":
            continue
        residue = int(cols["auth_seq_id"][i])
        if residue not in RESIDUES or cols["type_symbol"][i] in ("H", "D"):
            continue
        if cols["auth_comp_id"][i] != RESIDUES[residue]:
            raise ValueError(f"Unexpected residue at {residue} in {path}")
        chain = cols["auth_asym_id"][i]
        alt = cols["label_alt_id"][i]
        rank = (-float(cols["occupancy"][i]), "" if alt == "." else alt)
        key = (chain, residue, atom)
        if key in choices and choices[key] <= rank:
            continue
        choices[key] = rank
        point = np.array([float(cols[f"Cartn_{axis}"][i]) for axis in "xyz"])
        if not np.isfinite(point).all():
            raise ValueError("Nonfinite crystal coordinate")
        sites.setdefault(chain, {})[(residue, atom)] = (cols["type_symbol"][i], point)
        entry = metadata.setdefault(
            chain,
            {
                "label_chain": cols["label_asym_id"][i],
                "residue_map": {},
                "selected_alternates": {},
            },
        )
        entry["residue_map"][str(residue)] = cols["label_seq_id"][i]
        entry["selected_alternates"][f"{residue}:{atom}"] = alt
    for chain, entry in metadata.items():
        ca_rows = [
            i
            for i, atom in enumerate(cols["auth_atom_id"])
            if atom == "CA"
            and cols["auth_asym_id"][i] == chain
            and cols["group_PDB"][i] == "ATOM"
        ]
        entry["observed_sequence_auth_numbered"] = {
            cols["auth_seq_id"][i]: cols["auth_comp_id"][i] for i in ca_rows
        }
        entry["residue_232_present"] = "232" in entry["observed_sequence_auth_numbered"]
        entry["required_residues_present"] = all(
            any(key[0] == residue for key in sites[chain]) for residue in RESIDUES
        )
    if not sites:
        raise ValueError(f"No catalytic sites found in {path}")
    return sites, metadata


def gate_failures(row: dict[str, Any]) -> list[str]:
    """Return failed original criteria; equality passes every criterion."""
    values = {key: float(row[key]) for key in LIMITS}
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError("Nonfinite gate metric")
    if values["og1_shift"] < 0 or values["nterm_shift"] < 0:
        raise ValueError("Negative displacement")
    return [key for key, limit in LIMITS.items() if abs(values[key]) > limit]


def measure(
    reference: AtomMap, mobile: AtomMap, geometry: Any, fit: tuple[int, ...] = FIT
) -> dict[str, Any]:
    """Measure original frozen metrics after the specified paired local fit."""
    keys = [(residue, atom) for residue in fit for atom in BACKBONE]
    rotation, center, target, rmsd = geometry.fit_coordinates(
        np.array([mobile[key][1] for key in keys]),
        np.array([reference[key][1] for key in keys]),
    )
    aligned = {
        key: (element, (point - center) @ rotation + target)
        for key, (element, point) in mobile.items()
    }
    thr = geometry.state_metrics(reference, aligned)["thr233_held_out"]
    row = {
        "scaffold_rmsd": rmsd,
        "og1_shift": thr["heavy_atom_displacements_angstrom"]["OG1"],
        "nterm_shift": thr["heavy_atom_displacements_angstrom"]["N"],
        "psi_delta": thr["psi_delta_es_minus_reference_degrees"],
        "chi1_delta": thr["chi1_delta_es_minus_reference_degrees"],
    }
    row["failures"] = gate_failures(row)
    row["native_pass"] = not row["failures"]
    return row


def summarize_controls(path: Path) -> list[dict[str, Any]]:
    """Summarize existing occupied-site rows without duplicating site assignments."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    seen: set[tuple[str, ...]] = set()
    for row in rows:
        identity = tuple(row[key] for key in ("variant", "seed", "stage", "site"))
        if identity in seen:
            raise ValueError(f"Duplicated occupied site: {identity}")
        seen.add(identity)
    groups = sorted({(row["variant"], row["stage"]) for row in rows})
    summaries = []
    for variant, stage in groups:
        group = [
            row for row in rows if (row["variant"], row["stage"]) == (variant, stage)
        ]
        failures = [gate_failures(row) for row in group]
        if any(
            (not fail) != (row["native_pass"] == "True")
            for row, fail in zip(group, failures, strict=True)
        ):
            raise ValueError("Recomputed gate disagrees with stored native_pass")
        result: dict[str, Any] = {
            "variant": variant,
            "stage": stage,
            "sites": len(group),
            "models": len({row["seed"] for row in group}),
            "native_pass": sum(not fail for fail in failures),
            "joint_pass": sum(row["joint_pass"] == "True" for row in group),
        }
        for key in LIMITS:
            values = np.array([float(row[key]) for row in group])
            result[f"{key}_fail"] = sum(key in fail for fail in failures)
            for stat, value in zip(
                ("min", "median", "max"), np.quantile(values, [0, 0.5, 1]), strict=True
            ):
                result[f"{key}_{stat}"] = float(value)
        result["failure_combinations"] = {
            "+".join(fail) or "pass": failures.count(fail) for fail in failures
        }
        summaries.append(result)
    return summaries


def write_table(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a CSV projection; JSON remains the authoritative artifact."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(
            {
                key: json.dumps(value, sort_keys=True)
                if isinstance(value, (dict, list))
                else value
                for key, value in row.items()
            }
            for row in rows
        )


def summarize_md(paths: list[Path]) -> list[dict[str, Any]]:
    """Summarize existing MD scalar CSVs without reading trajectory coordinates."""
    groups: dict[str, list[dict[str, Any]]] = {}
    aliases = {
        "og1_shift": "native_og1_shift_A",
        "nterm_shift": "native_nalpha_shift_A",
        "psi_delta": "psi_delta_deg",
        "chi1_delta": "chi1_delta_deg",
    }
    for path in paths:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                failures = gate_failures(
                    {key: row[value] for key, value in aliases.items()}
                )
                if (not failures) != (row["native_pass"] == "True"):
                    raise ValueError("MD gate disagrees with stored flag")
                groups.setdefault(row["case"].split("__replica_")[0], []).append(
                    {
                        "failures": failures,
                        "fit": float(row["scaffold_rmsd_A"]),
                        **{key: float(row[value]) for key, value in aliases.items()},
                    }
                )
    summaries = []
    for state, rows in sorted(groups.items()):
        result: dict[str, Any] = {
            "state": state,
            "frames": len(rows),
            "native_pass": sum(not row["failures"] for row in rows),
        }
        for selection, selected in (
            ("all", rows),
            ("native_fail", [row for row in rows if row["failures"]]),
            ("native_pass", [row for row in rows if not row["failures"]]),
        ):
            result[f"fit_rmsd_{selection}_q05_q50_q95_A"] = (
                list(
                    map(
                        float,
                        np.quantile(
                            [row["fit"] for row in selected], [0.05, 0.5, 0.95]
                        ),
                    )
                )
                if selected
                else None
            )
        for key in LIMITS:
            result[f"{key}_fail"] = sum(key in row["failures"] for row in rows)
            result[f"{key}_q05_q50_q95"] = list(
                map(float, np.quantile([row[key] for row in rows], [0.05, 0.5, 0.95]))
            )
            result[f"fit_rmsd_vs_{key}_pearson_r"] = float(
                np.corrcoef(
                    [row["fit"] for row in rows], [abs(row[key]) for row in rows]
                )[0, 1]
            )
        summaries.append(result)
    if not summaries:
        raise ValueError("No MD scalar rows supplied")
    return summaries
