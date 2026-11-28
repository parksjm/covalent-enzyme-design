"""Extract the completed matched comparison without executing scientific jobs."""

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path("/data/jpc/covalent-enzyme-design")
COHORT = (
    ROOT
    / "design/nyl12_relay_ts_ensemble_20261126/broader_pilot/original_backbone_matched_20261127"
)
OUT = Path(__file__).resolve().parent
HASHES: dict[str, str] = {}
SITES = ["AE", "BF", "CG", "DH"]
GROUP = "A+B+C+D+E+F+G+H"


def read(path: Path) -> str:
    """Read and hash an immutable input."""
    data = path.read_bytes()
    HASHES[str(path)] = hashlib.sha256(data).hexdigest()
    return data.decode()


def load(path: Path) -> Any:
    """Read a JSON input while recording its hash."""
    return json.loads(read(path))


def save_csv(name: str, rows: list[dict[str, Any]]) -> None:
    """Write a new table, failing if its destination already exists."""
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with (OUT / name).open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Build six-case and per-site reports from existing evidence."""
    controller = load(COHORT / "controller/summary.json")["cases"]
    old = list(
        csv.DictReader(
            read(COHORT / "comparison/active_site_comparison.csv").splitlines()
        )
    )
    confidence = list(
        csv.DictReader(
            read(
                ROOT
                / "analysis/ipsae2_qualification_20261128/matched_ligand_confidence.csv"
            ).splitlines()
        )
    )
    cases, geometry, ligands, interfaces, all_amides = [], [], [], [], []
    details: dict[str, Any] = {}
    for case in sorted(controller):
        state = controller[case]
        base = COHORT / "relaxation" / case
        raw = load(base / "native_all48_positive_gate.json")
        mapping = load(base / "canonical/mapping_and_acceptance.json")
        relabel = load(base / "canonicalization_receipt.json")
        native_pairs = load(
            COHORT / "scoring/ipsae2/prepared_cases" / case / "physical_pairing.json"
        )["physical_pairs"]
        assembly = load(base / "assembly_precheck.json")
        manifest = load(
            COHORT / "scoring/ipsae2/prepared_cases" / case / "model_manifest.json"
        )
        native_confidence = load(Path(manifest["artifacts"]["confidence"]["path"]))
        pair_rows = list(
            csv.DictReader(
                read(COHORT / "scoring/ipsae2" / case / "chain_pairs.tsv").splitlines(),
                delimiter="\t",
            )
        )
        group_rows = list(
            csv.DictReader(
                read(
                    COHORT / "scoring/ipsae2" / case / "whole_protein.tsv"
                ).splitlines(),
                delimiter="\t",
            )
        )
        group_rows = [
            row for row in group_rows if row["chn1"] == GROUP and row["chn2"] in "IJKL"
        ]
        assert len(group_rows) == 4
        ligand_map = dict(
            sorted(
                set(
                    (row["native"][0], row["mature_local"][0])
                    for row in mapping["atom_mapping"]
                    if row["native"][0] in "IJKL"
                )
            )
        )
        assert ligand_map == {letter: letter for letter in "IJKL"}
        nearest = {row["site"]: row for row in mapping["nearest_d1_per_site"]}
        bijective = len({row["ligand"] for row in nearest.values()}) == 4
        site_by_ligand = (
            {row["ligand"]: site for site, row in nearest.items()} if bijective else {}
        )
        for row in group_rows:
            authoritative = next(
                item
                for item in confidence
                if item["case_id"] == case and item["native_ligand"] == row["chn2"]
            )
            assert (
                abs(
                    float(authoritative["protein_aligned_ipsae"])
                    - float(row["ipsae12"])
                )
                < 1e-10
            )
            ligands.append(
                {
                    "case": case,
                    "native_ligand": row["chn2"],
                    "canonical_ligand": ligand_map[row["chn2"]],
                    "canonical_site": site_by_ligand.get(row["chn2"]),
                    "site_join_status": "unique_nearest_d1_diagnostic"
                    if bijective
                    else "unresolved_nonbijective_invalid_assembly",
                    "protein_group": GROUP,
                    "protein_aligned_direction": "12",
                    "protein_aligned_ipsae2": float(row["ipsae12"]),
                    "ligand_aligned_ipsae21": float(row["ipsae21"]),
                    "native_group_iptm12": row["native_iptm12"],
                    "native_group_iptm21": row["native_iptm21"],
                    "assembly_valid": assembly["accepted"],
                }
            )
        for alpha, beta in native_pairs.items():
            row = next(
                row for row in pair_rows if row["chn1"] == alpha and row["chn2"] == beta
            )
            interfaces.append(
                {
                    "case": case,
                    "canonical_site": alpha
                    + relabel["actual_to_derived_beta_chain"][beta],
                    "native_pair": alpha + beta,
                    **{
                        key: row[key]
                        for key in [
                            "ipsae",
                            "ipsae12",
                            "ipsae21",
                            "native_iptm12",
                            "native_iptm21",
                            "n_scored_contact12",
                            "n_scored_contact21",
                        ]
                    },
                }
            )
        for row in raw["all_amide_site_rows"]:
            all_amides.append(
                {
                    "case": case,
                    "snapshot": "raw_native",
                    "assembly_valid": assembly["accepted"],
                    **row,
                }
            )
        before_tiers, after_tiers = [], []
        for snapshot in ["before", "after"]:
            for site in SITES:
                result: dict[str, Any] = {
                    "case": case,
                    "snapshot": snapshot,
                    "canonical_site": site,
                    "assembly_valid": assembly["accepted"],
                    "case_status": state["status"],
                }
                source = next(
                    row
                    for row in old
                    if row["case"] == case
                    and row["snapshot"] == snapshot
                    and row["site"] == site
                )
                if state["amber_launched"]:
                    prefix = (
                        "selected_"
                        if source["selected_ligand"]
                        else "diagnostic_nearest_d1_"
                    )
                    result.update(
                        {
                            "measurement_status": "measured_amber_input_or_endpoint",
                            "selection": "positive_ab_assignment"
                            if prefix == "selected_"
                            else "nearest_d1_diagnostic",
                            "positive_ab_assigned": source["positive_ab_assigned"]
                            == "True",
                            **{
                                key: source[prefix + key]
                                for key in [
                                    "ligand",
                                    "carbonyl_c",
                                    "carbonyl_o",
                                    "leaving_n",
                                    "tier",
                                    "d1",
                                    "d2",
                                    "d3",
                                    "a_bd",
                                    "herschlag_phi",
                                ]
                            },
                        }
                    )
                elif snapshot == "before":
                    result.update(
                        {
                            "measurement_status": "raw_native_invalid_assembly_diagnostic",
                            "selection": "nearest_d1_diagnostic",
                            "positive_ab_assigned": False,
                            **{
                                key: nearest[site][key]
                                for key in [
                                    "ligand",
                                    "carbonyl_c",
                                    "carbonyl_o",
                                    "leaving_n",
                                    "tier",
                                    "d1",
                                    "d2",
                                    "d3",
                                    "a_bd",
                                    "herschlag_phi",
                                ]
                            },
                        }
                    )
                else:
                    result["measurement_status"] = "unavailable_amber_not_launched"
                geometry.append(result)
                (before_tiers if snapshot == "before" else after_tiers).append(
                    result.get("tier", "NA")
                )
        case_row = {
            "case": case,
            "status": state["status"],
            "amber_launched": state["amber_launched"],
            "raw_positive_ab_sites": raw["assigned_ab_count"],
            "amber_input_positive_ab_sites": state.get("positive_ab_before"),
            "amber_endpoint_positive_ab_sites": state.get("positive_ab_after"),
            "before_tiers_AE_BF_CG_DH": "/".join(before_tiers),
            "after_tiers_AE_BF_CG_DH": "/".join(after_tiers),
            "final_rms_gradient": state.get("final_rms_gradient"),
            "force_converged": state.get("force_converged"),
            "raw_interchain_heavy_pairs_under_1A": assembly[
                "interchain_heavy_pairs_below_threshold"
            ],
            **{
                key: native_confidence[key]
                for key in ["iptm", "protein_iptm", "ligand_iptm", "confidence_score"]
            },
            "grouped_protein_aligned_ipsae2_I_J_K_L": "/".join(
                f"{float(next(row['ipsae12'] for row in group_rows if row['chn2'] == letter)):.6f}"
                for letter in "IJKL"
            ),
        }
        cases.append(case_row)
        details[case] = {
            "case_summary": case_row,
            "native_physical_alpha_beta_pairs": native_pairs,
            "actual_to_canonical_beta_chain": relabel["actual_to_derived_beta_chain"],
            "native_to_canonical_ligand": ligand_map,
            "nearest_d1_ligands_AE_BF_CG_DH": [
                nearest[site]["ligand"] for site in SITES
            ],
            "confidence_site_join_unique": bijective,
        }
    assert (
        len(cases) == 6
        and len(geometry) == 48
        and len(ligands) == 24
        and len(interfaces) == 24
        and len(all_amides) == 288
    )
    save_csv("case_summary.csv", cases)
    save_csv("site_geometry.csv", geometry)
    save_csv("ligand_confidence.csv", ligands)
    save_csv("protein_interface_confidence.csv", interfaces)
    save_csv("raw_all_amide_geometry.csv", all_amides)
    report = {
        "target_phi_sign": 1,
        "site_order": SITES,
        "cases": details,
        "site_geometry": geometry,
        "ligand_confidence": ligands,
        "protein_interface_confidence": interfaces,
        "input_sha256": HASHES,
        "limitations": [
            "Raw failed-assembly geometry is diagnostic only; no Amber endpoint exists.",
            "Nearest-d1 rows do not represent positive-phi Tier A/B acceptance.",
            "Before for accepted cases is serialized Amber input; raw native values are separately retained.",
            "Confidence is prediction confidence, not catalytic evidence; no catalytic confidence cutoff applied.",
            "Failed-case nearest ligand mapping is non-bijective; report all four ligand scores separately without fabricated site assignments.",
            "Force convergence is separate from geometry and assembly acceptance.",
        ],
    }
    (OUT / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# Matched WT / relay 03 / relay 06 comparison",
        "",
        "Sites ordered AE / BF / CG / DH in the canonical derived view. Failed-assembly before values are raw diagnostics; after values are unavailable. Accepted-case before values are the serialized Amber inputs. Signed positive-phi A/B assignment controls acceptance; nearest-d1 diagnostics can fail it.",
        "",
        "| Case | Tiers before → after (all four) | +phi A/B before → after | RMS gradient | <1 Å interchain pairs | Native protein ipTM | Group protein-aligned ipSAE2 I/J/K/L |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in cases:
        lines.append(
            f"| {row['case']} | {row['before_tiers_AE_BF_CG_DH']} → {row['after_tiers_AE_BF_CG_DH']} | {row['raw_positive_ab_sites']} → {row['amber_endpoint_positive_ab_sites'] if row['amber_launched'] else 'not run: invalid assembly'} | {row['final_rms_gradient']} | {row['raw_interchain_heavy_pairs_under_1A']} | {row['protein_iptm']:.4f} | {row['grouped_protein_aligned_ipsae2_I_J_K_L']} |"
        )
    lines += [
        "",
        "All three launched Amber calculations meet the saved force criterion. WT seed 20261126 retains 4/4 positive-phi A/B sites. Relay 03 seed 20261126 retains 3/4 (AE becomes C); seed 20261127 has 0/4 before and after. WT seed 20261127 and both relay 06 seeds fail assembly precheck and have no MM endpoint.",
        "",
        "Group scores use the exact A+B+C+D+E+F+G+H rows and explicit protein-aligned ipsae12 direction. Native ligand identifiers persist through canonicalization, which relabels beta chains only. WT26 site ligands are I/J/K/L; 03-26 I/J/L/K; 03-27 I/J/K/L. Failed assemblies have repeated nearest ligands, so no four-site confidence join is asserted. Per-native-ligand scores and both interface directions are retained separately.",
        "",
        "Detailed d1/d2/d3 (Å), alpha a_bd and Herschlag phi (degrees), selected carbonyl identity and diagnostic/accepted selection are in site_geometry.csv. Native raw 48-way measurements per case are in raw_all_amide_geometry.csv. Full input hashes and mapping evidence are in results.json. No predictions, minimizations or scientific jobs were rerun.",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
