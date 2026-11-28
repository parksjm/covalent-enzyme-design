"""Source-faithful protein-only templates for eight mature Nyl12 chains."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .ligandmpnn_pilot_audit import THREE_TO_ONE

ATOM_COLUMNS = (
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num",
)


def crystal_records(source: Path) -> list[dict[str, str]]:
    """Read original first-model protein atom records, retaining source values.

    Args:
        source: Original crystallographic mmCIF.

    Returns:
        Heavy protein records, using the first/A alternate location only.
    """
    columns = []
    rows = []
    active = False
    for line in source.read_text().splitlines():
        if line.startswith("_atom_site."):
            columns.append(line.strip().split(".")[1])
            active = True
        elif active and line.strip() == "#":
            break
        elif active and line.strip():
            d = dict(zip(columns, line.split(), strict=True))
            if (
                d["group_PDB"] == "ATOM"
                and d["label_comp_id"] in THREE_TO_ONE
                and d["type_symbol"] not in ("H", "D")
                and d["label_alt_id"] in (".", "?", "A")
                and d["pdbx_PDB_model_num"] == "1"
            ):
                rows.append(d)
    if not rows:
        raise ValueError("No original crystal protein atoms")
    return rows


def write_template(rows: list[dict[str, str]], output: Path) -> None:
    """Write a protein-only standard mmCIF with observed-residue sequences.

    Args:
        rows: Canonical mature-chain atom records from the original crystal.
        output: New template destination.

    Notes:
        Unresolved residues are absent rather than assigned coordinates.
        Explicit query/template maps preserve all corresponding query gaps.
    """
    if output.exists():
        raise FileExistsError(output)
    if not rows:
        raise ValueError("Empty template")
    chains = sorted({r["auth_asym_id"] for r in rows})
    residues = {
        c: dict(
            sorted(
                {
                    int(r["auth_seq_id"]): r["label_comp_id"]
                    for r in rows
                    if r["auth_asym_id"] == c
                }.items()
            )
        )
        for c in chains
    }
    lookup = {c: {r: i for i, r in enumerate(residues[c], 1)} for c in chains}
    lines = [
        "data_NYL12_PROTEIN_ONLY_TEMPLATE",
        "#",
        "_entry.id NYL12_PROTEIN_ONLY_TEMPLATE",
        "_exptl.method 'X-RAY DIFFRACTION'",
        "#",
        "loop_",
        "_entity.id",
        "_entity.type",
    ]
    lines += [f"{i} polymer" for i, c in enumerate(chains, 1)]
    lines += ["#", "loop_", "_struct_asym.id", "_struct_asym.entity_id"] + [
        f"{c} {i}" for i, c in enumerate(chains, 1)
    ]
    lines += [
        "#",
        "loop_",
        "_entity_poly.entity_id",
        "_entity_poly.type",
        "_entity_poly.pdbx_strand_id",
    ] + [f"{i} 'polypeptide(L)' {c}" for i, c in enumerate(chains, 1)]
    lines += [
        "#",
        "loop_",
        "_entity_poly_seq.entity_id",
        "_entity_poly_seq.num",
        "_entity_poly_seq.mon_id",
        "_entity_poly_seq.hetero",
    ]
    lines += [
        f"{i} {lookup[c][r]} {residue} n"
        for i, c in enumerate(chains, 1)
        for r, residue in residues[c].items()
    ]
    scheme_columns = (
        "asym_id",
        "entity_id",
        "seq_id",
        "mon_id",
        "ndb_seq_num",
        "pdb_seq_num",
        "auth_seq_num",
        "pdb_mon_id",
        "auth_mon_id",
        "pdb_strand_id",
        "pdb_ins_code",
    )
    lines += ["#", "loop_"] + ["_pdbx_poly_seq_scheme." + c for c in scheme_columns]
    lines += [
        f"{c} {i} {lookup[c][r]} {residue} {lookup[c][r]} {r} {r} {residue} {residue} {c} ."
        for i, c in enumerate(chains, 1)
        for r, residue in residues[c].items()
    ]
    lines += ["#", "loop_"] + ["_atom_site." + c for c in ATOM_COLUMNS]
    for serial, original in enumerate(rows, 1):
        d = dict(original)
        c = d["auth_asym_id"]
        d.update(
            id=str(serial),
            label_asym_id=c,
            label_entity_id=str(chains.index(c) + 1),
            label_seq_id=str(lookup[c][int(d["auth_seq_id"])]),
            label_alt_id=".",
            pdbx_PDB_ins_code="?",
        )
        lines.append(" ".join(d[k] for k in ATOM_COLUMNS))
    output.write_text("\n".join(lines + ["#"]) + "\n")


def prepare_crystal_templates(
    source: Path, sequences: dict[str, str], output: Path
) -> dict:
    """Split actual precursor copies into eight mature protein-only templates.

    Args:
        source: Unmodified n12c crystal with author chains A/F/H/V.
        sequences: Canonical full-query alpha and beta sequences.
        output: New artifact directory.

    Returns:
        Exact observed-residue alignments, gaps, substitutions and provenance.
    """
    output.mkdir(parents=True, exist_ok=False)
    source_rows = crystal_records(source)
    allrows = []
    alignments = {}
    for target, original in zip("ABCDEFGH", "AFHVAFHV", strict=True):
        role = "alpha" if target in "ABCD" else "beta"
        query = sequences[role]
        selected = []
        seen = set()
        for record in source_rows:
            if record["auth_asym_id"] != original:
                continue
            precursor = int(record["auth_seq_id"])
            if not (
                1 <= precursor <= 232 if role == "alpha" else 233 <= precursor <= 328
            ):
                continue
            d = dict(record)
            local = precursor if role == "alpha" else precursor - 232
            key = (local, d["label_atom_id"])
            if key in seen:
                raise ValueError(f"Duplicate source atom {original}:{key}")
            seen.add(key)
            d.update(auth_asym_id=target, auth_seq_id=str(local), label_asym_id=target)
            selected.append(d)
        residues = dict(
            sorted(
                {
                    int(r["auth_seq_id"]): r["label_comp_id"]
                    for r in selected
                    if r["label_atom_id"] == "CA"
                }.items()
            )
        )
        if not residues:
            raise ValueError(f"Missing crystal role {original}:{role}")
        if any(int(r["auth_seq_id"]) not in residues for r in selected):
            raise ValueError("Protein atom without observed CA")
        substitutions = [
            {
                "query_index0": r - 1,
                "query_residue": query[r - 1],
                "template_residue": THREE_TO_ONE[name],
                "source_precursor_residue": r if role == "alpha" else r + 232,
            }
            for r, name in residues.items()
            if query[r - 1] != THREE_TO_ONE[name]
        ]
        alignments[target] = {
            "source_auth_chain": original,
            "role": role,
            "query_sequence": query,
            "template_sequence": "".join(
                THREE_TO_ONE[name] for name in residues.values()
            ),
            "queryIndices": [r - 1 for r in residues],
            "templateIndices": list(range(len(residues))),
            "missing_query_indices0": [
                r - 1 for r in range(1, len(query) + 1) if r not in residues
            ],
            "substitutions": substitutions,
            "observed_atoms": len(selected),
            "template_file": f"chain_{target}.cif",
        }
        write_template(selected, output / f"chain_{target}.cif")
        allrows.extend(selected)
    write_template(allrows, output / "assembly_eight_mature_protein_chains.cif")
    manifest = {
        "status": "prepared_protein_only_templates_no_predictions",
        "source": str(source.resolve()),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_assembly_preserved": True,
        "chains": alignments,
        "coordinate_policy": "Only original observed protein coordinates, no transformation, interpolation, averaging or mutation; missing residues explicit. H228ALA remains sourceALA/queryPRO substitution.",
        "maturation": "Eight independent protein chains A-H; no struct_conn, no alpha-beta covalent links. beta1Thr is a free mature chain terminus.",
        "ligand_coordinates": 0,
        "waters": 0,
        "additives": 0,
        "templates_forced": False,
        "files_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(output.glob("*.cif"))
        },
    }
    (output / "template_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return manifest
