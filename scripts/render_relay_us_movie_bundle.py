"""Finalize a stripped relay-US export as a portable VMD movie bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

from parmed.amber import AmberParm

from nyl12_covalent.neb_us_movie_vmd import render_movie_vmd
from nyl12_covalent.neb_us_movie_scope import validated_sampling_scope
from nyl12_covalent.path_inspection import dcd_counts
from nyl12_covalent.vmd_movie_policy import POLICY_VERSION

LOGGER = logging.getLogger(__name__)


def _sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _lookup_atom_id(
    topology: Any, residue_number: int, residue_name: str, atom_name: str
) -> int:
    """Resolve one physical atom from the current Amber residue identity."""
    matches = [
        atom
        for atom in topology.atoms
        if atom.residue.idx + 1 == residue_number
        and atom.residue.name == residue_name
        and atom.name == atom_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {residue_name}{residue_number}.{atom_name}; "
            f"found {len(matches)}"
        )
    return matches[0].idx + 1


def _movie_roles(topology: Any, observed: dict[str, int]) -> dict[str, int]:
    """Map and verify reactive and active-site roles against this topology."""
    required = {"carbon", "leaving_n", "og1", "hg1", "nalpha"}
    if not required.issubset(observed):
        raise ValueError(
            f"Missing observable role mappings: {sorted(required - observed.keys())}"
        )
    carbon_id = observed["carbon"]
    carbon = topology.atoms[carbon_id - 1]
    oxygens = [
        atom
        for atom in carbon.bond_partners
        if atom.atomic_number == 8 and atom.residue.idx == carbon.residue.idx
    ]
    if len(oxygens) != 1:
        raise ValueError(f"Expected one bonded carbonyl oxygen at C{carbon_id}")
    roles = {
        "carbon": carbon_id,
        "leaving_n": observed["leaving_n"],
        "og1": observed["og1"],
        "hg1": observed["hg1"],
        "thr_n": observed["nalpha"],
        "oxygen": oxygens[0].idx + 1,
        "phe_n": _lookup_atom_id(topology, 113, "PHE", "N"),
        "asn_nd2": _lookup_atom_id(topology, 186, "ASN", "ND2"),
    }
    expected_names = {
        "carbon": "C",
        "leaving_n": "N",
        "og1": "OG1",
        "hg1": "HG1",
        "thr_n": "N",
        "oxygen": "O",
        "phe_n": "N",
        "asn_nd2": "ND2",
    }
    for role, atom_id in roles.items():
        if type(atom_id) is not int or atom_id < 1 or atom_id > len(topology.atoms):
            raise ValueError(f"Invalid physical atom ID for role {role}: {atom_id!r}")
        atom = topology.atoms[atom_id - 1]
        if atom.name != expected_names[role]:
            raise ValueError(
                f"Role {role} maps to {atom.name}, expected {expected_names[role]}"
            )
    for role in ("og1", "hg1", "thr_n"):
        atom = topology.atoms[roles[role] - 1]
        if atom.residue.name != "THR" or atom.residue.idx + 1 != 929:
            raise ValueError(f"Role {role} is not THR929 in the movie topology")
    return roles


def _run_cpptraj(script: Path, binary: Path, log_path: Path) -> None:
    """Run one deterministic cpptraj endpoint-coordinate probe."""
    environment = os.environ.copy()
    library_dir = str(binary.resolve().parent.parent / "lib")
    existing = environment.get("LD_LIBRARY_PATH")
    environment["LD_LIBRARY_PATH"] = (
        f"{library_dir}:{existing}" if existing else library_dir
    )
    result = subprocess.run(
        [str(binary.resolve()), "-i", str(script.resolve())],
        cwd=script.parent,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    log_path.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"cpptraj endpoint probe failed; see {log_path}")


def _validate_pdb_coordinates(
    path: Path, expected_atoms: int
) -> dict[str, float | int]:
    """Require real, finite, spatially extended coordinates in a probe PDB."""
    coordinates: list[tuple[float, float, float]] = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            coordinates.append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            )
    if len(coordinates) != expected_atoms:
        raise ValueError(f"Probe atom count mismatch in {path}: {len(coordinates)}")
    if not all(math.isfinite(value) for xyz in coordinates for value in xyz):
        raise ValueError(f"Nonfinite endpoint coordinates in {path}")
    spans = [
        max(xyz[axis] for xyz in coordinates) - min(xyz[axis] for xyz in coordinates)
        for axis in range(3)
    ]
    if (
        max(map(abs, (value for xyz in coordinates for value in xyz))) < 1.0
        or max(spans) < 1.0
    ):
        raise ValueError(f"Endpoint coordinates are zero or collapsed in {path}")
    return {"atoms": len(coordinates), "max_span_A": max(spans)}


def render_bundle(
    campaign: Path, bundle: Path, cpptraj: Path, prepared_manifest: Path | None = None
) -> dict[str, Any]:
    """Finalize a complete all-window NEB-US export for portable VMD viewing.

    Args:
        campaign: Campaign root containing ``prepared.json`` and topology.
        bundle: Export directory containing ``part.parm7``, ``part.dcd`` and
            ``manifest.json`` from :func:`export_part`.
        cpptraj: AmberTools cpptraj executable used to probe DCD endpoint frames.
        prepared_manifest: Explicit continuation manifest when exporting beyond
            the original declared production prefix.

    Returns:
        Final bundle manifest with viewer and coordinate-probe evidence.
    """
    campaign, bundle = campaign.resolve(), bundle.resolve()
    prepared_path = (prepared_manifest or campaign / "prepared.json").resolve(
        strict=True
    )
    prepared = json.loads(prepared_path.read_text())
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    labels = [window["label"] for window in prepared["windows"]]
    if manifest.get("windows") != labels:
        raise ValueError("Exported window order does not match prepared campaign")
    chart_counts = Counter(window["chart"] for window in prepared["windows"])
    if set(chart_counts) != {"activation", "delivery"}:
        raise ValueError(f"Expected both relay charts; found {dict(chart_counts)}")
    scope = validated_sampling_scope(campaign, prepared, manifest)

    source_topology = (campaign / prepared["topology"]).resolve(strict=True)
    topology = AmberParm(str(source_topology))
    raw_roles = manifest.get("roles_1based")
    if not isinstance(raw_roles, dict):
        raise ValueError("Export manifest lacks its native observable role mapping")
    roles = _movie_roles(topology, raw_roles)
    kept_ids = manifest["original_atom_ids_1based"]
    qm_ids = manifest["qm_atom_ids_1based"]
    if not set(qm_ids).issubset(kept_ids):
        raise ValueError("QM atom IDs are not retained by the stripped topology")
    if not set(roles.values()).issubset(kept_ids):
        raise ValueError(
            "One or more displayed roles are absent from the stripped topology"
        )
    if not set(raw_roles.values()).issubset(qm_ids):
        raise ValueError("Observable reactive roles do not match exported QM indices")

    old_names = (("part.parm7", "path.parm7"), ("part.dcd", "path.dcd"))
    for old_name, new_name in old_names:
        old_path, new_path = bundle / old_name, bundle / new_name
        if not old_path.is_file() or new_path.exists():
            raise FileNotFoundError(
                f"Expected one export file and no alias: {old_path}"
            )
        old_path.replace(new_path)
    stripped = AmberParm(str(bundle / "path.parm7"))
    frames, atoms = dcd_counts(bundle / "path.dcd")
    if frames != manifest["frames"] or atoms != len(stripped.atoms):
        raise ValueError("Final DCD dimensions do not match its stripped topology")

    mapping = {original: index for index, original in enumerate(kept_ids)}
    qm_local = [mapping[atom_id] for atom_id in qm_ids]
    index_map_path = bundle / "qm_indices.json"
    index_map = {
        "topology_sha256": _sha256(bundle / "path.parm7"),
        "original_qm_atom_ids_1based": qm_ids,
        "stripped_qm_atom_indices_0based": qm_local,
        "original_to_stripped_index_0based": mapping,
        "roles_1based": roles,
    }
    index_map_path.write_text(json.dumps(index_map, indent=2) + "\n")
    viewer_path = bundle / "view_path.vmd"
    viewer = render_movie_vmd(topology, kept_ids, qm_ids, roles, frames)
    viewer_path.write_text(viewer)

    endpoint_checks: dict[str, Any] = {}
    for label, frame in (("first", 1), ("last", frames)):
        probe_input = bundle / f"coordinate_probe_{label}.in"
        pdb_path = bundle / f"coordinate_probe_{label}.pdb"
        probe_input.write_text(
            "parm path.parm7\n"
            f"trajin path.dcd {frame} {frame} 1\n"
            f"trajout {pdb_path.name} pdb\nrun\nquit\n"
        )
        _run_cpptraj(probe_input, cpptraj, bundle / f"coordinate_probe_{label}.log")
        endpoint_checks[label] = {
            "frame_1based": frame,
            **_validate_pdb_coordinates(pdb_path, atoms),
            "pdb_sha256": _sha256(pdb_path),
        }

    manifest.update(
        {
            "topology_file": "path.parm7",
            "trajectory_file": "path.dcd",
            "topology_sha256": _sha256(bundle / "path.parm7"),
            "trajectory_sha256": _sha256(bundle / "path.dcd"),
            "viewer_file": "view_path.vmd",
            "viewer_sha256": _sha256(viewer_path),
            "qm_indices_file": index_map_path.name,
            "qm_indices_sha256": _sha256(index_map_path),
            "display_roles_1based": roles,
            "chart_window_counts": dict(chart_counts),
            **scope,
            "trajectory_semantics": "ordered biased windows; not a continuous reactive event",
            "sampling_converged": False,
            "coordinate_probe": endpoint_checks,
            "display_policy_version": POLICY_VERSION,
            "viewer_source_sha256": {
                str(Path(__file__).resolve()): _sha256(Path(__file__).resolve())
            },
        }
    )
    manifest.pop("part.dcd", None)
    manifest.pop("part.parm7", None)
    manifest["file_sha256"] = {
        name: _sha256(bundle / name)
        for name in ("path.parm7", "path.dcd", "view_path.vmd", "qm_indices.json")
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    (bundle / "README.md").write_text(
        "# N-terminal relay umbrella movie\n\n"
        f"The movie concatenates {len(labels)} ordered biased windows: "
        f"{chart_counts['activation']} activation and {chart_counts['delivery']} delivery. "
        f"Each window contributes {scope['sampled_frames_per_window']} spaced snapshots "
        f"from its first {scope['production_blocks_per_window']} accepted production blocks "
        f"({scope['sampled_horizon_per_window_ps']:g} ps/window). "
        "Window boundaries mark independent biased simulations; "
        "the concatenation is not a continuous reactive event and does not imply PMF "
        "sampling convergence.\n\n"
        f"Open `view_path.vmd` in VMD. Display policy: {POLICY_VERSION}. "
        "`path.parm7`, `path.dcd`, and `qm_indices.json` are the matching stripped "
        "viewer topology, trajectory, and atom-index map.\n\n"
        "The manifest records source hashes, window/frame order, topology mapping, "
        "and first/last DCD coordinate probes.\n"
    )
    checksum_files = sorted(
        path.name
        for path in bundle.iterdir()
        if path.is_file() and path.name != "SHA256SUMS.txt"
    )
    (bundle / "SHA256SUMS.txt").write_text(
        "".join(f"{_sha256(bundle / name)}  {name}\n" for name in checksum_files)
    )
    LOGGER.info("Finalized %d-frame movie bundle at %s", frames, bundle)
    return manifest


def main() -> None:
    """Parse paths and finalize the portable VMD movie bundle."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--cpptraj", type=Path, required=True)
    parser.add_argument("--prepared-manifest", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    render_bundle(args.campaign, args.bundle, args.cpptraj, args.prepared_manifest)


if __name__ == "__main__":
    main()
