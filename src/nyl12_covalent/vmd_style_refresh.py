"""Build reproducible, inline display corrections for existing VMD movies."""

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from parmed.amber import AmberParm

from .vmd_chain_groups import display_chain_groups
from .vmd_movie_policy import mapped_cc_bonds, render_movie_policy
from .vmd_style_manifests import digest

BEGIN = "# BEGIN NYL12 SHARED MOVIE DISPLAY POLICY"
END = "# END NYL12 SHARED MOVIE DISPLAY POLICY"


def topology_policy(
    topology: Any,
    qm_indices: Iterable[int],
    *,
    molecule: str = "[molinfo top]",
    original_topology: Any | None = None,
    original_ids_1based: list[int] | None = None,
) -> str:
    """Render the shared policy from validated physical atom identities.

    Args:
        topology: Displayed ParmEd topology, never modified.
        qm_indices: Explicit displayed zero-based physical QM atoms.
        molecule: Tcl expression identifying the movie's molecule.
        original_topology: Full source topology for compact displays.
        original_ids_1based: Full-source atom ID at each displayed index.

    Returns:
        Self-contained Tcl with explicit physical bonds and cartoon groups.

    Raises:
        ValueError: If supplied mappings or physical selections are invalid.
    """
    qm = list(qm_indices)
    if any(type(i) is not int or not 0 <= i < len(topology.atoms) for i in qm):
        raise ValueError("QM index is not a physical displayed atom")
    source = topology if original_topology is None else original_topology
    if (original_topology is None) != (original_ids_1based is None):
        raise ValueError("Original topology and atom mapping must be supplied together")
    ids = (
        list(range(1, len(topology.atoms) + 1))
        if original_ids_1based is None
        else original_ids_1based
    )
    if len(ids) != len(topology.atoms) or len(set(ids)) != len(ids):
        raise ValueError("Original mapping must cover each displayed atom exactly once")
    mapping = {}
    for displayed, original in enumerate(ids):
        if type(original) is not int or not 1 <= original <= len(source.atoms):
            raise ValueError("Original atom ID lies outside source topology")
        atom, native = topology.atoms[displayed], source.atoms[original - 1]
        if (atom.name, atom.atomic_number, atom.residue.name) != (
            native.name,
            native.atomic_number,
            native.residue.name,
        ):
            raise ValueError(
                f"Physical atom identity mismatch at display index {displayed}"
            )
        mapping[original] = displayed
    substrate = [
        displayed
        for original, displayed in mapping.items()
        if source.atoms[original - 1].residue.name in {"LIG", "P6N", "P6I", "P6C"}
    ]
    return render_movie_policy(
        qm,
        substrate,
        mapped_cc_bonds(source, {old - 1: new for old, new in mapping.items()}),
        molecule=molecule,
        chain_groups=display_chain_groups(source, mapping),
    )


def corrected_viewer(text: str, policy: str) -> str:
    """Return a viewer with one replaceable shared display-policy block.

    Args:
        text: Existing viewer, including original loading, labels and camera.
        policy: Inline policy generated from verified physical atom identities.

    Returns:
        Existing commands plus a single marked policy block.

    Raises:
        ValueError: If existing correction markers are malformed.
    """
    if text.count(BEGIN) != text.count(END) or text.count(BEGIN) > 1:
        raise ValueError("Malformed existing display-policy markers")
    if BEGIN in text:
        start, stop = text.index(BEGIN), text.index(END) + len(END)
        if text.index(END) < start:
            raise ValueError("Reversed existing display-policy markers")
        text = text[:start] + text[stop:].lstrip("\n")
    return text.rstrip() + "\n\n" + BEGIN + "\n" + policy.rstrip() + "\n" + END + "\n"


def append_topology_policy(
    viewer: Path, topology: Any, qm_indices: Iterable[int], **kwargs: Any
) -> None:
    """Append verified inline styling to a newly generated portable viewer.

    Args:
        viewer: Newly generated viewer before checksum publication.
        topology: Corresponding displayed topology.
        qm_indices: Explicit displayed physical QM indices.
        **kwargs: Additional arguments passed to ``topology_policy``.
    """
    viewer.write_text(
        corrected_viewer(
            viewer.read_text(), topology_policy(topology, qm_indices, **kwargs)
        )
    )


def prepare_correction(row: dict[str, Any]) -> tuple[Path, str]:
    """Validate a refresh-plan row and render its corrected viewer in memory.

    Args:
        row: Explicit viewer/topology/QM mapping and evidence information.

    Returns:
        Viewer path and corrected text; this function writes no files.

    Raises:
        ValueError: If paths touch immutable snapshots or evidence is missing.
    """
    viewer, topology = Path(row["viewer"]).resolve(), Path(row["topology"]).resolve()
    if any(
        part in {"source", "source_snapshot", "source_snapshots"}
        for part in viewer.parts
    ):
        raise ValueError("Immutable source viewers cannot be refreshed")
    if viewer.suffix not in {".vmd", ".tcl"} or not row["qm_evidence"]:
        raise ValueError("Explicit viewer and QM evidence are required")
    extra = {}
    if "source_topology" in row:
        extra = dict(
            original_topology=AmberParm(str(row["source_topology"])),
            original_ids_1based=row["original_ids_1based"],
        )
    before = digest(topology) if "topology_sha256" in row else None
    displayed = AmberParm(str(topology))
    if before is not None:
        if digest(topology) != before:
            raise ValueError("Topology changed while preparing viewer correction")
        if before != row["topology_sha256"]:
            from .vmd_style_watch import topology_identity

            if topology_identity(displayed) != row["topology_identity_sha256"]:
                raise ValueError(
                    "Pinned topology identity changed before viewer correction"
                )
    policy = topology_policy(
        displayed,
        row["qm_indices"],
        molecule=row.get("molecule", "[molinfo top]"),
        **extra,
    )
    text = viewer.read_text()
    if "display_caption" in row:
        old = 'puts "Frame0: exact sampled seed. Later frames: optimizer iterations, not MD time or reaction path."'
        if old not in text and row["display_caption"] not in text:
            raise ValueError("Expected original or corrected display caption")
        text = text.replace(old, row["display_caption"])
    return viewer, corrected_viewer(text, policy)


def append_file_policy(viewer: Path, topology: Path, qm_indices: Iterable[int]) -> None:
    """Style a newly exported viewer using its physical topology file.

    Args:
        viewer: Viewer being assembled before checksum publication.
        topology: Corresponding displayed Amber topology file.
        qm_indices: Explicit zero-based displayed physical QM atoms.
    """
    append_topology_policy(viewer, AmberParm(str(topology)), qm_indices)
