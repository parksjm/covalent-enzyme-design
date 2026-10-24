#!/usr/bin/env python3
"""Monitor Nyl12 PA6 active-site geometry from full or compact NetCDF files."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


LOGGER = logging.getLogger(__name__)
CAMPAIGN_DIR = Path(__file__).resolve().parent
SOLVATED_MANIFEST = CAMPAIGN_DIR / "solvated_systems.csv"
AMBERHOME = Path("/data/jpc/code/ambertools26/venv")
CPPTRAJ = AMBERHOME / "bin" / "cpptraj"
FULL_FRAME_SPACING_NS = 0.05
HERSCHLAG_COORDS_PDB = "herschlag_phi_atoms.pdb"

MEASUREMENTS = {
    "d1_d_OG1_scissile_C": "distance d1_d_OG1_scissile_C :929@OG1 :1313@C54",
    "d2_d_Phe113N_scissile_O": "distance d2_d_Phe113N_scissile_O :113@N :1313@O49",
    "d3_d_Asn186ND2_scissile_O": "distance d3_d_Asn186ND2_scissile_O :186@ND2 :1313@O49",
    "a1_alpha_Burgi_Dunitz": "angle a1_alpha_Burgi_Dunitz :929@OG1 :1313@C54 :1313@O49",
}

HERSCHLAG_ATOM_SELECTION = "(:929@OG1,:1313@C54,:1313@O49,:1313@N59,:1313@C62)"
HERSCHLAG_REQUIRED_ATOMS = {
    "thr_og1": (929, "OG1"),
    "carbonyl_c": (1313, "C54"),
    "carbonyl_o": (1313, "O49"),
    "substrate_n": (1313, "N59"),
    "substrate_c": (1313, "C62"),
}

Point = tuple[float, float, float]


@dataclass(frozen=True)
class ReplicaTarget:
    """A trajectory/topology pair to monitor."""

    variant: str
    model: str
    replica: str
    topology: Path
    trajectory: Path
    trajectory_kind: str
    frame_spacing_ns: float


def amber_env() -> dict[str, str]:
    """Return an environment suitable for AmberTools binaries."""
    env = os.environ.copy()
    env["AMBERHOME"] = str(AMBERHOME)
    amber_lib = str(AMBERHOME / "lib")
    current_path = env.get("LD_LIBRARY_PATH")
    env["LD_LIBRARY_PATH"] = (
        amber_lib if not current_path else f"{amber_lib}:{current_path}"
    )
    return env


def configure_logging() -> None:
    """Configure concise command-line logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--kind",
        choices=("full", "compact"),
        default="full",
        help="Monitor full solvated production trajectories or compact dry trajectories.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=10,
        help="Read every Nth frame from the selected trajectory.",
    )
    parser.add_argument(
        "--compact-stride",
        type=int,
        default=10,
        help="Compact trajectory stride encoded in 07_prod_100ns_dry_strideN.nc.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Skip missing trajectories instead of failing.",
    )
    parser.add_argument(
        "--allow-errors",
        action="store_true",
        help="Skip cpptraj failures, useful while NetCDF files are actively written.",
    )
    parser.add_argument(
        "--system",
        help="Only monitor system directories whose names contain this text.",
    )
    parser.add_argument(
        "--replica",
        help="Only monitor one replica name, for example replica_01.",
    )
    return parser.parse_args()


def manifest_rows() -> list[dict[str, str]]:
    """Return rows from the system manifest."""
    with SOLVATED_MANIFEST.open(newline="") as handle:
        return list(csv.DictReader(handle))


def compact_stride_from_path(path: Path, default: int) -> int:
    """Return the compact stride encoded in a trajectory filename."""
    match = re.search(r"stride(\d+)", path.name)
    return int(match.group(1)) if match else default


def build_targets(args: argparse.Namespace) -> list[ReplicaTarget]:
    """Build monitor targets from the campaign manifest."""
    targets: list[ReplicaTarget] = []
    for row in manifest_rows():
        system_dir = Path(row["system_dir"])
        if args.system and args.system not in system_dir.name:
            continue
        for replica_dir in sorted(system_dir.glob("replica_*")):
            if args.replica and replica_dir.name != args.replica:
                continue
            if args.kind == "full":
                topology = system_dir / "md_solvated.prmtop"
                trajectory = replica_dir / "07_prod_100ns.nc"
                frame_spacing = FULL_FRAME_SPACING_NS
            else:
                topology = system_dir / "md_dry.prmtop"
                trajectory = replica_dir / (
                    f"07_prod_100ns_dry_stride{args.compact_stride}.nc"
                )
                frame_spacing = FULL_FRAME_SPACING_NS * compact_stride_from_path(
                    trajectory, args.compact_stride
                )
            if not trajectory.is_file():
                if args.allow_missing:
                    LOGGER.info("Skipping missing trajectory %s", trajectory)
                    continue
                raise FileNotFoundError(trajectory)
            targets.append(
                ReplicaTarget(
                    variant=row["variant"],
                    model=row["model"],
                    replica=replica_dir.name,
                    topology=topology,
                    trajectory=trajectory,
                    trajectory_kind=args.kind,
                    frame_spacing_ns=frame_spacing,
                )
            )
    return targets


def cpptraj_input(target: ReplicaTarget, output_dir: Path, stride: int) -> str:
    """Return cpptraj input for one monitor target."""
    lines = [
        f"parm {target.topology}",
        f"trajin {target.trajectory} 1 last {stride}",
    ]
    if target.trajectory_kind == "full":
        lines.extend(["autoimage anchor :1313", "strip :WAT,Na+,Cl-"])
    for name, command in MEASUREMENTS.items():
        lines.append(f"{command} out {output_dir / f'{name}.dat'}")
    lines.extend(
        [
            f"strip !{HERSCHLAG_ATOM_SELECTION}",
            f"trajout {output_dir / HERSCHLAG_COORDS_PDB} pdb multi nobox",
            "run",
        ]
    )
    return "\n".join(lines) + "\n"


def run_cpptraj(target: ReplicaTarget, output_dir: Path, stride: int) -> None:
    """Run cpptraj for one monitor target."""
    output_dir.mkdir(parents=True, exist_ok=True)
    script_path = output_dir / "monitor_active_site.cpptraj"
    log_path = output_dir / "monitor_active_site.log"
    script_path.write_text(cpptraj_input(target, output_dir, stride))
    with log_path.open("w") as log_handle:
        subprocess.run(
            [str(CPPTRAJ), "-i", str(script_path)],
            check=True,
            env=amber_env(),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )


def read_dat(path: Path) -> dict[int, float]:
    """Read a cpptraj single-measurement data file keyed by frame number."""
    values: dict[int, float] = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")):
            continue
        parts = stripped.split()
        if len(parts) < 2 or not parts[0].isdigit():
            continue
        values[int(parts[0])] = float(parts[1])
    return values


def vector_subtract(first: Point, second: Point) -> Point:
    """Return first minus second for two 3D points."""
    return (first[0] - second[0], first[1] - second[1], first[2] - second[2])


def scale_vector(vector: Point, scalar: float) -> Point:
    """Scale a 3D vector by a scalar."""
    return (vector[0] * scalar, vector[1] * scalar, vector[2] * scalar)


def dot_product(first: Point, second: Point) -> float:
    """Return the dot product of two 3D vectors."""
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2]


def cross_product(first: Point, second: Point) -> Point:
    """Return the cross product of two 3D vectors."""
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def vector_norm(vector: Point) -> float:
    """Return the Euclidean norm of a 3D vector."""
    return math.sqrt(dot_product(vector, vector))


def project_point_to_plane(point: Point, plane_point: Point, plane_normal: Point) -> Point:
    """Project a point onto a plane defined by one point and a normal vector."""
    normal_norm = vector_norm(plane_normal)
    if normal_norm == 0:
        raise ValueError("Cannot project onto a plane with a zero normal vector.")
    unit_normal = scale_vector(plane_normal, 1.0 / normal_norm)
    point_delta = vector_subtract(point, plane_point)
    distance_to_plane = dot_product(point_delta, unit_normal)
    return vector_subtract(point, scale_vector(unit_normal, distance_to_plane))


def signed_dihedral(point_1: Point, point_2: Point, point_3: Point, point_4: Point) -> float:
    """Return a signed dihedral angle in degrees for four points."""
    vector_1 = vector_subtract(point_2, point_1)
    vector_2 = vector_subtract(point_3, point_2)
    vector_3 = vector_subtract(point_4, point_3)
    vector_2_norm = vector_norm(vector_2)
    if vector_2_norm == 0:
        raise ValueError("Cannot calculate a dihedral with coincident middle points.")
    normal_1 = cross_product(vector_1, vector_2)
    normal_2 = cross_product(vector_2, vector_3)
    unit_vector_2 = scale_vector(vector_2, 1.0 / vector_2_norm)
    x_value = dot_product(normal_1, normal_2)
    y_value = dot_product(cross_product(normal_1, unit_vector_2), normal_2)
    return math.degrees(math.atan2(y_value, x_value))


def herschlag_phi(atoms: dict[str, Point]) -> float:
    """Compute the projected Herschlag attack-plane dihedral."""
    nucleophile = atoms["thr_og1"]
    carbonyl_c = atoms["carbonyl_c"]
    carbonyl_o = atoms["carbonyl_o"]
    substrate_n = atoms["substrate_n"]
    substrate_c = atoms["substrate_c"]
    amide_plane_normal = cross_product(
        vector_subtract(substrate_n, carbonyl_o),
        vector_subtract(substrate_c, carbonyl_o),
    )
    projected_carbon = project_point_to_plane(
        carbonyl_c, carbonyl_o, amide_plane_normal
    )
    return signed_dihedral(nucleophile, carbonyl_o, projected_carbon, substrate_c)


def parse_pdb_atom(line: str) -> tuple[int, str, Point] | None:
    """Parse one PDB ATOM/HETATM line from the cpptraj coordinate export."""
    if not line.startswith(("ATOM", "HETATM")):
        return None
    atom_name = line[12:16].strip()
    residue_number = int(line[22:26])
    point = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return residue_number, atom_name, point


def read_herschlag_phi(path: Path) -> dict[int, float]:
    """Read selected-atom PDB frames and calculate projected Herschlag phi."""
    values: dict[int, float] = {}

    def atoms_from_lines(lines: list[str], frame: int) -> dict[str, Point]:
        atoms: dict[str, Point] = {}
        for line in lines:
            parsed = parse_pdb_atom(line)
            if parsed is None:
                continue
            residue_number, atom_name, point = parsed
            for key, required in HERSCHLAG_REQUIRED_ATOMS.items():
                if (residue_number, atom_name) == required:
                    atoms[key] = point
                    break
        missing = sorted(set(HERSCHLAG_REQUIRED_ATOMS) - set(atoms))
        if missing:
            raise ValueError(f"Missing Herschlag atoms in frame {frame}: {missing}")
        return atoms

    split_paths = sorted(
        path.parent.glob(f"{path.name}.*"),
        key=lambda item: int(item.name.rsplit(".", 1)[-1]),
    )
    if split_paths:
        for split_path in split_paths:
            frame = int(split_path.name.rsplit(".", 1)[-1])
            values[frame] = herschlag_phi(
                atoms_from_lines(split_path.read_text().splitlines(), frame)
            )
        return values

    frame = 0
    current_lines: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            current_lines = []
            frame += 1
            continue
        if line.startswith("ENDMDL"):
            values[frame] = herschlag_phi(atoms_from_lines(current_lines, frame))
            current_lines = []
            continue
        current_lines.append(line)
    if current_lines:
        if frame == 0:
            frame = 1
        values[frame] = herschlag_phi(atoms_from_lines(current_lines, frame))
    return values

def merge_timeseries(
    target: ReplicaTarget, output_dir: Path, stride: int
) -> list[dict[str, str]]:
    """Merge cpptraj measurement files into one per-frame table."""
    measurement_values = {
        name: read_dat(output_dir / f"{name}.dat") for name in MEASUREMENTS
    }
    measurement_values["Herschlag_phi"] = read_herschlag_phi(
        output_dir / HERSCHLAG_COORDS_PDB
    )
    frames = sorted(set().union(*(values.keys() for values in measurement_values.values())))
    rows = []
    for frame in frames:
        time_ns = (frame - 1) * stride * target.frame_spacing_ns
        row = {
            "variant": target.variant,
            "model": target.model,
            "replica": target.replica,
            "trajectory_kind": target.trajectory_kind,
            "trajectory": str(target.trajectory),
            "frame": str(frame),
            "time_ns_estimate": f"{time_ns:.3f}",
        }
        for name, values in measurement_values.items():
            row[name] = f"{values[frame]:.6f}" if frame in values else ""
        rows.append(row)
    return rows


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    """Write CSV rows, preserving discovered field order."""
    if not rows:
        return
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]]) -> dict[str, str]:
    """Summarize the last available frame for one target."""
    if not rows:
        return {}
    return rows[-1]


def main() -> None:
    """Run active-site geometry monitoring."""
    configure_logging()
    args = parse_args()
    output_root = CAMPAIGN_DIR / "active_site_monitor" / (
        f"{args.kind}_stride{args.stride}"
    )
    all_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    for target in build_targets(args):
        output_dir = output_root / target.model / target.replica
        LOGGER.info("Monitoring %s %s", target.model, target.replica)
        try:
            run_cpptraj(target, output_dir, args.stride)
        except subprocess.CalledProcessError as exc:
            if args.allow_errors:
                LOGGER.warning("Skipping cpptraj failure for %s: %s", target, exc)
                continue
            raise
        rows = merge_timeseries(target, output_dir, args.stride)
        write_rows(output_dir / "active_site_geometry_timeseries.csv", rows)
        all_rows.extend(rows)
        summary = summarize(rows)
        if summary:
            summary_rows.append(summary)
    write_rows(output_root / "active_site_geometry_timeseries.csv", all_rows)
    write_rows(output_root / "active_site_geometry_latest.csv", summary_rows)
    LOGGER.info("Wrote monitor output under %s", output_root)


if __name__ == "__main__":
    main()
