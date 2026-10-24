"""File-based and rigid-transform checks for seed inspection artifacts."""

from __future__ import annotations

import csv
import json
import runpy
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import parmed
import pytest
from parmed.topologyobjects import Atom, AtomType, Bond, BondType

from scripts import inspect_nyl12_pa6_es_seeds as module

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/inspect_nyl12_pa6_es_seeds.py"


@pytest.fixture
def inputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple:
    """Write canonical-layout synthetic topology and nine actual coordinate files."""
    structure = parmed.Structure()
    xyz, bonds = [], []

    def add(name: str, element: int, resname: str, number: int, point: tuple) -> int:
        mass = {1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999}[element]
        atom_type = AtomType(f"T{element}", element, mass, element)
        atom_type.set_lj_params(0.1, 1.5)
        atom = Atom(name=name, type=f"T{element}", atomic_number=element, mass=mass)
        atom.atom_type = atom_type
        structure.add_atom(atom, resname, number)
        xyz.append(point)
        return atom.idx

    mapping = {}
    specials = {
        113: ("PHE", [("N", 7, (0, 1, 1), "oxyanion_donor_backbone")]),
        186: ("ASN", [("ND2", 7, (1, 1, 0), "oxyanion_donor_sidechain")]),
        929: (
            "THR",
            [
                ("OG1", 8, (1, 0, 1), "nucleophile"),
                ("N", 7, (0, 0, 2), "thr_n_terminus"),
            ],
        ),
        968: ("ASP", [("OD1", 8, (4, 3, 0), ""), ("OD2", 8, (4, 4, 0), "")]),
        970: ("ASP", [("OD1", 8, (3, 4, 0), ""), ("OD2", 8, (3, 5, 0), "")]),
    }
    for number in range(1, 1313):
        resname, extras = specials.get(number, ("ALA", []))
        add(
            "CA",
            6,
            resname,
            number,
            (20 + number % 10, 20 + number % 17, 20 + number % 23),
        )
        for name, element, point, role in extras:
            index = add(name, element, resname, number, point)
            if role:
                mapping[role] = {"residue": number, "resname": resname, "atom": name}
            if name in {"OG1", "N"} and number == 929:
                hydrogen = add(
                    "HG1" if name == "OG1" else "H1",
                    1,
                    resname,
                    number,
                    tuple(np.asarray(point) + (0, 0.5, 0)),
                )
                bonds.append((index, hydrogen))
    for name, element, point, role in [
        ("C54", 6, (0, 0, 1), "scissile_carbon"),
        ("O49", 8, (0, 1, 0), "scissile_oxygen"),
        ("N59", 7, (0, 2, 0), "leaving_nitrogen"),
        ("C62", 6, (1, 2, 0), ""),
    ]:
        add(name, element, "LIG", 1313, point)
        if role:
            mapping[role] = {"residue": 1313, "resname": "LIG", "atom": name}
    oxygen = add("O", 8, "WAT", 1314, (2, 2, 2))
    for name, point in [("H1", (2.5, 2, 2)), ("H2", (2, 2.5, 2))]:
        bonds.append((oxygen, add(name, 1, "WAT", 1314, point)))
    bond_type = BondType(100, 1)
    structure.bond_types.append(bond_type)
    for left, right in bonds:
        structure.bonds.append(
            Bond(structure.atoms[left], structure.atoms[right], type=bond_type)
        )
    structure.coordinates = np.asarray(xyz, dtype=float)
    parm = parmed.amber.AmberParm.from_structure(structure)
    parm.save(str(tmp_path / "md_solvated.prmtop"))
    reference = tmp_path / "geometry_reference.py"
    reference.write_text("synthetic reference provenance")
    original = module.measure_active_site
    monkeypatch.setattr(
        module,
        "measure_active_site",
        lambda xyz, substrate, state: original(xyz, substrate, state, reference),
    )
    metrics = module._geometry(parm, mapping)["metrics"]
    frames, rows = [], []
    for index in range(9):
        parm.coordinates = np.asarray(xyz) + index * np.asarray((3, -2, 1))
        coordinates = tmp_path / f"frame{index}.rst7"
        parm.save(str(coordinates), format="rst7")
        row = {"model": "synthetic", "replica": "replica_01", "frame": index + 1}
        rows.append(
            {**row, **{column: metrics[key] for key, column in module.METRICS.items()}}
        )
        frames.append(
            {
                **row,
                "full_frame": index + 1,
                "solvated_coordinates": {
                    "path": str(coordinates),
                    "sha256": module._hash(coordinates),
                },
            }
        )
    legacy_dir = tmp_path / "synthetic/replica_01"
    legacy_dir.mkdir(parents=True)
    for index in range(9):
        parm.coordinates = np.asarray(xyz) + index * np.asarray((3, -2, 1))
        selected = [
            _atom
            for _atom in parm.atoms
            if (_atom.residue.idx + 1, _atom.name)
            in {
                (929, "OG1"),
                (1313, "C54"),
                (1313, "O49"),
                (1313, "N59"),
                (1313, "C62"),
            }
        ]
        lines = []
        coords = parm.coordinates
        for atom in selected:
            x, y, z = coords[atom.idx]
            lines.append(
                f"ATOM  {atom.idx + 1:5d} {atom.name:^4s} {atom.residue.name:>3s}  {atom.residue.idx + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  "
            )
        (legacy_dir / f"herschlag_phi_atoms.pdb.{index + 1}").write_text(
            "\n".join(lines) + "\n"
        )
    csv_path = tmp_path / "geometry.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    audit, manifest = tmp_path / "audit.json", tmp_path / "manifest.json"
    audit.write_text(
        json.dumps(
            {"selected_frames": frames, "acceptance": {"reuse_audit_complete": True}}
        )
    )
    manifest.write_text(
        json.dumps(
            {
                "system": str(tmp_path),
                "atom_map": mapping,
                "geometry_csv": str(csv_path),
            }
        )
    )
    data = json.loads(audit.read_text())
    data["sources"] = [
        {"path": str(path), "sha256": module._hash(path)}
        for path in [manifest, csv_path, tmp_path / "md_solvated.prmtop"]
    ]
    audit.write_text(json.dumps(data))
    return audit, manifest, tmp_path / "output"


def test_complete_artifacts_and_portable_viewer(
    inputs: tuple, monkeypatch: pytest.MonkeyPatch
) -> None:
    audit, manifest, output = inputs
    result = module.inspect_seeds(audit, manifest, output)
    assert len(result["frames"]) == 9
    assert len(list(output.glob("replica*.pdb"))) == 9
    assert result["frames"][8]["display_alignment"]["protein_ca_rmsd_A"] < 1e-5
    assert result["frames"][0]["water_inventory"]["local_water_count"] == 1
    first = (output / "replica_01_frame_0001.pdb").read_text()
    last = (output / "replica_01_frame_0009.pdb").read_text()

    def pdb_xyz(text: str) -> np.ndarray:
        return np.asarray(
            [
                [float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]]
                for line in text.splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
            ]
        )

    np.testing.assert_allclose(pdb_xyz(first), pdb_xyz(last), atol=0.001)
    assert "CONECT" in first
    mapping = json.loads((output / "replica_01_frame_0001_mapping.json").read_text())[
        "atoms"
    ]
    assert any(
        row["display_chain"] == "E" and row["display_residue"] == 233 for row in mapping
    )
    assert sum(row["display_chain"] == "W" for row in mapping) == 3
    assert module._display(328 + 232 * 3, {}) == ("E", 328)
    calls = []
    cmd = SimpleNamespace(
        **{
            name: lambda *args, _name=name: calls.append((_name, args))
            for name in ["load", "hide", "show", "color", "set", "disable", "zoom"]
        }
    )
    monkeypatch.setitem(
        sys.modules,
        "pymol",
        SimpleNamespace(cmd=cmd, util=SimpleNamespace(cnc=lambda: None)),
    )
    monkeypatch.chdir(output.parent)
    runpy.run_path(str(output / "view_es_seeds.py"))
    loads = [args for name, args in calls if name == "load"]
    assert len(loads) == 10 and all(Path(args[0]).parent == output for args in loads)
    assert len([1 for name, args in calls if name == "disable"]) == 8
    assert ("show", ("sticks", "(resi 113+186 and chain A) and not context")) in calls
    assert ("set", ("valence", 0)) in calls
    with pytest.raises(FileExistsError):
        module.inspect_seeds(audit, manifest, output)


@pytest.mark.parametrize("mutation", ["geometry", "hash", "count", "csv", "sourcehash"])
def test_rejects_inconsistent_sources(inputs: tuple, mutation: str) -> None:
    audit, manifest, output = inputs
    data = json.loads(audit.read_text())
    if mutation == "sourcehash":
        data["sources"][0]["sha256"] = "bad"
    elif mutation == "hash":
        data["selected_frames"][0]["solvated_coordinates"]["sha256"] = "bad"
    elif mutation == "count":
        data["selected_frames"].pop()
    else:
        csv_path = Path(json.loads(manifest.read_text())["geometry_csv"])
        rows = csv_path.read_text().splitlines()
        if mutation == "geometry":
            fields = rows[1].split(",")
            fields[3] = "99"
            rows[1] = ",".join(fields)
        else:
            rows.pop(1)
        csv_path.write_text("\n".join(rows) + "\n")
    for record in data["sources"]:
        if Path(record["path"]).suffix == ".csv":
            record["sha256"] = module._hash(Path(record["path"]))
    audit.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        module.inspect_seeds(audit, manifest, output)
    if mutation == "geometry":
        assert (output / "geometry_mismatch.json").is_file()


def test_wrapped_phi_and_nonfinite() -> None:
    metrics = {key: 1 for key in module.METRICS}
    reference = {column: 1 for column in module.METRICS.values()}
    metrics["herschlag_phi_projected"] = -179.999
    reference["Herschlag_phi"] = 180.001
    assert module._compare(metrics, reference, -179.999)["herschlag_phi_projected"][
        "passed"
    ]
    metrics["d1_thr929_og1_scissile_c"] = float("nan")
    with pytest.raises(ValueError):
        module._compare(metrics, reference, -179.999)


def test_legacy_precision_retains_raw_phi(inputs: tuple) -> None:
    _, manifest, _ = inputs
    config = json.loads(manifest.read_text())
    topology = Path(config["system"]) / "md_solvated.prmtop"
    parm = parmed.amber.AmberParm(
        str(topology), xyz=str(topology.parent / "frame0.rst7")
    )
    xyz = parm.coordinates.copy()
    xyz[module._atom(parm, config["atom_map"]["nucleophile"]).idx, 0] += 0.00049
    parm.coordinates = xyz
    geometry = module._geometry(parm, config["atom_map"])
    raw = geometry["metrics"]["herschlag_phi_projected"]
    assert abs(raw - geometry["NC_rounded_phi"]) > 0.0001
    parm.coordinates = np.asarray([[float(f"{v:.3f}") for v in point] for point in xyz])
    rounded = module._geometry(parm, config["atom_map"])
    assert geometry["NC_rounded_phi"] == pytest.approx(
        rounded["metrics"]["herschlag_phi_projected"]
    )
    reference = {
        column: geometry["metrics"][key] for key, column in module.METRICS.items()
    }
    reference["Herschlag_phi"] = geometry["NC_rounded_phi"]
    checks = module._compare(geometry["metrics"], reference, geometry["NC_rounded_phi"])
    assert all(row["passed"] for row in checks.values())
    assert checks["herschlag_phi_projected"]["raw_NC_observed"] == raw
    assert checks["herschlag_phi_projected"]["tolerance"] == 0.00001


def test_historical_float32_tie_bound(inputs: tuple) -> None:
    _, manifest, _ = inputs
    config = json.loads(manifest.read_text())
    directory = Path(config["system"])
    parm = parmed.amber.AmberParm(
        str(directory / "md_solvated.prmtop"), xyz=str(directory / "frame0.rst7")
    )
    original = directory / "synthetic/replica_01/herschlag_phi_atoms.pdb.1"
    shifted = directory / "tie.pdb"
    lines = []
    for line in original.read_text().splitlines():
        x = float(line[30:38]) + 53.348
        lines.append(line[:30] + f"{x:8.3f}" + line[38:])
    shifted.write_text("\n".join(lines) + "\n")
    xyz = parm.coordinates.copy() + (53.348, 0, 0)
    oxygen = module._atom(parm, config["atom_map"]["scissile_oxygen"]).idx
    xyz[oxygen, 0] = 53.3474998474
    parm.coordinates = xyz
    result = module.inspect_legacy_phi(parm, shifted)
    assert float(f"{xyz[oxygen, 0]:.3f}") == 53.347
    assert result["maximum_component_difference_A"] == pytest.approx(0.0005001526)
    assert result["phi_degrees"] == pytest.approx(
        module._geometry(parm, config["atom_map"])["NC_rounded_phi"], abs=0.1
    )
    xyz[oxygen, 0] -= 0.001
    parm.coordinates = xyz
    with pytest.raises(ValueError, match="bound exceeded"):
        module.inspect_legacy_phi(parm, shifted)
    with pytest.raises(FileNotFoundError):
        module.inspect_legacy_phi(parm, directory / "missing.pdb")
