"""Tests for matched four-ligand topology preparation and validation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pytest

from nyl12_covalent import design_relaxation_topology as topology_tools

RELEASE = Path(
    "/data/jpc/polyamide-amber-ff/parameters/releases/"
    "pa6_l4_modular_resp_gaff2_v1_20261031"
)


def test_leap_input_uses_neutral_thr_and_optional_salt(tmp_path: Path) -> None:
    """Build dry, unsalted and explicitly salted preparation scripts."""
    dry = topology_tools.leap_input(RELEASE, None)
    salted = topology_tools.leap_input(RELEASE, (114, 106))
    assert "source leaprc.protein.ff19SB" in dry
    assert "source leaprc.water.tip3p" in dry
    assert f"source {RELEASE}/leaprc.pa6_shared_p6i" in dry
    assert "loadpdb complex_heavy_modular.pdb" in dry
    assert dry.count("remove complex complex.") == 4
    assert dry.count("charge -0.4834") == 4
    assert "set complex.929.H1 charge 0.1224" in dry
    assert "solvateOct complex TIP3PBOX 12.0" in dry
    assert "noions.prmtop noions.rst7" in dry
    assert "addIonsRand complex Na+ 114 Cl- 106 6.0" in salted
    assert "complex.prmtop complex.rst7" in salted
    assert "-0.4844" not in salted


@dataclass
class Atom:
    """Minimal ParmEd-compatible atom for validator unit fixtures."""

    name: str
    type: str
    charge: float
    residue: "Residue | None" = None


@dataclass
class Residue:
    """Minimal ParmEd-compatible residue for validator unit fixtures."""

    name: str
    idx: int
    atoms: list[Atom] = field(default_factory=list)

    def add(self, atom: Atom) -> None:
        atom.residue = self
        self.atoms.append(atom)


@dataclass
class Bond:
    """Minimal ParmEd-compatible bond."""

    atom1: Atom
    atom2: Atom


class System:
    """Minimal topology collection accepted by the semantic validator."""

    def __init__(self, solvated: bool = False) -> None:
        self.residues: list[Residue] = []
        self.bonds: list[Bond] = []
        self.parm_data = {
            "RADIUS_SET": ["ArgH and AspGluO modified Bondi2 radii (mbondi3)"]
        }
        for number in range(1, 1329):
            if number in topology_tools.BETA_THR_RESIDUES:
                residue = self._neutral_thr(number)
            elif number in topology_tools.ASP_RESIDUES:
                residue = self._charged("ASP", number - 1, "ASP", -0.999998)
            elif number in (100, 332, 564, 796):
                residue = self._charged("LYS", number - 1, "NZ", 1.0)
            elif number in (1, 233, 465, 697):
                residue = self._charged("ALA", number - 1, "CA", 1.0)
            elif number in (232, 464, 696, 928, 1024, 1120, 1216, 1312):
                residue = self._charged("ALA", number - 1, "CA", -1.0)
            elif number >= topology_tools.FIRST_PA6_RESIDUE:
                template_name = topology_tools.PA6_SEQUENCE[(number - 1313) % 4]
                residue = self._pa6(template_name, number - 1)
            else:
                residue = self._charged("ALA", number - 1, "CA", 0.0)
            self.residues.append(residue)
        for start in (1313, 1317, 1321, 1325):
            for number in range(start, start + 3):
                carbon = next(
                    a for a in self.residues[number - 1].atoms if a.name == "C"
                )
                nitrogen = next(a for a in self.residues[number].atoms if a.name == "N")
                self.bonds.append(Bond(carbon, nitrogen))
        if solvated:
            water = self._charged("WAT", len(self.residues), "O", 0.0)
            self.residues.append(water)
            for index in range(8):
                self.residues.append(
                    self._charged("Na+", len(self.residues), "Na+", 1.0)
                )

    @staticmethod
    def _charged(resname: str, idx: int, atom_name: str, charge: float) -> Residue:
        residue = Residue(resname, idx)
        residue.add(Atom(atom_name, "c3", charge))
        return residue

    @staticmethod
    def _neutral_thr(number: int) -> Residue:
        residue = Residue("THR", number - 1)
        for name, charge in (
            ("N", -0.4834),
            ("H1", 0.1224),
            ("H2", 0.1224),
            ("HG1", 0.2386),
        ):
            residue.add(Atom(name, "n", charge))
        residue.add(Atom("CA", "c3", 0.0))
        return residue

    @staticmethod
    def _pa6(name: str, idx: int) -> Residue:
        residue = Residue(name, idx)
        template = topology_tools._mol2_atoms(RELEASE / f"{name}.mol2")
        for atom_name, (atom_type, charge) in template.items():
            residue.add(Atom(atom_name, atom_type, charge))
        return residue

    @property
    def atoms(self) -> list[Atom]:
        return [atom for residue in self.residues for atom in residue.atoms]


def test_validator_accepts_exact_four_ligand_dry_topology(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Validate released PA6 copies, twelve links, matched Thr and net charge."""
    path = tmp_path / "dry.prmtop"
    path.touch()
    system = System()
    monkeypatch.setattr(topology_tools.parmed, "load_file", lambda _path: system)
    result = topology_tools.validate_four_ligand_topology(path, RELEASE, solvated=False)
    assert result["pa6_residue_count"] == 16
    assert result["pa6_atom_count"] > result["pa6_residue_count"]
    assert result["pa6_interresidue_link_count"] == 12
    assert result["pa6_charge"] == pytest.approx(-16e-6, abs=1e-8)
    assert result["pa6_copy_charges"] == pytest.approx([-4e-6] * 4, abs=1e-8)
    assert result["total_charge"] == pytest.approx(-8.0)
    assert result["formal_solute_charge"] == -8
    assert len(result["neutral_thr_residue_charges"]) == 4
    assert result["radius_set"].endswith("mbondi3)")


def test_validator_rejects_bad_thr_charge_and_missing_pa6_link(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reject the historical N charge and an omitted polymer C–N connection."""
    path = tmp_path / "dry.prmtop"
    path.touch()
    system = System()
    system.residues[928].atoms[0].charge = -0.4844
    monkeypatch.setattr(topology_tools.parmed, "load_file", lambda _path: system)
    with pytest.raises(ValueError, match="Thr 929:N charge mismatch"):
        topology_tools.validate_four_ligand_topology(path, RELEASE, solvated=False)

    system = System()
    system.bonds.pop()
    monkeypatch.setattr(topology_tools.parmed, "load_file", lambda _path: system)
    with pytest.raises(ValueError, match="twelve within-copy PA6 C-N links"):
        topology_tools.validate_four_ligand_topology(path, RELEASE, solvated=False)


def test_validator_checks_solvated_neutrality(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Accept a neutralized solvated topology with the declared mbondi3 set."""
    path = tmp_path / "solvated.prmtop"
    path.touch()
    system = System(solvated=True)
    monkeypatch.setattr(topology_tools.parmed, "load_file", lambda _path: system)
    result = topology_tools.validate_four_ligand_topology(path, RELEASE, solvated=True)
    assert result["total_charge"] == pytest.approx(0.0)
    assert result["formal_solute_charge"] == -8
    assert result["solvated"] is True


def test_validator_rejects_counterions_that_do_not_neutralize_solute(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Require explicit sodium/chloride inventory to cancel formal solute charge."""
    path = tmp_path / "solvated.prmtop"
    path.touch()
    system = System(solvated=True)
    system.residues.pop()
    monkeypatch.setattr(topology_tools.parmed, "load_file", lambda _path: system)
    with pytest.raises(ValueError, match="Counterion inventory"):
        topology_tools.validate_four_ligand_topology(path, RELEASE, solvated=True)


def test_per_residue_charge_gate_remains_tighter_than_global_roundoff(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Accept the observed four-copy rounding residual but reject a larger error."""
    path = tmp_path / "dry.prmtop"
    path.touch()
    system = System()
    for number in topology_tools.ASP_RESIDUES:
        system.residues[number - 1].atoms[0].charge = -1.0
    monkeypatch.setattr(topology_tools.parmed, "load_file", lambda _path: system)
    result = topology_tools.validate_four_ligand_topology(path, RELEASE, solvated=False)
    assert result["total_charge"] == pytest.approx(-8.000016, abs=1e-8)
    system.residues[0].atoms[0].charge += 0.001
    with pytest.raises(ValueError, match="Protein residue 1 ALA charge"):
        topology_tools.validate_four_ligand_topology(path, RELEASE, solvated=False)
