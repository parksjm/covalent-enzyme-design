"""Portable PyMOL inspection views for canonical native controls."""

from __future__ import annotations

from pathlib import Path


def write_views(output: Path, models: list[str], fit_residues: list[int]) -> None:
    """Write relative-path PyMOL loader and one scene per prediction.

    Args:
        output: Curated inspection bundle directory.
        models: Relative model directory/object names.
        fit_residues: Verified alpha CA identities shared with crystal chain A.
    """
    script = '''"""Run from the bundle directory in PyMOL; no coordinate minimization."""
from pathlib import Path
import sys
import numpy as np
from pymol import cmd
HERE = Path(sys._getframe().f_code.co_filename).resolve().parent
cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("valence", 0)
cmd.set("label_color", "black")
cmd.set("stick_radius", 0.16)
cmd.load(str(HERE / "crystal_reference.pdb"), "crystal_reference")
cmd.hide("everything", "all")
cmd.show("cartoon", "crystal_reference")
cmd.color("gray70", "crystal_reference")
MODELS = __MODELS__
RESIDUES = "__RESIDUES__"
for name in MODELS:
    cmd.load(str(HERE / name / "canonical_precursor.pdb"), name)
    # One rigid protein-core fit only. No pocket, Thr or ligand atoms select it.
    cmd.pair_fit(f"{name} and chain A and name CA and resi {RESIDUES}",
                 f"crystal_reference and chain A and name CA and resi {RESIDUES}")
    cmd.show("cartoon", name)
    cmd.color("cyan", name)
    pocket = f"{name} and ((chain A+B+C+D and resi 113+186) or (chain E+F+G+H and resi 233+272+274) or chain I+J+K+L)"
    cmd.show("sticks", pocket)
    cmd.color("orange", f"{name} and chain I+J+K+L and elem C")
    cmd.color("blue", f"{name} and elem N")
    cmd.color("red", f"{name} and elem O")
    cmd.color("yellow", f"{name} and elem S")
cmd.show("sticks", "crystal_reference and chain A and resi 113+186+233+272+274")
phe = np.array(cmd.get_atom_coords("crystal_reference and chain A and resi 113 and name CA"))
asn = np.array(cmd.get_atom_coords("crystal_reference and chain A and resi 186 and name CG"))
thr = np.array(cmd.get_atom_coords("crystal_reference and chain A and resi 233 and name N"))
y = (phe + asn) / 2 - thr
y /= np.linalg.norm(y)
x = asn - phe
x -= np.dot(x, y) * y
x /= np.linalg.norm(x)
z = np.cross(x, y)
for name in MODELS:
    cmd.disable("all")
    cmd.enable("crystal_reference")
    cmd.enable(name)
    cmd.orient(f"{name} or crystal_reference")
    cmd.zoom(f"{name} or crystal_reference", 4)
    view = list(cmd.get_view())
    view[:9] = np.column_stack((x, y, z)).ravel().tolist()
    cmd.set_view(view)
    cmd.scene(name, "store")
cmd.scene(MODELS[0], "recall")
'''.replace("__MODELS__", repr(models)).replace(
        "__RESIDUES__", "+".join(map(str, fit_residues))
    )
    (output / "view_controls.py").write_text(script)
    (output / "load_all_controls.pml").write_text("run view_controls.py\n")
    for model in models:
        (output / f"view_{model}.pml").write_text(
            f"run view_controls.py\nscene {model}, recall\n"
        )
