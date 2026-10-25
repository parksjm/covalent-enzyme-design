"""Portable PyMOL view: optimizer evaluation states are not MD time points."""

from pathlib import Path

from pymol import cmd, util

base = Path(__file__).resolve().parent
cmd.load(str(base / "initial_fullprotein_ligand.pdb.gz"), "context")
cmd.load(str(base / "pocket_endpoint.pdb"), "endpoint")
cmd.load(str(base / "optimization_evaluations.pdb.gz"), "evaluations")
cmd.hide("everything", "all")
cmd.show("cartoon", "context")
cmd.color("gray80", "context")
cmd.set("cartoon_transparency", 0.7, "context")
focus = (
    "(chain A and resi 113+186) or (chain E and resi 233+272+274) or chain I or chain W"
)
cmd.show("sticks", "endpoint and (" + focus + ")")
cmd.show("sticks", "evaluations and (" + focus + ")")
cmd.show("spheres", "endpoint and chain W and elem H")
cmd.set("sphere_scale", 0.18, "chain W and elem H")
cmd.set("valence", 0)
util.cnc("endpoint or evaluations")
cmd.disable("evaluations")
cmd.orient("endpoint and (chain I or (chain E and resi 233))")
cmd.zoom("endpoint", 4)
