"""Render scene PNGs and save a portable PyMOL session."""
from pathlib import Path
import sys
from pymol import cmd

HERE = Path(sys._getframe().f_code.co_filename).resolve().parent
exec(compile((HERE / 'view_poses.py').read_text(), str(HERE / 'view_poses.py'), 'exec'))
for scene in ('flipped', 'positive', 'overlay'):
    cmd.scene(scene, 'recall')
    cmd.png(str(HERE / f'{scene}.png'), width=1200, height=900, ray=1)
cmd.scene('flipped', 'recall')
cmd.save(str(HERE / 'comparison.pse'))
