"""Load the two aligned WT sites in PyMOL, with individual and overlay scenes."""
from pathlib import Path
import sys
import numpy as np
from pymol import cmd

HERE = Path(sys._getframe().f_code.co_filename).resolve().parent
cmd.reinitialize()
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)
cmd.set('valence', 0)
cmd.set('label_color', 'black')
cmd.set('label_size', 14)
cmd.set('dash_width', 3)
cmd.set('dash_gap', 0.22)
cmd.set('sphere_scale', 0.38)
cmd.set('stick_radius', 0.16)
for obj, beta, ligand, color in [('flipped', 'E', 'I', 'black'), ('positive', 'F', 'J', 'yellow')]:
    cmd.load(str(HERE / f'{obj}.pdb'), obj)
    cmd.hide('everything', obj)
    cmd.show('sticks', obj)
    cmd.color('cyan', f'{obj} and elem C and not chain {ligand}')
    cmd.color('orange', f'{obj} and elem C and chain {ligand}')
    for element, shade in [('N', 'blue'), ('O', 'red'), ('S', 'yellow')]:
        cmd.color(shade, f'{obj} and elem {element}')
    nalpha = f'{obj} and chain {beta} and resi 1 and name N'
    leaving = f'{obj} and chain {ligand} and resi 1 and name N59'
    cmd.show('spheres', f'({nalpha}) or ({leaving})')
    cmd.distance(f'{obj}_N_distance', nalpha, leaving)
    cmd.color(color, f'{obj}_N_distance')
    cmd.label(nalpha, '"Thr233 Nα"')
    cmd.label(leaving, '"N59"')
    cmd.set('label_position', [1.8, -2.0, 0], nalpha)
    cmd.set('label_position', [-1.5, 2.0, 0], leaving)
    cmd.hide('labels', f'{obj}_N_distance')
    for chain, number, name in [(('A' if obj == 'flipped' else 'B'), 113, 'Phe113'), (('A' if obj == 'flipped' else 'B'), 186, 'Asn186'), (beta, 40, 'Asp272'), (beta, 42, 'Asp274')]:
        selection = f'{obj} and chain {chain} and resi {number} and name CA'
        cmd.label(selection, f'"{name}"')
        cmd.set('label_position', {113: [-1.5, 1.8, 0], 186: [1.5, 1.8, 0], 40: [-2.8, -1.0, 0], 42: [2.8, -1.5, 0]}[number], selection)
cmd.orient('flipped or positive')
cmd.zoom('flipped or positive', 4)
# Use one right-handed basis for both rigidly aligned poses.
phe = np.array(cmd.get_atom_coords('flipped and chain A and resi 113 and name CA'))
asn = np.array(cmd.get_atom_coords('flipped and chain A and resi 186 and name CG'))
thr = np.array(cmd.get_atom_coords('flipped and chain E and resi 1 and name N'))
y = (phe + asn) / 2 - thr
y /= np.linalg.norm(y)
x = asn - phe
x -= np.dot(x, y) * y
x /= np.linalg.norm(x)
z = np.cross(x, y)
view = list(cmd.get_view())
view[:9] = np.column_stack((x, y, z)).ravel().tolist()
cmd.set_view(view)
cmd.disable('positive')
cmd.disable('positive_N_distance')
cmd.scene('flipped', 'store')
cmd.enable('positive')
cmd.enable('positive_N_distance')
cmd.disable('flipped')
cmd.disable('flipped_N_distance')
cmd.scene('positive', 'store')
cmd.enable('flipped')
cmd.enable('flipped_N_distance')
cmd.scene('overlay', 'store')
cmd.scene('flipped', 'recall')
