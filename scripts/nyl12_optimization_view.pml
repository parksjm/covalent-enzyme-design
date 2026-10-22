# Full-system replay snapshots. Failed TI1 diagnostic, including trial points.
delete optimization
load optimization_trajectory.pdb.gz, optimization, discrete=0
hide everything
show cartoon, chain A+E
show sticks, chain I or (chain E and resi 1+2)
show sticks, resi 113+186 and chain A
show sticks, chain J
show spheres, chain J and name O
set sphere_scale, 0.25, chain J
set valence, 0
util.cnc
set movie_fps, 5
zoom chain I+J or (chain E and resi 1), 8
# Use the frame slider or run mplay.
