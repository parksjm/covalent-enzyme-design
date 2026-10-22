# Diagnostic only: this model failed TI1 bond and hydrogen-bond checks.
load af3_starting_with_hydrogens.pdb, af3_start
load pa6_l4_dftb3_diagnostic.pdb, dftb3_diagnostic
hide everything
show cartoon, chain A+E
show sticks, chain I or (chain E and resi 1+2)
show sticks, resi 113+186 and chain A
show spheres, chain J and name O
set sphere_scale, 0.25, chain J
set valence, 0
util.cnc
disable af3_start
distance reactive_bond, dftb3_diagnostic and chain E and resi 1 and name OG1, dftb3_diagnostic and chain I and name C54
distance water_to_N, dftb3_diagnostic and chain J and name O, dftb3_diagnostic and chain E and resi 1 and name N
distance water_to_OG1, dftb3_diagnostic and chain J and name O, dftb3_diagnostic and chain E and resi 1 and name OG1
zoom dftb3_diagnostic and (chain I+J or (chain E and resi 1)), 8
