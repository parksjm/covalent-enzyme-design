load full_pocket.pdb, pocket
load rejected_oxygen_probes.pdb, probes
hide everything
show sticks, pocket within 5 of probes
hide everything, pocket and chain W
show spheres, probes
set sphere_scale, 0.3, probes
color red, probes
util.cnc pocket
label probes and resi 1, "User O: protein clash"
label probes and resi 2, "Best grid O: still clashes"
set valence, 0
zoom probes, 6
