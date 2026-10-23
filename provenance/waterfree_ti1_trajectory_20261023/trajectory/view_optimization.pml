delete waterfree_ti1_opt
load optimization_trajectory.pdb, waterfree_ti1_opt, discrete=0
hide everything, waterfree_ti1_opt
show sticks, waterfree_ti1_opt
set valence, 0
util.cnc waterfree_ti1_opt
label (waterfree_ti1_opt and resn THR and name OG1), "OG1"
label (waterfree_ti1_opt and resn TIL and name C54), "C54"
label (waterfree_ti1_opt and resn TIL and name O49), "O49"
label (waterfree_ti1_opt and resn TIL and name N59), "N59"
set movie_fps, 3
mset 1 -14
zoom waterfree_ti1_opt
