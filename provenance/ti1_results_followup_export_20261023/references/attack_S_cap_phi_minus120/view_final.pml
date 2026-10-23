delete ti1_final_reference
load final_structure.sdf, ti1_final_reference
hide everything, ti1_final_reference
show sticks, ti1_final_reference
set valence, 0
util.cnc
set sphere_scale, 0.25
show spheres, ti1_final_reference and elem O+N
label ti1_final_reference and index 4, "OG1"
label ti1_final_reference and index 23, "C*"
label ti1_final_reference and index 24, "O-"
label ti1_final_reference and index 25, "leaving N"
orient ti1_final_reference
zoom ti1_final_reference, 3
