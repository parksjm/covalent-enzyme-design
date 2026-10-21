load pa6_l4_ti1_starting.pdb, pa6_l4
load maxa_maxa_ti1_starting.pdb, maxa_maxa
hide everything
set valence, 0
show cartoon, chain A+B+C+D+E+F+G+H
color gray80, chain A+B+C+D+F+G+H
color cyan, chain E
show sticks, chain I+J or (chain E and resi 1+2) or (chain A and resi 186)
color yellow, chain I
color red, chain J
set cartoon_transparency, 0.5
disable maxa_maxa
zoom pa6_l4 and (chain I+J or (chain E and resi 1)), 8
show sticks, resi 113+186 and chain A
util.cnc
