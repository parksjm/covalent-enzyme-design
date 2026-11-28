load /data/jpc/covalent-enzyme-design/design/nyl12_relay_ts_ensemble_20261126/templates/medoid_01_delivery_full_010_frame_068/complex_heavy.pdb, model
bg_color white
set ray_opaque_background, 1
set max_threads, 1
set valence, 0
set stick_radius, 0.2
hide everything
show cartoon, chain A+B+C+D+E+F+G+H
set_color protected_magenta, [0.83137,0,0.70980]
set_color mechanistic_teal, [0,0.61961,0.45098]
set_color variable_gray, [0.90196,0.90196,0.90196]
set ambient, 0.65
set direct, 0.35
set reflect, 0.15
set specular, 0
set ray_shadows, off
set depth_cue, off
color variable_gray, all
color protected_magenta, (chain A+B+C+D and resi 1+41+46+47+48+49+50+51+52+53+54+55+56+57+58+59+61+62+64+66+77+78+80+81+88+91+96+97+98+99+100+101+102+104+113+115+121+123+126+156+157+186+195+208+212+215+218+219+221+222+223+224+225+226+227+228+229+230+231+232) or (chain E+F+G+H and resi 233+234+235+246+247+250+251+254+255+257+258+261+262+265+268+269+270+271+272+273+274+275+276+282+287+292+293+294+295+296+298+299+300+301+302+304+305+306+308+309+312)
color mechanistic_teal, (chain A+B+C+D and resi 113+156+186) or (chain E+F+G+H and resi 233+272+274)
show sticks, chain I+J+K+L
color orange, chain I+J+K+L and elem C
color blue, chain I+J+K+L and elem N
color red, chain I+J+K+L and elem O
set orthoscopic, on
set antialias, 2
orient all
zoom all, 3
set_view (0.982273896268952,-0.1371227590979435,-0.1278097869725924,0.18411869489973598,0.8337668577543137,0.5205123755482823,0.035189471394581145,-0.5348178903588262,0.844234283392415,0.0,0.0,-286.6705627441406,64.26124572753906,34.18516540527344,56.89478302001953,226.0132598876953,347.327880859375,20.0)
