# Teaching-only SOC-mixing demonstration; not physical Fe parameters.
ase
Fe_wire.extxyz
spin = yes
soc = yes
stoner_fe_d = 0.88
soc_fe_d = 0.3
initial_moment = 2.5
mixing = 0.3
temperature_ev = 0.025
tolerance = 0.00001
max_iter = 140
kmesh = 24
band_path = G:0,0,0; X:0.5,0,0
divisions = 75
images = 3
cutoff = 5.6
quadrature = 20 18 12
workers = 2
output = ../results/Fe_wire_soc_exaggerated
plots = yes
