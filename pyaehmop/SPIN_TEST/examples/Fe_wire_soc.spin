# Atomic SOC without Stoner: nonmagnetic inversion+time-reversal state.
ase
Fe_wire.extxyz
spin = no
soc = yes
stoner_fe_d = 0.88
soc_fe_d = 0.06
kmesh = 24
band_path = G:0,0,0; X:0.5,0,0
divisions = 75
images = 3
cutoff = 5.6
quadrature = 20 18 12
workers = 2
output = ../results/Fe_wire_soc
plots = yes
