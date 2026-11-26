# Monatomic Fe wire: a=4.16 bohr = 2.2013 angstroms.
# Geometry based on Autes et al., J. Phys.: Condens. Matter 18, 6785 (2006).
ase
Fe_wire.extxyz
spin = yes
soc = no
stoner_fe_d = 0.88
soc_fe_d = 0.06
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
output = ../results/Fe_wire_spin
plots = yes
