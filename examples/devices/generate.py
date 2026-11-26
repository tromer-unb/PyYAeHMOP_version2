"""Generate illustrative geometries, not relaxed or interface-fitted structures."""
from pathlib import Path
import json
root=Path(__file__).parent
def xyz(name,atoms):
    (root/(name+".xyz")).write_text(str(len(atoms))+"\n"+name+"\n"+"\n".join("%s %.9f %.9f %.9f"%a for a in atoms)+"\n")
def job(name,extra="",left="left",right="right",device=None):
    (root/(name+".device.job")).write_text(f"""electrode_left = {left}.xyz
electrode_right = {right}.xyz
scattering_region = {device or name}.xyz
transverse_vectors = 2.459512 0 0
left_repeat = 0 0 8.52
right_repeat = 0 0 8.52
parameters = ../transport/carbon_pdz7_4B.dat
K = 2.0068527138311567
cutoff = 4.5
quadrature = 16 12 8
energy_reference_eV = -18.800896640571224
energy_min = -2
energy_max = 2
energy_points = 81
transverse_kpoints = 15
workers = 4
calculations = transmission, conductance, iv
temperature = 300
bias_min = -2
bias_max = 2
bias_points = 21
output = results/{name}
plots = true
"""+extra)
def graphene(cells):
    return [("C",x,0,z+4.26*i) for i in cells for x,z in [(0,0),(0,1.42),(1.229756,2.13),(1.229756,3.55)]]
L=graphene(range(-2,0)); D=graphene(range(0,2)); R=graphene(range(2,4))
xyz("left",L);xyz("right",R);xyz("pristine",D);job("pristine")
xyz("vacancy",D[:4]+D[5:]);job("vacancy")
xyz("adsorbate",D+[("C",0,3.3,3),("O",0,3.3,4.13)]);job("adsorbate")
mag={"device":{"exchange_eV":{"2":{"p":0.2}},"soc_eV":{"C":{"p":0.02}}}}
(root/"magnetic.json").write_text(json.dumps(mag,indent=2))
for name,opts in [("spin","spin = true\n"),("soc","soc = true\n"),("spin_soc","spin = true\nsoc = true\n")]:
    job(name,opts+"magnetic_parameters = magnetic.json\n",device="vacancy")
# Asymmetric example: different longitudinal repeats and orbital counts in a
# one-dimensional C--Mo--S chain; NOT a physical graphene/MoS2 heterojunction.
xyz("chain_left",[("C",0,0,-1.5)])
xyz("chain_device",[("C",0,0,0),("Mo",0,0,1.8)])
xyz("chain_right",[("S",0,0,3.6)])
(root/"asymmetric.device.job").write_text("""electrode_left = chain_left.xyz
scattering_region = chain_device.xyz
electrode_right = chain_right.xyz
left_repeat = 0 0 1.5
right_repeat = 0 0 1.8
cutoff = 2.0
transverse_kpoints =
energy_reference_eV = -11
energy_min = -2
energy_max = 2
energy_points = 41
workers = 2
calculations = transmission, conductance, iv
output = results/asymmetric
plots = true
""")
# Commensurate illustrative lateral 4x graphene / 3x MoS2 interface.
# MoS2 is strained to the common transverse period; no relaxation or band alignment.
import math
width=4*2.459512
a=width/3
b=math.sqrt(3)*a/2
def gwide(cells):
    return [(s,x+j*2.459512,y,z) for s,x,y,z in graphene(cells) for j in range(4)]
def mos(cells,origin):
    atoms=[]
    for i in cells:
        for j in range(3):
            for x,z in ((0,0),(a/2,b)):
                atoms.append(("Mo",x+j*a,0,z+2*b*i+origin))
                for y in (-1.56,1.56):
                    atoms.append(("S",(x+a/2)%a+j*a,y,(z+b/3)%(2*b)+2*b*i+origin))
    return atoms
# 2 graphene and 1 MoS2 rectangular cell in the device.
origin=8.52+0.6
xyz("hetero_left",gwide(range(-2,0)))
xyz("hetero_device",gwide(range(0,2))+mos(range(0,1),origin))
xyz("hetero_right",mos(range(1,2),origin))
(root/"graphene_mos2.device.job").write_text(f"""electrode_left = hetero_left.xyz
scattering_region = hetero_device.xyz
electrode_right = hetero_right.xyz
transverse_vectors = {width} 0 0
left_repeat = 0 0 8.52
right_repeat = 0 0 {2*b}
parameters = ../transport/carbon_pdz7_4B.dat
K = 2.0068527138311567
cutoff = 3.5
quadrature = 16 12 8
energy_reference_eV = -18.800896640571224
energy_min = -2
energy_max = 2
energy_points = 21
transverse_kpoints = 3
workers = 8
calculations = transmission, conductance, iv
output = results/graphene_mos2
plots = true
""")
