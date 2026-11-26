"""Build additional input-only tutorials from the shipped device geometries."""
from pathlib import Path
import json
root=Path(__file__).resolve().parent
src=root.parent/"devices"
def readxyz(name):
    return [(v[0],*map(float,v[1:4])) for line in (src/(name+".xyz")).read_text().splitlines()[2:] for v in [line.split()]]
def xyz(folder,name,atoms):
    (folder/(name+".xyz")).write_text(str(len(atoms))+"\n"+name+" (illustrative geometry)\n"+"\n".join("%s %.9f %.9f %.9f"%a for a in atoms)+"\n")
def make(name,L,D,R,v,l,r,extra=None):
    folder=root/name;folder.mkdir(exist_ok=True)
    for k,a in zip(("left","device","right"),(L,D,R)):xyz(folder,k,a)
    opts=dict(electrode_left="left.xyz",scattering_region="device.xyz",electrode_right="right.xyz",transverse_vectors=v,left_repeat=l,right_repeat=r,parameters="../../transport/carbon_pdz7_4B.dat",K=2.0068527138311567,cutoff=4.5,quadrature="16 12 8",energy_reference_eV=-18.800896640571224,energy_min=-2,energy_max=2,energy_points=81,transverse_kpoints=15,workers=4,calculations="transmission, conductance, iv",temperature=300,bias_min=-2,bias_max=2,bias_points=21,output="results",plots="true")
    opts.update(extra or {})
    (folder/"run.device.job").write_text("\n".join(f"{k} = {v}" for k,v in opts.items())+"\n")
    return folder
L,D,R=[readxyz(n) for n in ("left","pristine","right")]
def widen(a,n): return [(s,x+j*2.459512,y,z) for j in range(n) for s,x,y,z in a]
def rotate(a): return [(s,z,y,-x) for s,x,y,z in a]
# A true change of crystalline transport direction: build layers from x cells.
def along_x(cells):
    return [(s,x+2.459512*i,y,z) for i in cells for s,x,y,z in D[:4]]
make("graphene_x",along_x(range(-3,0)),along_x(range(0,3)),along_x(range(3,6)),"0 0 4.26","7.378536 0 0","7.378536 0 0")
# Finite width strip, not periodic in x. Edges intentionally unpassivated.
make("nanoribbon",widen(L,2),widen(D,2),widen(R,2),"","0 0 8.52","0 0 8.52",{"transverse_kpoints":""})
# A wider transverse supercell and one removed carbon only.
d=widen(D,4); del d[12]
make("wide_vacancy",widen(L,4),d,widen(R,4),"9.838048 0 0","0 0 8.52","0 0 8.52",{"transverse_kpoints":5,"energy_points":41})
# Graphene--MoS2--graphene, illustrative commensurate lateral channel.
hl,hd,hr=[readxyz(n) for n in ("hetero_left","hetero_device","hetero_right")]
channel=[a for a in hd if a[0]!="C"]
device=widen(D,4)+channel
start=max(a[3] for a in channel)+1.42
right=[(s,x,y,z-8.52+start) for s,x,y,z in widen(R,4)]
make("graphene_mos2_graphene",hl,device,right,"9.838048 0 0","0 0 8.52","0 0 8.52",{"cutoff":3.5,"energy_points":21,"transverse_kpoints":3})
# Identical leads with independent frozen electrostatic offset in the device.
f=make("potential_barrier",L,D,R,"2.459512 0 0","0 0 8.52","0 0 8.52",{"magnetic_parameters":"potentials.json"})
(f/"potentials.json").write_text(json.dumps({"device":{"offset_eV":.5}},indent=2))
# Parallel/antiparallel magnetizations in contacts and frozen central exchange.
for name,sign in (("spin_parallel",1),("spin_antiparallel",-1)):
    f=make(name,L,D,R,"2.459512 0 0","0 0 8.52","0 0 8.52",{"spin":"true","magnetic_parameters":"magnetic.json"})
    def exchange(a,scale):return {str(i):{"p":.2*scale} for i in range(len(a))}
    (f/"magnetic.json").write_text(json.dumps({"left":{"exchange_eV":exchange(L,1)},"device":{"exchange_eV":exchange(D,1)},"right":{"exchange_eV":exchange(R,sign)}},indent=2))
# Full-resolution original multiband jobs plus inexpensive runnable smoke copies.
folder=root/"bandfit";folder.mkdir(exist_ok=True)
params=json.loads((root.parent/"transport"/"multiband_parameters.json").read_text())
for name,p in params.items():
    zp=p["zeta_p"];c=p["coeff"]
    (folder/(name+".dat")).write_text(f"C 6 4 1 2 s {p['Es']:.17g} {p['zeta_s']:.17g} 0 1 0\nC 6 4 {len(zp)} 2 p {p['Ep']:.17g} {zp[0]:.17g} {(zp[1] if len(zp)>1 else 0):.17g} {c[0]:.17g} {(c[1] if len(c)>1 else 0):.17g}\n")
    raw=(root.parent/"transport"/"graphene.transport.job").read_text()
    opts=dict(line.split(" = ",1) for line in raw.splitlines() if " = " in line)
    opts.update(structure="../../transport/graphene.xyz",parameters=name+".dat",K=p["K"],energy_reference_eV=p["dirac_eV"],output="results/"+name)
    for smoke in (False,True):
        job=opts.copy()
        if smoke:job.update(transverse_kpoints=15,energy_points=81,output="results/"+name+"_smoke")
        (folder/(name+("_smoke" if smoke else "")+".transport.job")).write_text("\n".join(f"{k} = {v}" for k,v in job.items())+"\n")
