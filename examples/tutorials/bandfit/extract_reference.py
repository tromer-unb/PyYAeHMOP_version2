"""Extract the six frontier DFT bands from the original electrode EIG/KP."""
import argparse
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument("electrode",type=Path);a=p.parse_args()
ref=a.electrode; t=(ref/"elec.EIG").read_text().split()
ef=float(t[0]); nb,ns,nk=map(int,t[1:4])
assert ns==1
raw=np.array(t[4:],float).reshape(nk,nb+1)[:,1:]-ef
kr=np.loadtxt(ref/"elec.KP",skiprows=1)
k=kr[:,[1,3]]*np.array([2.459512,8.52])/0.529177210903/(2*np.pi)
k[:,1]/=2
out=Path(__file__).parent/"reference"/"frontier_bands.csv"
np.savetxt(out,np.c_[k,raw[:,13:19]],delimiter=",",fmt="%.17g",header="kx_reduced,kz_reduced,E_v3_eV,E_v2_eV,E_v1_eV,E_c1_eV,E_c2_eV,E_c3_eV")
print("SIESTA E_F eV",ef,"kpoints",nk,"bands",nb)
