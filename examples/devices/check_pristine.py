"""Compare the two NEGF paths at the SAME reduced numerical resolution."""
from pathlib import Path
import numpy as np
from pyaehmop.structure_io import read_xyz
from pyaehmop.parameter_files import load_parameter_selection
from pyaehmop.transport import realspace,spectrum
root=Path(__file__).resolve().parent
st=read_xyz(root.parent/"transport/graphene.xyz",load_parameter_selection("../transport/carbon_pdz7_4B.dat",root))
st.vectors=np.array([[2.459512,0,0],[0,0,4.26]]);st.bounds=(3,3);st.k_constant=2.0068527138311567
R,H,S=realspace(st,4.5,(16,12,8),8)
e=np.linspace(-2,2,81);k=((np.arange(15)+.5)/15-.5)[:,None]
old,_=spectrum(R,H,S,e-18.800896640571224,k,axis=1,workers=8)
new=np.loadtxt(root/"results/pristine/transmission.dat")[:,1]
err=np.abs(old-new)
print("graphene old/new MAE",err.mean(),"MAX",err.max())
assert err.max()<0.002
np.savetxt(root/"results/pristine_comparison.dat",np.c_[e,old,new],header="E_eV old_T new_T")
