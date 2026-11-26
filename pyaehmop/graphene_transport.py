"""Reproduce the frozen BandFit graphene benchmark locally."""
import argparse
from dataclasses import replace
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from .basis import Atom, read_parameter_table
from .input import Structure
from .transport import realspace, spectrum

MODELS={
 "pdz7":(-22.55867857634238,-19.536735435070756,2.0016866418080506,2.8360120667296433,(1.6792937388858635,2.7402308379800138),(1.,.781232649690947)),
 "sp5":(-13.190198983457268,-9.431812309188313,3.2240164367812105,2.634063113636377,(2.0384886384108434,),(1.,))
}
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model",choices=MODELS,default="pdz7")
    p.add_argument("--nk",type=int,default=61)
    p.add_argument("--ne",type=int,default=101)
    p.add_argument("--workers",type=int,default=1)
    p.add_argument("--eta",type=float,default=1e-5)
    p.add_argument("--output",type=Path,default=Path("graphene_transport"))
    p.add_argument("--reference",type=Path,help="Optional two-column DFT curve, already aligned")
    a=p.parse_args()
    if a.nk<1 or a.ne<2 or a.workers<1: p.error("Invalid grid or workers")
    es,ep,K,zs,zp,cp=MODELS[a.model]
    shells=tuple(replace(sh,energy=es if sh.l==0 else ep,zetas=(zs,) if sh.l==0 else zp,coeffs=(1.,) if sh.l==0 else cp) for sh in read_parameter_table()["C"]["shells"])
    pos=[[0,0,0],[0,0,1.42],[1.229756,0,2.13],[1.229756,0,3.55]]
    st=Structure("graphene",[Atom("C",np.array(r),4,shells) for r in pos],16,np.array([[2.459512,0,0],[0,0,4.26]]),(5,4),[],0,[])
    st.k_constant=K
    R,H,S=realspace(st,workers=a.workers)
    phase=np.exp(2j*np.pi*(R@np.array([1/3,0])))
    ev=eigh(np.einsum("r,rij->ij",phase,H),np.einsum("r,rij->ij",phase,S),eigvals_only=True)
    ef=(ev[7]+ev[8])/2
    x=np.linspace(-2,2,a.ne); ks=((np.arange(a.nk)+.5)/a.nk-.5)[:,None]
    y,tk=spectrum(R,H,S,x+ef,ks,eta=a.eta,workers=a.workers)
    a.output.mkdir(parents=True,exist_ok=True)
    np.savetxt(a.output/"transmission.dat",np.c_[x,y],header="E_minus_Dirac_eV T_per_spin")
    np.savez_compressed(a.output/"transmission_k.npz",k=ks,energy=x,T=tk,dirac_eV=ef)
    print("Dirac reference (eV):",ef)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig,ax=plt.subplots()
    if a.reference:
        e,t=np.loadtxt(a.reference,unpack=True)
        ax.plot(e,t,label="DFT reference")
    ax.plot(x,y,label=a.model)
    ax.set(xlabel="Energy relative to Dirac / reference (eV)",ylabel="Transmission per spin")
    ax.legend(); fig.tight_layout(); fig.savefig(a.output/"transmission.png",dpi=200)
    plt.close(fig)
if __name__=="__main__": main()
