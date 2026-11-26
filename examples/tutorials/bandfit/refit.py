"""Reproduce the partial BandFit calibration (fixed radial parameters)."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import least_squares
from dataclasses import replace
from pyaehmop.basis import make_basis
from pyaehmop.engine import coupling,_realspace_blocks
from pyaehmop.structure_io import read_xyz
p=argparse.ArgumentParser()
p.add_argument("--bandfit",required=True,type=Path,help="Path to yaehmop-bandfit checkout")
p.add_argument("--workers",type=int,default=4)
a=p.parse_args();sys.path.insert(0,str(a.bandfit.resolve()/"parameter_fitter"))
from global_optimizers import run_global
root=Path(__file__).resolve().parent
arr=np.loadtxt(root/"reference/frontier_bands.csv",delimiter=",")
k=arr[:,:2];raw=arr[:,2:];nk=len(k);k2=k.copy();k2[:,1]+=.5
train=np.arange(nk)%5!=0
out=root/"results/refit";out.mkdir(parents=True,exist_ok=True)
summary={}
for model in ("pdz7","shared"):
    stock=json.loads((root/"reference"/(model+"_initial.json")).read_text())
    st=read_xyz(root.parent.parent/"transport/graphene.xyz");st.vectors=np.array([[2.459512,0,0],[0,0,4.26]]);st.bounds=(5,4)
    for atom in st.atoms:
        atom.shells=tuple(replace(sh,energy=stock["Es"] if sh.l==0 else stock["Ep"],zetas=(stock["zeta_s"],) if sh.l==0 else tuple(stock["zeta_p"]),coeffs=(1.,) if sh.l==0 else tuple(stock["coeff"])) for sh in atom.shells)
    basis=make_basis(st.atoms);blocks=_realspace_blocks(st,basis,(24,20,12),11.,a.workers)
    R=np.array(list(blocks));S=np.array(list(blocks.values()))
    sk=np.einsum("kr,rij->kij",np.exp(2j*np.pi*np.vstack([k,k2])@R.T),S)
    sk=(sk+sk.conj().transpose(0,2,1))/2
    assert np.linalg.eigvalsh(sk).min()>0
    ep=stock["Ep"]
    def predict(v):
        es,K,shift=v;eps=np.tile([es,ep,ep,ep],4)
        hk=coupling(eps,K,True)*sk+np.diag(eps*(1-K))
        ev=np.array([eigh(h,s,eigvals_only=True,check_finite=False) for h,s in zip(hk,sk)])
        return np.sort(np.concatenate([ev[:nk],ev[nk:]],axis=1),axis=1)+shift
    for n in (2,4,6):
        sl=slice(16-n//2,16+n//2);target=raw[:,3-n//2:3+n//2]
        w=np.sqrt(1/(1+(target/2)**2));x0=[stock["Es"],stock["K"],-ep]
        bounds=[(-40,-5),(1.01,4.5),(0,35)]
        def res(v): return ((predict(v)[:,sl]-target)*w)[train].ravel()
        def obj(v,eid,stage):return float(np.mean(res(v)**2)),{}
        opt=run_global("differential-evolution",obj,x0,bounds,180,seed=42)
        fit=least_squares(res,opt["best_x"],bounds=np.array(bounds).T,max_nfev=80)
        pred=predict(fit.x)[:,sl];err=pred-target;name=f"{model}_{n}B"
        summary[name]=dict(stock,Es=float(fit.x[0]),K=float(fit.x[1]),alignment=float(fit.x[2]),RMSE_train_eV=float(np.sqrt(np.mean(err[train]**2))),RMSE_test_eV=float(np.sqrt(np.mean(err[~train]**2))))
        np.savetxt(out/(name+".csv"),np.c_[k,train,target,pred],delimiter=",",header="kx,kz,train,"+",".join([f"DFT_{i+1}" for i in range(n)]+[f"EHT_{i+1}" for i in range(n)]))
        print(name,summary[name],flush=True)
(out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
