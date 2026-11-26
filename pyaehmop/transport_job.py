"""Friendly opt-in transport jobs; periodic pristine coherent devices only."""
import argparse, json
from pathlib import Path
from itertools import product
import numpy as np
from .input import read_input
from .structure_io import read_xyz
from .parameter_files import apply_parameter_selection
from .transport import realspace,spectrum
from .transport_observables import conductance,current_voltage

def run_job(path,workers=None):
    path=Path(path).resolve(); base=path.parent
    opts={}
    for line in path.read_text().splitlines():
        line=line.split("#",1)[0].strip()
        if not line: continue
        key,sep,value=line.partition("=")
        if not sep or key.strip() in opts: raise ValueError("Expected unique key=value: "+line)
        opts[key.strip()]=value.strip()
    allowed=set("structure_format structure parameters lattice_vectors images transport_axis transverse_kpoints energy_min energy_max energy_points energy_reference_eV eta cutoff workers calculations temperature bias_min bias_max bias_points output plots K weighted_hij".split())
    if set(opts)-allowed: raise ValueError("Unknown options: "+str(set(opts)-allowed))
    def get(key,default=None): return opts.get(key,default)
    def ints(key,default): return tuple(int(x) for x in get(key,default).replace(","," ").split())
    kind=get("structure_format","xyz"); src=base/opts["structure"]
    cutoff=float(get("cutoff",11)); n=workers if workers is not None else int(get("workers",1))
    if kind=="xyz":
        st=read_xyz(src)
        st.vectors=np.array([[float(x) for x in row.split()] for row in opts["lattice_vectors"].split(";")])
        st.bounds=ints("images","5 4")
    elif kind=="yaehmop":
        st=read_input(src)
    else:
        from ase.io import read
        from .structure_io import from_ase_atoms
        atoms=read(src,format="cif" if kind=="cif" else None)
        if "lattice_vectors" in opts: raise ValueError("For ASE formats use the file cell and periodicity")
        st=from_ase_atoms(atoms,src.stem,cutoff,ints("images","5 4") if "images" in opts else "auto")
    if st.vectors is None or not 1<=len(st.vectors)<=3 or np.linalg.matrix_rank(st.vectors)!=len(st.vectors) or len(st.bounds)!=len(st.vectors) or min(st.bounds)<1:
        raise ValueError("Valid periodic lattice vectors and image bounds required")
    apply_parameter_selection(get("parameters","default"),base,st,kind)
    st.k_constant=float(get("K",st.k_constant)); st.weighted=get("weighted_hij","true").lower()=="true"
    dim=len(st.vectors); axis=int(get("transport_axis",dim-1))
    if not 0<=axis<dim: raise ValueError("transport_axis is a zero-based lattice index")
    mesh=ints("transverse_kpoints","" if dim==1 else " ".join(["61"]*(dim-1)))
    if len(mesh)!=dim-1 or any(x<1 for x in mesh): raise ValueError("One mesh size per transverse lattice direction")
    ks=np.array(list(product(*[((np.arange(m)+.5)/m-.5) for m in mesh]))).reshape(-1,dim-1) if mesh else np.empty((1,0))
    ne=int(get("energy_points",401)); emin=float(get("energy_min",-2)); emax=float(get("energy_max",2)); eta=float(get("eta",1e-5))
    if ne<2 or emin>=emax or n<1 or cutoff<=0 or eta<=0: raise ValueError("Invalid calculation grid")
    e=np.linspace(emin,emax,ne); ef=float(opts["energy_reference_eV"])
    tasks=set(get("calculations","transmission,conductance,iv").replace(","," ").split())
    if not tasks or tasks-{"transmission","conductance","iv"}: raise ValueError("Supported transport tasks: transmission, conductance, iv")
    R,H,S=realspace(st,cutoff=cutoff,workers=n)
    T,Tk=spectrum(R,H,S,e+ef,ks,axis=axis,eta=eta,workers=n)
    out=base/get("output","results"); out.mkdir(parents=True,exist_ok=True)
    data={"transmission":(e,T,"energy_relative_eV T_per_spin")}
    if "conductance" in tasks: data["conductance"]=(e,conductance(e,T),"energy_relative_eV G_S_zero_temperature")
    if "iv" in tasks:
        nb=int(get("bias_points",81)); vmin=float(get("bias_min",-2)); vmax=float(get("bias_max",2))
        if nb<2 or vmin>=vmax: raise ValueError("Invalid bias grid")
        v=np.linspace(vmin,vmax,nb); I=current_voltage(e,T,v,float(get("temperature",300)))
        data["current_voltage"]=(v,I,"bias_V current_A_frozen_transmission")
    for name,(x,y,header) in data.items(): np.savetxt(out/(name+".dat"),np.c_[x,y],header=header)
    np.savez_compressed(out/"transmission_k.npz",k=ks,energy=e,T=Tk)
    (out/"summary.json").write_text(json.dumps(dict(options=opts,workers=n,energy_reference_eV=ef,spin_degeneracy=2,iv_method="frozen equilibrium transmission"),indent=2))
    if get("plots","false").lower()=="true":
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for name,(x,y,header) in data.items():
            fig,ax=plt.subplots(); ax.plot(x,y); ax.set(xlabel=header.split()[0],ylabel=header.split()[1]); fig.tight_layout(); fig.savefig(out/(name+".png"),dpi=200); plt.close(fig)
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("job"); p.add_argument("--workers",type=int)
    a=p.parse_args(); print(run_job(a.job,a.workers))
if __name__=="__main__": main()
