"""Explicit two-terminal device input (XYZ, shared Cartesian coordinates)."""
import argparse, json
from pathlib import Path
from itertools import product
import numpy as np
from .structure_io import read_xyz
from .parameter_files import load_parameter_selection
from .device_transport import Region,build_device,spectrum_device
from .transport_observables import conductance,current_voltage

def run_job(path,workers=None):
    path=Path(path).resolve(); base=path.parent; o={}
    for line in path.read_text().splitlines():
        line=line.split("#",1)[0].strip()
        if not line: continue
        k,sep,v=line.partition("="); k=k.strip()
        if not sep or k in o: raise ValueError("Expected unique key=value: "+line)
        o[k]=v.strip()
    allowed=set("electrode_left electrode_right scattering_region transverse_vectors left_repeat right_repeat parameters K weighted_hij cutoff quadrature workers spin soc magnetic_parameters energy_reference_eV energy_min energy_max energy_points eta transverse_kpoints calculations temperature bias_min bias_max bias_points output plots".split())
    if set(o)-allowed: raise ValueError("Unknown options: "+str(set(o)-allowed))
    def yes(k): return o.get(k,"false").lower() in ("yes","true","1")
    def vec(k): return np.array([float(x) for x in o[k].split()])
    def num(k,d): return float(o.get(k,d))
    def integer(k,d): return int(o.get(k,d))
    v=np.array([[float(x) for x in row.split()] for row in o.get("transverse_vectors","").split(";") if row.strip()]).reshape(-1,3)
    params=load_parameter_selection(o.get("parameters","default"),base)
    spin,soc=yes("spin"),yes("soc")
    magnet=json.loads((base/o["magnetic_parameters"]).read_text()) if "magnetic_parameters" in o else {}
    if (spin or soc) and not magnet: raise ValueError("Explicit calibrated magnetic_parameters JSON required")
    regions=[]
    for label,key in (("left","electrode_left"),("device","scattering_region"),("right","electrode_right")):
        st=read_xyz(base/o[key],params); cfg=magnet.get(label,{})
        if set(cfg)-{"offset_eV","exchange_eV","soc_eV"}: raise ValueError("Unknown magnetic region key")
        regions.append(Region(st.atoms,float(cfg.get("offset_eV",0)),cfg.get("exchange_eV",{}) if spin else {},cfg.get("soc_eV",{}) if soc else {}))
    if spin and not any(np.any(r.shift) for r in regions): raise ValueError("spin enabled but no exchange")
    if soc and not any(r.soc for r in regions): raise ValueError("SOC enabled but no xi")
    mesh=tuple(int(x) for x in o.get("transverse_kpoints"," ".join(["31"]*len(v))).split())
    if len(mesh)!=len(v) or any(m<1 for m in mesh): raise ValueError("One k mesh per transverse vector")
    ks=np.array(list(product(*[((np.arange(m)+.5)/m-.5) for m in mesh]))) if mesh else np.empty((1,0))
    nw=workers if workers is not None else integer("workers",1)
    cutoff=num("cutoff",11); eta=num("eta",1e-5); ne=integer("energy_points",401)
    emin,emax=num("energy_min",-2),num("energy_max",2)
    if nw<1 or ne<2 or not emin<emax or cutoff<=0 or eta<=0: raise ValueError("Invalid numerical grid")
    e=np.linspace(emin,emax,ne); ef=float(o["energy_reference_eV"])
    tasks=set(o.get("calculations","transmission,conductance,iv").replace(","," ").split())
    if not tasks or tasks-{"transmission","conductance","iv"}: raise ValueError("Unknown calculation")
    M=build_device(*regions,v,vec("left_repeat"),vec("right_repeat"),cutoff,num("K",1.75),o.get("weighted_hij","true").lower() in ("true","yes"),nw,tuple(int(x) for x in o.get("quadrature","24 20 12").split()),spin or soc)
    t,tk=spectrum_device(M,e+ef,ks,eta,nw,spin and not soc)
    factor=.5 if spin or soc else 1. # existing Landauer helper assumes per-spin
    out=base/o.get("output","device_results"); out.mkdir(parents=True,exist_ok=True)
    names="T_total_spinor" if spin or soc else "T_per_spin"
    np.savetxt(out/"transmission.dat",np.c_[e,t],header="E_relative_eV "+names+(" T_up T_down" if spin and not soc else ""))
    np.savez_compressed(out/"transmission_k.npz",energy=e,k=ks,T=tk)
    data={"transmission":(e,t[:,0],"E (eV)",names)}
    if "conductance" in tasks: data["conductance"]=(e,conductance(e,t[:,0])*factor,"E (eV)","G (S)")
    if "iv" in tasks:
        nb=integer("bias_points",81); lo,hi=num("bias_min",-2),num("bias_max",2)
        if nb<2 or lo>=hi: raise ValueError("Invalid bias grid")
        bias=np.linspace(lo,hi,nb); data["current_voltage"]=(bias,current_voltage(e,t[:,0],bias,num("temperature",300))*factor,"V (V)","I (A)")
        if spin and not soc:
            currents=np.column_stack([current_voltage(e,t[:,i],bias,num("temperature",300))*.5 for i in (1,2)])
            np.savetxt(out/"current_spin.dat",np.c_[bias,currents],header="bias_V I_up_A I_down_A")
    for name,(x,y,xlabel,ylabel) in data.items():
        if name!="transmission": np.savetxt(out/(name+".dat"),np.c_[x,y],header=xlabel+" "+ylabel)
        if yes("plots"):
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig,ax=plt.subplots(); ax.plot(x,y); ax.set(xlabel=xlabel,ylabel=ylabel); fig.tight_layout(); fig.savefig(out/(name+".png"),dpi=180); plt.close(fig)
    (out/"summary.json").write_text(json.dumps(dict(options=o,spin=spin,soc=soc,spin_degeneracy=1 if spin or soc else 2,magnetism="frozen",iv="equilibrium transmission",atom_counts=[len(r.atoms) for r in regions]),indent=2))
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("job");p.add_argument("--workers",type=int)
    a=p.parse_args();print(run_job(a.job,a.workers))
if __name__=="__main__": main()
