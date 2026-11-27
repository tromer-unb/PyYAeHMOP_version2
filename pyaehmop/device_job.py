"""Explicit two-terminal device input (XYZ or ordered P1 CIF, common frame)."""
import argparse, json
from pathlib import Path
from itertools import product
import numpy as np
from .structure_io import read_xyz,from_ase_atoms
from .parameter_files import load_parameter_selection
from .device_transport import Region,build_device,spectrum_device
from .transport_observables import conductance,current_voltage

def read_device_cif(path,params,metadata=None,*,return_cell=False):
    """Read a fully occupied P1 interchange CIF in a common Cartesian box.

    The artificial CIF cell is not a lead repeat. All three CIFs must share
    the same box and coordinate frame. A global translation is irrelevant;
    legacy origin metadata is optional.
    """
    from ase.io import read
    atoms=read(path,format="cif",index=0)
    group=atoms.info.get("spacegroup")
    if getattr(group,"no",None)!=1:
        raise ValueError("Device CIF requires explicit P1 sites; expand symmetry before building the device")
    for site,comps in atoms.info.get("occupancy",{}).items():
        if len(comps)!=1 or abs(float(next(iter(comps.values())))-1)>1e-6:
            raise ValueError(f"Device CIF requires fully occupied ordered sites: {site}")
    cell=np.asarray(atoms.cell,float)
    lengths=np.diag(cell)
    if not np.isfinite(cell).all() or np.any(lengths<=0) or not np.allclose(cell,np.diag(lengths),atol=1e-8,rtol=0):
        raise ValueError("Device CIF requires a positive orthogonal common box")
    origin=np.zeros(3)
    if metadata is not None:
        origin=np.asarray(metadata["cartesian_origin_A"],float)
        legacy_lengths=np.asarray(metadata["box_lengths_A"],float)
        if origin.shape!=(3,) or legacy_lengths.shape!=(3,) or not np.isfinite(origin).all() or not np.isfinite(legacy_lengths).all() or np.any(legacy_lengths<=0):
            raise ValueError("Invalid legacy CIF origin or box lengths")
        if not np.allclose(cell,np.diag(legacy_lengths),atol=1e-8,rtol=0):
            raise ValueError("CIF region cell does not match legacy origin box")
    atoms.positions=np.round(atoms.positions+origin,12)
    atoms.pbc=False
    structure=from_ase_atoms(atoms,Path(path).stem,parameter_table=params)
    return (structure,cell) if return_cell else structure

def run_job(path,workers=None):
    path=Path(path).resolve(); base=path.parent; o={}
    for line in path.read_text().splitlines():
        line=line.split("#",1)[0].strip()
        if not line: continue
        k,sep,v=line.partition("="); k=k.strip()
        if not sep or k in o: raise ValueError("Expected unique key=value: "+line)
        o[k]=v.strip()
    allowed=set("structure_format cif_origin electrode_left electrode_right scattering_region transverse_vectors left_repeat right_repeat parameters K weighted_hij cutoff quadrature workers spin soc magnetic_parameters energy_reference_eV energy_min energy_max energy_points eta transverse_kpoints calculations temperature bias_min bias_max bias_points output plots".split())
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
    kind=o.get("structure_format","xyz").lower()
    if kind not in ("xyz","cif"): raise ValueError("Device structure_format must be xyz or cif")
    if kind=="xyz" and "cif_origin" in o: raise ValueError("cif_origin is only valid for device CIF")
    metadata=json.loads((base/o["cif_origin"]).read_text()) if "cif_origin" in o else None
    regions=[]
    common_cell=None
    for label,key in (("left","electrode_left"),("device","scattering_region"),("right","electrode_right")):
        source=base/o[key]
        if kind=="cif":
            st,cell=read_device_cif(source,params,metadata,return_cell=True)
            if common_cell is None: common_cell=cell
            elif not np.allclose(cell,common_cell,atol=1e-8,rtol=0):
                raise ValueError("All three device CIFs must share the same cell and Cartesian frame")
        else:
            st=read_xyz(source,params)
        cfg=magnet.get(label,{})
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
