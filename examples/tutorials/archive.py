"""Archive compact, reproducible results after PBS suites finish."""
from pathlib import Path
import json,re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
repo=root.parent.parent
out=root/"reference";out.mkdir(exist_ok=True)
for status in (root/"results").glob("*_status.txt"):
    if status.read_text(): (out/status.name).write_text(status.read_text())
def table(src,dst):
    a=np.loadtxt(src); np.savetxt(dst,a,fmt="%.10g",header=src.read_text().splitlines()[0].lstrip("# "))
    return a
summary={}
for f in sorted(root.glob("*/run.device.job")):
    r=f.parent/"results"; dest=f.parent/"reference";dest.mkdir(exist_ok=True)
    if not (r/"summary.json").exists():continue
    t=table(r/"transmission.dat",dest/"transmission.dat")
    g=table(r/"conductance.dat",dest/"conductance.dat")
    i=table(r/"current_voltage.dat",dest/"current_voltage.dat")
    if (r/"current_spin.dat").exists(): table(r/"current_spin.dat",dest/"current_spin.dat")
    md=json.loads((r/"summary.json").read_text())
    summary[f.parent.name]={"T_at_reference":float(t[len(t)//2,1]),"T_max":float(t[:,1].max()),"I_at_2V_uA":float(i[-1,1]*1e6),"G_at_reference_uS":float(g[len(g)//2,1]*1e6),"spin":md["spin"],"soc":md["soc"],"atom_counts":md["atom_counts"]}
    fig,axs=plt.subplots(1,3,figsize=(10,3))
    axs[0].plot(t[:,0],t[:,1],label="total" if md["spin"] or md["soc"] else "per spin")
    if t.shape[1]==4:
        axs[0].plot(t[:,0],t[:,2],label="up");axs[0].plot(t[:,0],t[:,3],label="down");axs[0].legend(fontsize=8)
    axs[0].set(xlabel="E - Eref (eV)",ylabel="Transmission")
    axs[1].plot(g[:,0],g[:,1]*1e6);axs[1].set(xlabel="E - Eref (eV)",ylabel="Conductance (µS)")
    axs[2].plot(i[:,0],i[:,1]*1e6);axs[2].set(xlabel="Bias (V)",ylabel="Current (µA)")
    fig.suptitle(f.parent.name.replace("_"," ")+" — illustrative reduced-grid model");fig.tight_layout()
    fig.savefig(dest/"spectra.svg");plt.close(fig)
(out/"devices_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
# Curated scalar reference from every executed molecular/periodic and spin example.
electronic={}
for f in sorted((repo/"pyaehmop/results").glob("*/summary.json")):
    d=json.loads(f.read_text())
    keys=("type","title","atoms","orbitals","electrons","electronic_energy_eV","homo_eV","lumo_eV","gap_homo_lumo_eV","band_path_gap_eV","minimum_direct_gap_eV","vbm_eV","cbm_eV","mulliken_charges","warnings")
    electronic[f.parent.name]={k:d[k] for k in keys if k in d}
for f in sorted((repo/"pyaehmop/SPIN_TEST/results").glob("*/spin_soc_result.json")):
    d=json.loads(f.read_text())
    keys=("mode","converged","iterations","final_residual","fermi_eV","total_moment_muB","soc_max_element_eV","spin","soc","n_orbitals")
    electronic[f.parent.name]={k:d[k] for k in keys if k in d}
(out/"electronic_summary.json").write_text(json.dumps(electronic,indent=2)+"\n")
# All six transport spectra per resolution, as compact shared-grid tables.
b=root/"bandfit";dest=b/"reference"; names=sorted(json.loads((root.parent/"transport"/"multiband_parameters.json").read_text()))
for mode in ("smoke","full"):
    cols={key:[] for key in ("transmission","conductance","current_voltage")}
    for n in names:
        d=b/"results"/(n+("_smoke" if mode=="smoke" else ""))
        if not (d/"summary.json").exists():break
        for key in cols:cols[key].append(np.loadtxt(d/(key+".dat")))
    else:
        for key,arr in cols.items():
            assert all(np.allclose(a[:,0],arr[0][:,0]) for a in arr)
            np.savetxt(dest/(mode+"_"+key+".csv"),np.column_stack([arr[0][:,0]]+[a[:,1] for a in arr]),delimiter=",",fmt="%.10g",header=("bias_V" if key=="current_voltage" else "relative_energy_eV")+","+",".join(names))
        fig,ax=plt.subplots(figsize=(7,4))
        for n,a in zip(names,cols["transmission"]):ax.plot(a[:,0],a[:,1],label=n)
        ax.set(xlabel="E - Eref (eV)",ylabel="T per spin",title="Six frozen BandFit parameter sets — "+mode);ax.legend(ncol=2,fontsize=8);fig.tight_layout()
        fig.savefig(dest/(mode+"_transmission.svg"));plt.close(fig)
        if mode=="full":
            old=np.genfromtxt(root.parent/"transport"/"multiband.csv",delimiter=",",names=True)
            # Record reproducibility against archived curves by exact column name, when available.
            print("Archived columns",old.dtype.names)
            metric={}
            for n,a in zip(names,cols["transmission"]):
                if n in old.dtype.names:
                    e=old[old.dtype.names[0]];err=np.abs(np.interp(e,a[:,0],a[:,1])-old[n])
                    metric[n]={"MAE":float(err.mean()),"max_abs":float(err.max())}
            (dest/"reproduction_metrics.json").write_text(json.dumps(metric,indent=2)+"\n")
fit=b/"results/refit/summary.json"
if fit.exists():
    d=json.loads(fit.read_text());(dest/"refit_summary.json").write_text(json.dumps(d,indent=2)+"\n")
    archived=json.loads((root.parent/"transport"/"multiband_parameters.json").read_text())
    diff={n:{k:float(d[n][k]-archived[n][k]) for k in ("Es","K","alignment")} for n in names}
    (dest/"refit_parameter_differences.json").write_text(json.dumps(diff,indent=2)+"\n")
print(json.dumps(summary,indent=2))
