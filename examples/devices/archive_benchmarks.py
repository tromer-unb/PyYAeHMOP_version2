"""Archive compact postprocessed spectra from completed PBS device examples."""
from pathlib import Path
import shutil
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
out=root/"benchmarks";out.mkdir(exist_ok=True)
names=["pristine","vacancy","adsorbate","asymmetric","spin","soc","spin_soc","graphene_mos2"]
data={}
summary={}
for name in names:
    d=root/"results"/name
    t=np.loadtxt(d/"transmission.dat"); data[name]=t
    i=np.loadtxt(d/"current_voltage.dat")
    summary[name]={"T_max":float(t[:,1].max()),"T_at_reference":float(t[len(t)//2,1]),"I_at_max_bias_uA":float(i[-1,1]*1e6),"n_energy":len(t),"spin_convention":"total spinor" if name in ("spin","soc","spin_soc") else "per spin"}
    for f in ("transmission.dat","conductance.dat","current_voltage.dat","current_spin.dat"):
        if (d/f).exists(): shutil.copyfile(d/f,out/(name+"_"+f))
a=np.loadtxt(root/"results/pristine_comparison.dat")
shutil.copyfile(root/"results/pristine_comparison.dat",out/"pristine_comparison.dat")
summary["pristine_regression"]={"MAE":float(abs(a[:,1]-a[:,2]).mean()),"max_abs":float(abs(a[:,1]-a[:,2]).max())}
(out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
def plot(filename,curves,ylabel,title):
    fig,ax=plt.subplots(figsize=(6.6,4))
    for x,y,label in curves: ax.plot(x,y,label=label,lw=1.4)
    ax.set(xlabel="E - Eref (eV)",ylabel=ylabel,title=title)
    ax.legend();fig.tight_layout();fig.savefig(out/(filename+".svg"));plt.close(fig)
plot("pristine_regression",[(a[:,0],a[:,1],"Original NEGF"),(a[:,0],a[:,2],"Device NEGF")],"T per spin","Graphene: matched reduced grid")
plot("graphene_defects",[(data[n][:,0],data[n][:,1],n) for n in ("pristine","vacancy","adsorbate")],"T per spin","Illustrative periodic vacancy and CO")
plot("magnetic",[(data["spin"][:,0],data["spin"][:,j],label) for j,label in [(1,"spin total"),(2,"up"),(3,"down")]]+[(data[n][:,0],data[n][:,1],n) for n in ("soc","spin_soc")],"T (explicit spin)","Frozen spin and SOC demonstration")
plot("graphene_mos2",[(data["graphene_mos2"][:,0],data["graphene_mos2"][:,1],"unfitted interface")],"T per spin","Graphene/MoS2: numerical smoke test only")
print(json.dumps(summary))
