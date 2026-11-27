"""Generate paired inputs and compare native CIF device jobs to XYZ."""
import json,sys,subprocess
from pathlib import Path
import numpy as np
from examples.transport.prepare_cif_examples import main as prepare,options,xyz,ROOT,TR
from examples.transport.validate_cif import compare
from pyaehmop.device_job import read_device_cif
from pyaehmop.parameter_files import load_parameter_selection
def run(p):
    subprocess.run([sys.executable,"-m","pyaehmop",str(p),"--workers","8"],check=True)
    return p.parent/options(p)["output"]
def main():
    prepare()
    out=TR/"cif_reference";out.mkdir(exist_ok=True)
    summary={}
    jobs=[p for p in sorted((ROOT/"examples/devices").glob("*.device.job")) if not p.name.endswith("_cif.device.job")]+sorted((ROOT/"examples/tutorials").glob("*/run.device.job"))
    for p in jobs:
        q=p.with_name(p.name.replace(".device.job","_cif.device.job"))
        o,c=options(p),options(q)
        metadata=json.loads((q.parent/c["cif_origin"]).read_text())
        params=load_parameter_selection(c.get("parameters","default"),q.parent)
        err=0.
        for k in ("electrode_left","scattering_region","electrode_right"):
            s,x=xyz(p.parent/o[k])
            st=read_device_cif(q.parent/c[k],params,metadata)
            y=np.array([a.position for a in st.atoms])
            assert [v.upper() for v in s]==[a.symbol.upper() for a in st.atoms],(p,k,s,[a.symbol for a in st.atoms])
            err=max(err,float(np.max(abs(x-y))))
            assert np.allclose(x,y,atol=1e-10,rtol=0)
        a=p.parent/o["output"]
        if not (a/"transmission.dat").exists():a=run(p)
        b=run(q)
        summary[str(p.relative_to(ROOT))]={"max_coordinate_error_A":err,**compare(a,b)}
        if (a/"current_spin.dat").exists():
            x,y=np.loadtxt(a/"current_spin.dat"),np.loadtxt(b/"current_spin.dat")
            assert np.allclose(x,y,atol=2e-10,rtol=2e-6)
            summary[str(p.relative_to(ROOT))]["current_spin.dat"]=float(np.max(abs(x-y)))
        (out/"direct_device_validation.json").write_text(json.dumps(summary,indent=2)+"\n")
        print(p.name,summary[str(p.relative_to(ROOT))],flush=True)
    print("PASS",len(summary))
if __name__=="__main__":main()
