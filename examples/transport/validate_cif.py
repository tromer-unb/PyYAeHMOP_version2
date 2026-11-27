"""Validate the six quick CIF jobs and 15 explicit device CIF conversions."""
from pathlib import Path
import json, argparse, subprocess, sys
import numpy as np
from ase.io import read
from examples.transport.device_from_cif import convert
from examples.transport.prepare_cif_examples import options,xyz,ROOT,TR,BF

def compare(a,b):
    result={}
    for f in ("transmission.dat","conductance.dat","current_voltage.dat"):
        x,y=np.loadtxt(a/f),np.loadtxt(b/f)
        assert x.shape==y.shape and np.allclose(x[:,0],y[:,0],atol=1e-12)
        err=np.max(abs(x[:,1:]-y[:,1:]))
        result[f]=float(err)
        assert np.allclose(x[:,1:],y[:,1:],atol=2e-6 if f=="transmission.dat" else 2e-10,rtol=2e-6), (f,err)
    return result

def run(path,workers=8):
    # Fresh interpreter for every job, avoiding accumulated native library state.
    subprocess.run([sys.executable,"-m","pyaehmop",str(path),"--workers",str(workers)],check=True)
    return path.parent/options(path)["output"]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume",action="store_true",help="Reuse completed records from a previous interrupted validation")
    args=parser.parse_args()
    out=TR/"cif_reference";out.mkdir(exist_ok=True)
    summary=json.loads((out/"validation.json").read_text()) if args.resume and (out/"validation.json").exists() else {"uniform":{},"devices":{}}
    for p in sorted(BF.glob("*B_smoke.transport.job")):
        if args.resume and p.stem in summary["uniform"]: continue
        q=p.with_name(p.name.replace(".transport.job","_cif.transport.job"))
        a=p.parent/options(p)["output"]
        if not (a/"transmission.dat").exists(): a=run(p)
        b=run(q)
        summary["uniform"][p.stem]=compare(a,b)
        print(p.name,summary["uniform"][p.stem],flush=True)
        if "pdz7_4B" in p.name:
            x,y=np.loadtxt(a/"transmission.dat"),np.loadtxt(b/"transmission.dat")
            np.savetxt(out/"graphene_xyz_cif.dat",np.c_[x,y[:,1]],header="E_relative_eV T_XYZ T_CIF")
        (out/"validation.json").write_text(json.dumps(summary,indent=2)+"\n")
    jobs=[*sorted((ROOT/"examples/devices").glob("*.device.job")),*sorted((ROOT/"examples/tutorials").glob("*/run.device.job"))]
    for p in jobs:
        if p.name.endswith("_cif.device.job"): continue
        if args.resume and str(p.relative_to(ROOT)) in summary["devices"]: continue
        q=convert(p);o=options(p);new=options(q)
        err=0.
        shift=None
        for k in ("electrode_left","scattering_region","electrode_right"):
            s,x=xyz(p.parent/o[k]);t,y=xyz(q.parent/new[k])
            assert s==t,(p,k,"atom order")
            if shift is None: shift=y[0]-x[0]
            err=max(err,float(np.max(abs(x+shift-y))))
            assert np.allclose(x+shift,y,atol=1e-8,rtol=0),(p,k,err)
        a=p.parent/o["output"]
        if not (a/"transmission.dat").exists(): a=run(p)
        b=run(q)
        summary["devices"][str(p.relative_to(ROOT))]={"max_coordinate_error_A":err,**compare(a,b)}
        print(p.name,summary["devices"][str(p.relative_to(ROOT))],flush=True)
        (out/"validation.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("PASS",len(summary["uniform"]),len(summary["devices"]),flush=True)
if __name__=="__main__":main()
