"""Convert the shipped three-region CIF sources into a runnable XYZ device job."""
import argparse,os
from pathlib import Path
from ase.io import read

def convert(job):
    job=Path(job).resolve()
    opts=dict(tuple(x.strip() for x in l.split("#",1)[0].split("=",1)) for l in job.read_text().splitlines() if "=" in l.split("#",1)[0])
    name=job.name.removesuffix(".device.job")
    src=job.parent/"cif_sources"/name
    dest=job.parent/"results"/"cif_roundtrip"/name
    dest.mkdir(parents=True,exist_ok=True)
    for key,label in zip(("electrode_left","scattering_region","electrode_right"),("left","device","right")):
        a=read(src/(label+".cif"),format="cif")
        # Preserve the common CIF frame; no origin restoration is needed.
        a.pbc=False
        rows=["%s %.12f %.12f %.12f"%(s,*r) for s,r in zip(a.get_chemical_symbols(),a.positions)]
        (dest/(label+".xyz")).write_text(str(len(a))+"\nCartesian coordinates from CIF\n"+"\n".join(rows)+"\n")
        opts[key]=label+".xyz"
    for key in ("parameters","magnetic_parameters"):
        if key in opts and opts[key]!="default":
            opts[key]=os.path.relpath(job.parent/opts[key],dest)
    opts["output"]="calculated"
    new=dest/"converted.device.job"
    new.write_text("# Generated from the three CIF sources in a common frame.\n"+"\n".join(f"{k} = {v}" for k,v in opts.items())+"\n")
    return new

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("job");a=p.parse_args();print(convert(a.job))
