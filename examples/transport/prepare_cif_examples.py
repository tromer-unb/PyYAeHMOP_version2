"""Create CIF alternatives from the distributed geometries; no calculations.

Uniform graphene uses a physical periodic cell with vacuum. Explicit device
CIFs are P1 coordinate interchange boxes, NOT the bulk electrode lattice.
"""
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
TR=ROOT/"examples/transport"
BF=ROOT/"examples/tutorials/bandfit"

def xyz(path):
    rows=[l.split() for l in path.read_text().splitlines()[2:] if l.strip()]
    return [r[0] for r in rows],np.array([list(map(float,r[1:4])) for r in rows])

def write_cif(path,symbols,coords,lengths):
    path.parent.mkdir(parents=True,exist_ok=True)
    header=["data_"+path.stem,"_symmetry_space_group_name_H-M 'P 1'","_symmetry_Int_Tables_number 1"]
    header += [f"_cell_length_{a} {v:.15g}" for a,v in zip("abc",lengths)]
    header += [f"_cell_angle_{a} 90" for a in ("alpha","beta","gamma")]
    header += ["loop_","_atom_site_label","_atom_site_type_symbol","_atom_site_fract_x","_atom_site_fract_y","_atom_site_fract_z","_atom_site_occupancy"]
    for i,(s,r) in enumerate(zip(symbols,coords/lengths),1):
        header.append(f"{s}{i} {s} "+" ".join(f"{x:.16f}" for x in r)+" 1")
    path.write_text("\n".join(header)+"\n")

def options(path):
    return dict(tuple(x.strip() for x in line.split("#",1)[0].split("=",1)) for line in path.read_text().splitlines() if "=" in line.split("#",1)[0])

def dump_job(path,opts,header):
    path.write_text("# "+header+"\n"+"\n".join(f"{k} = {v}" for k,v in opts.items())+"\n")

def main():
    s,p=xyz(TR/"graphene.xyz")
    write_cif(TR/"graphene.cif",s,p,np.array([2.459512,30.,4.26]))
    # Keep all existing XYZ jobs unchanged; create separate output paths.
    for path in [TR/"graphene.transport.job",*sorted(BF.glob("*B.transport.job")),*sorted(BF.glob("*B_smoke.transport.job"))]:
        o=options(path)
        o.pop("lattice_vectors")
        o.update(structure_format="cif",structure="graphene.cif" if path.parent==TR else "../../transport/graphene.cif",images="5 1 4",transport_axis="2",transverse_kpoints=o["transverse_kpoints"]+" 1",output=o["output"]+"_cif")
        dump_job(path.with_name(path.name.replace(".transport.job","_cif.transport.job")),o,"CIF has a 30 angstrom vacuum cell along b (y); transport along c (z).")
    # Explicit devices use independent regions in one Cartesian frame.
    # A common padded CIF box prevents wraparound; a global translation is harmless.
    jobs=[*sorted((ROOT/"examples/devices").glob("*.device.job")),*sorted((ROOT/"examples/tutorials").glob("*/run.device.job"))]
    for path in jobs:
        if path.name.endswith("_cif.device.job"): continue
        o=options(path);name=path.name.removesuffix(".device.job")
        folder=path.parent/"cif_sources"/name
        keys=("electrode_left","scattering_region","electrode_right")
        data=[xyz(path.parent/o[k]) for k in keys]
        coords=np.concatenate([p for s,p in data])
        origin=coords.min(axis=0)-5
        lengths=coords.max(axis=0)-origin+5
        for label,(s,p) in zip(("left","device","right"),data):
            write_cif(folder/(label+".cif"),s,p-origin,lengths)

        direct={"structure_format":"cif",**o}
        for key,label in zip(keys,("left","device","right")):
            direct[key]=str((folder/(label+".cif")).relative_to(path.parent))
        direct["output"]="results/cif_direct/"+name
        dump_job(path.with_name(name+"_cif.device.job"),direct,"Direct CIF input: common P1 coordinate box, physical repeats specified below.")

if __name__=="__main__":main()
