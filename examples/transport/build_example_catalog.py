"""Generate the paired-input index from checked-in job files."""
from pathlib import Path
from examples.transport.prepare_cif_examples import ROOT,TR,BF,options
import os
def link(p,base,label=None):
    return f"[{label or p.name}]({os.path.relpath(p,base)})"
def main():
    dest=ROOT/"examples/README.md";base=dest.parent
    lines=["# All runnable examples — start here","","Clone the repository, install with pip or Conda ([installation](../docs/installation.md)), then run the following commands **from the repository root**. No PBS or external transport software is required. Each transport row has two real executable job files, not just geometry files.","","## Uniform periodic transport (XYZ and CIF)","","Run either linked job with `python -m pyaehmop PATH --workers 4`. The CIF reads a 3D cell with 30 Å vacuum; the XYZ uses two explicit periodic vectors. Same material, parameters and energy reference. See [format comparison](../docs/transport-cif.md).","","| Model | XYZ full | CIF full | XYZ quick | CIF quick |","|---|---|---|---|---|"]
    for p in sorted(BF.glob("*B.transport.job")):
        name=p.name.removesuffix(".transport.job")
        q=BF/(name+"_cif.transport.job")
        a=BF/(name+"_smoke.transport.job");b=BF/(name+"_smoke_cif.transport.job")
        lines.append("| "+name+" | "+" | ".join(link(x,base) for x in (p,q,a,b))+" |")
    lines += ["","Shared geometries: "+link(TR/"graphene.xyz",base)+" · "+link(TR/"graphene.cif",base)+". Full p-DZ7 4B shortcut: "+link(TR/"graphene.transport.job",base)+" · "+link(TR/"graphene_cif.transport.job",base)+".","","## Independent electrodes and scattering region (15 paired cases)","","Run **either** job directly with `python -m pyaehmop PATH --workers 4`. No manual conversion is required. The CIF jobs use ordered P1 region files with a shared padded box and coordinate frame (no auxiliary JSON), not three independent bulk CIF cells. The physical repeats remain in the job. The XYZ jobs use the original Cartesian geometries.","","| Case | XYZ input | CIF input | CIF region files | Explanation / results |","|---|---|---|---|---|"]
    dev=ROOT/"examples/devices";tut=ROOT/"examples/tutorials"
    names=["pristine","vacancy","adsorbate","asymmetric","graphene_mos2","spin","soc","spin_soc","graphene_x","nanoribbon","wide_vacancy","potential_barrier","spin_parallel","spin_antiparallel","graphene_mos2_graphene"]
    for n in names:
        old=(dev/(n+".device.job")).exists();p=dev/(n+".device.job") if old else tut/n/"run.device.job"
        q=p.with_name(p.name.replace(".device.job","_cif.device.job"))
        src=p.parent/"cif_sources"/p.name.removesuffix(".device.job")
        expl=dev/"README.md" if old else p.parent/"README.md"
        o=options(p)
        lines.append("| "+n+" | "+link(p,base)+" | "+link(q,base)+" | "+link(src,base,"L/D/R")+" | "+link(expl,base,"Guide")+" |")
    lines += ["","## Electronic structure and spin/SOC","","These are the original molecular/crystalline inputs, not transport jobs. Standard XYZ is Cartesian and nonperiodic for molecules, CIF provides crystal cell and fractional sites, ASE/POSCAR retain their lattice. All remain in the full tutorial.","","| Case | Executable input |","|---|---|"]
    for p in sorted((ROOT/"pyaehmop/examples_friendly").glob("*.job")):
        lines.append(f"| {p.stem} | {link(p,base)} |")
    for p in sorted((ROOT/"pyaehmop/SPIN_TEST/examples").glob("*.spin")):
        lines.append(f"| {p.stem} | {link(p,base)} |")
    lines += ["","Run ordinary jobs with `python -m pyaehmop PATH`; run .spin jobs with `python -m pyaehmop.SPIN_TEST PATH`.","","## Full illustrated manual and validation","","- [Complete PDF with all input files and device plots](../docs/PyYAeHMOP_v2_tutorial.pdf)","- [Original tutorial index and numerical reference results](tutorials/README.md)","- [Direct CIF device validation](transport/cif_reference/direct_device_validation.json)","- [Uniform XYZ/CIF validation](transport/cif_reference/validation.json)","- [Parameters](../docs/parameters.md) and [electrode construction](../docs/device-geometry.md)","","All original examples are retained. New devices are illustrative frozen-Hamiltonian models; I(V) is not self-consistent finite-bias NEGF."]
    dest.write_text("\n".join(lines)+"\n")
if __name__=="__main__":main()
