"""Standalone command line: python -m pyaehmop input [options]."""
import argparse
import json
from pathlib import Path
import numpy as np
from .input import read_input
from .engine import calculate
from .export import export_csv

def _tolist(x):
    if isinstance(x,np.ndarray):
        return x.tolist()
    if isinstance(x,np.generic):
        return x.item()
    raise TypeError(f"Cannot serialize: {type(x)}")

def _main():
    p=argparse.ArgumentParser(description="Extended Hückel in Python/NumPy/SciPy, independent of the YAeHMOP executable")
    p.add_argument("input", type=Path,help="Traditional YAeHMOP input or friendly .job file")
    p.add_argument("-o","--output",type=Path,help="Destination JSON file")
    p.add_argument("-j","--workers",type=int,default=None,help="CPU processes (overrides friendly configuration)")
    p.add_argument("--quadrature",nargs=3,type=int,default=None,
                   metavar=("LAGUERRE","LEGENDRE","AZIMUTH"),help="Overlap quadrature nodes")
    p.add_argument("--cutoff",type=float,default=None,help="Maximum periodic overlap radius (Å)")
    p.add_argument("--no-matrices",action="store_true",help="Omit matrices and coefficients from JSON")
    p.add_argument("--csv-dir",type=Path,help="Export bands, DOS, charges and orbitals as CSV")
    p.add_argument("--png-dir",type=Path,help="Optional PNG plots in friendly mode")
    p.add_argument("--easy",action="store_true",help="Interpret input as friendly configuration")
    args=p.parse_args()
    if args.input.name.endswith('.device.job'):
        from .device_job import run_job
        print(run_job(args.input,args.workers))
        return 0
    if args.input.name.endswith('.transport.job'):
        from .transport_job import run_job
        print(run_job(args.input,args.workers))
        return 0
    from .friendly import is_easy_input
    if args.easy or is_easy_input(args.input):
        return _run_easy(args)
    s=read_input(args.input)
    out=calculate(s,tuple(args.quadrature or (24,20,12)),
                  args.workers or 1,args.cutoff if args.cutoff is not None else 11.)
    if args.no_matrices:
        out={k:v for k,v in out.items() if k not in
             ("overlap_matrix","hamiltonian_matrix","coefficients")}
        if "fmo" in out:
            fmo=out["fmo"]
            out["fmo"]={k:v for k,v in fmo.items() if k not in
                        ("ao_to_fmo_transform","fmo_overlap_matrix","fmo_hamiltonian_matrix")}
            out["fmo"]["fragments"]=[{k:v for k,v in frag.items()
                                     if k!="local_ao_coefficients"}
                                    for frag in fmo["fragments"]]
    text=json.dumps(out,default=_tolist,indent=2,ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+"\n",encoding="utf-8")
    compact={k:v for k,v in out.items() if k not in
             ("overlap_matrix","hamiltonian_matrix","coefficients","orbital_energies_eV",
              "occupations","mulliken_charges","orbital_populations",
              "orbital_names","kpoints","band_energies_eV","bond_analyses","fmo")}
    if "fmo" in out:
        compact["fmo"]=[{"fragment":r["fragment"],"atoms":r["atoms"],
                         "energy_fragment_eV":r["fragment_electronic_energy_eV"]}
                        for r in out["fmo"]["fragments"]]
    if "bond_analyses" in out:
        compact["computed_bonds"]=[{"atoms":r["atoms"],
                                      "integrated_coop":r["integrated_coop"],
                                      "integrated_cohp":r["integrated_cohp"]}
                                     for r in out["bond_analyses"]]
    if "sampled_properties" in compact:
        av=compact["sampled_properties"]
        compact["sampled_properties"]={k:v for k,v in av.items()
                                           if k not in ("dos","orbital_populations")}
    print(json.dumps(compact,default=_tolist,indent=2,ensure_ascii=False))
    if args.output:
        print("Full result:",args.output)
    if args.csv_dir:
        paths=export_csv(out,args.csv_dir)
        print("CSV files:",", ".join(paths))

def _run_easy(args):
    """Human-friendly interface; does not alter the legacy YAeHMOP reader."""
    from .friendly import read_easy_job, calculate_easy_job
    from .friendly import _normalize
    job=read_easy_job(args.input)
    if args.workers is not None:
        if args.workers<1:
            raise ValueError("--workers must be >=1")
        job.workers=args.workers
    if args.quadrature is not None:
        job.quadrature=tuple(args.quadrature)
    if args.cutoff is not None:
        job.cutoff=args.cutoff
    report=calculate_easy_job(job)
    if args.no_matrices:
        report={k:v for k,v in report.items() if k not in
                ("overlap_matrix","hamiltonian_matrix","coefficients")}
        if "fmo" in report:
            fmo=report["fmo"]
            report["fmo"]={k:v for k,v in fmo.items() if k not in
                           ("ao_to_fmo_transform","fmo_overlap_matrix","fmo_hamiltonian_matrix")}
            report["fmo"]["fragments"]=[{k:v for k,v in frag.items()
                                         if k!="local_ao_coefficients"}
                                        for frag in fmo["fragments"]]
    target=args.output or (job.output_dir/"summary.json")
    write_json=args.output is not None or "json" in job.outputs
    if write_json:
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps(report,default=_tolist,indent=2,ensure_ascii=False)+"\n",
                          encoding="utf-8")
    csv_path=args.csv_dir if args.csv_dir else job.output_dir
    csv_files=[]
    if "csv" in job.outputs or args.csv_dir:
        csv_files=export_csv(report,csv_path)
    print(f"Format: {job.source_type.upper()} | {report['type']} | "
          f"{report['atoms']} atoms, {report['orbitals']} orbitals")
    print("Calculations:",", ".join(sorted(job.tasks)),"| workers:",job.workers)
    if "electronic_energy_eV" in report:
        print("Electronic energy (eV):",round(report["electronic_energy_eV"],8))
    if "mulliken_charges" in report:
        print("Mulliken charges:",[round(float(c),6) for c in report["mulliken_charges"]])
    if "sampled_properties" in report:
        props=report["sampled_properties"]
        if "average_electronic_energy_eV" in props:
            print("Average electronic energy (eV):",round(props["average_electronic_energy_eV"],8))
        if "mulliken_charges" in props:
            print("Mulliken charges:",[round(float(c),6) for c in props["mulliken_charges"]])
    if "kpoints" in report:
        print("Band structure k points:",len(report["kpoints"]))
    if "fmo" in report:
        print("FMO fragments:",len(report["fmo"]["fragments"]))
    if "bond_analyses" in report:
        print("Analyzed bonds:",len(report["bond_analyses"]))
    for warning in report.get("warnings",[]):
        print("WARNING:",warning)
    if write_json:
        print("JSON:",target)
    if csv_files:
        print("CSV:",", ".join(csv_files))
    if "png" in job.outputs or args.png_dir:
        from .plots import export_plots
        plots=export_plots(report,args.png_dir or job.output_dir)
        if plots:
            print("Plots:",", ".join(plots))
    return 0



def main():
    """Command line entry point with readable error messages."""
    try:
        return _main()
    except (ValueError,FileNotFoundError,ImportError,NotImplementedError) as exc:
        raise SystemExit("Input/calculation error: "+str(exc)) from None


if __name__=="__main__":
    main()
