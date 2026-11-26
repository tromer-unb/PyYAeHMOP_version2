"""Isolated experimental interface: never modifies pyaehmop outside SPIN_TEST.

Usage:
  python3 -m pyaehmop.SPIN_TEST examples/Fe_wire_spin.spin
  python3 -m pyaehmop.SPIN_TEST --demo
"""
import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path
import numpy as np

from .config import read_spin_job
from .solver import scf, evaluate_bands, spin_dos
from .plotting import plot_single, comparison_plot, plot_soc_zoom


def jsonify(x):
    if isinstance(x,np.ndarray):
        return x.tolist()
    if isinstance(x,np.generic):
        return x.item()
    raise TypeError(f"JSON serialization not available for {type(x)}")


def _csv(path,header,rows):
    with Path(path).open("w",newline="",encoding="utf-8") as stream:
        writer=csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def run_case(filename,verbose=False):
    job=read_spin_job(filename)
    job.outdir.mkdir(parents=True,exist_ok=True)
    result,reference=scf(
        job.structure,job.kpoints,job.weights,job.settings,
        quadrature=job.quadrature,cutoff=job.cutoff,verbose=verbose)
    if not result["converged"]:
        raise RuntimeError(
            "SCF did not converge within the configured iterations. "
            f"Final residual {result['final_residual']:.3e}; "
            "check mixing, initial_moment, kmesh and max_iter.")
    bands=np.empty((0,0))
    projected=np.empty((0,0))
    if job.band_path is not None:
        bands,projected=evaluate_bands(
            job.structure,job.band_path,reference,job.quadrature,job.cutoff)
    grid,up,down=spin_dos(
        result["integration_energies_eV"],result["sz_integration"],
        job.weights,sigma=job.dos_sigma,points=job.dos_points)
    out={
        "project":"Experimental SPIN_TEST; not validated for predictive use",
        "input":str(job.path),
        "structure":str(job.source),
        "species":[atom.symbol for atom in job.structure.atoms],
        "mode":("spin+soc" if job.settings.spin and job.settings.soc else
                "spin" if job.settings.spin else "soc"),
        "parameters":asdict(job.settings),
        "periodic":job.structure.vectors is not None,
        "k_path":job.band_path,
        "special_kpoints":job.special_labels,
        "bands_eV":bands,
        "spin_z_bands":projected,
        "dos":{"energy_eV":grid,"up":up,"down":down,
               "sigma_eV":job.dos_sigma},
        **result,
        "warnings":[
            "Not DFT; EHT has not been reparametrized for Stoner/SOC.",
            "Onsite atomic L.S SOC, without two-center SOC terms.",
            "Self-consistent collinear-z Stoner, not noncollinear spin SCF.",
            "With SOC, up/down are not exact channels; use <sigma_z> projections.",
            "Check convergence in k mesh, cutoff, quadrature and temperature."
        ]}
    outfile=job.outdir/"spin_soc_result.json"
    outfile.write_text(json.dumps(out,default=jsonify,indent=2,
                                  ensure_ascii=False)+"\n",encoding="utf-8")
    if len(bands):
        _csv(job.outdir/"bands_spin.csv",
             ["k_index","k_x","k_y","k_z","band","energy_eV","spin_z"],
             ((i+1,*k,j+1,float(E),float(proj))
              for i,(k,energies,polar) in enumerate(zip(job.band_path,bands,projected))
              for j,(E,proj) in enumerate(zip(energies,polar))))
    _csv(job.outdir/"spin_dos.csv",["energy_eV","dos_up","dos_down","dos_total"],
         zip(grid,up,down,up+down))
    _csv(job.outdir/"scf_history.csv",
         ["iteration","residuo_m_muB","total_moment_muB","fermi_eV"],
         ((v["iter"],v["spin_residual"],v["total_moment_muB"],v["fermi_eV"])
          for v in result["history"]))
    if len(bands) and job.image:
        path=plot_single(job,out,job.outdir)
        print("Bands and DOS:",path)
    print(f"SCF {out['mode']:8s}: converged={result['converged']}, "
          f"iterations={result['iterations']}, "
          f"m_z={result['total_moment_muB']:.6f} muB, "
          f"E_model={result['stoner_model_energy_eV']:.6f} eV")
    print("Result:",outfile)
    return out,job


def demo():
    root=Path(__file__).resolve().parent
    cases={}
    for kind in ("spin","soc","both"):
        report,job=run_case(root/"examples"/f"Fe_wire_{kind}.spin")
        cases[kind]=report
    out=root/"results"/"Fe_wire_literature_comparison.png"
    comparison_plot(cases,
        "Autes et al. (2006) DOI: 10.1088/0953-8984/18/29/018",out)
    zoom=root/"results"/"Fe_wire_SOC_spin_mixing_detail.png"
    plot_soc_zoom(cases["spin"],cases["both"],zoom)
    print("SOC SPIN-MIXING DETAIL:",zoom)
    summary={
        kind:{
            "moment_muB":float(case["total_moment_muB"]),
            "converged":case["converged"],
            "iterations":case["iterations"],
            "EF_eV":float(case["fermi_eV"]),
            "model_energy_eV":float(case["stoner_model_energy_eV"]),
            "mixed_spin_states_on_path":int(np.sum(
                np.abs(case["spin_z_bands"])<0.8)),
        } for kind,case in cases.items()
    }
    diff=(cases["both"]["bands_eV"]-cases["both"]["fermi_eV"] -
          (cases["spin"]["bands_eV"]-cases["spin"]["fermi_eV"]))
    summary["ordered_band_difference_meV"]={
        "max_abs":float(np.max(abs(diff))*1000),
        "rms":float(np.sqrt(np.mean(diff**2))*1000)}
    path=root/"results"/"Fe_wire_comparison_summary.json"
    path.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",
                    encoding="utf-8")
    print("\nCOMPARISON FIGURE:",out)
    print("SUMMARY:",path)
    return out


def main():
    parser=argparse.ArgumentParser(
        description="Experimental EHT Stoner SCF and spin-orbit, isolated in SPIN_TEST")
    parser.add_argument("input",nargs="?",help="Spin input file; first lines specify format and structure")
    parser.add_argument("--demo",action="store_true",
                        help="Generate three Fe-wire calculations and a comparison figure")
    parser.add_argument("--verbose",action="store_true",
                        help="Print the magnetic SCF residual at each iteration")
    args=parser.parse_args()
    try:
        if args.demo:
            demo()
        elif args.input:
            run_case(args.input,args.verbose)
        else:
            parser.error("Provide a .spin file or --demo")
    except (ValueError,FileNotFoundError,NotImplementedError,RuntimeError) as exc:
        parser.exit(2,"Erro SPIN_TEST: "+str(exc)+"\n")


if __name__=="__main__":
    main()
