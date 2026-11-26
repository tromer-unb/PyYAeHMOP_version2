"""CSV exporters that respect selected tasks, without pandas."""
import csv
from pathlib import Path


def _write(directory,name,header,rows):
    p=directory/name
    with p.open("w",newline="",encoding="utf-8") as out:
        writer=csv.writer(out)
        writer.writerow(header)
        writer.writerows(rows)
    return str(p)


def export_csv(report,directory):
    """Export only datasets present in the calculation report."""
    outdir=Path(directory)
    outdir.mkdir(parents=True,exist_ok=True)
    files=[]
    if report["type"]=="molecule":
        if "orbital_energies_eV" in report and "occupations" in report:
            files.append(_write(outdir,"molecular_levels.csv",
                         ["orbitals","energy_eV","electron_occupation"],
                         ((i+1,float(E),float(o)) for i,(E,o) in
                          enumerate(zip(report["orbital_energies_eV"],report["occupations"])))))
        if "mulliken_charges" in report:
            files.append(_write(outdir,"mulliken.csv",
                         ["atom","electronic_charge"],
                         ((i+1,float(c)) for i,c in enumerate(report["mulliken_charges"]))))
    else:
        if "kpoints" in report and "band_energies_eV" in report:
            files.append(_write(outdir,"bands.csv",
                         ["k_index","k_x","k_y","k_z","band","energy_eV"],
                         ((i,float(k[0]),float(k[1]),float(k[2]),j+1,float(E))
                          for i,(k,energies) in
                          enumerate(zip(report["kpoints"],report["band_energies_eV"]))
                          for j,E in enumerate(energies))))
        props=report.get("sampled_properties",{})
        if "mulliken_charges" in props:
            files.append(_write(outdir,"mulliken.csv",
                         ["atom","electronic_charge"],
                         ((i+1,float(c)) for i,c in enumerate(props["mulliken_charges"]))))
    props=report.get("sampled_properties",{})
    dos=report.get("dos",props.get("dos"))
    if dos is not None:
        files.append(_write(outdir,"dos.csv",["energy_eV","states_per_eV"],
                            zip(dos["energy_eV"],dos["states_per_eV"])))
    energies={}
    for name in ("electronic_energy_eV","homo_eV","lumo_eV",
                 "gap_homo_lumo_eV","band_path_gap_eV","vbm_eV","cbm_eV"):
        if name in report:
            energies[name]=report[name]
    for name in ("fermi_eV","average_electronic_energy_eV"):
        if name in props:
            energies[name]=props[name]
    if energies:
        files.append(_write(outdir,"energies.csv",["property","value_eV"],
                            energies.items()))
    fmo=report.get("fmo")
    if fmo:
        files.append(_write(outdir,"fmo_levels.csv",
                     ["fragment","orbital_index","energy_eV","electron_occupation"],
                     ((f["fragment"],i+1,float(E),float(o))
                      for f in fmo["fragments"] for i,(E,o) in
                      enumerate(zip(f["orbital_energies_eV"],f["occupations"])))))
    bonds=report.get("bond_analyses")
    if bonds:
        files.append(_write(outdir,"bonds.csv",
                     ["type","curva","atom_a","atom_b","R_a","R_b","R_c",
                      "integrated_coop","integrated_cohp"],
                     ((bond["requested_type"],bond["curva"],*bond["atoms"],
                       *(bond.get("translation",(0,0,0))),
                       bond["integrated_coop"],bond["integrated_cohp"])
                      for bond in bonds)))
    return files
