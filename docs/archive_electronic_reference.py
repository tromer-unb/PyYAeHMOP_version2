"""Compact numeric reference files used by build_tutorial.py.
Run after the electronic example suite; no new calculations.
"""
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
out=root/"examples/tutorials/reference/electronic";out.mkdir(parents=True,exist_ok=True)
sets=[("water",root/"pyaehmop/results/friendly_water_full",["molecular_levels.csv","dos.csv","mulliken.csv"]),("diamond",root/"pyaehmop/results/friendly_diamond_cif",["bands.csv","dos.csv"]),("fe_both",root/"pyaehmop/SPIN_TEST/results/Fe_wire_both",["bands_spin.csv","spin_dos.csv","scf_history.csv"])]
for name,src,files in sets:
    for f in files:
        p=src/f
        a=np.loadtxt(p,skiprows=1,delimiter=",")
        np.savetxt(out/(name+"_"+f),a,delimiter=",",fmt="%.9g",header=p.read_text().splitlines()[0],comments="")
