"""Optional Matplotlib plots without changing the scientific engine."""
from pathlib import Path

import numpy as np


def export_plots(report,directory):
    """Generate PNG plots for bands, DOS and charges in the filtered report."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("PNG plots require matplotlib. "
                          "Install with: python3 -m pip install matplotlib") from exc
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    saved=[]

    def finish(fig,name):
        path=directory/name
        fig.savefig(path,format="png",dpi=170,bbox_inches="tight")
        plt.close(fig)
        saved.append(str(path))

    if "band_energies_eV" in report:
        bands=np.asarray(report["band_energies_eV"],dtype=float)
        if bands.ndim==2 and len(bands):
            fig,ax=plt.subplots(figsize=(8.5,5))
            x=np.arange(len(bands))
            for ib in range(bands.shape[1]):
                ax.plot(x,bands[:,ib],lw=1.0)
            special=report.get("special_kpoints",[])
            for index,label in special:
                ax.axvline(index,color="0.7",ls="--",lw=0.7)
            if special:
                ax.set_xticks([i for i,label in special],
                              [label for i,label in special])
            ax.set_xlabel("Reciprocal-space path (k points)")
            ax.set_ylabel("Energy (eV)")
            ax.set_title("Band structure — EHT")
            ax.grid(True,alpha=.2)
            finish(fig,"bands.png")

    props=report.get("sampled_properties",{})
    dos=report.get("dos",props.get("dos"))
    if dos is not None:
        fig,ax=plt.subplots(figsize=(8,4.6))
        ax.plot(dos["energy_eV"],dos["states_per_eV"],lw=1.4)
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("DOS (states/eV)")
        ax.set_title("Density of states — EHT")
        ax.grid(True,alpha=.2)
        finish(fig,"dos.png")

    charges=report.get("mulliken_charges",props.get("mulliken_charges"))
    if charges is not None:
        values=np.array(charges,dtype=float)
        fig,ax=plt.subplots(figsize=(max(7,min(20,2+len(values)*.38)),4.2))
        xpos=np.arange(len(values))
        ax.bar(xpos,values)
        symbols=report.get("species",["A"]*len(values))
        labels=[f"{s}{i+1}" for i,s in enumerate(symbols)]
        ax.set_xticks(xpos,labels,rotation=90 if len(values)>16 else 0)
        ax.axhline(0,color="0.3",lw=0.8)
        ax.set_xlabel("Atom")
        ax.set_ylabel("Net Mulliken charge (e)")
        ax.set_title("Atomic charges — EHT")
        ax.grid(True,axis="y",alpha=.2)
        finish(fig,"charges.png")
    return saved
