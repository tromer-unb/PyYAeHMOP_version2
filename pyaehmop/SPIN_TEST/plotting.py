"""Plots of spin-up/down bands and SOC with spin projection."""
from pathlib import Path
import numpy as np

def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _color_bands(ax,bands,sz,fermi,alpha=1.0,lw=1.05):
    """One line per band, colored according to <sigma_z>."""
    from matplotlib.collections import LineCollection
    import matplotlib as mpl
    nk,nb=bands.shape
    x=np.arange(nk,dtype=float)
    colors=mpl.colormaps["coolwarm"]
    for ib in range(nb):
        pts=np.column_stack((x,bands[:,ib]-fermi))
        seg=np.stack((pts[:-1],pts[1:]),axis=1)
        pol=(sz[:-1,ib]+sz[1:,ib])/2
        lc=LineCollection(seg,array=pol,cmap=colors,norm=mpl.colors.Normalize(-1,1),
                          lw=lw,alpha=alpha)
        ax.add_collection(lc)
    ax.set_xlim(0,nk-1)
    ax.autoscale_view(scaley=True)
    return mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1,1),cmap=colors)


def _exact_up_down(ax,bands,sz,fermi,**kwargs):
    """Without SOC, up/down are good quantum numbers."""
    n=len(bands[0])//2
    up=[];dn=[]
    for en,pol in zip(bands,sz):
        iu=np.where(pol>0.5)[0];idn=np.where(pol<-.5)[0]
        if len(iu)!=n or len(idn)!=n:
            raise ValueError("Esperados canais up/down exatos quando SOC=0")
        up.append(np.sort(en[iu]));dn.append(np.sort(en[idn]))
    ax.plot(np.array(up)-fermi,color="#cf3f3f",lw=1.0,alpha=.8)
    ax.plot(np.array(dn)-fermi,color="#356abb",lw=1.0,alpha=.8)
    return np.array(up),np.array(dn)


def _style(ax,labels,title):
    ax.axhline(0,c="0.25",ls="--",lw=.8)
    ax.set_title(title,fontsize=11)
    ax.set_ylabel("$E-E_F$ (eV)")
    ax.grid(alpha=.17)
    if labels:
        ax.set_xticks([i for i,name in labels],
                      [name.replace("G","Γ") if name=="G" else name for i,name in labels])


def plot_single(job,result,out_dir):
    plt=_plt()
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),layout="constrained")
    bands=result["bands_eV"]; sz=result["spin_z_bands"]; mu=result["fermi_eV"]
    if job.settings.spin and not job.settings.soc:
        _exact_up_down(axes[0],bands,sz,mu)
        axes[0].plot([],[],color="#cf3f3f",label="Spin ↑ (majority)")
        axes[0].plot([],[],color="#356abb",label="Spin ↓ (minority)")
        axes[0].legend(fontsize=9,loc="upper right")
    else:
        sm=_color_bands(axes[0],bands,sz,mu)
        if job.settings.spin:
            fig.colorbar(sm,ax=axes[0],label="⟨σz⟩ (−1 a +1)",shrink=.8)
    _style(axes[0],job.special_labels,
           "Spin-resolved bands" + (" e SOC" if job.settings.soc else ""))
    E=np.array(result["dos"]["energy_eV"])-mu
    axes[1].plot(E,result["dos"]["up"],c="#cf3f3f",label="Projected DOS ↑")
    axes[1].plot(E,result["dos"]["down"],c="#356abb",label="Projected DOS ↓")
    axes[1].set_xlabel("$E-E_F$ (eV)")
    axes[1].set_ylabel("DOS (states/eV/cell)")
    axes[1].set_title("Spin-resolved density of states")
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=.18)
    fig.suptitle("PyYAeHMOP / SPIN_TEST — experimental EHT + Stoner + SOC",fontsize=13)
    path=out_dir/"bands_spin_soc.png"
    fig.savefig(path,dpi=175)
    plt.close(fig)
    return str(path)


def comparison_plot(cases,source,outfile):
    """Spin-only, SOC-only, spin+SOC and eigenvalue shifts.

    Based on published Fe-wire geometry; does not reproduce published bands quantitatively.
    """
    plt=_plt()
    outfile=Path(outfile);outfile.parent.mkdir(parents=True,exist_ok=True)
    spin=cases["spin"];soc=cases["soc"];both=cases["both"]
    labels=spin["special_kpoints"]
    fig,axs=plt.subplots(2,2,figsize=(14,10))
    fig.subplots_adjust(left=.08,right=.90,top=.88,bottom=.16,
                        hspace=.40,wspace=.27)
    band_up,band_dn=_exact_up_down(
        axs[0,0],spin["bands_eV"],spin["spin_z_bands"],spin["fermi_eV"])
    axs[0,0].plot([],[],c="#cf3f3f",label="Spin ↑ (majority)")
    axs[0,0].plot([],[],c="#356abb",label="Spin ↓ (minority)")
    axs[0,0].legend(fontsize=9)
    _style(axs[0,0],labels,"(a) Stoner SCF, no SOC: exact ↑ and ↓ channels")

    sm=_color_bands(axs[0,1],both["bands_eV"],both["spin_z_bands"],both["fermi_eV"])
    _style(axs[0,1],labels,"(b) Stoner SCF + SOC: color denotes ⟨σz⟩ projection")
    fig.colorbar(sm,ax=axs[0,1],label="⟨σz⟩ (projected spin)",shrink=.85)

    axs[1,0].plot(soc["bands_eV"]-soc["fermi_eV"],c="0.45",lw=.75,alpha=.75)
    _style(axs[1,0],labels,"(c) SOC only: no spontaneous magnetization")

    # Pairwise differences are based only on energy-sorted bands,
    # NOT on adiabatic tracking of individual states.
    de=(both["bands_eV"]-both["fermi_eV"] -
        (spin["bands_eV"]-spin["fermi_eV"]))*1000
    axs[1,1].plot(de,c="#5a349e",lw=.8,alpha=.6)
    axs[1,1].axhline(0,c="0.3",lw=.8)
    axs[1,1].set_ylabel("ΔE of energy-sorted bands (meV)")
    axs[1,1].grid(alpha=.17)
    axs[1,1].set_title("(d) SOC shifts (energy-sorted bands, not state-tracked)",fontsize=10)
    if labels:
        axs[1,1].set_xticks([i for i,n in labels],
                          ["Γ" if n=="G" else n for i,n in labels])
    for ax in (axs[0,0],axs[0,1],axs[1,0]):
        ax.set_ylim(-5,7)
    for ax in axs.flat:ax.set_xlabel("Points along the k path")
    fig.suptitle("Monatomic Fe wire · EHT + Stoner + SOC (independent simulation)",
                 fontsize=15,fontweight="semibold",y=.96)
    fig.text(.5,.052,
        "Geometry inspired by Autès et al., J. Phys.: Condens. Matter 18, 6785 (2006), DOI 10.1088/0953-8984/18/29/018.\n"
        "EHT parameters and SPIN_TEST approximations do not quantitatively reproduce published bands.",
        ha="center",va="bottom",fontsize=9)
    fig.savefig(outfile,dpi=180)
    plt.close(fig)
    return str(outfile)



def plot_soc_zoom(spin_case,both_case,outfile,krange=(54,66),
                  yrange=(-1.35,-.78)):
    """Zoomed region in which SOC changes <sigma_z>.

    Plot limits are set for the Fe_wire demonstration case,
    not an automatic selection of a topological gap.
    """
    plt=_plt()
    outfile=Path(outfile);outfile.parent.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.5),layout="constrained")
    spin_energies=spin_case["bands_eV"]
    both_energies=both_case["bands_eV"]
    _exact_up_down(axes[0],spin_energies,spin_case["spin_z_bands"],
                   spin_case["fermi_eV"])
    sm=_color_bands(axes[1],both_energies,both_case["spin_z_bands"],
                    both_case["fermi_eV"],lw=2.4)
    for ax,title in zip(axes,("Stoner SCF without SOC — ↑ and ↓",
                              "Stoner SCF + SOC — ⟨σz⟩ projection")):
        ax.set_xlim(*krange)
        ax.set_ylim(*yrange)
        ax.set_xlabel("K-point index along Γ–X")
        ax.set_ylabel("$E-E_F$ (eV)")
        ax.set_title(title)
        ax.grid(alpha=.24)
    fig.colorbar(sm,ax=axes[1],label="⟨σz⟩")
    fig.suptitle("Monatomic Fe: detail of SOC-induced spin mixing",
                 fontsize=14,fontweight="semibold")
    fig.savefig(outfile,dpi=185)
    plt.close(fig)
    return str(outfile)
