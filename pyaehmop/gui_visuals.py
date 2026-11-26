"""Matplotlib visualizations for the optional Tk desktop interface.

Figures are presentation only and do not alter or recompute the
EHT/Stoner/SOC Hamiltonian or electronic populations.
"""
from itertools import product
import numpy as np


def _note(ax,message):
    ax.text(.5,.52,message,ha="center",va="center",transform=ax.transAxes,
            fontsize=11,color="#687487",wrap=True)
    ax.set_axis_off()


def _style(fig,ax,title,xlabel="",ylabel=""):
    fig.patch.set_facecolor("#ffffff")
    ax.set_title(title,fontsize=12,fontweight="medium",pad=14)
    if xlabel:ax.set_xlabel(xlabel)
    if ylabel:ax.set_ylabel(ylabel)
    ax.grid(alpha=.15)


def draw_structure(fig,view):
    """Atomistic 3D preview and simple distance-based bonds, not energy data."""
    fig.clear()
    ax=fig.add_subplot(111,projection="3d")
    coords=np.asarray(view.positions,dtype=float)
    if not len(coords):
        _note(ax,"No atoms found")
        return
    from matplotlib import colormaps
    symbols=list(view.symbols)
    types=sorted(set(symbols))
    ctab=colormaps["tab20"]
    colors={sym:ctab(i/max(1,len(types)-1)) for i,sym in enumerate(types)}
    if len(coords)<=250:
        try:
            from ase.data import covalent_radii,atomic_numbers
            radii=[float(covalent_radii[atomic_numbers.get(s,0)]) for s in symbols]
        except Exception:
            radii=[.75]*len(symbols)
        for i in range(len(coords)):
            for j in range(i+1,len(coords)):
                dist=np.linalg.norm(coords[i]-coords[j])
                if .38<dist<1.24*(radii[i]+radii[j]):
                    p,q=coords[i],coords[j]
                    ax.plot([p[0],q[0]],[p[1],q[1]],[p[2],q[2]],
                            color="#9ca7b5",lw=1.45,alpha=.7,zorder=1)
    for sym in types:
        indices=[i for i,s in enumerate(symbols) if s==sym]
        pts=coords[indices]
        ax.scatter(pts[:,0],pts[:,1],pts[:,2],label=sym,
                   c=[colors[sym]],s=120,edgecolors="#273448",
                   linewidths=.45,depthshade=True,zorder=3)
    if len(coords)<=36:
        for i,(sym,p) in enumerate(zip(symbols,coords),1):
            ax.text(*p,f"  {sym}{i}",fontsize=7.5,color="#26354a")
    vectors=np.asarray(view.vectors,dtype=float)
    all_points=[*coords]
    if vectors.ndim==2 and 0<len(vectors)<=3:
        origin=np.zeros(3)
        for codes in product((0.,1.),repeat=len(vectors)):
            point=origin+np.sum(np.asarray(codes)[:,None]*vectors,axis=0)
            all_points.append(point)
            for axis,coord in enumerate(codes):
                if coord==0:
                    q=point+vectors[axis]
                    ax.plot([point[0],q[0]],[point[1],q[1]],
                            [point[2],q[2]],lw=1.3,c="#0096a6",alpha=.55)
    bounds=np.asarray(all_points,dtype=float)
    center=np.mean([bounds.min(axis=0),bounds.max(axis=0)],axis=0)
    span=max(float(np.max(bounds.max(axis=0)-bounds.min(axis=0))),1.7)
    for fn,c in zip((ax.set_xlim,ax.set_ylim,ax.set_zlim),center):
        fn(c-span*.57,c+span*.57)
    ax.set_box_aspect((1,1,1))
    ax.set_xlabel("x (Å)");ax.set_ylabel("y (Å)");ax.set_zlabel("z (Å)")
    ax.view_init(elev=20,azim=34)
    ax.set_title(f"{view.title}  ·  {len(coords)} atoms",pad=12,fontsize=12)
    if len(types)<12:
        ax.legend(loc="upper right",fontsize=8)


def draw_bands(fig,report):
    fig.clear()
    ax=fig.add_subplot(111)
    energies=report.get("bands_eV",report.get("band_energies_eV"))
    if energies is None or not np.asarray(energies).size:
        _note(ax,"No band structure was requested or available.\n"
                 "Select Bands and a periodic structure.")
        return
    E=np.asarray(energies,dtype=float)
    if E.ndim!=2:
        _note(ax,"Band data has an unexpected shape.")
        return
    mu=float(report.get("fermi_eV",
             report.get("sampled_properties",{}).get("fermi_eV",0.0)))
    x=np.arange(E.shape[0])
    spin=report.get("spin_z_bands")
    if spin is not None and np.asarray(spin).shape==E.shape:
        from matplotlib.collections import LineCollection
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        polar=np.asarray(spin)
        norm=Normalize(-1,1)
        cmap="coolwarm"
        for j in range(E.shape[1]):
            pts=np.column_stack([x,E[:,j]-mu])
            if len(pts)<2:continue
            segments=np.stack([pts[:-1],pts[1:]],axis=1)
            values=(polar[:-1,j]+polar[1:,j])/2
            lc=LineCollection(segments,cmap=cmap,norm=norm,linewidths=1.5)
            lc.set_array(values)
            ax.add_collection(lc)
        ax.autoscale()
        cbar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=ax,shrink=.84,pad=.025)
        cbar.set_label("Projected spin ⟨σz⟩")
    else:
        for j in range(E.shape[1]):
            ax.plot(x,E[:,j]-mu,color="#2466b4",lw=1.2,alpha=.8)
    ax.axhline(0,color="#d06f28",lw=1.1,ls="--",alpha=.85)
    labels=report.get("special_kpoints",[])
    if labels:
        ax.set_xticks([int(i) for i,_ in labels],
                      ["Γ" if name=="G" else name for _,name in labels])
        for i,name in labels:
            ax.axvline(int(i),color="#b0b7c5",lw=.7,ls=":")
    ax.set_xlim(0,max(E.shape[0]-1,1))
    _style(fig,ax,"Electronic band structure","Reciprocal-space k path",
           "E − EF (eV)" if mu else "Energy (eV)")


def draw_dos(fig,report):
    fig.clear()
    ax=fig.add_subplot(111)
    d=report.get("dos") or report.get("sampled_properties",{}).get("dos")
    if not d or "energy_eV" not in d:
        _note(ax,"No density of states available.\nSelect DOS to calculate it.")
        return
    E=np.asarray(d["energy_eV"],dtype=float)
    mu=report.get("fermi_eV",report.get("sampled_properties",{}).get("fermi_eV"))
    shift=float(mu) if mu is not None else 0.
    if "up" in d and "down" in d:
        ax.plot(E-shift,d["up"],color="#ce4949",label="Spin ↑",lw=1.5)
        ax.plot(E-shift,d["down"],color="#3269ae",label="Spin ↓",lw=1.5)
        ax.plot(E-shift,np.asarray(d["up"])+np.asarray(d["down"]),
                color="#536171",label="Total",ls="--",lw=1.2,alpha=.6)
        ax.legend()
    else:
        ax.plot(E-shift,d["states_per_eV"],color="#147f89",lw=1.8)
        ax.fill_between(E-shift,d["states_per_eV"],color="#147f89",alpha=.14)
    if shift:
        ax.axvline(0,ls="--",c="#d06f28",alpha=.8,lw=1.1)
    _style(fig,ax,"Density of states","E − EF (eV)" if shift else "Energy (eV)",
           "States / (eV · cell)")


def draw_charges(fig,report):
    fig.clear()
    ax=fig.add_subplot(111)
    props=report.get("sampled_properties",{})
    charges=report.get("mulliken_charges",report.get(
        "atomic_charges",props.get("mulliken_charges")))
    if charges is None:
        _note(ax,"No Mulliken charges available.\nSelect Charges.")
        return
    c=np.asarray(charges,dtype=float)
    species=report.get("species",["A"]*len(c))
    labels=[f"{species[i] if i<len(species) else 'A'}{i+1}" for i in range(len(c))]
    ax.bar(np.arange(len(c)),c,color=["#d26759" if x>0 else "#2c8b9a" for x in c])
    ax.axhline(0,color="#526173",lw=.9)
    ax.set_xticks(np.arange(len(c)),labels,rotation=90 if len(c)>20 else 0)
    _style(fig,ax,"Mulliken atomic charges","Atom index","Charge (e)")


def draw_scf(fig,report):
    fig.clear()
    ax=fig.add_subplot(111)
    history=report.get("history",[])
    if not history:
        _note(ax,"SCF history is shown for spin-polarized calculations.")
        return
    iterations=[int(r.get("iter",i+1)) for i,r in enumerate(history)]
    moments=[float(r.get("total_moment_muB",0)) for r in history]
    errors=[float(r.get("spin_residual",1e-15)) for r in history]
    ax.plot(iterations,moments,color="#167d87",lw=1.8,label="Total moment")
    ax.set_ylabel("Moment (μB)")
    ax2=ax.twinx()
    ax2.semilogy(iterations,np.maximum(errors,1e-15),c="#d06f28",
                 ls="--",lw=1.4,label="Spin residual")
    ax2.set_ylabel("SCF residual (μB)")
    ax.set_xlabel("SCF iteration")
    ax.grid(alpha=.16)
    ax.set_title("Magnetic self-consistency",fontsize=12,pad=13)


def format_summary(report):
    values=[]
    for key,label in (
        ("project","Model"),
        ("mode","Magnetic mode"),
        ("type","Geometry"),
        ("atoms","Atoms"),
        ("orbitals","Orbitals"),
        ("electrons","Electrons"),
        ("electronic_energy_eV","Electronic energy (eV)"),
        ("average_electronic_energy_eV","Average electronic energy (eV)"),
        ("fermi_eV","Fermi level (eV)"),
        ("homo_eV","HOMO (eV)"),
        ("lumo_eV","LUMO (eV)"),
        ("band_path_gap_eV","Band-path gap (eV)"),
        ("gap_homo_lumo_eV","HOMO–LUMO gap (eV)"),
        ("total_moment_muB","Magnetic moment (μB)"),
        ("stoner_model_energy_eV","Stoner model energy (eV)"),
        ("iterations","SCF iterations"),
        ("converged","SCF converged")):
        v=report.get(key,report.get("sampled_properties",{}).get(key))
        if v is not None:
            values.append(f"{label:<32} {v}")
    warnings=report.get("warnings",report.get("alertas",[]))
    if warnings:
        values.append("\nWarnings / scientific qualifications:")
        values.extend(f"  • {v}" for v in warnings)
    return "\n".join(values) or "Open a file and run a calculation."
