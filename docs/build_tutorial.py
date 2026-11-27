"""Create the illustrated manual from the shipped inputs and archived results.
Run from any directory with numpy and matplotlib. No scientific run is needed.
"""
from pathlib import Path
import textwrap, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"; DEV=ROOT/"examples/devices"; TUT=ROOT/"examples/tutorials"; BF=TUT/"bandfit"; TR=ROOT/"examples/transport"
NAVY="#17334a"; BLUE="#245d7b"; TEAL="#007f83"; ORANGE="#bb641e"; GRAY="#56636c"; LIGHT="#e9f0f3"
COLORS={"C":"#334b5d","O":"#be4938","Mo":"#247e8a","S":"#d39a19","H":"#8c9da8","Fe":"#965d36"}
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"axes.labelsize":8,"axes.titlesize":9,"xtick.labelsize":8,"ytick.labelsize":8,"pdf.fonttype":42,"pdf.compression":9,"svg.fonttype":"none"})
URL="https://github.com/tromer-unb/PyYAeHMOP_version2"
def xyz(p):
    lines=p.read_text().splitlines()
    names=[];coords=[]
    for l in lines[2:]:
        w=l.split()
        if w:names.append(w[0]);coords.append([float(x) for x in w[1:4]])
    return np.array(names),np.array(coords)
def opts(p):return dict(tuple(x.strip() for x in l.split("=",1)) for l in p.read_text().splitlines() if "=" in l and not l.lstrip().startswith("#"))
def jobpaths(p):
    o=opts(p);return [p.parent/o[k] for k in ("electrode_left","scattering_region","electrode_right")]
def geo(ax,p,side=False,labels=True):
    o=opts(p); paths=jobpaths(p); rr=np.array([float(x) for x in o["right_repeat"].split()]); axis=int(np.argmax(abs(rr))); others=[i for i in range(3) if i!=axis]
    a=others[0] if axis!=0 else 2
    if side:a=1
    names=["Left electrode","Scattering region","Right electrode"]; colors=["#dce8f1","#d9eee6","#f6e6d6"]; allpos=[]
    for path in paths:allpos.append(xyz(path))
    points=[q[:,axis] for _,q in allpos]; boundaries=[(points[i].max()+points[i+1].min())/2 for i in (0,1)]
    lo=points[0].min()-.5;hi=points[-1].max()+.5
    limits=[lo,*boundaries,hi]
    for i,(sym,q) in enumerate(allpos):
        ax.axvspan(limits[i],limits[i+1],color=colors[i],alpha=.9,zorder=0)
        for j in range(len(q)):
            for k in range(j):
                d=np.linalg.norm(q[j]-q[k])
                if d<1.65 and sym[j]==sym[k]=="C" or 0.8<d<1.8 and {"C","O"}=={sym[j],sym[k]}:
                    ax.plot(q[[j,k],axis],q[[j,k],a],color="#9aa7ac",lw=.7,zorder=1)
        for el in sorted(set(sym)):
            m=sym==el;ax.scatter(q[m,axis],q[m,a],s=15 if len(q)>24 else 28,c=COLORS.get(el,"#805b8d"),edgecolors="white",linewidths=.4,zorder=3)
        if labels:ax.text((limits[i]+limits[i+1])/2,1.04,names[i]+f" ({len(q)})",transform=ax.get_xaxis_transform(),ha="center",va="bottom",fontsize=7,color=NAVY)
    ax.set_xlim(lo,hi);ax.set_xlabel("xyz"[axis]+" (Å) — transport direction")
    ax.set_ylabel("xyz"[a]+" (Å)");ax.set_ylim(min(q[:,a].min() for _,q in allpos)-.65,max(q[:,a].max() for _,q in allpos)+.65)
    el=sorted({s for sy,_ in allpos for s in sy})
    ax.legend(handles=[Line2D([],[],marker="o",linestyle="none",color=COLORS.get(e,"#805b8d"),label=e,markersize=5) for e in el],loc="upper left",frameon=False,ncol=len(el),fontsize=7)
    ax.spines[["top","right"]].set_visible(False)
    return axis
def geofile(name,p,side=False):
    f,ax=plt.subplots(figsize=(8.1,2.2));geo(ax,p,side);f.tight_layout(pad=1.1)
    f.savefig(DOC/"figures"/(name+"_regions.svg"));plt.close(f)
for name,p,side in [("pristine",DEV/"pristine.device.job",False),("vacancy",DEV/"vacancy.device.job",False),("adsorbate",DEV/"adsorbate.device.job",True),("nanoribbon",TUT/"nanoribbon/run.device.job",False),("graphene_mos2",DEV/"graphene_mos2.device.job",False),("graphene_mos2_graphene",TUT/"graphene_mos2_graphene/run.device.job",False)]:
    geofile(name,p,side)
out=DOC/"PyYAeHMOP_v2_tutorial.pdf"
P=0
def page(pdf,section,title,sub=""):
    global P;P+=1
    f=plt.figure(figsize=(8.27,11.69),facecolor="white")
    f.patches.append(Rectangle((0,.965),1,.035,transform=f.transFigure,color=NAVY))
    f.text(.075,.978,"PyYAeHMOP  /  VERSION 2",color="white",fontsize=9,fontweight="bold",va="center")
    f.text(.925,.978,section.upper(),color="white",fontsize=8,ha="right",va="center")
    f.text(.075,.928,title,color=NAVY,fontsize=22,fontweight="bold",va="top")
    if sub:f.text(.075,.884,sub,color=GRAY,fontsize=10,va="top")
    f.add_artist(Line2D([.075,.925],[.047,.047],transform=f.transFigure,color="#c5d2d9",lw=.7))
    f.text(.075,.029,"USER TUTORIAL  •  0.5.0.dev0  •  26 November 2026",fontsize=7,color=GRAY)
    f.text(.925,.029,f"{P:02d}",fontsize=8,color=GRAY,ha="right")
    return f
def para(f,y,text,width=94,size=10,indent=.075):
    for p in text.split("\n"):
        for line in textwrap.wrap(p,width=width,break_long_words=False,break_on_hyphens=False) or [""]:
            f.text(indent,y,line,fontsize=size,color=NAVY,va="top");y-=.0185 if size>=10 else .016
        y-=.010
    return y
def head(f,y,t):
    f.text(.075,y,t,fontsize=12,color=TEAL,fontweight="bold",va="top");return y-.032
def code(f,y,text,size=8.2):
    lines=text.strip("\n").splitlines()
    h=.019*len(lines)+.025
    f.patches.append(Rectangle((.075,y-h),.85,h,transform=f.transFigure,facecolor=LIGHT,edgecolor="#d4e0e6",lw=.5))
    for i,l in enumerate(lines):f.text(.09,y-.013-i*.019,l,fontfamily="DejaVu Sans Mono",fontsize=size,color=NAVY,va="top")
    return y-h-.022
def compact_code(f,y,text,size=8):
    lines=text.strip("\n").splitlines()
    h=.0158*len(lines)+.019
    f.patches.append(Rectangle((.075,y-h),.85,h,transform=f.transFigure,facecolor=LIGHT,edgecolor="#d4e0e6",lw=.5))
    for i,line in enumerate(lines):f.text(.09,y-.009-i*.0158,line,fontfamily="DejaVu Sans Mono",fontsize=size,color=NAVY,va="top")
    return y-h-.016
def save(pdf,f):pdf.savefig(f);plt.close(f)
def tblock(f,y,title,body):
    return para(f,head(f,y,title),body)
def plot_curves(ax,t,g,i,which):
    if which=="T":
        ax.plot(t[:,0],t[:,1],color=BLUE,lw=1.4,label="total" if t.shape[1]==4 else None)
        if t.shape[1]==4:
            ax.plot(t[:,0],t[:,2],color=TEAL,lw=1,label="up");ax.plot(t[:,0],t[:,3],color=ORANGE,lw=1,label="down");ax.legend(fontsize=6,frameon=False)
        ax.set(xlabel="E - Eref (eV)",ylabel="T(E)")
    elif which=="G":ax.plot(g[:,0],g[:,1]*1e6,color=TEAL,lw=1.4);ax.set(xlabel="E - Eref (eV)",ylabel="G (µS)")
    else:ax.plot(i[:,0],i[:,1]*1e6,color=ORANGE,lw=1.4);ax.set(xlabel="Bias (V)",ylabel="I (µA)")
    ax.grid(alpha=.18)
def spectra(name,old):
    if old:
        return [np.loadtxt(DEV/"benchmarks"/(name+"_"+s+".dat")) for s in ("transmission","conductance","current_voltage")]
    return [np.loadtxt(TUT/name/"reference"/(s+".dat")) for s in ("transmission","conductance","current_voltage")]
def case(pdf,name,title,desc,old=True,side=False):
    p=DEV/(name+".device.job") if old else TUT/name/"run.device.job"
    rel=p.relative_to(ROOT).as_posix();o=opts(p)
    f=page(pdf,"Device examples",title,"Reproduce this case directly from the source checkout.")
    y=para(f,.851,desc)
    y=code(f,y,"# XYZ and CIF alternatives (choose one)\npython -m pyaehmop "+rel+" --workers 4\npython -m pyaehmop "+rel.replace(".device.job","_cif.device.job")+" --workers 4",7.0)
    # fixed geometry and plots below intro
    ax=f.add_axes([.12,.515,.78,.155]);geo(ax,p,side)
    f.text(.075,.449,"Actual XYZ coordinates, one transverse repeat. Shaded regions: L / D / R.",fontsize=8,color=GRAY)
    t,g,i=spectra(name,old)
    for j,what in enumerate(("T","G","I")):plot_curves(f.add_axes([.12+j*.282,.262,.226,.151]),t,g,i,what)
    conv="explicit spinor total" if o.get("spin")=="true" or o.get("soc")=="true" else "per spin"
    f.text(.075,.212,f"T(E) convention: {conv}; G at 0 K; I(V) at {o.get('temperature','300')} K (frozen T).",fontsize=8,color=GRAY)
    y=.183
    counts=[len(xyz(q)[0]) for q in jobpaths(p)]
    y=para(f,y,f"Atoms L / D / R: {counts[0]} / {counts[1]} / {counts[2]}. Cutoff: {o.get('cutoff','11')} Å; transverse mesh: {o.get('transverse_kpoints','') or 'none (1D)'}; {o.get('energy_points','401')} energies. At +2 V: {i[-1,1]*1e6:.6g} µA. Diagnostic model only; not a converged prediction.",width=94,size=9)
    f.text(.075,.083,"Input paths and numerical reference files: "+str(p.parent.relative_to(ROOT)),fontsize=8,color=BLUE)
    save(pdf,f)

with PdfPages(out,metadata={"Title":"PyYAeHMOP Version 2: Reproducible Electronic Structure and Transport","Author":"PyYAeHMOP contributors","Subject":"Installation, input conventions, parameter fitting and device tutorials"}) as pdf:
    f=page(pdf,"Getting started","Electronic structure &\nquantum transport","")
    y=para(f,.827,"A practical guide to the independent Python extended Hückel implementation: molecules, crystals, parameter fitting, spin, spin–orbit coupling and two-terminal nonorthogonal NEGF.",width=86,size=11)
    f.text(.075,.705,"START WITH THE SOURCE CHECKOUT",fontsize=11,fontweight="bold",color=TEAL)
    f.text(.075,.674,URL,fontsize=10,color=BLUE,url=URL)
    f.text(.075,.636,"No access to a specific research cluster is required.",fontsize=11,color=NAVY)
    ax=f.add_axes([.12,.420,.78,.15]);geo(ax,DEV/"pristine.device.job")
    f.text(.075,.373,"What this manual covers",fontsize=12,fontweight="bold",color=TEAL)
    for j,l in enumerate(["02  Installation: pip, virtual environments and Conda","03  Electronic structure and magnetic calculations","04–05  Parameters and reproducible BandFit calibration","06–08  Geometry, device input and interpreting output","09–23  Fifteen individual device walkthroughs with T, G and I(V)","24–29  Benchmarks, electronic plots, contact setup and PBS","30–31  XYZ/CIF: direct uniform and explicit-device inputs","32–46  Appendix A: all 15 paired device inputs","47–52  Appendix B: all six BandFit pairs","53–67  Appendix C: all original electronic and magnetic inputs"]):
        f.text(.075,.34-j*.025,l,fontsize=10,color=NAVY)
    f.text(.075,.074,"Experimental code: numerical demonstrations are not automatically physical predictions.",fontsize=9,color=ORANGE)
    save(pdf,f)
    f=page(pdf,"Getting started","Install and run","Choose either a Python virtual environment or Conda.")
    y=head(f,.85,"A. pip and venv (recommended)")
    y=code(f,y,"git clone "+URL+".git\ncd PyYAeHMOP_version2\npython -m venv .venv\nsource .venv/bin/activate\npython -m pip install --upgrade pip\npython -m pip install -e \".[all]\"")
    y=para(f,y,"Windows PowerShell: use .venv\\Scripts\\Activate.ps1. Python 3.10+ is supported; use Python 3.11 for a straightforward setup. Git and an internet connection are needed to install dependencies.")
    y=head(f,y,"B. Conda (alternative)")
    y=code(f,y,"git clone "+URL+".git\ncd PyYAeHMOP_version2\nconda env create -f environment.yml\nconda activate pyaehmop")
    y=para(f,y,"Conda installs dependencies and then the checkout via pip. No separate PyPI or conda-forge release is claimed.")
    y=head(f,y,"Verify the installation")
    y=code(f,y,"python -c \"import pyaehmop; print(pyaehmop.__version__)\"\npython -m pyaehmop --help\npython -m pyaehmop pyaehmop/examples_friendly/water_full.job")
    para(f,y,"Run commands from the repository root. Use the full source checkout because the examples and reference datasets are not all packaged in a wheel. For an editable installation, source changes take effect immediately.")
    save(pdf,f)
    f=page(pdf,"Electronic structure","First electronic calculations","The original electronic-structure command-line workflows are retained.")
    y=head(f,.85,"Ordinary .job files")
    y=code(f,y,"python -m pyaehmop pyaehmop/examples_friendly/water_full.job\npython -m pyaehmop pyaehmop/examples_friendly/diamond_cif.job\npython -m pyaehmop pyaehmop/examples_friendly/xef4_original.job",7.6)
    y=para(f,y,"Water: molecular levels, HOMO/LUMO, occupations, energy, Mulliken charges and DOS. Diamond: bands, DOS and charges. XeF4: native YAeHMOP input, FMO and bond analysis. Other jobs demonstrate XYZ, CIF, POSCAR, ASE and external parameters. See the example index.")
    y=head(f,y,"Spin and SOC (electronic structure)")
    y=code(f,y,"python -m pyaehmop.SPIN_TEST pyaehmop/SPIN_TEST/examples/Fe_wire_spin.spin\npython -m pyaehmop.SPIN_TEST pyaehmop/SPIN_TEST/examples/Fe_wire_soc.spin\npython -m pyaehmop.SPIN_TEST pyaehmop/SPIN_TEST/examples/Fe_wire_both.spin",6.9)
    y=para(f,y,"These use the existing constrained collinear Stoner SCF / atomic L·S model. Transport spin uses frozen exchange shifts instead; the periodic SCF is not an open-device SCF.")
    f.text(.075,.329,"Archived numerical checks",fontsize=12,fontweight="bold",color=TEAL)
    d=json.loads((TUT/"reference/electronic_summary.json").read_text())
    vals=[("Water HOMO (eV)",d["friendly_water_full"]["homo_eV"]),("Diamond sampled band-path gap (eV)",d["friendly_diamond_cif"]["band_path_gap_eV"]),("Fe wire spin + SOC moment (µB)",d["Fe_wire_both"]["total_moment_muB"])]
    for j,(k,v) in enumerate(vals):
        f.text(.075,.286-j*.047,k,fontsize=10,color=NAVY);f.text(.925,.286-j*.047,f"{v:.6f}",fontsize=11,color=TEAL,ha="right")
    para(f,.120,"These are model-specific EHT values, not experimental observables guaranteed by the method. The archived scalar results and all input files are in examples/tutorials and pyaehmop/examples_friendly.",size=9)
    save(pdf,f)
    f=page(pdf,"Parameters","Choose an EHT parameter set","Three alternatives: bundled defaults, custom YAeHMOP table, or BandFit.")
    y=head(f,.85,"1. Bundled defaults")
    y=code(f,y,"parameters = default\nK = 1.75\nweighted_hij = true")
    y=para(f,y,"The native table pyaehmop/eht_parms.dat defines atomic Slater shells. K=1.75 is a customary starting value, not a universal fit.")
    y=head(f,y,"2. Supply a native table")
    y=code(f,y,"parameters = my_parameters.dat   # in the job\n# Contents of my_parameters.dat:\n# element Z valence n_zeta n shell Hii_eV zeta1 zeta2 c1 c2\nC 6 4 1 2 s -21.4 1.625 0 1 0\nC 6 4 1 2 p -11.4 1.625 0 1 0")
    y=para(f,y,"The two C rows are illustrative conventional values. The classic PARAMETERS record is also accepted: C 6 4 2 1.71 -21.4 2 1.625 -11.4. Use native 11-column rows for double-zeta contractions (both exponents and coefficients supplied). Energies are eV, exponents inverse Bohr; the contraction is normalized.")
    y=head(f,y,"3. Use a calibrated BandFit table")
    y=code(f,y,"parameters = ../transport/carbon_pdz7_4B.dat\nK = 2.0068527138311567\nweighted_hij = true\nenergy_reference_eV = -18.800896640571224")
    para(f,y,"Paths are relative to the job. External tables replace all shells for listed species; absent species retain bundled defaults. K is a separate job option. Do not combine independently fitted material offsets or K values without checking interface transferability. See docs/parameters.md.",size=9)
    save(pdf,f)
    f=page(pdf,"Parameters","Refit 2, 4 and 6 bands","A self-contained DFT frontier-band reference is shipped.")
    y=para(f,.85,"The original PBE/DZP SIESTA electrode reference has 338 k points. Every fifth point is held out. The 2-, 4-, and 6-band targets are frontier valence/conduction energies of the folded cell. This partial recalibration freezes p onsite and radial data; it optimizes s onsite, K and an alignment shift. No transmission is fitted.")
    y=code(f,y,"git clone https://github.com/tromer-unb/yaehmop-bandfit.git ../yaehmop-bandfit\ngit -C ../yaehmop-bandfit checkout c4a7fa0326b04409bd5ba2b0556985b4e5e0e428\npython -m examples.tutorials.bandfit.refit \\\n  --bandfit ../yaehmop-bandfit --workers 4",7.1)
    y=para(f,y,"Run from repository root. The included reference/frontier_bands.csv and initial JSONs are enough; no SIESTA executable is required. The BandFit optimizer uses differential evolution (seed 42, 180 evaluations), followed by least squares. All six fits reproduced exactly in the archived run.")
    f.text(.075,.49,"Validation RMSE (eV)",fontsize=12,fontweight="bold",color=TEAL)
    a=f.add_axes([.13,.263,.76,.19]);x=np.arange(3);a.bar(x-.16,[.1807107,.34115645,.48365584],.32,label="p-DZ7",color=BLUE);a.bar(x+.16,[.20616767,.28224699,.33004110],.32,label="sp5 (shared)",color=TEAL);a.set_xticks(x,["2 bands","4 bands","6 bands"]);a.set_ylabel("RMSE (eV)");a.legend(frameon=False);a.grid(axis="y",alpha=.2)
    y=para(f,.218,"The seven/five parameter names do not denote the number of bands. Each carbon still has four angular orbitals. Compare RMSE only with awareness of different target counts. For transport use the separately archived Dirac reference, not the fitted band-alignment shift.",size=9)
    para(f,y,"Quick job: python -m pyaehmop examples/tutorials/bandfit/pdz7_4B_smoke.transport.job --workers 4. For all six full and quick jobs, see examples/tutorials/bandfit/README.md.",size=9)
    save(pdf,f)
    f=page(pdf,"Transport geometry","How electrodes are defined","Actual pristine graphene positions, not a schematic device.")
    y=para(f,.85,"The rectangular cell contains four C atoms, with x period 2.459512 Å and z period 4.26 Å. The explicit device uses two cells per principal layer. Left and right have identical internal structure and parameters, but different ABSOLUTE positions.")
    ax=f.add_axes([.12,.565,.78,.17]);geo(ax,DEV/"pristine.device.job")
    f.text(.075,.505,"Same material; translated copies of the central eight-atom region.",fontsize=9,color=GRAY)
    y=code(f,.477,"electrode_left = left.xyz       # z = -8.52 ... -0.71 Å\nscattering_region = pristine.xyz # z = 0 ... 7.81 Å\nelectrode_right = right.xyz     # z = 8.52 ... 16.33 Å\ntransverse_vectors = 2.459512 0 0\nleft_repeat = 0 0 8.52\nright_repeat = 0 0 8.52")
    y=para(f,y,"Both repeat vectors point towards positive transport. Left is repeated toward negative multiples, right toward positive multiples. The files give the NEAREST layers, not the entire semi-infinite contacts. Never supply the same unshifted XYZ for both leads.")
    y=para(f,y,"The common x period repeats the sheet indefinitely; this example is infinite 2D graphene, not a nanoribbon. For x transport construct the layer along x and set the transverse repeat along z. For a ribbon, use finite width and omit transverse vectors.")
    para(f,y,"Increase lead thickness and the central buffer before increasing cutoff. Long-range couplings must not bypass the nearest principal layers. See docs/device-geometry.md.",size=9)
    save(pdf,f)
    f=page(pdf,"Transport input","Two distinct transport solvers","Uniform periodic materials versus explicit two-terminal devices.")
    y=head(f,.85,"Uniform .transport.job")
    y=code(f,y,"structure_format = xyz\nstructure = graphene.xyz\nlattice_vectors = 2.459512 0 0; 0 0 4.26\nimages = 5 4\ntransport_axis = 1\ntransverse_kpoints = 301")
    y=para(f,y,"The zero-based index 1 chooses the second lattice vector (z). The solver folds interactions into principal layers automatically. It describes pristine uniform materials, not an isolated central defect.")
    y=head(f,y,"Explicit .device.job (XYZ or CIF)")
    y=code(f,y,"electrode_left = left.xyz\nscattering_region = pristine.xyz\nelectrode_right = right.xyz\ntransverse_vectors = 2.459512 0 0\nleft_repeat = 0 0 8.52\nright_repeat = 0 0 8.52\ncutoff = 4.5\nquadrature = 16 12 8")
    y=para(f,y,"All XYZs must share a Cartesian origin and coordinate axes. CIF alternatives use a common P1 box and coordinate frame (see page 31). The nearest layers can have different atom counts, species and longitudinal periods. A common transverse lattice and a shared energy convention are required. For 1D leave transverse_kpoints blank.")
    para(f,y,"These excerpts are not complete inputs. Use the supplied runnable .job files, which also specify parameters, energy grid and output options. The device solver checks projected layer separation conservatively; enlarge layers if the check fails.",size=9)
    save(pdf,f)
    f=page(pdf,"Transport output","What the calculated curves mean","Coherent, elastic NEGF in a nonorthogonal basis.")
    y=code(f,.85,"energy_reference_eV = -18.800896640571224\nenergy_min = -2\nenergy_max = 2\nenergy_points = 81\neta = 0.00001\ncalculations = transmission, conductance, iv\ntemperature = 300\nbias_min = -2\nbias_max = 2\nbias_points = 21\nworkers = 4\nplots = true")
    y=para(f,y,"The energy axis is E minus the explicit reference. For another material determine its reference; the code does not find a common Fermi level. Retarded Green functions use zS-H and lead self-energies; the Caroli trace gives T(E,k), averaged over transverse k.")
    y=head(f,y,"Output files")
    y=para(f,y,"transmission.dat: energy (eV) and T; transmission_k.npz: k-resolved data; conductance.dat: G in siemens at 0 K; current_voltage.dat: V and I in amperes; summary.json: options and conventions. plots=true writes PNG files. In collinear spin mode, current_spin.dat separates spin channels.")
    f.text(.075,y-.019,r"$G(E)=g_s e^2 T(E)/h$",fontsize=14,color=NAVY)
    f.text(.075,y-.067,r"$I(V)=g_s e/h\,\int T(E,0)[f_L(E)-f_R(E)]\,dE$",fontsize=14,color=NAVY)
    para(f,y-.118,"Use g_s=2 for per-spin T and g_s=1 for explicit spinors. Symmetric bias, fixed equilibrium T and thermal tails (20 kBT). No charge or electrostatic self-consistency under bias. G(E) is not a finite-temperature differential conductance.",size=9)
    save(pdf,f)
    cases=[
      ("pristine","Pristine graphene — z","Identical semi-infinite graphene leads, with a pristine central layer. This is the explicit counterpart of uniform transport.",True,False),
      ("vacancy","Graphene with a vacancy","One carbon is removed only from the device. Leads stay pristine. The vacancy repeats across the transverse cell; the structure is not relaxed.",True,False),
      ("adsorbate","CO adsorbed on graphene","CO is added only to the central sheet. The side projection shows the molecule above the sheet; it is an illustrative, unrelaxed geometry.",True,True),
      ("asymmetric","Different lead basis sizes","The C–Mo–S chain tests rectangular couplings and different longitudinal periods. It is a numerical model, not a realistic chain prediction.",True,False),
      ("graphene_mos2","Lateral graphene–MoS2","Dissimilar semi-infinite leads share a commensurate transverse repeat. Mo/S use defaults; no interface-specific fit or alignment is assumed.",True,False),
      ("spin","Collinear spin transport","A frozen exchange shift is applied to a selected central carbon p shell. See magnetic.json. T_up, T_down and their currents are saved.",True,False),
      ("soc","Spin–orbit transport","SOC is a local atomic spinor L·S term, with a prescribed illustrative carbon p strength. T includes both spin components; no spin-flip decomposition.",True,False),
      ("spin_soc","Combined spin and SOC","The same vacancy geometry has frozen exchange and atomic SOC together. The solver uses complex spinor matrices and total T.",True,False),
      ("graphene_x","Pristine graphene — x","The transport direction is changed to the original x lattice vector; the transverse period is now z. Compare per unit transverse width, not raw currents.",False,False),
      ("nanoribbon","Finite-width nanoribbon","The strip has no transverse periodicity and unpassivated edges. No edge relaxation or magnetic SCF is included; the reference is illustrative.",False,False),
      ("wide_vacancy","Vacancy in a wider supercell","Four transverse repeats, with one missing central C. This reduces periodic-defect concentration; isolated-vacancy convergence is still required.",False,False),
      ("potential_barrier","Central potential barrier","Identical graphene contacts surround a device with a frozen +0.5 eV offset in potentials.json. The offset is not a Poisson solution.",False,False),
      ("spin_parallel","Parallel magnetic contacts","Prescribed + exchange shifts in both leads and the device. These potentials are numerical examples, not spontaneous graphene magnetism.",False,False),
      ("spin_antiparallel","Antiparallel magnetic contacts","The right-lead exchange sign is reversed; left and device are unchanged. Compare T_up/T_down and I_up/I_down with the parallel case.",False,False),
      ("graphene_mos2_graphene","Graphene–MoS2–graphene","A strained commensurate MoS2 channel between graphene leads. The geometry, default Mo/S parameters and energy alignment are not calibrated.",False,False)]
    for name,title,desc,old,side in cases:case(pdf,name,title,desc,old,side)
    f=page(pdf,"Validation","Graphene versus TranSIESTA","Archived matched transmission benchmark, not a finite-bias DFT simulation.")
    a=np.genfromtxt(TR/"multiband.csv",delimiter=",",names=True)
    ax=f.add_axes([.13,.55,.76,.26]);e=a["energy_eV"];ax.plot(e,a["DFT"],color=NAVY,lw=2,label="TranSIESTA PBE/DZP")
    for n,c in [("pdz7_2B",BLUE),("pdz7_4B",TEAL),("pdz7_6B",ORANGE)]:ax.plot(e,a[n],color=c,label=n)
    ax.axvspan(-1,1,color="#dbeee6",alpha=.6);ax.set(xlabel="E - reference (eV)",ylabel="T per spin");ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.2)
    ax=f.add_axes([.13,.32,.76,.16]);mask=abs(e)<=1;ax.plot(e[mask],a["DFT"][mask],color=NAVY,lw=2,label="DFT");ax.plot(e[mask],a["pdz7_4B"][mask],color=TEAL,ls="--",label="p-DZ7 4B");ax.set(xlabel="E - reference (eV)",ylabel="T per spin");ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.2)
    y=para(f,.265,"The p-DZ7 four-band fit gives MAE 0.003914 over ±1 eV (41 shared points over ±2 eV; 21 in the narrow window). The two-band model gives the lowest MAE over ±2 eV. Larger differences occur above +1.5 eV. The independent DFT and EHT energy references must still be rigorously aligned.",size=9)
    para(f,y,"All six 301-k / 401-energy EHT reruns reproduced archived EHT spectra to within 5e-11 at the 41 shared samples. This is a numerical reproducibility result, not proof of DFT accuracy. No T(E) values were used in the band fits.",size=9)
    save(pdf,f)
    f=page(pdf,"Validation","All six graphene transport fits","Full reruns, 301 transverse k points and 401 energies; frozen transmission.")
    full=np.loadtxt(BF/"reference/full_transmission.csv",delimiter=",")
    iv=np.loadtxt(BF/"reference/full_current_voltage.csv",delimiter=",")
    cond=np.loadtxt(BF/"reference/full_conductance.csv",delimiter=",")
    colors=[BLUE,TEAL,ORANGE,"#6e6a97","#a0556d","#72884c"]
    names=["p-DZ7 2B","p-DZ7 4B","p-DZ7 6B","sp5 2B","sp5 4B","sp5 6B"]
    ax=f.add_axes([.13,.57,.76,.24])
    ax.plot(e,a["DFT"],color=NAVY,lw=2,label="TranSIESTA")
    for j in range(6):ax.plot(full[:,0],full[:,j+1],color=colors[j],lw=1.1,label=names[j])
    ax.axvspan(-1,1,color="#dbeee6",alpha=.45);ax.set(xlabel="E - reference (eV)",ylabel="T per spin");ax.legend(ncol=2,frameon=False,fontsize=7);ax.grid(alpha=.2)
    ax=f.add_axes([.13,.30,.35,.19]);bx=f.add_axes([.57,.30,.32,.19])
    for j in range(6):
        ax.plot(cond[:,0],cond[:,j+1]*1e6,color=colors[j],lw=1.2)
        bx.plot(iv[:,0],iv[:,j+1]*1e6,color=colors[j],lw=1.2)
    ax.set(xlabel="E - reference (eV)",ylabel="G (µS)");bx.set(xlabel="Bias (V)",ylabel="I (µA)");ax.grid(alpha=.2);bx.grid(alpha=.2)
    y=para(f,.24,"The six native parameter files, full/quick jobs and the compact DFT reference are in examples/tutorials/bandfit. The complete EHT transport rerun matched the archived EHT curves; the I(V) curves are Landauer postprocessing at 300 K, not finite-bias DFT.",size=9)
    para(f,y,"The shaded region marks ±1 eV around the reported reference. Agreement is best there for the fitted four-band p-DZ7 set. Different reference alignments and transverse meshes remain a limitation of the DFT comparison.",size=9)
    save(pdf,f)
    f=page(pdf,"Electronic structure","Molecules and periodic bands","Numerical plots from the shipped water and diamond inputs.")
    ed=TUT/"reference/electronic"
    w=np.loadtxt(ed/"water_molecular_levels.csv",delimiter=",",skiprows=1)
    wd=np.loadtxt(ed/"water_dos.csv",delimiter=",",skiprows=1)
    ax=f.add_axes([.13,.61,.33,.22]);bx=f.add_axes([.56,.61,.33,.22])
    for j,energy,occ in w:
        ax.hlines(energy,j-.32,j+.32,color=BLUE if occ>0 else ORANGE,lw=2)
    ax.set(xlabel="Molecular orbital index",ylabel="Energy (eV)",title="Water orbital levels");ax.grid(axis="y",alpha=.2)
    bx.plot(wd[:,0],wd[:,1],color=TEAL);bx.set(xlabel="Energy (eV)",ylabel="DOS (states/eV)",title="Water DOS");bx.grid(alpha=.2)
    d=np.loadtxt(ed/"diamond_bands.csv",delimiter=",",skiprows=1);dos=np.loadtxt(ed/"diamond_dos.csv",delimiter=",",skiprows=1)
    ax=f.add_axes([.13,.29,.33,.23]);bx=f.add_axes([.56,.29,.33,.23])
    for b in np.unique(d[:,4]):
        q=d[d[:,4]==b];ax.plot(q[:,0],q[:,5],color=BLUE,lw=.8)
    ax.set(xlabel="Sampled k-path index",ylabel="Energy (eV)",title="Diamond bands");ax.grid(alpha=.2)
    bx.plot(dos[:,0],dos[:,1],color=TEAL);bx.set(xlabel="Energy (eV)",ylabel="DOS (states/eV)",title="Diamond DOS");bx.grid(alpha=.2)
    y=para(f,.23,"Run water_full.job and diamond_cif.job in pyaehmop/examples_friendly. The DOS depends on the chosen k sampling and broadening; the reported gap is on the sampled path, not necessarily a global Brillouin-zone gap.",size=9)
    para(f,y,"These numerical references are archived in examples/tutorials/reference/electronic. Full CSVs and optional plots are regenerated in the job output directories.",size=9)
    save(pdf,f)
    f=page(pdf,"Magnetism","Spin and spin–orbit coupling","Electronic SCF versus frozen magnetic transport.")
    ed=TUT/"reference/electronic"
    b=np.loadtxt(ed/"fe_both_bands_spin.csv",delimiter=",",skiprows=1)
    d=np.loadtxt(ed/"fe_both_spin_dos.csv",delimiter=",",skiprows=1)
    h=np.loadtxt(ed/"fe_both_scf_history.csv",delimiter=",",skiprows=1)
    ax=f.add_axes([.13,.60,.46,.23]);bx=f.add_axes([.68,.60,.21,.23])
    q=ax.scatter(b[:,0],b[:,5],c=b[:,6],vmin=-1,vmax=1,cmap="coolwarm",s=3)
    ax.set(xlabel="k-path index",ylabel="Energy (eV)",title="Fe wire: spin + SOC");f.colorbar(q,ax=ax,label="spin-z",fraction=.045,pad=.025)
    bx.plot(d[:,0],d[:,1],label="up",color=BLUE);bx.plot(d[:,0],d[:,2],label="down",color=ORANGE);bx.set(xlabel="Energy (eV)",ylabel="DOS");bx.legend(fontsize=7,frameon=False)
    ax=f.add_axes([.13,.375,.76,.15]);ax.semilogy(h[:,0],h[:,1],color=TEAL);ax.set(xlabel="SCF iteration",ylabel="Moment residual",title="Constrained Stoner convergence");ax.grid(alpha=.2)
    y=para(f,.32,"The Fe-wire .spin calculation solves a constrained collinear-z SCF with atomic SOC. In the device job, spin=true instead reads frozen exchange I*m/2 from magnetic_parameters. Index 2 means the third atom in that region's XYZ. SOC strengths are in eV, keyed by species and shell.",size=9)
    y=code(f,y,'{"device": {"exchange_eV": {"2": {"p": 0.2}},\n            "soc_eV": {"C": {"p": 0.02}}}}',8)
    para(f,y,"The example values are illustrative. Add left/right JSON sections to assign contact exchange and SOC; offsets may be used without spin. The code does not solve open-boundary magnetic SCF. No spin-flip-resolved transmission is reported with SOC.",size=9)
    save(pdf,f)
    f=page(pdf,"Transport geometry","From bulk to contacts","Build each principal layer in one global coordinate frame.")
    ax=f.add_axes([.12,.61,.78,.19])
    for label,path,shift,color in [("L−1",DEV/"left.xyz",-8.52,"#bdd4e3"),("L",DEV/"left.xyz",0,"#dce8f1"),("D",DEV/"pristine.xyz",0,"#d9eee6"),("R",DEV/"right.xyz",0,"#f6e6d6"),("R+1",DEV/"right.xyz",8.52,"#e9c6a4")]:
        sym,pos=xyz(path);z=pos[:,2]+shift
        ax.axvspan(z.min()-.355,z.max()+.355,color=color)
        ax.scatter(z,pos[:,0],s=20,color=NAVY,edgecolors="white",linewidths=.3)
        ax.text((z.min()+z.max())/2,1.04,label,transform=ax.get_xaxis_transform(),ha="center",fontsize=8)
    ax.annotate("repeat to −∞",xy=(-17,.4),xytext=(-12,.4),arrowprops=dict(arrowstyle="->",color=GRAY),ha="left",fontsize=8)
    ax.annotate("repeat to +∞",xy=(25,.4),xytext=(20,.4),arrowprops=dict(arrowstyle="->",color=GRAY),ha="right",fontsize=8)
    ax.set(xlabel="z (Å)",ylabel="x (Å)",ylim=(-.5,1.85))
    y=para(f,.55,"Only L, D and R are supplied. The repeated L−1 and R+1 layers are generated implicitly by the lead repeat vectors. D is finite along transport. In this pristine example all three supplied layers have the same local geometry. For a defect, change D only, preserving the contact buffer.")
    y=head(f,y,"Layer thickness and transverse periodicity")
    y=para(f,y,"The cutoff is a distance in Å, not a count of cells. A principal layer must be thick enough that no interactions reach beyond adjacent lead layers. Also prevent direct L-to-R coupling or coupling from D into deeper lead layers. If rejected, lengthen the contact layers and device; do not simply lower cutoff without convergence tests.")
    y=para(f,y,"For 2D graphene transverse_vectors = 2.459512 0 0 and transverse_kpoints = 15 in the diagnostic job. For a ribbon there is no transverse vector; for 3D use two independent transverse vectors and two mesh integers. In a heterogeneous contact, commensurate transverse periods and validated interface distances are required.")
    para(f,y,"A single vacancy in a periodic transverse cell is an array of defects, not an isolated defect. The wider example is a first step; check supercell size, lead buffers, k sampling, cutoff, eta and parameter transferability.",size=9)
    save(pdf,f)
    f=page(pdf,"Reproduction","Run locally or on your cluster","The examples need no access to any particular institution.")
    y=head(f,.85,"Workstation (default)")
    y=code(f,y,"export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\npython -m pyaehmop examples/devices/pristine.device.job --workers 4\npython -m pyaehmop examples/tutorials/nanoribbon/run.device.job --workers 4",7.2)
    y=para(f,y,"Set workers to the CPU count you can use. Avoid nested BLAS threads. A single quick case is the recommended first step; full runs require more resources.")
    y=head(f,y,"Optional PBS scheduler")
    y=code(f,y,"qsub -v SUITE=devices,PYTHON=\"$(command -v python)\" \\\n  examples/tutorials/run_suite.pbs\nqsub -v SUITE=full_bandfit,PYTHON=\"$(command -v python)\" \\\n  examples/tutorials/run_suite.pbs")
    y=para(f,y,"Other SUITE values: electronic, bandfit (quick). Submit from repository root. Adjust queue, account, walltime and memory to your scheduler. An absolute Python executable avoids site-specific paths. No GitHub workflow runs calculations.")
    y=head(f,y,"Validation and references")
    y=para(f,y,"15 device examples, 15 electronic examples, six quick/full transport models and six refits completed. Original device/core/spin tests: 81 passed, two skipped. Numerical validation is not DFT validation of defects or interfaces.",size=9)
    y=para(f,y,"Troubleshooting: check geometry and cutoff if layers are rejected; check overlap positivity, quadrature, k mesh and eta. For I(V) extend the energy window to include thermal tails. For SOC use a justified xi; for frozen spin specify exchange_eV by zero-based atom index.",size=9)
    para(f,y,"References: Hoffmann, JCP 39, 1397 (1963), doi:10.1063/1.1734456; López Sancho et al., J. Phys. F 15, 851 (1985); Caroli et al., J. Phys. C 4, 916 (1971); Brandbyge et al., PRB 65, 165401 (2002); Papior et al., CPC 212, 8 (2017). See examples/transport/README.md for links.",size=8)
    save(pdf,f)
    f=page(pdf,"CIF and XYZ","Uniform transport from CIF","Same physical graphene model, different lattice representation.")
    y=para(f,.85,"The XYZ example describes a true 2D sheet in x-z. The CIF has three lattice vectors a=x, b=y (30 angstrom vacuum), c=z. Its vacuum exceeds the 11 angstrom cutoff, so periodic images across y do not interact. ASE reads the cell from the CIF.")
    y=code(f,y,"# XYZ                              # CIF\nstructure_format = xyz             structure_format = cif\nstructure = graphene.xyz           structure = graphene.cif\nlattice_vectors = 2.459512 0 0; 0 0 4.26  # omit for CIF\nimages = 5 4                       images = 5 1 4\ntransport_axis = 1                 transport_axis = 2\ntransverse_kpoints = 301            transverse_kpoints = 301 1",7.8)
    y=para(f,y,"The index is a zero-based lattice-vector index, not a Cartesian component. The second XYZ vector is z; the third CIF vector is z. All remaining settings (including carbon parameters, K and energy reference) stay the same.")
    y=code(f,y,"python -m pyaehmop examples/transport/graphene_cif.transport.job --workers 4",7.0)
    f.text(.075,.425,"Numerical XYZ / CIF equivalence (full p-DZ7 4B)",fontsize=12,fontweight="bold",color=TEAL)
    data=np.loadtxt(TR/"cif_reference/graphene_xyz_cif_full.dat")
    ax=f.add_axes([.13,.225,.76,.16])
    ax.plot(data[:,0],data[:,1],color=BLUE,lw=2,label="XYZ")
    ax.plot(data[:,0],data[:,2],color=ORANGE,ls="--",lw=1.1,label="CIF")
    ax.set(xlabel="E - reference (eV)",ylabel="T per spin");ax.legend(frameon=False);ax.grid(alpha=.2)
    checks=json.loads((TR/"cif_reference/validation.json").read_text())
    m=max(x["transmission.dat"] for x in checks["uniform"].values())
    fullerr=json.loads((TR/"cif_reference/full_validation.json").read_text())["transmission.dat"]
    y=para(f,.171,f"Full 4B (301 k, 401 energies): maximum absolute T difference {fullerr:.3g}. All six quick pairs: {m:.3g}. The lines overlap. Every model has paired full/quick CIF and XYZ jobs.",size=9)
    para(f,y,"For a quick run use examples/tutorials/bandfit/pdz7_4B_smoke_cif.transport.job. See docs/transport-cif.md for all paths, direction changes and vacuum convergence cautions.",size=9)
    save(pdf,f)
    f=page(pdf,"CIF and XYZ","Independent devices from CIF","Run either input directly; no separate conversion command.")
    y=para(f,.85,"All 15 devices have paired .device.job files. XYZ gives Cartesian coordinates; the CIF variant reads three ordered P1 region files sharing a padded orthogonal box and coordinate frame. The physical repeats remain in the job.")
    y=code(f,y,"python -m pyaehmop examples/devices/pristine.device.job --workers 4\npython -m pyaehmop examples/devices/pristine_cif.device.job --workers 4\npython -m pyaehmop examples/devices/graphene_mos2.device.job --workers 4\npython -m pyaehmop examples/devices/graphene_mos2_cif.device.job --workers 4",7.2)
    y=para(f,y,"The reader converts fractional positions to Cartesian coordinates in the common box, preserving atom order. A global shift does not change transport. No origin JSON is needed. The input contains structure_format=cif and three CIF paths. All physical settings stay the same. See Appendix A for all 15 complete input pairs.")
    ax=f.add_axes([.12,.41,.78,.115]);geo(ax,DEV/"graphene_mos2.device.job")
    f.text(.075,.368,"Actual graphene-MoS2 coordinates. Both input formats preserve these regions.",fontsize=8,color=GRAY)
    checks=json.loads((TR/"cif_reference/direct_device_validation.json").read_text())
    m=max(x["transmission.dat"] for x in checks.values())
    c=max(x["max_coordinate_error_A"] for x in checks.values())
    y=para(f,.331,f"Direct CIF validation: {len(checks)} devices; maximum coordinate difference after common translation {c:.3g} angstrom; maximum absolute T difference {m:.3g}. This establishes format equivalence, not physical calibration.",size=9)
    y=para(f,y,"For new devices, read real bulk CIFs with ASE, expand symmetry, build commensurate principal layers and place L/D/R in one frame. Export XYZ or matching ordered P1 region CIFs in one common box. Never infer lead repeats from the padded interchange box. No automatic interface assembly is performed.",size=9)
    para(f,y,"The optional converter remains available for exporting XYZ. For exact commands and geometry paths see the catalog at examples/README.md and docs/transport-cif.md.",size=9)
    save(pdf,f)
    # Appendix A: retain every original plotted case and provide its full job.
    for name,title,desc,old,side in cases:
        p=DEV/(name+".device.job") if old else TUT/name/"run.device.job"
        q=p.with_name(p.name.replace(".device.job","_cif.device.job"))
        x,c=opts(p),opts(q)
        f=page(pdf,"Appendix A / devices",title,"Complete paired inputs. Choose either format; same Hamiltonian and geometry.")
        y=para(f,.85,desc,width=96,size=9)
        y=compact_code(f,y,"# From repository root\npython -m pyaehmop "+p.relative_to(ROOT).as_posix()+" --workers 4\npython -m pyaehmop "+q.relative_to(ROOT).as_posix()+" --workers 4",7.0)
        y=head(f,y,"Geometry block: XYZ input")
        y=compact_code(f,y,"structure_format = xyz  # default\n"+"\n".join(k+" = "+x[k] for k in ("electrode_left","scattering_region","electrode_right")),7.6)
        y=head(f,y,"Geometry block: CIF input")
        y=compact_code(f,y,"\n".join(k+" = "+c[k] for k in ("structure_format","electrode_left","scattering_region","electrode_right")),7.2)
        y=head(f,y,"Shared job settings (append to either block)")
        keys=[k for k in x if k not in ("structure_format","cif_origin","electrode_left","scattering_region","electrode_right","output")]
        # Two columns, actual values from shipped file, with an explicit output row.
        lines=[k+" = "+x[k] for k in keys]
        split=(len(lines)+1)//2
        left,right=lines[:split],lines[split:]
        h=.0155*split+.022
        f.patches.append(Rectangle((.075,y-h),.85,h,transform=f.transFigure,facecolor=LIGHT,edgecolor="#d4e0e6",lw=.5))
        for j,col in enumerate((left,right)):
            for i,line in enumerate(col):f.text(.09+j*.425,y-.009-i*.0155,line,fontfamily="DejaVu Sans Mono",fontsize=6.6,color=NAVY,va="top")
        y-=h+.018
        y=para(f,y,"For XYZ add output = "+x["output"]+"; for CIF add output = "+c["output"]+". Paths are relative to the job directory. The CIF box is coordinate interchange only; the listed repeat/transverse vectors determine the physical periodicity.",width=105,size=8)
        data=json.loads((TR/"cif_reference/direct_device_validation.json").read_text())
        err=data[p.relative_to(ROOT).as_posix()]["transmission.dat"]
        para(f,y,f"Direct CIF versus XYZ max absolute transmission difference: {err:.3g}. For plots and physical limitations see the corresponding Device examples page (09–23).",width=105,size=8)
        save(pdf,f)
    # Appendix B: all six complete uniform transport pairs.
    for p in sorted(BF.glob("*B.transport.job")):
        name=p.name.removesuffix(".transport.job")
        q=p.with_name(name+"_cif.transport.job")
        x,c=opts(p),opts(q)
        f=page(pdf,"Appendix B / BandFit",name+" transport","Full 301-k / 401-energy jobs, plus quick alternatives.")
        y=compact_code(f,.85,"python -m pyaehmop "+p.relative_to(ROOT).as_posix()+" --workers 4\npython -m pyaehmop "+q.relative_to(ROOT).as_posix()+" --workers 4",7.0)
        y=head(f,y,"XYZ geometry and sampling")
        y=compact_code(f,y,"\n".join(k+" = "+x[k] for k in ("structure_format","structure","lattice_vectors","images","transport_axis","transverse_kpoints")),7.5)
        y=head(f,y,"CIF geometry and sampling")
        y=compact_code(f,y,"\n".join(k+" = "+c[k] for k in ("structure_format","structure","images","transport_axis","transverse_kpoints")),7.5)
        y=head(f,y,"Shared settings")
        keys=[k for k in x if k not in ("structure_format","structure","lattice_vectors","images","transport_axis","transverse_kpoints","output")]
        lines=[k+" = "+x[k] for k in keys];split=(len(lines)+1)//2
        h=.016*split+.022
        f.patches.append(Rectangle((.075,y-h),.85,h,transform=f.transFigure,facecolor=LIGHT,edgecolor="#d4e0e6",lw=.5))
        for j,col in enumerate((lines[:split],lines[split:])):
            for i,line in enumerate(col):f.text(.09+j*.425,y-.009-i*.016,line,fontfamily="DejaVu Sans Mono",fontsize=6.8,color=NAVY,va="top")
        y-=h+.015
        y=para(f,y,"Add output = "+x["output"]+" (XYZ) or output = "+c["output"]+" (CIF). Quick jobs: "+name+"_smoke.transport.job and "+name+"_smoke_cif.transport.job (15 k, 81 energies). The fitted parameter file and K are shared.",size=8,width=105)
        y=compact_code(f,y,(BF/x["parameters"]).read_text().strip(),6.9)
        para(f,y,"Native EHT table above; Hii in eV and Slater exponents in inverse Bohr. A 30 Å CIF vacuum removes coupling across y for the 11 Å cutoff. The six quick pairs and full 4B CIF were validated; other five full CIF jobs are supplied but not rerun.",size=8,width=105)
        save(pdf,f)
    # Appendix C: every previously distributed friendly and magnetic input.
    explanations={
        "water_full":"Water: Cartesian XYZ, default EHT parameters, molecular levels, charges and DOS. The original DOS and HOMO/LUMO plots appear earlier.",
        "water_charges":"Molecular water with only Mulliken charge calculation.",
        "water_default_parameters":"Water with explicit selection of the bundled YAeHMOP parameter table.",
        "water_custom_parameters":"Water with a user-provided native parameter file; the remaining options define the molecular calculation.",
        "diamond_cif":"Primitive diamond: CIF stores cell and fractional sites. Custom carbon parameters, sampled band path, DOS and Mulliken charges.",
        "diamond_ase":"Diamond in ASE extended XYZ. Unlike ordinary XYZ, extended XYZ includes lattice and periodic boundary flags.",
        "diamond_poscar":"Diamond using POSCAR through the structure reader; the cell and atomic sites come from the file.",
        "diamond_reference":"Native YAeHMOP diamond input; reproduces the legacy electronic-structure example.",
        "xef4_original":"Native molecular YAeHMOP example including fragment molecular orbital / bond analysis.",
        "Fe_wire_spin":"Periodic Fe wire with constrained collinear Stoner self-consistency, without SOC.",
        "Fe_wire_soc":"Periodic Fe wire with atomic spin-orbit coupling only.",
        "Fe_wire_both":"Periodic Fe wire with Stoner exchange and atomic SOC; spin DOS and SCF history plotted earlier.",
        "Fe_wire_external_parameters":"Fe wire with external EHT parameters and its prescribed spin/SOC settings.",
        "Fe_wire_soc_exaggerated":"Artificially large SOC to demonstrate the numerical effect; not a physical Fe value.",
        "C_atom_both":"Isolated carbon atom with spin and atomic SOC; electronic, not transport, calculation."
    }
    for p in [*sorted((ROOT/"pyaehmop/examples_friendly").glob("*.job")),*sorted((ROOT/"pyaehmop/SPIN_TEST/examples").glob("*.spin"))]:
        f=page(pdf,"Appendix C / electronic",p.stem.replace("_"," "),"Original molecular, crystal or magnetic example — complete input.")
        y=para(f,.85,explanations.get(p.stem,"Original electronic structure example."),width=94,size=10)
        mod="pyaehmop.SPIN_TEST" if p.suffix==".spin" else "pyaehmop"
        y=compact_code(f,y,"python -m "+mod+" "+p.relative_to(ROOT).as_posix(),7.0)
        y=head(f,y,"Input file (comments and blank lines omitted)")
        lines=[l.strip() for l in p.read_text().splitlines() if l.strip() and not l.lstrip().startswith("#")]
        y=compact_code(f,y,"\n".join(lines),7.1)
        y=para(f,y,"Run from the checkout root. Structure and parameter paths are relative to the input directory. Read the output directory specified in the job; compare with examples/tutorials/reference/electronic_summary.json. See the plots and discussion in the electronic-structure and magnetism chapters.",size=9)
        save(pdf,f)
print(out,P,"pages",out.stat().st_size,"bytes")
