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
from lead_repeat_chapter import CHAPTER_PAGES, add_chapter

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

def repo_path(p):
    return Path(p).resolve().relative_to(ROOT).as_posix()
def linkline(f,y,p,label=None,size=8.5,x=.075):
    path=Path(p).resolve()
    assert path.exists(),path
    rel=repo_path(path);kind="tree" if path.is_dir() else "blob"
    f.text(x,y,label or rel,fontsize=size,color=BLUE,va="top",url=URL+"/"+kind+"/main/"+rel)
    return y-.021
def column_code(f,x,y,w,text,size=7.1,step=.0145):
    lines=text.strip("\n").splitlines()
    h=step*len(lines)+.018
    f.patches.append(Rectangle((x,y-h),w,h,transform=f.transFigure,facecolor=LIGHT,edgecolor="#d4e0e6",lw=.5))
    for i,l in enumerate(lines):
        f.text(x+.012,y-.009-i*step,l,fontfamily="DejaVu Sans Mono",fontsize=size,color=NAVY,va="top")
    return y-h-.014
def cif_content(p):
    # The shipped region CIFs are ordered P1 with one explicit atom-site loop.
    lines=Path(p).read_text().splitlines()
    tags={l.split()[0]:l.split()[1] for l in lines if l.startswith(("_cell_length_","_cell_angle_"))}
    rows=[l.split() for l in lines if len(l.split())==6 and l.split()[0][0].isalpha() and l.split()[0][-1].isdigit()]
    return tags,rows
def xyz_display(p):
    names,coords=xyz(p)
    return str(len(names))+"\n"+p.read_text().splitlines()[1]+"\n"+"\n".join(f"{a:2s} {x:9.6f} {y:9.6f} {z:10.6f}" for a,(x,y,z) in zip(names,coords))
def cif_display(p):
    tags,rows=cif_content(p)
    return "\n".join(f"{r[0]:4s} {r[1]:2s} "+ " ".join(f"{float(v):10.8f}" for v in r[2:5])+" "+r[5] for r in rows)
def job_folder(f,y,p,command=True):
    y=linkline(f,y,p.parent,"Folder: "+repo_path(p.parent))
    if command:
        y=column_code(f,.075,y,.85,"cd "+repo_path(p.parent),7.4)
    return y

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
    f=page(pdf,"Device results",title,"Archived T(E), G(E) and I(V) for the preceding XYZ and CIF jobs.")
    y=para(f,.851,desc)
    y=linkline(f,y,p.parent,"Folder: "+repo_path(p.parent),8.7)
    y=linkline(f,y,p,"XYZ job: "+p.name,8.1)
    y=linkline(f,y,p.with_name(p.name.replace(".device.job","_cif.device.job")),"CIF job: "+p.name.replace(".device.job","_cif.device.job"),8.1)
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
    err=json.loads((TR/"cif_reference/direct_device_validation.json").read_text())[rel]["transmission.dat"]
    f.text(.075,.105,f"XYZ/CIF max absolute T difference: {err:.3g}",fontsize=8,color=GRAY)
    linkline(f,.083,DEV/"benchmarks" if old else p.parent/"reference","Reference data: "+repo_path(DEV/"benchmarks" if old else p.parent/"reference"),8.0)
    save(pdf,f)

# Each case is a continuous chapter: XYZ, CIF, and the archived result.
# Large geometries continue immediately within their format, never in a distant atlas.
def device_groups(p):
    paths=jobpaths(p)
    counts=[len(xyz(q)[0]) for q in paths]
    groups=[]; pending=[]; total=0
    for label,path,n in zip(("L — Left electrode","D — Scattering region","R — Right electrode"),paths,counts):
        # Keep small regions together; use balanced chunks for larger regions.
        nchunk=(n+27)//28
        for j in range(nchunk):
            lo=j*n//nchunk;hi=(j+1)*n//nchunk
            if pending and (total+hi-lo>28 or len(pending)==3):
                groups.append(pending);pending=[];total=0
            pending.append((label,path,lo,hi,n));total+=hi-lo
    if pending:groups.append(pending)
    return groups

def source_code(f,x,y,w,text,size=7.8,step=.0135):
    # Bound monospace lines inside the available column, retaining exact text.
    lines=text.strip("\n").splitlines()
    longest=max(map(len,lines),default=1)
    size=min(size,(w-.026)*8.27*72/(.605*max(longest,1)))
    return column_code(f,x,y,w,text,size,step)

def settings_text(p):
    return "\n".join(l.strip() for l in p.read_text().splitlines()
                     if l.strip() and not l.lstrip().startswith("#"))

def heading_at(f,x,y,text,size=10):
    f.text(x,y,text,fontsize=size,color=TEAL,fontweight="bold",va="top")
    return y-.024

def path_ref(f,y,p,base,x=.075,size=8.1):
    return linkline(f,y,p,str(Path(p).relative_to(base)),size,x)

def data_refs(f,y,p,keys=("parameters","magnetic_parameters")):
    o=opts(p)
    for k in keys:
        if k in o and o[k]!="default":
            path=(p.parent/o[k]).resolve()
            if path.is_file():
                y=linkline(f,y,path,k+": "+repo_path(path),7.1)
    return y

def device_format(pdf,p,q,title,desc,kind,number):
    job=p if kind=="XYZ" else q
    o=opts(job); paths=jobpaths(job); groups=device_groups(p)
    if kind=="CIF":
        tags,_=cif_content(paths[0])
        for path in paths:
            t,_=cif_content(path)
            assert tags==t, path
        box=", ".join(tags["_cell_length_"+a] for a in "abc")
    for j,group in enumerate(groups):
        f=page(pdf,f"Device {number:02d} / {kind}",title+" — "+kind,
               "Full input, local command and corresponding structure files." if j==0 else f"Geometry continuation {j+1}/{len(groups)} — same {kind} input.")
        y=job_folder(f,.85,job,False)
        y=linkline(f,y,job,"Input: "+job.name,9.0)
        if j==0:
            y=source_code(f,.075,y,.85,"cd "+repo_path(job.parent)+"\npython -m pyaehmop "+job.name+" --workers 4",8.0,.014)
            y=para(f,y,desc,width=104,size=8.7)
        else:
            y=para(f,y,"Continue reading the geometry of this same job. All coordinates remain in their original file order; no additional calculation is needed.",width=104,size=9)
        top=y
        if j==0:
            left=heading_at(f,.075,top,"Complete job (without comments)")
            left=source_code(f,.075,left,.415,settings_text(job),7.65,.0129)
            if kind=="XYZ":
                left=para(f,left,"XYZ: atom count, comment, element and x y z in Å. All three regions use one Cartesian frame.",width=49,size=8)
                if number==1:
                    left=para(f,left,"Both repeats point along +z. Left repeats toward negative multiples; right toward positive multiples. Each nearest layer has 8 C (two cells). The transverse x repeat makes this an infinite sheet.",width=49,size=8)
            else:
                left=para(f,left,"CIF: common ordered P1 box; a,b,c = "+box+" Å, all angles 90°. Physical repeats are defined only in the job. Do not independently center L, D or R; no origin JSON is needed.",width=49,size=8)
            if left>.11:
                # File paths are shown in full in the job; explicit local output
                # path makes it easy to find the results after changing directory.
                left=para(f,left,"Output folder: "+o["output"],width=49,size=8)
            if left<.062:raise RuntimeError(("job overflow",job,left))
        else:
            left=heading_at(f,.075,top,"Same input and physical geometry")
            left=para(f,left,"This page continues "+job.name+". No changes to the input, parameters, energy grid or physical electrode repeat vectors.",width=49,size=9)
            if kind=="CIF":
                left=para(f,left,"Common CIF cell (Å): "+box+". P1, angles 90°. The column labels fx, fy, fz are the fractional atomic coordinates; occupancy is 1. The physical contact repeats are the job vectors.",width=49,size=9)
            else:
                left=para(f,left,"The first line of each XYZ file gives the number of atoms. All x, y, z values are in Å. Do not translate the regions separately.",width=49,size=9)
            left=linkline(f,left,job,job.name,8.4)
        right=heading_at(f,.52,top,"Structure files — "+kind)
        for label,xp,lo,hi,n in group:
            # Associate the original XYZ region with the matching CIF path.
            idx=jobpaths(p).index(xp);path=paths[idx]
            right=heading_at(f,.52,right,label+f" ({lo+1}–{hi} of {n})",9)+.006
            right=path_ref(f,right,path,job.parent,.52,7.6)+.005
            if kind=="XYZ":
                sym,coords=xyz(path)
                lines=xyz_display(path).splitlines()
                if lo==0 and hi==n:
                    text="\n".join(lines)
                else:
                    text=f"# {n} atoms; rows {lo+1}–{hi}\n"+"\n".join(lines[2+lo:2+hi])
                right=source_code(f,.52,right,.405,text,7.3,.0103)+.006
            else:
                rows=cif_display(path).splitlines()
                text="label el     fx         fy         fz    occ\n"+"\n".join(rows[lo:hi])
                right=source_code(f,.52,right,.405,text,7.1,.0103)+.006
        right=para(f,right,"Actual sites, rounded for display. Click the blue filenames for the exact file contents.",width=50,size=7.4,indent=.52)
        if right<.062:raise RuntimeError(("geometry overflow",job,j,right))
        if j==0 and min(left,right)>.13:
            yy=min(left,right)-.01
            yy=data_refs(f,yy,job)
        save(pdf,f)

def cif_full_text(path):
    # Preserve all CIF tags and sites, with eight decimal places for narrow
    # printed columns. The linked source has its original precision.
    lines=path.read_text().splitlines()
    out=[]
    for line in lines:
        w=line.split()
        if len(w)==6 and w[0][0].isalpha() and w[0][-1].isdigit():
            line=f"{w[0]:4s} {w[1]:2s} "+" ".join(f"{float(v):.8f}" for v in w[2:5])+" "+w[5]
        out.append(line)
    return "\n".join(out)

def cif_page_count(q):
    paths=jobpaths(q)
    # The first page shows the full L file next to the job. D and R
    # are printed together if both fit; otherwise each gets a page.
    n=[len(cif_content(p)[1]) for p in paths]
    return 2 if max(n[1:])<=34 else 3

def device_cif_format(pdf,p,q,title,desc,number):
    o=opts(q);paths=jobpaths(q)
    tags,_=cif_content(paths[0])
    for path in paths:
        t,_=cif_content(path);assert tags==t,path
    box=", ".join(tags["_cell_length_"+a] for a in "abc")
    f=page(pdf,f"Device {number:02d} / CIF",title+" — CIF",
           "Same device as the preceding XYZ example; full CIF input and structures.")
    y=job_folder(f,.85,q,False)
    y=linkline(f,y,q,"Input: "+q.name,9.0)
    y=source_code(f,.075,y,.85,"cd "+repo_path(q.parent)+"\npython -m pyaehmop "+q.name+" --workers 4",8.0,.014)
    y=para(f,y,desc,width=104,size=8.7)
    top=y
    left=heading_at(f,.075,top,"Complete job (without comments)")
    left=source_code(f,.075,left,.415,settings_text(q),7.65,.0129)
    left=para(f,left,"CIF: common ordered P1 box; a,b,c = "+box+" Å, all angles 90°. Physical repeats are defined only in the job. Do not independently center L, D or R. No origin JSON is needed.",width=49,size=8)
    left=para(f,left,"Output: "+o["output"],width=49,size=8)
    left=heading_at(f,.075,left,"Other structure files (next page"+("s)" if cif_page_count(q)==3 else ")"),9)
    left=path_ref(f,left,paths[1],q.parent,.075,7.2)
    left=path_ref(f,left,paths[2],q.parent,.075,7.2)
    right=heading_at(f,.52,top,"Left electrode — complete CIF",10)
    right=path_ref(f,right,paths[0],q.parent,.52,7.3)
    right=source_code(f,.52,right,.405,cif_full_text(paths[0]),7.1,.0103)
    right=para(f,right,"All tags and atom sites shown. Fractions rounded for print; click the blue file link for original precision.",width=49,size=7.4,indent=.52)
    if min(left,right)<.063:raise RuntimeError(("CIF first page overflow",q,left,right))
    if min(left,right)>.13:data_refs(f,min(left,right)-.01,q)
    save(pdf,f)
    groups=[[(1,.075,.415),(2,.52,.405)]] if cif_page_count(q)==2 else [[(1,.075,.85)],[(2,.075,.85)]]
    for j,group in enumerate(groups,1):
        f=page(pdf,f"Device {number:02d} / CIF",title+" — CIF",
               "Complete CIF geometry files (continued); same job, no additional command.")
        y=job_folder(f,.85,q,False)
        y=linkline(f,y,q,"Input: "+q.name,8.5)
        y=para(f,y,"The physical repeats are defined in the job; all CIF files share one coordinate frame.",width=110,size=8)
        ys=[]
        for i,x,w in group:
            yy=heading_at(f,x,y,("Scattering region" if i==1 else "Right electrode")+" — complete CIF",10)
            yy=path_ref(f,yy,paths[i],q.parent,x,7.5)
            # For a large file, use the full width and the exact original text.
            body=paths[i].read_text() if w>.5 else cif_full_text(paths[i])
            step=.0096 if len(body.splitlines())>55 else .0126
            yy=source_code(f,x,yy,w,body,7.7 if w>.5 else 7.1,step)
            ys.append(yy)
        end=min(ys)
        if end>.10:
            para(f,end,"The linked CIF files retain full precision. Any rounded decimal rows above are for printed illustration only.",size=8,width=104)
        if end<.062:raise RuntimeError(("CIF structure overflow",q,j,end))
        save(pdf,f)


def uniform_format(pdf,p,kind,title,quick=None):
    o=opts(p);g=(p.parent/o["structure"]).resolve()
    f=page(pdf,f"Periodic transport / {kind}",title+" — "+kind,
           "Complete job and its actual geometry; use the same installed CLI.")
    y=job_folder(f,.85,p,False)
    y=linkline(f,y,p,"Input: "+p.name,9.0)
    y=source_code(f,.075,y,.85,"cd "+repo_path(p.parent)+"\npython -m pyaehmop "+p.name+" --workers 4",8.0,.014)
    if quick:
        y=linkline(f,y,quick,"Quick alternative (15 k, 81 energies): "+quick.name,8.1)
    top=y
    left=heading_at(f,.075,top,"Complete job (without comments)")
    left=source_code(f,.075,left,.415,settings_text(p),7.65,.0134)
    if kind=="XYZ":
        left=para(f,left,"Two periodic vectors in x and z. transport_axis=1 selects the second lattice vector (z). The cell and periodic images are given in the job.",width=49,size=8)
    else:
        left=para(f,left,"Three cell vectors; y has 30 Å vacuum, larger than cutoff 11 Å. transport_axis=2 selects c (z); one transverse k point along the vacuum is sufficient in this model.",width=49,size=8)
    left=para(f,left,"Output folder: "+o["output"],width=49,size=8)
    right=heading_at(f,.52,top,"Structure file — "+kind)
    right=linkline(f,right,g,g.name,8.5,.52)
    # The tiny uniform cells fit in full, including all CIF tags and atom sites.
    right=source_code(f,.52,right,.405,g.read_text(),7.3,.0123)
    right=para(f,right,"Geometry location: "+repo_path(g),width=48,size=8,indent=.52)
    right=heading_at(f,.52,right,"EHT parameter file")
    par=(p.parent/o["parameters"]).resolve()
    right=linkline(f,right,par,par.name,8.1,.52)
    # Preserve all 11 values, split long native records into two display rows.
    lines=[]
    for line in par.read_text().splitlines():
        if len(line)>47:
            w=line.split()
            lines.append(" ".join(w[:6]))
            lines.append("  "+" ".join(w[6:]))
        else:lines.append(line)
    right=source_code(f,.52,right,.405,"\n".join(lines),6.8,.012)
    right=para(f,right,"Native table: element, Z, valence, n_zeta, n, shell, Hii, zeta1, zeta2, c1, c2. Long records wrapped for display only. Use the exact linked table with the job's K.",width=49,size=8,indent=.52)
    if quick:
        right=para(f,right,"Quick run from this same folder: python -m pyaehmop "+quick.name,width=48,size=8,indent=.52)
    if min(left,right)<.062:raise RuntimeError(("uniform overflow",p,left,right))
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

device_index=[]
next_page=6+CHAPTER_PAGES
for n,(name,title,desc,old,side) in enumerate(cases,1):
    p=DEV/(name+".device.job") if old else TUT/name/"run.device.job"
    k=len(device_groups(p))
    q=p.with_name(p.name.replace('.device.job','_cif.device.job'))
    nc=cif_page_count(q)
    device_index.append((name,p,next_page,next_page+k,next_page+k+nc))
    next_page+=k+nc+1
periodic_start=next_page
periodic_end=periodic_start+14
guide_start=periodic_end+1
electronic_start=guide_start+8
last_page=electronic_start+18
with PdfPages(out,metadata={"Title":"PyYAeHMOP Version 2: Reproducible Electronic Structure and Transport","Author":"PyYAeHMOP contributors","Subject":"Installation, input conventions, parameter fitting and device tutorials"}) as pdf:
    f=page(pdf,"Getting started","Electronic structure &\nquantum transport","")
    y=para(f,.827,"A practical guide to the independent Python extended Hückel implementation: molecules, crystals, parameter fitting, spin, spin–orbit coupling and two-terminal nonorthogonal NEGF.",width=86,size=11)
    f.text(.075,.705,"START WITH THE SOURCE CHECKOUT",fontsize=11,fontweight="bold",color=TEAL)
    f.text(.075,.674,URL,fontsize=10,color=BLUE,url=URL)
    f.text(.075,.636,"No access to a specific research cluster is required.",fontsize=11,color=NAVY)
    linkline(f,.612,ROOT/"examples","Open examples/ — complete clickable example catalog",9.5)
    ax=f.add_axes([.12,.420,.78,.15]);geo(ax,DEV/"pristine.device.job")
    f.text(.075,.373,"What this manual covers",fontsize=12,fontweight="bold",color=TEAL)
    for j,l in enumerate(["02  Installation: pip, virtual environments and Conda",
        "03–05  Electronic structure and EHT / BandFit parameters",
        f"06–{5+CHAPTER_PAGES}  Essential transport: electrode repeats and principal layers",
        f"{6+CHAPTER_PAGES:02d}–{periodic_start-1}  Device-by-device: XYZ, CIF, results (15 cases)",
        f"{periodic_start}–{periodic_end}  Uniform graphene and six BandFit models: XYZ then CIF",
        f"{guide_start}–{electronic_start-1}  Transport concepts, benchmarks, electronic plots and PBS",
        f"{electronic_start}–{last_page}  Molecular, crystalline, magnetic and native input examples",
        f"{last_page+1}  Example locator: exact directories and XYZ / CIF page numbers"]):
        f.text(.075,.34-j*.028,l,fontsize=10,color=NAVY)
    linkline(f,.111,ROOT/"examples","Each example includes its local directory and clickable source files.",9)
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
    y=para(f,y,"Run from repository root. Folder: examples/tutorials/bandfit. The included reference/frontier_bands.csv and initial JSONs are enough; no SIESTA executable is required. The BandFit optimizer uses differential evolution (seed 42, 180 evaluations), followed by least squares. All six fits reproduced exactly in the archived run.")
    f.text(.075,.49,"Validation RMSE (eV)",fontsize=12,fontweight="bold",color=TEAL)
    a=f.add_axes([.13,.263,.76,.19]);x=np.arange(3);a.bar(x-.16,[.1807107,.34115645,.48365584],.32,label="p-DZ7",color=BLUE);a.bar(x+.16,[.20616767,.28224699,.33004110],.32,label="sp5 (shared)",color=TEAL);a.set_xticks(x,["2 bands","4 bands","6 bands"]);a.set_ylabel("RMSE (eV)");a.legend(frameon=False);a.grid(axis="y",alpha=.2)
    y=para(f,.218,"The seven/five parameter names do not denote the number of bands. Each carbon still has four angular orbitals. Compare RMSE only with awareness of different target counts. For transport use the separately archived Dirac reference, not the fitted band-alignment shift.",size=9)
    para(f,y,"Quick job: python -m pyaehmop examples/tutorials/bandfit/pdz7_4B_smoke.transport.job --workers 4. For all six full and quick jobs, see examples/tutorials/bandfit/README.md.",size=9)
    save(pdf,f)
    add_chapter(pdf,page,para,head,code,linkline,save,DEV,DOC)
    # All device examples are continuous and complete, XYZ immediately followed by CIF.
    for n,(name,title,desc,old,side) in enumerate(cases,1):
        p=DEV/(name+".device.job") if old else TUT/name/"run.device.job"
        q=p.with_name(p.name.replace(".device.job","_cif.device.job"))
        device_format(pdf,p,q,title,desc,"XYZ",n)
        device_cif_format(pdf,p,q,title,desc,n)
        case(pdf,name,title+" — results",desc,old,side)
    # Uniform periodic calculations follow the same XYZ -> CIF order.
    for p in [TR/"graphene.transport.job",*sorted(BF.glob("*B.transport.job"))]:
        q=p.with_name(p.name.replace(".transport.job","_cif.transport.job"))
        name="Uniform graphene" if p.parent==TR else p.name.removesuffix(".transport.job")
        quick=p.with_name(p.name.replace(".transport.job","_smoke.transport.job")) if p.parent==BF else None
        quick_cif=quick.with_name(quick.name.replace(".transport.job","_cif.transport.job")) if quick else None
        uniform_format(pdf,p,"XYZ",name,quick)
        uniform_format(pdf,q,"CIF",name,quick_cif)
    f=page(pdf,"Periodic transport / validation","XYZ and CIF equivalence","Same physical graphene model, independent geometrical representation.")
    y=para(f,.85,"Both uniform p-DZ7 4B jobs use the same parameter table, K, reference energy and cutoff. The CIF has a 30 Å vacuum in y. The lines below overlap at the full 301-k / 401-energy resolution.")
    y=linkline(f,y,TR,"Folder: examples/transport")
    data=np.loadtxt(TR/"cif_reference/graphene_xyz_cif_full.dat")
    ax=f.add_axes([.13,.36,.76,.34])
    ax.plot(data[:,0],data[:,1],color=BLUE,lw=2,label="XYZ")
    ax.plot(data[:,0],data[:,2],color=ORANGE,ls="--",lw=1.1,label="CIF")
    ax.set(xlabel="E - reference (eV)",ylabel="T per spin");ax.legend(frameon=False);ax.grid(alpha=.2)
    checks=json.loads((TR/"cif_reference/validation.json").read_text())
    m=max(x["transmission.dat"] for x in checks["uniform"].values())
    fullerr=json.loads((TR/"cif_reference/full_validation.json").read_text())["transmission.dat"]
    y=para(f,.30,f"Full 4B: maximum absolute T difference {fullerr:.3g}; all six quick pairs: {m:.3g}. Other five full-resolution CIF jobs are supplied but not rerun. This establishes file-format equivalence, not physical calibration.",size=9)
    y=linkline(f,y,TR/"cif_reference/graphene_xyz_cif_full.dat")
    y=linkline(f,y,BF,"All BandFit inputs: examples/tutorials/bandfit")
    y=para(f,y,"For device CIFs, the common P1 padded box is coordinate interchange only. The physical repeats are specified in the job; never center the three regions independently.",size=9)
    save(pdf,f)
    f=page(pdf,"Transport input","Two distinct transport solvers","Uniform periodic materials versus explicit two-terminal devices.")
    y=head(f,.85,"Uniform .transport.job")
    y=linkline(f,y,TR,"Folder: examples/transport")
    y=code(f,y,"structure_format = xyz\nstructure = graphene.xyz\nlattice_vectors = 2.459512 0 0; 0 0 4.26\nimages = 5 4\ntransport_axis = 1\ntransverse_kpoints = 301")
    y=para(f,y,"The zero-based index 1 chooses the second lattice vector (z). The solver folds interactions into principal layers automatically. It describes pristine uniform materials, not an isolated central defect.")
    y=head(f,y,"Explicit .device.job (XYZ or CIF)")
    y=linkline(f,y,DEV,"Folder: examples/devices")
    y=code(f,y,"electrode_left = left.xyz\nscattering_region = pristine.xyz\nelectrode_right = right.xyz\ntransverse_vectors = 2.459512 0 0\nleft_repeat = 0 0 8.52\nright_repeat = 0 0 8.52\ncutoff = 4.5\nquadrature = 16 12 8")
    y=para(f,y,"All XYZs must share a Cartesian origin and coordinate axes. CIF alternatives use a common P1 box and coordinate frame (see the preceding CIF examples). The nearest layers can have different atom counts, species and longitudinal periods. A common transverse lattice and a shared energy convention are required. For 1D leave transverse_kpoints blank.")
    para(f,y,"These excerpts are not complete inputs; the complete jobs are in the individual preceding examples. Use the supplied runnable .job files, which also specify parameters, energy grid and output options. The device solver checks projected layer separation conservatively; enlarge layers if the check fails.",size=9)
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
    y=para(f,y,"15 device examples, 15 electronic examples, six quick/full transport models and six refits completed. Current CLI core/spin regression: 74 passed after GUI removal. Numerical validation is not DFT validation of defects or interfaces.",size=9)
    y=para(f,y,"Troubleshooting: check geometry and cutoff if layers are rejected; check overlap positivity, quadrature, k mesh and eta. For I(V) extend the energy window to include thermal tails. For SOC use a justified xi; for frozen spin specify exchange_eV by zero-based atom index.",size=9)
    para(f,y,"References: Hoffmann, JCP 39, 1397 (1963), doi:10.1063/1.1734456; López Sancho et al., J. Phys. F 15, 851 (1985); Caroli et al., J. Phys. C 4, 916 (1971); Brandbyge et al., PRB 65, 165401 (2002); Papior et al., CPC 212, 8 (2017). See examples/transport/README.md for links.",size=8)
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
        f=page(pdf,"Electronic examples",p.stem.replace("_"," "),"Original molecular, crystal or magnetic example — complete input.")
        y=para(f,.85,explanations.get(p.stem,"Original electronic structure example."),width=94,size=10)
        mod="pyaehmop.SPIN_TEST" if p.suffix==".spin" else "pyaehmop"
        y=job_folder(f,y,p,False)
        y=compact_code(f,y,"cd "+repo_path(p.parent)+"\npython -m "+mod+" "+p.name,7.3)
        lines=[l.strip() for l in p.read_text().splitlines() if l.strip() and not l.lstrip().startswith("#")]
        geom=(p.parent/lines[1]).resolve()
        f.text(.075,y,"Input file",fontsize=11,color=TEAL,fontweight="bold",va="top")
        f.text(.52,y,"Structure file",fontsize=11,color=TEAL,fontweight="bold",va="top")
        yy=y-.025
        y1=column_code(f,.075,yy,.415,"\n".join(lines),6.9,.016)
        y2=linkline(f,yy,geom,geom.name,8.2,.52)
        gl=geom.read_text().splitlines()
        # Long native inputs are reproduced in Appendix D. Wrap very long
        # lines for display; source file retains its exact line structure.
        preview=gl if len(gl)<=20 else gl[:12]
        shown=[]
        for l in preview:shown.extend(textwrap.wrap(l,width=45,replace_whitespace=False,drop_whitespace=False,break_long_words=True,break_on_hyphens=False) or [""])
        y2=column_code(f,.52,y2,.405,"\n".join(shown),6.4,.015)
        if len(gl)>20:y2=para(f,y2,"First 12 lines shown. Full native input reproduced at the end of this manual and linked above.",width=46,size=8,indent=.52)
        elif any(len(l)>45 for l in preview):y2=para(f,y2,"Long lines wrapped for display only; use the linked file for the exact text.",width=46,size=8,indent=.52)
        y=min(y1,y2)-.012
        y=linkline(f,y,geom,size=8)
        for l in lines:
            if l.startswith("parameters") and "=" in l:
                v=l.split("=",1)[1].strip()
                if (p.parent/v).is_file(): y=linkline(f,y,p.parent/v,"Parameter file: "+repo_path(p.parent/v),7.8)
        para(f,y,"Run after activating the installed environment. Paths inside the job are relative to its directory. Compare the output with examples/tutorials/reference/electronic_summary.json.",size=8,width=105)
        save(pdf,f)

    for geom in (ROOT/"pyaehmop/inputs/diamond",ROOT/"pyaehmop/inputs/XeF4"):
        lines=geom.read_text().splitlines()
        for lo in range(0,len(lines),44):
            f=page(pdf,"Native input examples",geom.name+" — YAeHMOP source",f"Original file lines {lo+1}–{min(lo+44,len(lines))} of {len(lines)}.")
            y=linkline(f,.85,geom,size=9)
            y=column_code(f,.075,y,.85,"\n".join(lines[lo:lo+44]),7.3,.0155)
            save(pdf,f)
    f=page(pdf,"Example locator","Where to find each example","Use these folders from the cloned repository root; each chapter is continuous.")
    y=.85
    for name,p,x,c,r in device_index:
        y=linkline(f,y,p.parent,name+" — "+repo_path(p.parent),8.1)
        f.text(.925,y+.021,f"XYZ {x} / CIF {c} / results {r}",fontsize=8.1,color=NAVY,ha="right",va="top")
        y-=.004
    y-=.013
    y=linkline(f,y,TR,f"Uniform graphene — examples/transport (pages {periodic_start}–{periodic_start+1})",8.4)
    y=linkline(f,y,BF,f"BandFit — examples/tutorials/bandfit (pages {periodic_start+2}–{periodic_end-1})",8.4)
    y=linkline(f,y,ROOT/"pyaehmop/examples_friendly",f"Molecules / crystals — pyaehmop/examples_friendly (pages {electronic_start}–{electronic_start+8})",8.1)
    y=linkline(f,y,ROOT/"pyaehmop/SPIN_TEST/examples",f"Spin / SOC — pyaehmop/SPIN_TEST/examples (pages {electronic_start+9}–{electronic_start+14})",8.1)
    save(pdf,f)
    assert P==last_page+1,(P,last_page)
print(out,P,"pages",out.stat().st_size,"bytes")
