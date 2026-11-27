"""Illustrated lead-repeat chapter; documentation only, no scientific solver runs."""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

CHAPTER_PAGES=6
NAVY="#17334a"; BLUE="#245d7b"; TEAL="#007f83"; ORANGE="#bb641e"; GRAY="#56636c"
C={"C":"#334b5d","Mo":"#247e8a","S":"#d39a19"}
BG=["#dce8f1","#d9eee6","#f6e6d6"]

def xyz(path):
    rows=[l.split() for l in path.read_text().splitlines()[2:] if l.strip()]
    return np.array([r[0] for r in rows]),np.array([list(map(float,r[1:4])) for r in rows])

def points(ax,s,pos,shift=0,vertical=0,size=26):
    for el in sorted(set(s)):
        m=s==el
        ax.scatter(pos[m,2]+shift,pos[m,vertical],s=size,c=C[el],edgecolors="white",linewidths=.35,zorder=3)

def arrow(ax,x1,x2,y,label,color=BLUE):
    ax.annotate("",(x2,y),(x1,y),arrowprops=dict(arrowstyle="<->",color=color,lw=1.25))
    ax.text((x1+x2)/2,y+.075,label,ha="center",va="bottom",color=color,fontsize=9)

def convention(ax,dev):
    for i,(x0,w,label,kind) in enumerate([(-3,1,"L2",0),(-2,1,"L1",0),(-1,1,"L0",0),(0,1.4,"D",1),(1.4,.7,"R0",2),(2.1,.7,"R1",2),(2.8,.7,"R2",2)]):
        ax.add_patch(Rectangle((x0,.15),w,.75,facecolor=BG[kind],edgecolor="white"))
        for x in np.linspace(x0+w*.25,x0+w*.75,2):
            ax.scatter([x,x],[.35,.7],s=18,c=[BLUE,TEAL,ORANGE][kind])
        ax.text(x0+w/2,.98,label,ha="center",fontsize=9)
    arrow(ax,-2.5,-1.5,1.25,r"$\ell_L$")
    arrow(ax,1.75,2.45,1.25,r"$\ell_R$",ORANGE)
    ax.annotate("−∞",(-3.12,.53),ha="right");ax.annotate("+∞",(3.62,.53))
    ax.annotate("positive transport direction",(3.3,-.05),(-2.9,-.05),arrowprops=dict(arrowstyle="->",color=GRAY),ha="left",va="center",fontsize=9)
    ax.set(xlim=(-3.55,3.95),ylim=(-.17,1.55));ax.axis("off")

def chain(ax,dev):
    s,p=xyz(dev/"chain_left.xyz");t,q=xyz(dev/"chain_device.xyz");u,r=xyz(dev/"chain_right.xyz")
    ax.axvspan(-7.7,-.75,color=BG[0]);ax.axvspan(-.75,2.7,color=BG[1]);ax.axvspan(2.7,11,color=BG[2])
    ax.axhline(0,color="#b3bdc3",lw=.8)
    for i in range(5):
        points(ax,s,p,-1.5*i,0,65);points(ax,u,r,1.8*i,0,65)
        ax.text(-1.5*(i+1),-.16,"C",ha="center",fontsize=8);ax.text(3.6+1.8*i,-.16,"S",ha="center",fontsize=8)
    points(ax,t,q,0,0,72);ax.text(0,-.16,"C",ha="center",fontsize=8);ax.text(1.8,-.16,"Mo",ha="center",fontsize=8)
    arrow(ax,-3,-1.5,.29,"1.5 Å");arrow(ax,3.6,5.4,.29,"1.8 Å",ORANGE)
    ax.set(xlim=(-8,11.3),ylim=(-.32,.67),xlabel="z (Å)");ax.set_yticks([])
    ax.set_xticks([-7.5,-6,-4.5,-3,-1.5,0,1.8,3.6,5.4,7.2,9,10.8]);ax.tick_params(axis="x",labelsize=7)
    for x,t in [(-4.2,"C lead"),(.9,"D"),(6.7,"S lead")]:ax.text(x,.53,t,ha="center",fontsize=9)
    ax.spines["left"].set_visible(False)

def graphene(ax,dev):
    for i,(name,shift,bg) in enumerate([("left.xyz",-8.52,BG[0]),("left.xyz",0,BG[0]),("pristine.xyz",0,BG[1]),("right.xyz",0,BG[2]),("right.xyz",8.52,BG[2])]):
        s,p=xyz(dev/name);z0=p[:,2].min()+shift
        ax.axvspan(z0,z0+8.52,color=bg,alpha=.65)
        points(ax,s,p,shift,0,22)
        ax.axvline(z0+4.26,color=GRAY,lw=.7,ls=":")
        ax.text(z0+4.26,1.68,["L1","L0","D","R0","R1"][i],ha="center",fontsize=8)
    arrow(ax,-8.52,-4.26,2.05,"4.26 Å")
    arrow(ax,8.52,17.04,2.05,"8.52 Å",ORANGE)
    ax.set(xlabel="z (Å)",ylabel="x (Å)",xlim=(-17.6,26),ylim=(-.35,2.48))
    ax.set_xticks([-17.04,-8.52,0,8.52,17.04,25.56]);ax.tick_params(axis="x",labelsize=8)

def hetero(ax,dev):
    for name,shift,bg,label in [("hetero_left.xyz",-8.52,BG[0],"L1"),("hetero_left.xyz",0,BG[0],"L0"),("hetero_device.xyz",0,BG[1],"D"),("hetero_right.xyz",0,BG[2],"R0"),("hetero_right.xyz",5.679999661100459,BG[2],"R1")]:
        s,p=xyz(dev/name);lo=p[:,2].min()+shift;hi=p[:,2].max()+shift
        ax.axvspan(lo-.25,hi+.25,color=bg,alpha=.7);points(ax,s,p,shift,0,13)
        ax.text((lo+hi)/2,9.25,label,ha="center",fontsize=8)
    ax.annotate("",(-8.52,10.5),(-17.04,10.5),arrowprops=dict(arrowstyle="<->",color=BLUE));ax.text(-12.78,10.8,"8.52 Å",ha="center",color=BLUE)
    ax.annotate("",(20.48,10.5),(14.8,10.5),arrowprops=dict(arrowstyle="<->",color=ORANGE));ax.text(17.64,10.8,"5.68 Å",ha="center",color=ORANGE)
    ax.set(xlim=(-18,25.3),ylim=(-.5,12),xlabel="z (Å)",ylabel="x (Å)")
    ax.legend(handles=[Line2D([],[],ls="none",marker="o",color=C[e],label=e,markersize=5) for e in C],loc="upper center",bbox_to_anchor=(.5,1.18),ncol=3,frameon=False,fontsize=8)

def cutoff(ax,dev):
    for row,y in enumerate([1.15,0.15]):
        for i in range(3):
            ax.add_patch(Rectangle((i*2,y),2,.65,fc=BG[0],ec="white"))
            ax.scatter([i*2+.2,i*2+.75,i*2+1.3,i*2+1.8],[y+.25]*4,s=21,c=NAVY)
            ax.text(i*2+1,y+.48,f"Layer {i}",ha="center",fontsize=8)
        end=4.2 if row==0 else 2.2
        ax.annotate("",(end,y+.3),(1.8,y+.3),arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=-.55",color="#ac3535" if row==0 else TEAL,lw=1.7))
        ax.text(6.15,y+.32,"Rejected:\nsecond-neighbor coupling" if row==0 else "Allowed:\nnearest layers only",va="center",fontsize=8,color="#ac3535" if row==0 else TEAL)
    ax.set(xlim=(-.1,9),ylim=(0,2.1));ax.axis("off")

def crystal(ax,dev):
    # oblique drawing of orthogonal parallelepipeds, c parallel to z.
    for base,n in [(0,1),(4.5,2)]:
        for j in range(n):
            x=base+j*1.7
            a=np.array([[x,.1],[x+1.7,.1],[x+2.1,.6],[x+.4,.6],[x,.1]])
            b=a+[0,1.1]
            ax.plot(*a.T,color=GRAY);ax.plot(*b.T,color=GRAY)
            for k in range(4):ax.plot([a[k,0],b[k,0]],[a[k,1],b[k,1]],color=GRAY)
        ax.annotate("",(base+n*1.7,-.15),(base,-.15),arrowprops=dict(arrowstyle="<->",color=TEAL,lw=1.5))
        ax.text(base+n*.85,-.3,"6.8 Å = c" if n==1 else "13.6 Å = 2c",ha="center",color=TEAL,fontsize=10)
        ax.text(base+n*.85,1.95,"One cell" if n==1 else "Two cells",ha="center",fontsize=10)
    ax.set(xlim=(-.4,8.5),ylim=(-.5,2.15));ax.axis("off")

FIGURES=[("convention",convention),("chain",chain),("graphene",graphene),("hetero",hetero),("cutoff",cutoff),("crystal",crystal)]

def write_figures(doc,dev):
    for name,fn in FIGURES:
        f,ax=plt.subplots(figsize=(8.1,2.25))
        fn(ax,dev);f.tight_layout(pad=1.0)
        target=doc/"figures"/("lead_repeat_"+name+".svg")
        f.savefig(target,metadata={"Date":None})
        target.write_text("\n".join(line.rstrip() for line in target.read_text().splitlines())+"\n")
        plt.close(f)

def add_chapter(pdf,page,para,head,code,linkline,save,dev,doc):
    write_figures(doc,dev)
    def draw(f,fn,y=.56,h=.23):
        fn(f.add_axes([.11,y,.79,h]),dev)
    def eq(f,y,s,size=15):
        f.text(.09,y,s,fontsize=size,color=NAVY,va="top")
        return y-.052
    def note(f,y,s):
        return para(f,y,s,width=100,size=9)
    def links(f,y,name):
        y=linkline(f,y,dev,"Folder: examples/devices",8.8)
        y=linkline(f,y,dev/(name+".device.job"),"XYZ input: "+name+".device.job",8.5)
        return linkline(f,y,dev/(name+"_cif.device.job"),"CIF input: "+name+"_cif.device.job",8.5)

    f=page(pdf,"Lead repeats / 1","What is an electrode repeat?","Essential reading before the explicit two-terminal transport examples.")
    y=para(f,.85,"The supplied L and R files contain the nearest principal layers of two semi-infinite electrodes. D is the finite scattering region. The repeat vectors reproduce the entire electrode motif; they do not specify an interface gap.")
    draw(f,convention,.565,.20)
    y=eq(f,.535,r"$\mathbf{r}_L^{(n)}=\mathbf{r}_L^{(0)}-n\mathbf{a}_L,\quad \mathbf{r}_R^{(n)}=\mathbf{r}_R^{(0)}+n\mathbf{a}_R$")
    y=note(f,y,"n = 0, 1, 2, ... . Both repeat vectors point along the same positive transport normal. The code applies the negative translation for the left electrode. D is not periodically repeated longitudinally.")
    y=head(f,y,"The physical rule")
    y=eq(f,y,r"$\mathbf{a}_{\mathrm{repeat}}=N\mathbf{c}$",19)
    y=note(f,y,"Choose a physical bulk longitudinal lattice vector c and build N consecutive cells in the principal layer. Select N for each electrode independently. Identical leads with identical layer construction have equal repeats; different materials need not.")
    y=note(f,y,"A repeat is neither the first-to-last atom extent nor a nearest-neighbor bond length in general. An ordinary XYZ file does not uniquely determine a bulk period. Use the physical lattice or an explicitly defined periodic motif.")
    linkline(f,y,doc/"lead-repeats.md","Read this chapter on GitHub: docs/lead-repeats.md",8.7)
    save(pdf,f)

    f=page(pdf,"Lead repeats / 2","Start with a one-atom chain","Actual asymmetric example: independent C and S lead periods.")
    y=links(f,.85,"asymmetric")
    draw(f,chain,.53,.23)
    y=code(f,.49,"# Actual supplied files (z in angstrom)\nchain_left.xyz:   C   0 0 -1.5\nchain_device.xyz: C   0 0  0.0\n                  Mo  0 0  1.8\nchain_right.xyz:  S   0 0  3.6\nleft_repeat  = 0 0 1.5\nright_repeat = 0 0 1.8\ncutoff = 2.0",8.5)
    y=eq(f,y,r"$z_L(n)=-1.5-1.5n,\qquad z_R(n)=3.6+1.8n$")
    y=note(f,y,"Left C positions: -1.5, -3.0, -4.5, -6.0, -7.5 Å. Right S positions: 3.6, 5.4, 7.2, 9.0, 10.8 Å. With one atom per layer the repeat is the spacing between neighboring atoms. These spacings are illustrative, not optimized. No transverse vectors.")
    code(f,y,"cd examples/devices\npython -m pyaehmop asymmetric.device.job --workers 2\npython -m pyaehmop asymmetric_cif.device.job --workers 2",7.7)
    save(pdf,f)

    f=page(pdf,"Lead repeats / 3","Graphene: why 8.52 Å?","Two rectangular bulk cells form the eight-carbon principal layer.")
    y=links(f,.85,"pristine")
    draw(f,graphene,.555,.22)
    y=eq(f,.515,r"$\mathbf{a}_x=(2.459512,0,0),\quad \mathbf{c}=(0,0,4.26)\ \mathrm{\AA}$")
    y=eq(f,y,r"$\mathbf{a}_L=\mathbf{a}_R=2\mathbf{c}=(0,0,8.52)\ \mathrm{\AA}$")
    y=code(f,y,"# Rectangular cell: C (x,y,z), angstrom\nC  0.000000  0  0.00\nC  0.000000  0  1.42\nC  1.229756  0  2.13\nC  1.229756  0  3.55\n# Repeat two cells along z to form L and R\ntransverse_vectors = 2.459512 0 0\nleft_repeat = 0 0 8.52\nright_repeat = 0 0 8.52",7.8)
    y=note(f,y,"Left z: -8.52 to -0.71 Å; D: 0 to 7.81 Å; right: 8.52 to 16.33 Å. The left atom extent is 7.81 Å, NOT the translation. The additional 0.71 Å is a projected gap, not a C-C bond length. The full XYZ and CIF files appear in the following pristine example.")
    code(f,y,"cd examples/devices\npython -m pyaehmop pristine.device.job --workers 4\npython -m pyaehmop pristine_cif.device.job --workers 4",7.7)
    save(pdf,f)

    f=page(pdf,"Lead repeats / 4","Graphene–MoS2: unequal repeats","Derivation of 8.52 Å and 5.679999661100459 Å from the actual generator.")
    y=links(f,.85,"graphene_mos2")
    draw(f,hetero,.565,.18)
    y=note(f,.525,"Four graphene cells along x are matched to three MoS2 cells in the illustrative junction. The MoS2 model is strained isotropically in-plane to this width; no structural relaxation or interface calibration was performed.")
    y=eq(f,y,r"$W=4(2.459512)=9.838048\ \mathrm{\AA}$")
    y=eq(f,y,r"$a_{\mathrm{MoS_2}}=W/3=3.279349333\ldots\ \mathrm{\AA}$")
    y=eq(f,y,r"$L_L=2(4.26)=8.52\ \mathrm{\AA},\quad L_R=\sqrt{3}\,a_{\mathrm{MoS_2}}=5.679999661\ldots\ \mathrm{\AA}$",12)
    y=code(f,y,"transverse_vectors = 9.838048 0 0\nleft_repeat = 0 0 8.52\nright_repeat = 0 0 5.679999661100459",8.2)
    y=note(f,y,"L: 32 C; D: 32 C + 6 Mo + 12 S; R: 6 Mo + 12 S. Equivalent Mo sites at x=y=0 have z=9.12 in D and 14.799999661 in R. Their difference confirms the right period because the D buffer is a copy of the lead.")
    y=linkline(f,y,dev/"generate.py","Exact construction: examples/devices/generate.py",8.5)
    code(f,y,"cd examples/devices\npython -m pyaehmop graphene_mos2.device.job --workers 8\npython -m pyaehmop graphene_mos2_cif.device.job --workers 8",7.4)
    save(pdf,f)

    f=page(pdf,"Lead repeats / 5","Thickness versus cutoff","The actual conservative geometry checks in the device solver.")
    y=note(f,.85,"Only neighboring principal-layer blocks are retained. Let n be the unit normal to the shared transverse lattice (or the transport direction for a 1D device). Project positions and lead vectors onto it.")
    y=eq(f,y,r"$\ell_i=\mathbf{a}_i\cdot\hat{\mathbf{n}},\quad p_i=\mathbf{r}_i\cdot\hat{\mathbf{n}},\quad w_i=p_i^{\max}-p_i^{\min}$",13)
    y=eq(f,y,r"$g_L=p_D^{\min}-p_L^{\max},\quad g_R=p_R^{\min}-p_D^{\max}$",13)
    draw(f,cutoff,.49,.175)
    y=note(f,.465,"For cutoff rc, the following five strict inequalities must hold (and the regions must be ordered L, D, R):")
    y=eq(f,y,r"$2\ell_L-w_L>r_c,\quad 2\ell_R-w_R>r_c$",14)
    y=eq(f,y,r"$g_L+\ell_L>r_c,\quad g_R+\ell_R>r_c,\quad p_R^{\min}-p_L^{\max}>r_c$",12)
    y=note(f,y,"Pristine graphene: repeat 8.52 Å, atom extent 7.81 Å and projected contact gap 0.71 Å. All five bounds are 9.23 Å: cutoff 4.5 Å passes, 11 Å fails.")
    y=code(f,y,"N cells    atoms/layer    repeat (Å)    2 repeat - extent (Å)\n  1             4            4.26              4.97\n  2             8            8.52              9.23\n  3            12           12.78             13.49\n  4            16           17.04             17.75\n  5            20           21.30             22.01",7.6)
    y=note(f,y,"Three cells pass the lead-only check for cutoff 11 Å; also check D length and buffers. Rebuild coordinates with every change of N. The tests are sufficient, not necessary; passing does not prove numerical convergence.")
    save(pdf,f)

    f=page(pdf,"Lead repeats / 6","From a bulk CIF or XYZ to L/D/R","A general workflow for building a new explicit device.")
    y=note(f,.85,"For an orthogonal bulk crystal with a=3.20 Å, b=4.50 Å, c=6.80 Å and transport along c, the longitudinal vector is (0,0,6.80) Å.")
    draw(f,crystal,.59,.19)
    y=code(f,.56,"# One cell (N=1)             # Two cells (N=2)\nleft_repeat = 0 0 6.8        left_repeat = 0 0 13.6",8.1)
    y=note(f,y,"For nonorthogonal lattices use the actual Cartesian components of the longitudinal lattice vector or an integer supercell combination. A lattice length is not generally its z component.")
    y=head(f,y,"Construction checklist")
    for t in ["1. Start from the physical bulk cells of the leads. Choose transport direction and a common transverse lattice; no transverse vector for a chain or finite-width ribbon.",
              "2. Build N cells for each principal layer; place L, D and R in one common frame and include contact buffers in D.",
              "3. Export three XYZ or three ordered P1 CIF files in the same orthogonal auxiliary box. Do not independently center or wrap the regions. The auxiliary box does not determine the physical repeats.",
              "4. Set the repeats from the original bulk vectors. Check all five cutoff inequalities; converge cutoff, k grid, eta, quadrature and device length. Validate a pristine reference before a defect or interface."]:
        y=note(f,y,t)
    y=note(f,y,"Example along x: graphene_x uses 3 × 2.459512 = 7.378536 Å along x and the transverse vector (0,0,4.26) Å. The following walkthroughs give exact input paths, structures and commands.")
    y=linkline(f,y,dev.parent/"tutorials/graphene_x","Folder: examples/tutorials/graphene_x",8.4)
    linkline(f,y,doc/"lead-repeats.md","Full chapter and figures: docs/lead-repeats.md",8.4)
    save(pdf,f)
