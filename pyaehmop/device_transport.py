"""Two-terminal nonorthogonal EHT with independent, explicit principal layers.

All XYZ coordinates share a Cartesian frame. Left/right repeat vectors point
from left to right. The supplied lead layers must contain the nearest contact
cells and be thick enough for the chosen real-space cutoff.
"""
from dataclasses import dataclass
from itertools import product
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from .basis import make_basis
from .engine import coupling, _atom_slices
from .overlap import overlap_block
from .transport import surface

@dataclass
class Region:
    atoms: list
    offset: float = 0.
    exchange: dict | None = None  # atom-index -> shell -> frozen I*m/2 (eV)
    soc: dict | None = None       # element -> shell -> xi (eV)
    def __post_init__(self):
        self.basis=make_basis(self.atoms)
        self.pos=np.array([a.position for a in self.atoms],float)
        self.eps=np.array([o.energy for o in self.basis])
        self.slices=_atom_slices(self.basis,len(self.atoms))
        self.shift=np.array([float((self.exchange or {}).get(str(o.atom),{}).get("spdf"[o.l],0.)) for o in self.basis])
        if not np.isfinite(self.pos).all() or not np.isfinite(self.shift).all() or not np.isfinite(self.offset):
            raise ValueError("Nonfinite region geometry/potential")
        for i, shells in (self.exchange or {}).items():
            if not 0<=int(i)<len(self.atoms) or any(l not in "spdf" for l in shells):
                raise ValueError("Invalid exchange atom/shell")
    def onsite_soc(self):
        from .SPIN_TEST.operators import soc_matrix
        # Original routine defaults to d but internally uses each actual shell.
        return soc_matrix(self.basis,[a.symbol for a in self.atoms],self.soc or {})

def _pair_task(arg):
    t,i,j,oa,ob,r,q=arg
    return t,i,j,overlap_block(oa,ob,r,q)

def blocks(a,b,transverse,translation,cutoff,K,weighted,workers,quadrature,spinor=False,onsite=False):
    """Real-space H/S from a to translated b, all relevant transverse images."""
    v=np.asarray(transverse,float).reshape(-1,3)
    delta=b.pos[None,:,:]+translation-a.pos[:,None,:]
    if len(v):
        dual=np.linalg.pinv(v.T)
        extent=np.ceil(np.max(abs(delta@dual.T),axis=(0,1))+cutoff*np.linalg.norm(dual,axis=1)).astype(int)
    else: extent=[]
    tasks=[]
    for t in product(*(range(-n,n+1) for n in extent)):
        shift=np.asarray(t,float)@v
        for i,j in np.argwhere(np.linalg.norm(delta+shift,axis=2)<=cutoff):
            tasks.append((t,i,j,tuple(a.basis[x] for x in a.slices[i]),tuple(b.basis[x] for x in b.slices[j]),delta[i,j]+shift,quadrature))
    mats={}
    def collect(results):
        for t,i,j,s in results:
            mat=mats.setdefault(t,np.zeros((len(a.basis),len(b.basis))))
            mat[np.ix_(a.slices[i],b.slices[j])]=s
    if workers>1 and len(tasks)>8:
        with ProcessPoolExecutor(max_workers=workers) as pool: collect(pool.map(_pair_task,tasks,chunksize=16))
    else: collect(map(_pair_task,tasks))
    n,m=len(a.basis),len(b.basis)
    F=coupling(np.r_[a.eps,b.eps],K,weighted)[:n,n:]
    out={}
    for t,S in mats.items():
        H=F*S
        if onsite and not any(t):
            H+=np.diag(a.eps*(1-K))
            S=(S+S.T)/2; H=(H+H.T)/2
        H=H+(a.offset+b.offset)/2*S
        if spinor:
            ex=(a.shift[:,None]+b.shift[None,:])/2*S
            zero=np.zeros_like(S)
            H=np.block([[H-ex,zero],[zero,H+ex]]).astype(complex)
            S=np.block([[S,zero],[zero,S]])
            if onsite and not any(t) and a.soc: H+=a.onsite_soc()
        out[t]=(H,S)
    shape=(n*(2 if spinor else 1),m*(2 if spinor else 1))
    return out,shape

def bloch(block,k):
    mats,shape=block
    H=np.zeros(shape,complex); S=H.copy()
    for t,(h,s) in mats.items():
        phase=np.exp(2j*np.pi*np.dot(t,k))
        H+=phase*h; S+=phase*s
    return H,S

def _normal(transverse,left,right):
    v=np.asarray(transverse,float).reshape(-1,3)
    if len(v)>2 or (len(v) and np.linalg.matrix_rank(v)!=len(v)):
        raise ValueError("Transverse vectors must be independent")
    l=np.array(left,float); r=np.array(right,float)
    if len(v):
        l-=v.T@np.linalg.solve(v@v.T,v@l)
    if np.linalg.norm(l)<1e-9: raise ValueError("Left repeat parallel to transverse plane")
    n=l/np.linalg.norm(l)
    rp=r-v.T@np.linalg.solve(v@v.T,v@r) if len(v) else r
    if np.dot(rp,n)<=0 or np.linalg.norm(rp-n*np.dot(rp,n))>1e-7:
        raise ValueError("Lead repeats must have a common positive longitudinal direction")
    return n

def build_device(L,D,R,transverse,left_repeat,right_repeat,cutoff=11.,K=1.75,weighted=True,workers=1,quadrature=(24,20,12),spinor=False):
    """Build rectangular interface matrices, rejecting omitted long-range bonds."""
    v=np.asarray(transverse,float).reshape(-1,3); l=np.asarray(left_repeat,float); r=np.asarray(right_repeat,float)
    n=_normal(v,l,r); pl,pd,pr=(x.pos@n for x in (L,D,R))
    if not (pl.max()<pd.min()-1e-7 and pd.max()<pr.min()-1e-7):
        raise ValueError("Contact layers must be ordered L, D, R in one Cartesian frame")
    # Projected separation gives a conservative proof of no omitted interactions.
    if 2*(l@n)-(pl.max()-pl.min())<=cutoff or 2*(r@n)-(pr.max()-pr.min())<=cutoff:
        raise ValueError("Lead principal layer too thin: enlarge it until second-neighbor couplings vanish")
    if pd.min()-pl.max()+l@n<=cutoff or pr.min()-pd.max()+r@n<=cutoff:
        raise ValueError("Device couples to a deeper lead layer: enlarge contact principal layers")
    if pr.min()-pl.max()<=cutoff:
        raise ValueError("Direct lead-to-lead coupling possible: lengthen device")
    def b(a,z,t,onsite=False):
        return blocks(a,z,v,t,cutoff,K,weighted,workers,quadrature,spinor,onsite)
    zero=np.zeros(3)
    return { "L":b(L,L,zero,True),"Lm":b(L,L,-l),"D":b(D,D,zero,True),"R":b(R,R,zero,True),"Rp":b(R,R,r),"DL":b(D,L,zero),"DR":b(D,R,zero)}

def _adj(pair): return tuple(x.conj().T for x in pair)
def _a(pair,z): return z*pair[1]-pair[0]

def solve_k(matrices,energies,k,eta=1e-5,channels=False):
    B={name:bloch(block,k) for name,block in matrices.items()}
    for name in ("L","D","R"):
        H,S=B[name]
        if not np.allclose(H,H.conj().T,atol=1e-8) or not np.allclose(S,S.conj().T,atol=1e-8):
            raise ValueError("Non-Hermitian onsite block")
        if np.linalg.eigvalsh(S).min()<=1e-9: raise ValueError("Overlap not positive definite")
    values=[]
    for e in energies:
        z=e+1j*eta
        dl,ld=_a(B["DL"],z),_a(_adj(B["DL"]),z)
        dr,rd=_a(B["DR"],z),_a(_adj(B["DR"]),z)
        lm,lp=_a(B["Lm"],z),_a(_adj(B["Lm"]),z)
        rp,rm=_a(B["Rp"],z),_a(_adj(B["Rp"]),z)
        sl=dl@surface(_a(B["L"],z),lm,lp)@ld
        sr=dr@surface(_a(B["R"],z),rp,rm)@rd
        gl=1j*(sl-sl.conj().T); gr=1j*(sr-sr.conj().T)
        G=np.linalg.solve(_a(B["D"],z)-sl-sr,np.eye(len(sl)))
        A=gl@G@gr@G.conj().T
        vals=[np.trace(A).real]
        if channels:
            q=len(A)//2
            vals.extend([np.trace(A[:q,:q]).real,np.trace(A[q:,q:]).real])
        if not np.isfinite(vals).all() or min(vals)<-1e-7: raise RuntimeError("Unphysical transmission")
        values.append(np.maximum(vals,0.))
    return np.array(values)

def _solve_task(args): return solve_k(*args)
def spectrum_device(matrices,energies,kpoints,eta=1e-5,workers=1,channels=False):
    if eta<=0 or workers<1 or not len(kpoints): raise ValueError("Invalid grid")
    args=[(matrices,energies,k,eta,channels) for k in kpoints]
    if workers==1: tk=np.array(list(map(_solve_task,args)))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool: tk=np.array(list(pool.map(_solve_task,args)))
    return tk.mean(axis=0),tk

