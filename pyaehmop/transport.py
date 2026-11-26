"""Coherent, non-orthogonal EHT transport for identical periodic leads."""
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from .basis import make_basis
from .engine import _realspace_blocks, coupling

def realspace(structure, cutoff=11., quadrature=(24,20,12), workers=1):
    basis=make_basis(structure.atoms)
    eps=np.array([o.energy for o in basis])
    blocks=_realspace_blocks(structure,basis,quadrature,cutoff,workers)
    R=np.array(list(blocks))
    S=np.array(list(blocks.values()))
    H=coupling(eps,structure.k_constant,structure.weighted)*S
    H[np.all(R==0,axis=1)]+=np.diag(eps*(1-structure.k_constant))
    return R,H,S

def surface(a,b,c,tol=1e-10,maxiter=100):
    """Retarded surface inverse by Lopez Sancho decimation."""
    a0=a.copy(); a=a.copy(); b=b.copy(); c=c.copy()
    for _ in range(maxiter):
        gb=np.linalg.solve(a,b); gc=np.linalg.solve(a,c)
        bc=b@gc; cb=c@gb
        a0-=bc; a-=bc+cb
        b=-b@gb; c=-c@gc
        if max(np.max(abs(b)),np.max(abs(c)))<tol:
            return np.linalg.inv(a0)
    raise RuntimeError("Surface Green function did not converge")

def principal_layers(R,H,S,k,axis=-1):
    """Group long-range blocks into nearest-neighbor principal layers.
    k contains reduced reciprocal coordinates for the remaining axes.
    Returns blocks in order (-1,0,+1).
    """
    R=np.asarray(R); axis=axis%R.shape[1]
    transverse=np.delete(R,axis,axis=1)
    k=np.asarray(k,dtype=float)
    if k.shape!=(transverse.shape[1],):
        raise ValueError("Wrong number of transverse k coordinates")
    m=max(1,int(np.max(abs(R[:,axis])))); n=H.shape[1]
    hh=np.zeros((3,m*n,m*n),complex); ss=hh.copy()
    for r,h,s,phase in zip(R,H,S,np.exp(2j*np.pi*(transverse@k))):
        for i in range(m):
            layer,j=divmod(i+int(r[axis]),m)
            ii=slice(i*n,(i+1)*n); jj=slice(j*n,(j+1)*n)
            hh[layer+1,ii,jj]+=phase*h; ss[layer+1,ii,jj]+=phase*s
    for a in (hh,ss):
        if not np.allclose(a[0],a[2].conj().T,atol=1e-8) or not np.allclose(a[1],a[1].conj().T,atol=1e-8):
            raise ValueError("Non-Hermitian real-space blocks")
    return hh,ss

def transmission(energies,hh,ss,eta=1e-5):
    """One pristine principal layer between identical semi-infinite leads.
    Returns transmission per spin; no spin degeneracy factor is included.
    """
    if eta<=0: raise ValueError("eta must be positive")
    out=[]
    for e in energies:
        z=e+1j*eta
        a=z*ss[1]-hh[1]; b=z*ss[2]-hh[2]; c=z*ss[0]-hh[0]
        L=c@surface(a,c,b)@b; R=b@surface(a,b,c)@c
        gl=1j*(L-L.conj().T); gr=1j*(R-R.conj().T)
        G=np.linalg.inv(a-L-R)
        t=float(np.trace(gl@G@gr@G.conj().T).real)
        if not np.isfinite(t) or t < -1e-7:
            raise RuntimeError(f"Unphysical transmission at {e}: {t}")
        out.append(max(0.,t))
    return np.array(out)

def _task(arg):
    R,H,S,k,axis,energies,eta=arg
    return transmission(energies,*principal_layers(R,H,S,k,axis),eta)

def spectrum(R,H,S,energies,kpoints,axis=-1,eta=1e-5,workers=1,weights=None):
    """Return weighted T(E) and T(k,E); energies are absolute eV."""
    if workers<1 or len(kpoints)==0: raise ValueError("Empty grid or invalid workers")
    args=[(R,H,S,k,axis,energies,eta) for k in kpoints]
    if workers==1: values=np.array(list(map(_task,args)))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            values=np.array(list(pool.map(_task,args)))
    w=np.ones(len(values)) if weights is None else np.asarray(weights,float)
    if w.shape!=(len(values),) or not np.isfinite(w).all() or np.any(w<0) or w.sum()<=0:
        raise ValueError("Invalid k weights")
    return np.average(values,axis=0,weights=w),values
