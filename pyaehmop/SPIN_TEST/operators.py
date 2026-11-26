"""Atomic SOC operators for the real s/p/d/f PyYAeHMOP basis.

Convention: H_SOC = xi * L.S with hbar=1 and S=sigma/2.
Spinor block ordering: all AO-up followed by all AO-down.
"""
from functools import lru_cache
import numpy as np
try:
    from scipy.special import sph_harm_y
except ImportError:  # SciPy < 1.15
    from scipy.special import sph_harm
    def sph_harm_y(l,m,theta,phi):
        return sph_harm(m,l,phi,theta)
from pyaehmop.basis import angular_functions


@lru_cache(maxsize=4)
def angular_momentum(l):
    """Hermitian Lx,Ly,Lz in the real-harmonic basis ordering.

    The transformation is obtained through exact angular projection (for spherical
    harmonics up to l=3) onto complex spherical harmonics Y_lm.
    """
    if l < 0 or l > 3:
        raise NotImplementedError("SOC is tested only for s/p/d/f harmonics")
    n=2*l+1
    if l==0:
        return tuple(np.zeros((1,1),complex) for _ in range(3))
    theta_nodes,weights=np.polynomial.legendre.leggauss(24)
    theta=np.arccos(theta_nodes)
    phi=np.arange(64)*2*np.pi/64
    xyz=np.stack(np.broadcast_arrays(
        np.sqrt(1-theta_nodes[:,None]**2)*np.cos(phi)[None,:],
        np.sqrt(1-theta_nodes[:,None]**2)*np.sin(phi)[None,:],
        theta_nodes[:,None]+0*phi[None,:]),axis=-1)
    real=np.stack(angular_functions(xyz,l),axis=-1)
    complex_y=np.stack([sph_harm_y(l,m,theta[:,None],phi[None,:])
                        for m in range(-l,l+1)],axis=-1)
    U=np.einsum("ipm,ipj,i->mj",complex_y.conj(),real,weights)*2*np.pi/len(phi)
    if not np.allclose(U.conj().T@U,np.eye(n),atol=1e-11):
        raise AssertionError(f"Base angular l={l} no unitaria")
    m_vals=np.arange(-l,l+1)
    Lz=np.diag(m_vals).astype(complex)
    Lplus=np.zeros((n,n),complex)
    for idx,m in enumerate(m_vals[:-1]):
        Lplus[idx+1,idx]=np.sqrt(l*(l+1)-m*(m+1))
    Lminus=Lplus.conj().T
    Lx=(Lplus+Lminus)/2
    Ly=(Lplus-Lminus)/(2j)
    result=tuple(U.conj().T@L@U for L in (Lx,Ly,Lz))
    return tuple((L+L.conj().T)/2 for L in result)


def soc_matrix(basis, atom_symbols, xi_by_species, orbital_labels=("d",)):
    """SOC onsite 2N x 2N.

    xi_by_species: {'FE': {'d': .06, 'p':0}} in eV.
    The EHT orbitals are nonorthogonal between atomic sites; SOC is approximated
    as an onsite atomic term, without two-center corrections.
    """
    n=len(basis)
    L_tot=[np.zeros((n,n),complex) for _ in range(3)]
    used=set()
    for i,orb in enumerate(basis):
        if (orb.atom,orb.n,orb.l) in used:
            continue
        used.add((orb.atom,orb.n,orb.l))
        indices=[j for j,o in enumerate(basis)
                 if (o.atom,o.n,o.l)==(orb.atom,orb.n,orb.l)]
        if len(indices)!=2*orb.l+1:
            raise ValueError("Shell angular incompleta no SOC")
        symbol=atom_symbols[orb.atom].upper()
        strength=xi_by_species.get(symbol,{}).get("spdf"[orb.l],0.)
        if strength==0:
            continue
        for cart,L in enumerate(angular_momentum(orb.l)):
            L_tot[cart][np.ix_(indices,indices)]+=strength*L
    lx,ly,lz=L_tot
    h=np.zeros((2*n,2*n),complex)
    h[:n,:n]=lz/2
    h[n:,n:]=-lz/2
    h[:n,n:]=(lx-1j*ly)/2
    h[n:,:n]=(lx+1j*ly)/2
    h=(h+h.conj().T)/2
    return h


def spin_operators(overlap):
    """Metricas sigma_x,y,z no espaco espinorial generalizado."""
    S=np.asarray(overlap)
    z=np.zeros_like(S)
    Mx=np.block([[z,S],[S,z]])
    My=np.block([[z,-1j*S],[1j*S,z]])
    Mz=np.block([[S,z],[z,-S]])
    return Mx,My,Mz
