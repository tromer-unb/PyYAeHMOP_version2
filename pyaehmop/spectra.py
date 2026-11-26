"""Electronic properties with k weights, Gaussian DOS and bond analyses.

Independent of the C program. Discrete bands/DOS use only
the provided k points: physical convergence requires a suitable mesh.
"""
import numpy as np


def weighted_occupations(energies, weights, electrons, tolerance=1e-5):
    """Occupy spin-degenerate states (2 electrons/band) with normalized k weights.

    energies: shape=(nk,nb)
    weights: shape=(nk,) nonnegative integration weight
    returns (state occupations, Fermi energy), with 0<=f<=2.
    """
    E=np.asarray(energies,dtype=float)
    w=np.asarray(weights,dtype=float)
    if E.ndim!=2 or w.shape!=(E.shape[0],):
        raise ValueError("Incompatible band/weight shapes")
    if np.any(w<0) or w.sum()<=0:
        raise ValueError("K weights must be nonnegative and sum to a positive value")
    wn=w/w.sum()
    if electrons<0 or electrons>2*E.shape[1]+1e-8:
        raise ValueError("Electron count exceeds the cell orbital capacity")
    ids=np.argsort(E.ravel(),kind="stable")
    energies_flat=E.ravel()
    weights_flat=np.repeat(wn,E.shape[1])
    filled=np.zeros(len(ids),dtype=float)
    rest=float(electrons)
    i=0
    fermi=None
    while i<len(ids):
        j=i+1
        while j<len(ids) and abs(energies_flat[ids[j]]-energies_flat[ids[i]])<=tolerance:
            j+=1
        group=ids[i:j]
        capacity=2*weights_flat[group].sum()
        taken=min(rest,capacity)
        if capacity>0 and taken>0:
            filled[group]=2*taken/capacity
            fermi=energies_flat[group[-1]]
            rest-=taken
        i=j
    return filled.reshape(E.shape),float(fermi) if fermi is not None else None


def gaussian_dos(energies,weights,points=600,sigma=0.15,energy_range=None,
                 include_spin=True):
    """Total Gaussian-broadened DOS: states/eV per cell."""
    E=np.asarray(energies,dtype=float)
    w=np.asarray(weights,dtype=float)
    if E.ndim!=2 or w.shape!=(E.shape[0],):
        raise ValueError("Energy and weight array sizes are inconsistent")
    if sigma<=0 or points<2 or w.sum()<=0 or np.any(w<0):
        raise ValueError("Check sigma, point count and weights")
    if energy_range is None:
        emin=float(E.min()-5*sigma)
        emax=float(E.max()+5*sigma)
    else:
        emin,emax=map(float,energy_range)
        if emax<=emin:
            raise ValueError("Invalid energy interval")
    grid=np.linspace(emin,emax,points)
    dos=np.zeros(points)
    spin_factor=2.0 if include_spin else 1.0
    w=w/w.sum()
    for k,energies_k in enumerate(E):
        dd=(grid[:,None]-energies_k[None,:])/sigma
        dos+=spin_factor*w[k]*np.exp(-0.5*dd*dd).sum(axis=1)
    dos/=sigma*np.sqrt(2*np.pi)
    return grid,dos


def average_mulliken(Es,Cs,Ss,weights,occ,basis,atoms):
    """Mulliken charges per atom with discrete k-point integration."""
    w=np.asarray(weights,dtype=float)
    w=w/w.sum()
    pop=np.zeros(len(basis))
    for k,(C,S) in enumerate(zip(Cs,Ss)):
        density=(C*occ[k][None,:])@C.conj().T
        pop+=w[k]*np.real(np.diag(density@S))
    per_atom=np.bincount([orb.atom for orb in basis],weights=pop,minlength=len(atoms))
    charge=np.array([a.valence for a in atoms],dtype=float)-per_atom
    return charge,pop


def integrated_bond_populations(C,S,H,occupation,basis,first,second):
    """Total COOP and COHP for one atomic pair at a given k point.

    COOP = 2*Re[sum_(mu in A,nu in B,n) f_n C*_{mu,n} S_{mu,nu} C_{nu,n}]
    COHP uses H instead of S. For A!=B both Hermitian directions are included.
    """
    ai=np.array([i for i,o in enumerate(basis) if o.atom==first],dtype=int)
    bj=np.array([j for j,o in enumerate(basis) if o.atom==second],dtype=int)
    if first==second or not len(ai) or not len(bj):
        raise ValueError("Provide two distinct atoms with basis functions")
    weighted_Ca=np.conj(C[ai,:])*occupation[None,:]
    Cb=C[bj,:]
    overlap_bond=2*np.real(np.einsum("in,ij,jn->", weighted_Ca,S[np.ix_(ai,bj)],Cb))
    hamilton_bond=2*np.real(np.einsum("in,ij,jn->",weighted_Ca,H[np.ix_(ai,bj)],Cb))
    return float(overlap_bond),float(hamilton_bond)


def bond_curves(C,S,H,basis,first,second):
    """COOP/COHP contributions per molecular orbital (including unoccupied orbitals).

    Returns two length-nb arrays (states not yet weighted by occupation).
    This supports smoothed spectra or occupation-weighted integrals.
    """
    ai=np.array([i for i,o in enumerate(basis) if o.atom==first],dtype=int)
    bj=np.array([j for j,o in enumerate(basis) if o.atom==second],dtype=int)
    if first==second or not len(ai) or not len(bj):
        raise ValueError("COOP/COHP requires two distinct valid atoms")
    Ca=np.conj(C[ai,:]); Cb=C[bj,:]
    coop=2*np.real(np.einsum("in,ij,jn->n",Ca,S[np.ix_(ai,bj)],Cb))
    cohp=2*np.real(np.einsum("in,ij,jn->n",Ca,H[np.ix_(ai,bj)],Cb))
    return coop,cohp


def periodic_bond_curves(C, overlap_R, coupling_matrix, basis,
                         first, second, translation, k_reduced):
    """COOP/COHP per band between A in cell 0 and B in cell T.

    S(T)_{mu,nu} as in the Fourier sum S(k)=sum_T S(T)*exp(i2pi k.T).
    Includes the conjugate Hermitian bond contribution.
    Returns COOP and COHP arrays (one value per band, before occupation).
    """
    first=int(first); second=int(second)
    if first==second and all(t==0 for t in translation):
        raise ValueError("Same-atom same-cell bond requires a different decomposition")
    T=np.asarray(translation,dtype=int)
    k=np.asarray(k_reduced,dtype=float)
    if T.shape!=k.shape:
        raise ValueError("Translation and k point must have equal dimensionality")
    ai=np.flatnonzero([o.atom==first for o in basis])
    bj=np.flatnonzero([o.atom==second for o in basis])
    if not len(ai) or not len(bj):
        raise ValueError("Invalid bond atom indices")
    phase=np.exp(2j*np.pi*np.dot(T,k))
    matrix=overlap_R[np.ix_(ai,bj)]*phase
    ca=np.conj(C[ai,:]); cb=C[bj,:]
    coop=2*np.real(np.einsum("in,ij,jn->n",ca,matrix,cb))
    hmatrix=matrix*coupling_matrix[np.ix_(ai,bj)]
    cohp=2*np.real(np.einsum("in,ij,jn->n",ca,hmatrix,cb))
    return coop,cohp
