"""Collinear Stoner SCF plus atomic spinor SOC in a nonorthogonal EHT basis.

Independent, explicitly semiempirical prototype. The current SCF updates
only z magnetization in selected shells (no noncollinear spin SCF).
"""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial
import numpy as np
from scipy.linalg import eigh
from scipy.special import expit
from pyaehmop.basis import make_basis
from pyaehmop.engine import coupling, _realspace_blocks, molecular
from .operators import soc_matrix


@dataclass
class SpinSettings:
    spin: bool = True
    soc: bool = False
    stoner: dict | None = None
    xi: dict | None = None
    initial_moment: float = 2.0
    mixing: float = 0.3
    temperature: float = 0.025
    tolerance: float = 1e-5
    max_iter: int = 130
    workers: int = 1

    def validate(self):
        if not (self.spin or self.soc):
            raise ValueError("Choose spin = yes, soc = yes, or both")
        if not 0 < self.mixing <= 1:
            raise ValueError("mixing must belong to (0,1]")
        if self.temperature <= 0:
            raise ValueError("temperature_eV must be positive")
        if self.tolerance <= 0 or self.max_iter < 1:
            raise ValueError("tolerance positiva e max_iter>=1 necessarios")
        if self.workers<1:
            raise ValueError("workers must be at least 1")
        if self.spin and not any(v > 0 for shells in (self.stoner or {}).values()
                                 for v in shells.values()):
            raise ValueError("Self-consistent spin requires positive Stoner I for at least one shell")
        if self.soc and not any(abs(v)>0 for shells in (self.xi or {}).values()
                                for v in shells.values()):
            raise ValueError("SOC requires nonzero xi for at least one shell")


def fermi_dirac(energies,weights,ne,temp):
    """Fermi-Dirac occupations: maximum 1 per spinor; normalized k weights.

    SCF uses a single Fermi energy for both channels (N_up/N_down free).
    """
    energies=np.asarray(energies,dtype=float)
    weights=np.asarray(weights,dtype=float)
    if energies.ndim!=2 or weights.shape!=(len(energies),):
        raise ValueError("Energias/weights k invalidos")
    if np.any(weights<0) or weights.sum()<=0:
        raise ValueError("Pesos k precisam ser no negativos")
    nw=weights/weights.sum()
    n=energies.shape[1]
    if ne<0 or ne>n+1e-8:
        raise ValueError("Electron number exceeds spinor basis capacity")
    if ne==0:
        return np.zeros_like(energies),float(energies.min()-100*temp)
    if ne>=n:
        return np.ones_like(energies),float(energies.max()+100*temp)
    lo=energies.min()-150*temp
    hi=energies.max()+150*temp
    for _ in range(160):
        mu=(lo+hi)/2
        occ=expit((mu-energies)/temp)
        got=np.sum(nw[:,None]*occ)
        if got < ne:lo=mu
        else:hi=mu
    mu=(lo+hi)/2
    return expit((mu-energies)/temp),float(mu)


def bloch_matrices(k,blocks,F,eps,kconst):
    """S(k),H0(k) EHT conservando Hermiticidade."""
    n=len(eps)
    S=np.zeros((n,n),complex)
    for R,mat in blocks.items():
        kr=np.dot(np.asarray(k)[:len(R)],np.asarray(R))
        S+=np.exp(2j*np.pi*kr)*mat
    S=(S+S.conj().T)/2
    H=F*S
    np.fill_diagonal(H,eps*(1-kconst+kconst*np.real(np.diag(S))))
    H=(H+H.conj().T)/2
    return H,S


def assemble_reference(structure,kpoints,quadrature=(24,20,12),
                       cutoff=5.6,workers=1):
    """Reuses PyYAeHMOP Python integrals without touching original files."""
    basis=make_basis(structure.atoms)
    kpts=np.asarray(kpoints,dtype=float)
    if kpts.ndim!=2 or kpts.shape[1]!=3:
        raise ValueError("k points must have shape (nk,3)")
    if structure.vectors is None:
        if len(kpts)!=1 or np.linalg.norm(kpts[0])>1e-12:
            raise ValueError("Molecula admite only Gamma")
        initial=molecular(structure,quadrature,workers)
        H0=[np.asarray(initial["hamiltonian_matrix"],complex)]
        S0=[np.asarray(initial["overlap_matrix"],complex)]
        return basis,H0,S0,None
    eps=np.array([x.energy for x in basis],float)
    F=coupling(eps,structure.k_constant,structure.weighted)
    blocks=_realspace_blocks(structure,basis,quadrature,cutoff,workers)
    if not blocks:
        raise ValueError("No overlap integrals within the selected range")
    pairs=[bloch_matrices(k,blocks,F,eps,structure.k_constant) for k in kpts]
    H0=[x[0] for x in pairs]
    S0=[x[1] for x in pairs]
    for k,S in enumerate(S0):
        eigmin=np.linalg.eigvalsh(S)[0]
        if eigmin<1e-8:
            raise ValueError(f"S(k) matrix is not positive definite at {k}: {eigmin:g}; "
                             "adjust cutoff/quadrature/geometry")
    return basis,H0,S0,(blocks,F,eps)


def stoner_orbital_I(basis,atoms,species_shell_I):
    """Stoner vector per AO, grouped by (atom,l)."""
    strength=np.zeros(len(basis))
    groups=[]
    for i,orb in enumerate(basis):
        strength[i]=(species_shell_I.get(atoms[orb.atom].symbol,{}).
                     get("spdf"[orb.l],0.))
    for atom in range(len(atoms)):
        for l in range(4):
            ids=np.array([i for i,o in enumerate(basis) if
                          o.atom==atom and o.l==l],dtype=int)
            if len(ids)>0 and np.any(strength[ids]>0):
                groups.append((atom,l,ids,float(strength[ids[0]])))
    return strength,groups


def spinor_problem(H0,S,shifts,soc):
    """Stoner potential symmetrized in the nonorthogonal metric:
       D_mu = I_(A,l) * m_(A,l) / 2;
       H_up = H0 - (S D + D S)/2, H_down = H0 + (S D + D S)/2.
       Acrescenta SOC local em espaco spinorial.
    """
    n=len(S)
    diag=np.asarray(shifts,dtype=float)
    sym=(S*diag[None,:]+diag[:,None]*S)/2
    H=np.zeros((2*n,2*n),complex)
    # Ferromagnetic Stoner exchange: majority spin up has LOWER energy.
    H[:n,:n]=H0-sym
    H[n:,n:]=H0+sym
    H+=soc
    Sb=np.zeros((2*n,2*n),complex)
    Sb[:n,:n]=S
    Sb[n:,n:]=S
    return (H+H.conj().T)/2,Sb


def solve_one(args):
    H0,S,shifts,soc=args
    H,Sb=spinor_problem(H0,S,shifts,soc)
    values,C=eigh(H,Sb,check_finite=False)
    n=len(S)
    up=np.real(np.einsum("in,in->n",C[:n].conj(),S@C[:n]))
    dn=np.real(np.einsum("in,in->n",C[n:].conj(),S@C[n:]))
    sz=up-dn
    return values,C,np.clip(sz,-1,1)


def spin_populations(solution,occs,weights,overlaps,basis,atoms):
    """Mulliken spin up/down, incluindo espinores SOC."""
    w=np.asarray(weights,dtype=float);w=w/w.sum()
    n=len(basis)
    nup=np.zeros(n); ndn=np.zeros(n)
    for kw,(values,C,sz),occ,S in zip(w,solution,occs,overlaps):
        Puu=(C[:n,:]*occ[None,:])@C[:n,:].conj().T
        Pdd=(C[n:,:]*occ[None,:])@C[n:,:].conj().T
        nup+=kw*np.real(np.diag(Puu@S))
        ndn+=kw*np.real(np.diag(Pdd@S))
    natom=len(atoms)
    total=np.bincount([o.atom for o in basis],weights=nup+ndn,minlength=natom)
    mz=np.bincount([o.atom for o in basis],weights=nup-ndn,minlength=natom)
    charges=np.array([a.valence for a in atoms])-total
    return nup,ndn,mz,charges


def scf(structure,kpoints,weights,settings,quadrature=(24,20,12),
        cutoff=5.6,verbose=False):
    """Self-consistent SCF for m_{A,l}=n_up-n_down with linear mixing.

    With SOC: diagonalize spinors and update only the local
    z component. The transverse exchange is not optimized (constrained collinear).
    """
    settings.validate()
    basis,h0list,slist,data=assemble_reference(
        structure,kpoints,quadrature,cutoff,settings.workers)
    weights=np.asarray(weights,float)
    if np.any(weights<0) or weights.sum()<=0 or len(weights)!=len(kpoints):
        raise ValueError("Pesos k invalidos")
    weights=weights/weights.sum()
    ne=float(structure.electrons)
    if ne<0 or ne>2*len(basis):
        raise ValueError("Electron count is incompatible with the basis")
    strengths,groups=stoner_orbital_I(basis,structure.atoms,settings.stoner or {})
    symbols=[atom.symbol for atom in structure.atoms]
    soc=soc_matrix(basis,symbols,settings.xi or {}) if settings.soc else np.zeros((2*len(basis),)*2,complex)
    m=np.array([settings.initial_moment if settings.spin else 0.
                for _ in groups],float)
    last_m=np.array(m)
    history=[]
    converged=not settings.spin
    iters=settings.max_iter if settings.spin else 1

    def shifts_for(moment):
        shifts=np.zeros(len(basis))
        for val,(_,_,ids,I) in zip(moment,groups):
            shifts[ids]+=I*val/2
        return shifts

    with ProcessPoolExecutor(max_workers=settings.workers) if settings.workers>1 else _NullPool() as pool:
        for iteration in range(1,iters+1):
            shifts=shifts_for(m)
            args=[(H,S,shifts,soc) for H,S in zip(h0list,slist)]
            solutions=list(pool.map(solve_one,args)) if settings.workers>1 else list(map(solve_one,args))
            eig=np.array([x[0] for x in solutions])
            occ,mu=fermi_dirac(eig,weights,ne,settings.temperature)
            nu,nd,moments,charges=spin_populations(solutions,occ,weights,
                                                   slist,basis,structure.atoms)
            output_m=np.array([np.sum((nu-nd)[ids]) for _,_,ids,I in groups])
            err=float(np.max(abs(output_m-m))) if len(m) else 0.
            history.append({"iter":iteration,"spin_residual":err,
                            "total_moment_muB":float(sum(moments)),
                            "fermi_eV":mu})
            if verbose and (iteration<4 or iteration%10==0 or err<settings.tolerance):
                print(f"SCF iter={iteration:3d} spin_error={err:.3e} total_moment={sum(moments):.6f} muB")
            if err<settings.tolerance:
                converged=True
                m=output_m
                break
            m=m+settings.mixing*(output_m-m)
        # Recompute at final input m (ensures consistent Hamiltonian and populations).
        shifts=shifts_for(m)
        args=[(H,S,shifts,soc) for H,S in zip(h0list,slist)]
        solutions=list(pool.map(solve_one,args)) if settings.workers>1 else list(map(solve_one,args))
    eig=np.array([x[0] for x in solutions])
    occ,mu=fermi_dirac(eig,weights,ne,settings.temperature)
    nu,nd,moments,charges=spin_populations(solutions,occ,weights,slist,
                                           basis,structure.atoms)
    actual_m=np.array([np.sum((nu-nd)[ids]) for _,_,ids,I in groups])
    residual=float(max(abs(actual_m-m))) if len(m) else 0.
    converged=bool(converged and residual<settings.tolerance*2)
    band_sum=float(np.sum(weights[:,None]*eig*occ))
    doublecount=float(sum(I*x*x/4 for x,(_,_,ids,I) in zip(actual_m,groups)))
    result={
        "converged":converged, "iterations":len(history), "final_residual":residual,
        "history":history, "fermi_eV":mu,
        "electrons":ne, "band_energy_eV":band_sum,
        "stoner_double_counting_correction_eV":doublecount,
        "stoner_model_energy_eV":band_sum+doublecount,
        "atomic_moments_muB":moments,
        "total_moment_muB":float(sum(moments)),
        "stoner_moments_muB":actual_m,
        "atomic_charges":charges,
        "orbital_populations_up":nu, "orbital_populations_down":nd,
        "integration_energies_eV":eig,
        "sz_integration":np.array([x[2] for x in solutions]),
        "occupations":occ,"integration_kpoints":np.array(kpoints),
        "k_weights":weights,"n_orbitals":len(basis),
        "soc_max_element_eV":float(np.max(abs(soc))),
        "spin":settings.spin,"soc":settings.soc
    }
    reference=(basis,data,m,soc,settings)
    return result,reference


class _NullPool:
    def __enter__(self):return self
    def __exit__(self,*args):return False


def evaluate_bands(structure,kpoints,reference,quadrature=(24,20,12),cutoff=5.6):
    """K-resolved bands at frozen self-consistent moments, without mixed occupations."""
    basis,data,m,soc,settings=reference
    if structure.vectors is None:
        raise ValueError("Molecular systems do not have a band path")
    blocks,F,eps=data
    strengths,groups=stoner_orbital_I(basis,structure.atoms,settings.stoner or {})
    shifts=np.zeros(len(basis))
    for val,(_,_,ids,I) in zip(m,groups):
        shifts[ids]+=I*val/2
    args=[]
    for k in kpoints:
        H,S=bloch_matrices(k,blocks,F,eps,structure.k_constant)
        args.append((H,S,shifts,soc))
    if settings.workers>1:
        with ProcessPoolExecutor(max_workers=settings.workers) as pool:
            eigvecs=list(pool.map(solve_one,args))
    else:
        eigvecs=list(map(solve_one,args))
    return np.array([x[0] for x in eigvecs]),np.array([x[2] for x in eigvecs])


def spin_dos(eig,sz,weights,sigma=0.14,points=550):
    """Spin-z projected DOS with up weight (1+<sig_z>)/2."""
    eig=np.asarray(eig);sz=np.asarray(sz)
    weights=np.asarray(weights)/np.sum(weights)
    emin=float(np.min(eig)-5*sigma);emax=float(np.max(eig)+5*sigma)
    grid=np.linspace(emin,emax,points)
    up=np.zeros(points);dn=np.zeros(points)
    for k,(vals,polar) in enumerate(zip(eig,sz)):
        kernel=np.exp(-.5*((grid[:,None]-vals[None,:])/sigma)**2)
        up+=weights[k]*np.sum(kernel*(1+polar)[None,:]/2,axis=1)
        dn+=weights[k]*np.sum(kernel*(1-polar)[None,:]/2,axis=1)
    prefactor=1/(np.sqrt(2*np.pi)*sigma)
    return grid,up*prefactor,dn*prefactor
