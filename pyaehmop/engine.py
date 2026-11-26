"""Independent scientific engine: extended Hückel for molecules and bands."""
from concurrent.futures import ProcessPoolExecutor
from itertools import product
import numpy as np
from scipy.linalg import eigh
from scipy.integrate import trapezoid
from .basis import make_basis
from .overlap import overlap_block
from .spectra import (weighted_occupations, gaussian_dos,
                      average_mulliken, bond_curves, periodic_bond_curves)


def _atom_slices(basis, natoms):
    return [np.flatnonzero([o.atom==i for o in basis]) for i in range(natoms)]


def coupling(epsilon, k=1.75, weighted=True):
    """F_ij: H_ij=F_ij*S_ij; YAeHMOP weighted formula."""
    eps=np.asarray(epsilon,dtype=float)
    total=eps[:,None]+eps[None,:]
    if weighted:
        ratio=np.divide(eps[:,None]-eps[None,:],total,
                        out=np.zeros_like(total),where=np.abs(total)>1e-12)
        if np.any((np.abs(total)<=1e-12) & (~np.eye(len(eps),dtype=bool))):
            raise ValueError("Weighted formula is singular because epsilon_i + epsilon_j = 0")
        d2=ratio**2
        return .5*total*(k+d2+(1-k)*d2**2)
    return .5*k*total


def occupancies(energies, electrons, tolerance=1e-5):
    """Fill orbital levels with capacity 2 and redistribute over degenerate levels."""
    eig=np.asarray(energies)
    n=len(eig)
    if electrons<0 or electrons>2*n+1e-6:
        raise ValueError("Electron count exceeds the available orbital capacity")
    occ=np.zeros(n)
    rest=float(electrons)
    i=0
    while i<n:
        j=i+1
        while j<n and abs(eig[j]-eig[i])<=tolerance:
            j+=1
        take=min(rest,2*(j-i))
        if take>0:
            occ[i:j]=take/(j-i)
            rest-=take
        i=j
    return occ


def mulliken_populations(C,S,occ,basis,atoms):
    density=(C*occ[None,:])@C.conj().T
    pop=np.real(np.diag(density@S))
    atom_pops=np.bincount([o.atom for o in basis],weights=pop,minlength=len(atoms))
    return np.array([atom.valence for atom in atoms])-atom_pops, pop


def _eigensolver(H,S):
    H=(H+H.conj().T)/2
    S=(S+S.conj().T)/2
    # Generalized Hermitian eigenproblem: H c = E S c; no C code.
    v=eigh(S,eigvals_only=True)
    if v[0]<=1e-9:
        raise ValueError(f"Overlap matrix is not positive definite: minimum={v[0]:g}")
    E,C=eigh(H,S)
    return E,C


def _molecular_overlap_task(task):
    i,j,oa,ob,d,quadrature=task
    return i,j,overlap_block(oa,ob,d,quadrature)


def molecular(structure, quadrature=(24,20,12), workers=1):
    """Independent interatomic integrals enable multiprocessing."""
    basis=make_basis(structure.atoms)
    slices=_atom_slices(basis,len(structure.atoms))
    n=len(basis)
    S=np.eye(n)
    tasks=[]
    for i in range(len(slices)):
        for j in range(i+1,len(slices)):
            oa=tuple(basis[a] for a in slices[i])
            ob=tuple(basis[b] for b in slices[j])
            d=structure.atoms[j].position-structure.atoms[i].position
            tasks.append((i,j,oa,ob,d,quadrature))
    if workers<1:
        raise ValueError("workers must be at least 1")
    if workers>1 and len(tasks)>1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results=list(pool.map(_molecular_overlap_task,tasks))
    else:
        results=list(map(_molecular_overlap_task,tasks))
    for i,j,block in results:
        S[np.ix_(slices[i],slices[j])]=block
        S[np.ix_(slices[j],slices[i])]=block.T
    epsilon=np.array([o.energy for o in basis])
    F=coupling(epsilon,structure.k_constant,structure.weighted)
    H=F*S
    np.fill_diagonal(H,epsilon)
    eigenvalues,C=_eigensolver(H,S)
    occ=occupancies(eigenvalues,structure.electrons)
    charges,atom_orb_pops=mulliken_populations(C,S,occ,basis,structure.atoms)
    occ_ids=np.flatnonzero(occ>1e-8)
    empty_ids=np.flatnonzero(occ<1e-8)
    result={
        "type":"molecule","title":structure.title,"atoms":len(structure.atoms),
        "warnings":list(structure.warnings),
        "orbitals":n,"electrons":structure.electrons,
        "orbital_energies_eV":eigenvalues,"occupations":occ,
        "electronic_energy_eV":float(eigenvalues@occ),
        "mulliken_charges":charges,
        "orbital_populations":atom_orb_pops,
        "homo_eV":float(eigenvalues[occ_ids[-1]]) if len(occ_ids) else None,
        "lumo_eV":float(eigenvalues[empty_ids[0]]) if len(empty_ids) else None,
        "orbital_names":[f"{o.atom+1}:{o.label}" for o in basis],
        "overlap_matrix":S,"hamiltonian_matrix":H,"coefficients":C
    }
    if result["homo_eV"] is not None and result["lumo_eV"] is not None:
        result["gap_homo_lumo_eV"]=result["lumo_eV"]-result["homo_eV"]
    if structure.fmo_fragments:
        result["fmo"]=_molecular_fmo_analysis(structure,basis,S,H)
    if structure.bond_requests:
        requested=[]
        for kind,curve,a,b,translation in structure.bond_requests:
            if any(translation):
                raise NotImplementedError("Molecular COOP/COHP requires zero lattice translation")
            per_coop,per_cohp=bond_curves(C,S,H,basis,a,b)
            requested.append({"requested_type":kind,"curva":curve,"atoms":[a+1,b+1],
                              "coop_contributions":per_coop,
                              "cohp_contributions":per_cohp,
                              "integrated_coop":float(per_coop@occ),
                              "integrated_cohp":float(per_cohp@occ)})
        result["bond_analyses"]=requested
    return result


def _eigensystem_k(k,blocks,F,epsilon,k_const):
    S=np.zeros_like(next(iter(blocks.values())),dtype=np.complex128)
    for t,mat in blocks.items():
        S+=np.exp(2j*np.pi*np.dot(np.asarray(k)[:len(t)],np.array(t)))*mat
    S=(S+S.conj().T)/2
    H=F*S
    np.fill_diagonal(H,epsilon*(1-k_const+k_const*np.diag(S)))
    vals,C=_eigensolver(H,S)
    return vals,C,S,H


def _fourier_k(args):
    """Separate function to distribute k-point computations via ProcessPoolExecutor."""
    k,blocks,F,epsilon,k_const=args
    vals,_,_,_=_eigensystem_k(k,blocks,F,epsilon,k_const)
    return vals


def _overlap_task(task):
    """Independent task for parallel atom-pair integrals."""
    t,i,j,oa,ob,r,quadrature=task
    return t,i,j,overlap_block(oa,ob,r,quadrature)


def _realspace_blocks(structure,basis,quadrature,cutoff,workers=1):
    atoms=structure.atoms
    slices=_atom_slices(basis,len(atoms))
    n=len(basis)
    tasks=[]
    for t in product(*(range(-b,b+1) for b in structure.bounds)):
        shift=np.array(t,dtype=float)@structure.vectors
        for i in range(len(atoms)):
            for j in range(len(atoms)):
                r=atoms[j].position+shift-atoms[i].position
                if np.linalg.norm(r)>cutoff:
                    continue
                oa=tuple(basis[x] for x in slices[i])
                ob=tuple(basis[x] for x in slices[j])
                tasks.append((t,i,j,oa,ob,r,quadrature))
    if workers>1 and len(tasks)>4:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results=pool.map(_overlap_task,tasks,chunksize=8)
            computed=list(results)
    else:
        computed=list(map(_overlap_task,tasks))
    blocks={}
    for t,i,j,block in computed:
        mat=blocks.setdefault(t,np.zeros((n,n),dtype=float))
        mat[np.ix_(slices[i],slices[j])]=block
    blocks={t:mat for t,mat in blocks.items() if np.max(np.abs(mat))>1.e-14}
    # Numerical quadrature may break S(-T)=S(T)^T by tiny rounding errors.
    for t in list(blocks):
        opp=tuple(-x for x in t)
        m=blocks.get(opp,np.zeros((n,n)))
        avg=(blocks[t]+m.T)/2
        blocks[t]=avg
        blocks[opp]=avg.T
    return blocks


def interpolate_path(path, divisions):
    if divisions<=0:
        raise ValueError("Number of intervals per segment must be positive")
    if len(path)<2:
        raise ValueError("BAND requires at least two special k points")
    points=[]; labels=[]
    for (name1,start),(name2,end) in zip(path[:-1],path[1:]):
        for ix in range(divisions):
            p=start+(end-start)*ix/divisions
            points.append(np.asarray(p,dtype=float))
            if ix==0:
                labels.append((len(points)-1,name1))
    points.append(np.asarray(path[-1][1],dtype=float))
    labels.append((len(points)-1,path[-1][0]))
    return np.array(points),labels


def periodic(structure, quadrature=(24,20,12), workers=1, cutoff=11.0,
             kpoints=None, bands=True, dos_sigma=0.15, dos_points=600,
             compute_dos=True):
    if structure.vectors is None:
        raise ValueError("Periodic systems require lattice vectors")
    if structure.bond_requests and not structure.kpoints:
        raise ValueError("Periodic COOP/COHP requires K POINTS for property integration")
    basis=make_basis(structure.atoms)
    epsilon=np.array([o.energy for o in basis])
    F=coupling(epsilon,structure.k_constant,structure.weighted)
    blocks=_realspace_blocks(structure,basis,quadrature,cutoff,workers)
    if not blocks:
        raise ValueError("No overlap blocks could be computed")
    labels=[]
    if not bands:
        kpoints=np.zeros((0,3),dtype=float)
    elif kpoints is None:
        if structure.band_path:
            kpoints,labels=interpolate_path(structure.band_path,structure.band_divisions)
        elif structure.kpoints:
            kpoints=np.array([p for p,w in structure.kpoints])
        else:
            kpoints=np.zeros((1,3))
    kpoints=np.asarray(kpoints,dtype=float)
    if kpoints.ndim!=2 or kpoints.shape[1]!=3:
        raise ValueError("kpoints precisam ter forma (n,3)")
    args=((k,blocks,F,epsilon,structure.k_constant) for k in kpoints)
    if workers==1:
        energies=np.array([_fourier_k(arg) for arg in args])
    else:
        if workers<1:
            raise ValueError("workers must be >=1")
        with ProcessPoolExecutor(max_workers=workers) as pool:
            energies=np.array(list(pool.map(_fourier_k,args)))
    result={"type":"periodic","title":structure.title,"atoms":len(structure.atoms),
            "warnings":list(structure.warnings),
            "orbitals":len(basis),"electrons":structure.electrons,
            "kpoints":kpoints,"band_energies_eV":energies,
            "special_kpoints":labels,"total_blocos_R":len(blocks),
            "cutoff_angstrom":cutoff,"workers":workers}
    # Gap restricted to sampled k path; may exceed the fundamental gap.
    if bands and len(energies)>0 and structure.electrons>0 and structure.electrons%2==0:
        occupied=int(structure.electrons//2)
        if occupied<len(basis):
            vbm=np.max(energies[:,occupied-1])
            cbm=np.min(energies[:,occupied])
            result["band_path_gap_eV"]=float(max(0,cbm-vbm))
            result["minimum_direct_gap_eV"]=float(np.min(energies[:,occupied]-energies[:,occupied-1]))
            result["vbm_eV"]=float(vbm);result["cbm_eV"]=float(cbm)
    if structure.kpoints:
        sampled=[]
        weights=[]
        for p,w in structure.kpoints:
            sampled.append(_eigensystem_k(p,blocks,F,epsilon,structure.k_constant))
            weights.append(w)
        bands_avg=np.array([x[0] for x in sampled])
        occ,fermi=weighted_occupations(bands_avg,weights,structure.electrons)
        charges,pops=average_mulliken(
            bands_avg,[x[1] for x in sampled],[x[2] for x in sampled],
            weights,occ,basis,structure.atoms)
        wn=np.array(weights,dtype=float)/sum(weights)
        energy=float(np.sum(bands_avg*occ*wn[:,None]))
        if compute_dos:
            grid,dos=gaussian_dos(bands_avg,weights,points=dos_points,sigma=dos_sigma)
        result["sampled_properties"]={
            "total_kpoints":len(structure.kpoints),
            "fermi_eV":fermi,
            "average_electronic_energy_eV":energy,
            "mulliken_charges":charges,
            "orbital_populations":pops,
            "charge_sum":float(sum(charges))}
        if compute_dos:
            result["sampled_properties"]["dos"]={
                "energy_eV":grid,"states_per_eV":dos,
                "sigma_eV":dos_sigma,"integrated_states":float(trapezoid(dos,grid))}
        if structure.bond_requests:
            bond_results=[]
            groups={}
            for kind,curve,a,b,translation in structure.bond_requests:
                if not 0<=a<len(structure.atoms) or not 0<=b<len(structure.atoms):
                    raise ValueError("COOP/COHP atom index is outside the geometry")
                dim=len(structure.bounds)
                if any(translation[dim:]):
                    raise ValueError("COOP/COHP translation exceeds lattice dimensionality")
                t=tuple(translation[:dim])
                if any(abs(x)>bound for x,bound in zip(t,structure.bounds)):
                    raise ValueError("COOP/COHP translation exceeds LATTICE bounds")
                m=blocks.get(t,np.zeros((len(basis),len(basis))))
                per_k_coop=[];per_k_cohp=[]
                for (p,_),eigvecs in zip(structure.kpoints,sampled):
                    C=eigvecs[1]
                    c1,c2=periodic_bond_curves(C,m,F,basis,a,b,t,p[:dim])
                    per_k_coop.append(c1)
                    per_k_cohp.append(c2)
                co=np.asarray(per_k_coop)
                hp=np.asarray(per_k_cohp)
                integrated_coop=float(np.sum(co*occ*wn[:,None]))
                integrated_cohp=float(np.sum(hp*occ*wn[:,None]))
                bond_results.append({
                    "requested_type":kind,"curva":curve,
                    "atoms":[a+1,b+1],"translation":translation,
                    "coop_contributions_per_k_and_band":co,
                    "cohp_contributions_per_k_and_band":hp,
                    "integrated_coop":integrated_coop,
                    "integrated_cohp":integrated_cohp})
                group=groups.setdefault((kind,curve),[])
                group.append((integrated_coop,integrated_cohp))
            result["bond_analyses"]=bond_results
            result["bond_groups"]=[
                {"requested_type":kind,"curva":curve,"count":len(vals),
                 "average_integrated_coop":float(np.mean([v[0] for v in vals])),
                 "average_integrated_cohp":float(np.mean([v[1] for v in vals]))}
                for (kind,curve),vals in sorted(groups.items())]
    return result


def calculate(structure,quadrature=(24,20,12),workers=1,cutoff=11.0,
              kpoints=None, bands=True, dos_sigma=0.15, dos_points=600,
              compute_dos=True):
    if structure.vectors is None:
        return molecular(structure,quadrature,workers)
    return periodic(structure,quadrature,workers,cutoff,kpoints,bands,
                    dos_sigma,dos_points,compute_dos)


def _molecular_fmo_analysis(structure,basis,S,H):
    """Independent fragment diagonalization and FMO basis transformation.

    Each fragment is diagonalized using its own generalized eigenproblem.
    The transformation columns represent fragment orbitals expressed in the
    molecular AO basis. This never invokes YAeHMOP C.
    """
    n=len(basis)
    transformation=np.zeros((n,n))
    fragments=[]
    cursor=0
    for number,(atom_numbers,fragment_electrons) in enumerate(structure.fmo_fragments,1):
        included={x-1 for x in atom_numbers}
        ids=np.array([ix for ix,orb in enumerate(basis) if orb.atom in included])
        if len(ids)==0:
            raise ValueError("Empty FMO fragment or fragment without orbitals")
        subS=S[np.ix_(ids,ids)]
        subH=H[np.ix_(ids,ids)]
        eigenvals,coeff=_eigensolver(subH,subS)
        occ=occupancies(eigenvals,fragment_electrons)
        transformation[np.ix_(ids,np.arange(cursor,cursor+len(ids)))]=coeff
        cursor+=len(ids)
        fragments.append({
            "fragment":number,"atoms":atom_numbers,
            "electrons":fragment_electrons,"orbital_count":len(ids),
            "orbital_energies_eV":eigenvals,"occupations":occ,
            "fragment_electronic_energy_eV":float(eigenvals@occ),
            "global_basis_orbitals":(ids+1).tolist(),
            "local_ao_coefficients":coeff
        })
    if cursor!=n:
        raise ValueError("FMO fragments do not cover the full orbital basis")
    return {"fragments":fragments,
            "fragment_energy_sum_eV":float(sum(x["fragment_electronic_energy_eV"] for x in fragments)),
            "ao_to_fmo_transform":transformation,
            "fmo_overlap_matrix":transformation.T@S@transformation,
            "fmo_hamiltonian_matrix":transformation.T@H@transformation}
