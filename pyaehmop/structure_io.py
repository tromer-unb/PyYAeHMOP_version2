"""XYZ, CIF and ASE structure readers for an independent EHT engine.

Standard XYZ files are parsed natively. CIF and other formats use the
optional ASE library, never YAeHMOP C.
"""
from itertools import product
from math import ceil
from pathlib import Path

import numpy as np

from .basis import Atom, read_parameter_table
from .input import Structure, read_input


def _species_atoms(symbols, positions, parameter_table=None):
    table = read_parameter_table()
    if parameter_table:
        table={**table,**parameter_table}
    elements = []
    for name, pos in zip(symbols, positions):
        symbol = str(name).strip().upper()
        if symbol not in table or not table[symbol]["shells"]:
            raise ValueError(
                f"Element {symbol!r} lacks orbital parameters in eht_parms.dat. "
                "Provide valid EHT parameters for this element."
            )
        params = table[symbol]
        elements.append(Atom(symbol, np.array(pos,dtype=float),
                             float(params["valence"]),tuple(params["shells"])))
    if not elements:
        raise ValueError("The structure contains no atoms")
    return elements


def read_xyz(path,parameter_table=None):
    """Read standard Cartesian XYZ, one frame; lengths in angstroms."""
    source=Path(path)
    lines=source.read_text(encoding="utf-8").splitlines()
    if len(lines)<2:
        raise ValueError(f"XYZ incompleto: {source}")
    try:
        count=int(lines[0].strip())
    except ValueError as exc:
        raise ValueError("XYZ file must start with the number of atoms") from exc
    if count<=0 or len(lines)<count+2:
        raise ValueError("XYZ coordinate count does not match the header")
    title=lines[1].strip() or source.stem
    symbols=[]; positions=[]
    for entry in lines[2:2+count]:
        v=entry.split()
        if len(v)<4:
            raise ValueError("XYZ: expected atom format 'C x y z': "+entry)
        symbols.append(v[0])
        try:
            positions.append([float(x) for x in v[1:4]])
        except ValueError as exc:
            raise ValueError("XYZ: invalid Cartesian coordinate: "+entry) from exc
    if any(x.strip() for x in lines[2+count:]):
        raise ValueError("Multiple XYZ frames: use 'ase' and select a frame")
    atoms=_species_atoms(symbols,positions,parameter_table)
    electrons=sum(a.valence for a in atoms)
    return Structure(title,atoms,electrons,None,None,[],0,[])


def _vector_bounds(vectors, cutoff, *, max_auto=14):
    """Cover all lattice images that may have Cartesian distances <= cutoff.

    The dual-basis norm gives a safe upper bound after accounting for
    fractional position differences in [0,1). For materials with very
    elongated/skewed cells, require explicit images instead of truncating.
    """
    vectors=np.asarray(vectors,dtype=float)
    dual=np.linalg.pinv(vectors.T)
    extents=np.ceil(cutoff*np.linalg.norm(dual,axis=1)).astype(int)+1
    if any(extents>max_auto):
        raise ValueError("Cell is too skewed/small for images=auto; "
                         f"requires at least {tuple(extents)}. "
                         "Specify images explicitly and check convergence")
    return tuple(int(x) for x in extents)


def from_ase_atoms(ase_atoms, title, cutoff=11.0, images="auto",parameter_table=None):
    """Convert ase.Atoms into Structure, preserving cell and atomic sites."""
    pbc=np.asarray(ase_atoms.pbc,dtype=bool)
    coords=ase_atoms.get_positions()
    symbols=ase_atoms.get_chemical_symbols()
    atoms=_species_atoms(symbols,coords,parameter_table)
    electrons=sum(a.valence for a in atoms)
    if not np.any(pbc):
        return Structure(title,atoms,electrons,None,None,[],0,[])
    vectors=np.asarray(ase_atoms.cell,dtype=float)[pbc]
    if np.linalg.matrix_rank(vectors)<len(vectors):
        raise ValueError("Periodic cell is singular or rank deficient")
    if images=="auto":
        bounds=_vector_bounds(vectors,cutoff)
    else:
        bounds=tuple(int(x) for x in images)
        if len(bounds)!=len(vectors) or any(x<1 for x in bounds):
            raise ValueError("images: specify one positive integer per periodic direction")
    translations=int(np.prod([2*b+1 for b in bounds]))
    if translations>120000:
        raise ValueError(f"Too many periodic images ({translations}); adjust cutoff or images")
    return Structure(title,atoms,electrons,vectors,bounds,[],0,[])


def read_ase_structure(path, kind, cutoff=11.0, images="auto", ase_format=None,
                       index=0,parameter_table=None):
    try:
        from ase.io import read
    except ImportError as exc:
        raise ImportError("CIF/ASE exige a biblioteca ASE. "
                          "Install with: python3 -m pip install ase") from exc
    path=Path(path)
    fmt="cif" if kind=="cif" else ase_format
    try:
        atoms=read(str(path),format=fmt,index=index)
    except Exception as exc:
        raise ValueError(f"ASE could not read {path}: {exc}") from exc
    if isinstance(atoms,list):
        raise ValueError("Select exactly one frame for the calculation")
    # Reject fractional crystallographic occupancies, which EHT does not handle.
    occupancies=atoms.info.get("occupancy",{})
    if isinstance(occupancies,dict):
        for site,comps in occupancies.items():
            if isinstance(comps,dict):
                values=[float(v) for v in comps.values()]
                if len(values)!=1 or abs(values[0]-1)>1e-6:
                    raise ValueError("Disordered/partially occupied CIF is not supported: "
                                     f"site {site!r}: {comps}")
    return from_ase_atoms(atoms,path.stem,cutoff,images,parameter_table)


def build_kmesh(mesh,dim,shift=None):
    """Simplified Gamma-inclusive uniform k mesh with equal weights.

    k in reduced reciprocal-lattice coordinates.
    """
    mesh=tuple(int(x) for x in mesh)
    if len(mesh)!=dim or any(x<=0 for x in mesh):
        raise ValueError(f"kmesh must contain {dim} inteiros positivos")
    if np.prod(mesh)>8000:
        raise ValueError("K mesh is too large for interactive execution (>8000 points)")
    if shift is None:
        shift=(0.,)*dim
    if len(shift)!=dim:
        raise ValueError("k_shift has the wrong dimension")
    results=[]
    for parts in product(*(range(n) for n in mesh)):
        k=np.zeros(3,dtype=float)
        k[:dim]=[(i+float(s))/n for i,n,s in zip(parts,mesh,shift)]
        results.append((k,1.0))
    return results


def bandpath_automatic(vectors, npoints=80):
    """High-symmetry band path suggested by ASE for the cell."""
    try:
        from ase.cell import Cell
    except ImportError as exc:
        raise ImportError("band_path=auto requires ASE") from exc
    cell=np.zeros((3,3))
    cell[:len(vectors)]=vectors
    path=Cell(cell).bandpath(npoints=npoints)
    # ASE uses the cell's reciprocal coordinates.
    kpts=np.array(path.kpts,dtype=float)
    dist, ticks, labels=path.get_linear_kpoint_axis()
    indices=[int(np.argmin(abs(dist-t))) for t in ticks]
    return kpts,list(zip(indices,labels))
