"""Parser for the documented YAeHMOP input subset; unsupported modes fail explicitly."""
from dataclasses import dataclass, field
from pathlib import Path
import re
import numpy as np
from .basis import Atom, read_parameter_table, shells_from_custom
from .geometry import zmat_to_cartesian
from .fragments import parse_integer_ranges

@dataclass
class Structure:
    title: str
    atoms: list[Atom]
    electrons: float
    vectors: np.ndarray | None
    bounds: tuple[int, ...] | None
    band_path: list[tuple[str,np.ndarray]]
    band_divisions: int
    kpoints: list[tuple[np.ndarray, float]]
    weighted: bool = True
    k_constant: float = 1.75
    bond_requests: list = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    fmo_fragments: list = field(default_factory=list)

def _crystal_basis(lengths, angles):
    a,b,c = map(float,lengths)
    alpha,beta,gamma = np.deg2rad(angles)
    cg,cb,ca = np.cos(gamma),np.cos(beta),np.cos(alpha)
    sg = np.sin(gamma)
    if abs(sg)<1e-12:
        raise ValueError("Singular gamma angle")
    vz2 = 1-cb**2-((ca-cb*cg)/sg)**2
    if vz2 < -1e-10:
        raise ValueError("Invalid crystal cell parameters")
    return np.array([[a,0,0],[b*cg,b*sg,0],
                     [c*cb,c*(ca-cb*cg)/sg,c*np.sqrt(max(0,vz2))]])

def read_input(path):
    """Read Cartesian/crystallographic geometry, parameters, lattice, bands and electrons.

    Does not implement Walsh diagrams or charge iteration.
    """
    raw = Path(path).read_text(encoding="latin-1")
    lines = [line.split(";",1)[0].strip() for line in raw.splitlines()]
    lines = [line for line in lines if line]
    if not lines:
        raise ValueError("Empty input file")
    title = lines[0]
    table = {key:dict(value) for key,value in read_parameter_table().items()}
    custom = {}
    atoms_data = []
    coordinate_type = "cartesian"
    crystal_lengths = crystal_angles = None
    lattice_data = None
    is_molecular = False
    electrons = None
    charge = None
    band_path = []
    band_divisions = 0
    kpoints = []
    weighted = True; k_constant = 1.75
    bond_requests=[]
    warnings=[]
    fmo_fragments=[]
    i = 1
    while i < len(lines):
        line = lines[i]
        words = line.lower().replace("_"," ").split()
        if not words:
            i += 1; continue
        key = words[0]
        if key == "geometry":
            zmat = ("z" in words[1:] and "matrix" in words[1:]) or (
                "z-matrix" in line.lower() or "zmatrix" in line.lower())
            coordinate_type = ("zmatrix" if zmat else
                ("crystallographic" if len(words)>1 and words[1].startswith("cryst") else "cartesian"))
            count = int(lines[i+1].split()[0])
            records=lines[i+2:i+2+count]
            if len(records)!=count:
                raise ValueError("Number of atoms is smaller than specified by GEOMETRY")
            if zmat:
                atoms_data=zmat_to_cartesian(records)
            else:
                atoms_data=[]
                for entry in records:
                    v=entry.split()
                    if len(v)<5:
                        raise ValueError("Invalid atom record: "+entry)
                    atoms_data.append((v[1].upper(),np.array([float(t) for t in v[2:5]])))
            i += 2+count;continue
        if key == "parameters":
            i+=1
            while i<len(lines) and re.match(r"^[A-Za-z]{1,2}\s+\d+\s+[-+]?\d",lines[i]):
                sym,info=shells_from_custom(lines[i])
                custom[sym]=info
                i+=1
            continue
        if key=="molecular":
            is_molecular = True
        elif key=="electrons" and len(words)==1:
            electrons=float(lines[i+1].split()[0]);i+=1
        elif key=="charge" and len(words)==1:
            charge=float(lines[i+1].split()[0]);i+=1
        elif key=="crystal" and len(words)>1 and words[1].startswith("spec"):
            crystal_lengths = list(map(float,lines[i+1].split()[:3]))
            crystal_angles = list(map(float,lines[i+2].split()[:3]))
            i+=2
        elif key=="lattice":
            dim = int(lines[i+1].split()[0])
            if not 1<=dim<=3:
                raise ValueError("Lattice dimension must be 1, 2 or 3")
            bounds=tuple(map(int,lines[i+2].split()[:dim]))
            if len(bounds)!=dim:
                raise ValueError("Insufficient translation bounds")
            pairs=[]
            for entry in lines[i+3:i+3+dim]:
                v=entry.split()
                pairs.append((int(v[0])-1,int(v[1])-1))
            lattice_data=(dim,bounds,pairs)
            i+=2+dim
        elif key in ("band","bands"):
            band_divisions=int(lines[i+1].split()[0])
            count=int(lines[i+2].split()[0])
            band_path=[]
            for entry in lines[i+3:i+3+count]:
                v=entry.split()
                band_path.append((v[0],np.array(list(map(float,v[1:4])))))
            i+=2+count
        elif key=="k" and len(words)>1 and words[1]=="points":
            count=int(lines[i+1].split()[0])
            kpoints=[]
            for entry in lines[i+2:i+2+count]:
                v=list(map(float,entry.split()))
                if len(v)<3:
                    raise ValueError("Invalid k point")
                kpoints.append((np.array(v[:3]),v[3] if len(v)>3 else 1.))
            i+=1+count
        elif key in ("coop","cohp") and len(words)==1:
            count=int(lines[i+1].split()[0])
            for entry in lines[i+2:i+2+count]:
                v=entry.lower().split()
                if len(v)<7 or v[0] not in ("atom","h-atom"):
                    raise NotImplementedError("COOP/COHP: only atom and h-atom types are supported")
                curve=int(v[1]);a=int(v[2])-1;b=int(v[3])-1
                translation=tuple(map(int,v[4:7]))
                bond_kind="cohp" if key=="cohp" or v[0]=="h-atom" else "coop"
                bond_requests.append((bond_kind,curve,a,b,translation))
            i+=1+count
        elif key=="nonweighted":
            weighted=False
        elif key=="the" and len(words)>1 and words[1]=="const":
            k_constant=float(words[2]) if len(words)>2 else float(lines[i+1].split()[0])
            if len(words)<=2:
                i+=1
        elif key=="k" and len(words)>1 and words[1]=="value":
            k_constant=float(lines[i+1].split()[0]);i+=1
        elif key in ("walsh","electrostatics","fco") or (key=="charge" and len(words)>1 and words[1]=="iteration"):
            raise NotImplementedError("Unsupported feature that would affect the calculation: "+line)
        elif key=="fmo" and len(words)==1:
            nfrag=int(lines[i+1].split()[0])
            if nfrag<1:
                raise ValueError("FMO requires at least one fragment")
            electron_numbers=[float(x) for x in lines[i+2].replace(" ","").split(",")]
            if len(electron_numbers)!=nfrag:
                raise ValueError("FMO electron counts do not match the fragments")
            fragment_lines=lines[i+3:i+3+nfrag]
            if len(fragment_lines)!=nfrag:
                raise ValueError("Insufficient FMO fragment lines")
            fmo_fragments=[(parse_integer_ranges(x),electron_numbers[k])
                           for k,x in enumerate(fragment_lines)]
            i+=2+nfrag
        elif key=="mo" and len(words)>1 and words[1]=="print":
            warnings.append("MO PRINT is not implemented; full coefficients are available in JSON")
        i+=1

    table.update(custom)
    crystal = _crystal_basis(crystal_lengths,crystal_angles) if crystal_lengths is not None else None
    if coordinate_type=="crystallographic" and crystal is None:
        raise ValueError("CRYSTAL SPEC is required for crystallographic geometry")
    positions = [(sym,np.asarray(x)@crystal if coordinate_type=="crystallographic" else x)
                 for sym,x in atoms_data]
    vectors = None
    bounds = None
    if lattice_data is not None:
        dim,bounds,pairs=lattice_data
        if len(positions)<dim+1:
            raise ValueError("The geometry contains fewer atoms than the lattice specification")
        vectors=[]
        for src,dst in pairs:
            if min(src,dst)<0 or max(src,dst)>=len(positions):
                raise ValueError("Invalid lattice vector index")
            vectors.append(positions[dst][1]-positions[src][1])
        vectors=np.asarray(vectors)
        positions=positions[:-dim]  # Last 'dim' atom records are lattice translation endpoints
    elif not is_molecular:
        raise NotImplementedError("No mode defined; expected 'molecular' or 'lattice'")
    atoms=[]
    custom_order=iter(custom)
    for symbol,position in positions:
        if symbol=="&":
            continue  # Dummy atom: used for geometry references, not for electronic basis functions
        if symbol=="*":
            try:
                # Each * consumes the next PARAMETERS entry.
                symbol=next(custom_order)
            except StopIteration as exc:
                raise ValueError("Number of * symbols exceeds the supplied parameter entries") from exc
        else:
            symbol=re.match(r"[A-Za-z]+",symbol).group().upper()
            if symbol not in table:
                # Symbols such as Te2 and Br1: letters only
                raise ValueError("Missing parameters for "+symbol)
        param=table[symbol]
        atoms.append(Atom(symbol,np.asarray(position,dtype=float),
                          param["valence"],tuple(param["shells"])))
    if not atoms:
        raise ValueError("No physical atoms were defined")
    if electrons is None:
        electrons=sum(a.valence for a in atoms) - (0 if charge is None else charge)
    if electrons<0:
        raise ValueError("Negative electron count")
    if fmo_fragments:
        if lattice_data is not None:
            raise NotImplementedError("Periodic FMO is not supported; use FCO in the original software")
        atom_ids=[i for frag,n in fmo_fragments for i in frag]
        if sorted(atom_ids)!=list(range(1,len(atoms)+1)):
            raise ValueError("FMO must partition the physical atoms without duplication")
    return Structure(title, atoms, electrons,vectors,bounds,band_path,
                     band_divisions,kpoints,weighted,k_constant,bond_requests,warnings,fmo_fragments)
