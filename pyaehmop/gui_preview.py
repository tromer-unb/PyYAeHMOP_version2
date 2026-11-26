"""Read structures for 3D GUI preview without running an EHT calculation."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np


@dataclass
class StructureView:
    symbols:list
    positions:np.ndarray
    vectors:np.ndarray
    title:str


def load_structure_preview(source,kind,ase_format="",frame=0):
    path=Path(source).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"Structure file not found: {path}")
    if kind=="yaehmop":
        from .input import read_input
        structure=read_input(path)
        symbols=[a.symbol.capitalize() for a in structure.atoms]
        positions=np.asarray([a.position for a in structure.atoms],dtype=float)
        vectors=(np.asarray(structure.vectors,dtype=float)
                 if structure.vectors is not None else np.empty((0,3)))
        return StructureView(symbols,positions,vectors,structure.title)
    if kind not in ("xyz","cif","ase"):
        raise ValueError("Unsupported structure format")
    try:
        from ase.io import read
    except ImportError as exc:
        raise ImportError("Structure preview requires ASE (pip install ase)") from exc
    atoms=read(str(path),format=("cif" if kind=="cif" else (ase_format or None)),
               index=int(frame))
    symbols=atoms.get_chemical_symbols()
    cell=np.asarray(atoms.cell,dtype=float)
    periodic=np.asarray(atoms.pbc,dtype=bool)
    vectors=cell[periodic] if periodic.any() else np.empty((0,3))
    return StructureView(symbols,np.asarray(atoms.positions,dtype=float),
                         vectors,path.stem)
