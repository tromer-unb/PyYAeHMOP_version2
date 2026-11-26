"""Friendly SPIN_TEST input parser; preserves the PyYAeHMOP API."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from pyaehmop.input import read_input
from pyaehmop.structure_io import read_xyz,read_ase_structure,build_kmesh
from pyaehmop.engine import interpolate_path
from pyaehmop.friendly import _k_path, _numbers, _yesno
from .solver import SpinSettings


@dataclass
class SpinJob:
    path:Path
    source:Path
    structure:object
    settings:SpinSettings
    kpoints:np.ndarray
    weights:np.ndarray
    band_path:np.ndarray | None
    special_labels:list
    quadrature:tuple
    cutoff:float
    outdir:Path
    dos_sigma:float
    dos_points:int
    image:bool


def read_spin_job(path):
    """First lines: format, file; remaining lines are 'key = value'.

    Example:
        ase
        Fe_wire.extxyz
        spin = yes
        soc = yes
        stoner_fe_d = 0.88
        soc_fe_d = 0.06
        kmesh = 16
        band_path = G:0,0,0; X:0.5,0,0
    """
    path=Path(path).expanduser().resolve()
    lines=[]
    for item in path.read_text(encoding="utf-8").splitlines():
        item=item.strip()
        if item and not item.startswith(("#",";")):
            lines.append(item)
    if len(lines)<2:
        raise ValueError("Spin file requires format and geometry on the first two lines")
    fmt=lines[0].lower().split("=")[-1].strip()
    if fmt not in ("xyz","cif","ase","yaehmop"):
        raise ValueError("Expected format: xyz, cif, ase or yaehmop")
    loc=lines[1].split("=",1)[-1].strip().strip("'\"")
    source=Path(loc).expanduser()
    if not source.is_absolute():
        source=path.parent/source
    source=source.resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Geometry not found: {source}")
    opts={}
    for entry in lines[2:]:
        if "=" not in entry:
            raise ValueError("Invalid parameter (use name = value): "+entry)
        key,val=entry.split("=",1)
        key=key.strip().lower().replace("-","_")
        if key in opts:
            raise ValueError("Duplicate parameter: "+key)
        opts[key]=val.strip()
    known={"spin","soc","kmesh","band_path","divisions","cutoff",
           "images","workers","quadrature","temperature_ev","tolerance",
           "mixing","max_iter","initial_moment","output","dos_sigma",
           "dos_points","plots","electrons","charge","ase_format","parameters"}
    for k in opts:
        if k in known or k.startswith("stoner_") or k.startswith("soc_"):
            continue
        raise ValueError("Unsupported option: "+k)
    spin=_yesno(opts.get("spin","yes"))
    soc=_yesno(opts.get("soc","no"))
    image=_yesno(opts.get("plots","yes"))
    cutoff=float(opts.get("cutoff",5.6))
    if cutoff<=0:raise ValueError("cutoff must be positive")
    images=opts.get("images","auto")
    images=("auto" if images=="auto" else _numbers(images,int))
    from pyaehmop.parameter_files import load_parameter_selection
    parameter_table=load_parameter_selection(opts.get("parameters","default"),path.parent)
    if fmt=="yaehmop":
        structure=read_input(source)
        if images!="auto" and structure.vectors is not None:
            if len(images)!=len(structure.vectors):
                raise ValueError("Image dimension does not match periodicity")
            structure.bounds=images
    elif fmt=="xyz":
        structure=read_xyz(source,parameter_table=parameter_table)
    else:
        structure=read_ase_structure(source,fmt,cutoff,images,opts.get("ase_format"),parameter_table=parameter_table)
    if "parameters" in opts:
        from pyaehmop.parameter_files import apply_parameter_selection
        apply_parameter_selection(opts["parameters"],path.parent,structure,fmt)
    if "electrons" in opts and "charge" in opts:
        raise ValueError("Choose electrons OR charge, not both")
    if "electrons" in opts:
        structure.electrons=float(opts["electrons"])
    elif "charge" in opts:
        structure.electrons=sum(x.valence for x in structure.atoms)-float(opts["charge"])
    stoner={};xi={}
    for key,val in opts.items():
        if key.startswith("stoner_"):
            kind="stoner";parts=key.split("_");dict_target=stoner
        elif key.startswith("soc_") and key not in {"soc"}:
            kind="soc";parts=key.split("_");dict_target=xi
        else:
            continue
        # soc_fe_d, stoner_fe_d; default 'd' for soc_fe.
        if len(parts) not in (2,3) or (len(parts)==3 and parts[2] not in "spdf"):
            raise ValueError("Stoner/SOC parameters: use stoner_fe_d and soc_fe_d")
        species=parts[1].upper()
        shell=parts[2] if len(parts)==3 else "d"
        strength=float(val)
        if kind=="stoner" and strength<0:
            raise ValueError("Stoner constant must be nonnegative")
        dict_target.setdefault(species,{})[shell]=strength
    settings=SpinSettings(
        spin=spin,soc=soc,stoner=stoner,xi=xi,
        initial_moment=float(opts.get("initial_moment",2.0)),
        mixing=float(opts.get("mixing",0.3)),
        temperature=float(opts.get("temperature_ev",0.025)),
        tolerance=float(opts.get("tolerance",1e-5)),
        max_iter=int(opts.get("max_iter",130)),
        workers=int(opts.get("workers",1)))
    settings.validate()
    quad=_numbers(opts.get("quadrature","20 18 12"),int)
    if len(quad)!=3 or min(quad)<4:
        raise ValueError("quadrature requires three integers >=4")
    if structure.vectors is None:
        kpts=np.zeros((1,3))
        weights=np.array([1.])
        bandpath=None; labels=[]
    else:
        dim=len(structure.vectors)
        mesh=_numbers(opts.get("kmesh"," ".join(["12"]*dim)),int)
        if len(mesh)!=dim:
            raise ValueError(f"kmesh requires {dim} integers")
        points=build_kmesh(mesh,dim)
        kpts=np.array([k for k,w in points])
        weights=np.array([w for k,w in points])
        bandpath=None;labels=[]
        if "band_path" in opts:
            path_value=opts["band_path"]
            if path_value.lower()=="auto":
                from pyaehmop.structure_io import bandpath_automatic
                bandpath,labels=bandpath_automatic(
                    structure.vectors,int(opts.get("divisions",80)))
            else:
                special=_k_path(path_value)
                bandpath,labels=interpolate_path(
                    special,int(opts.get("divisions",50)))
    outdir=Path(opts.get("output","results/"+path.stem)).expanduser()
    if not outdir.is_absolute():
        outdir=path.parent/outdir
    return SpinJob(path,source,structure,settings,kpts,weights,
                   bandpath,labels,tuple(quad),cutoff,outdir.resolve(),
                   float(opts.get("dos_sigma",0.14)),
                   int(opts.get("dos_points",550)),image)
