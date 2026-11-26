"""Human-friendly input: first two lines specify format and file.

Example:
    xyz
    water.xyz
    calculations = energy, charges, dos
    charge = 0
    workers = 2
    output = results/water
"""
from dataclasses import dataclass
from pathlib import Path
import re
import unicodedata

import numpy as np

from .input import read_input
from .basis import shells_from_custom
from .structure_io import read_xyz,read_ase_structure,build_kmesh,bandpath_automatic


KINDS = {"yaehmop","xyz","cif","ase"}
TASKS = {"energy","charges","dos","bands","bonds","fmo","orbitals","matrices"}
SUPPORTED = {
    "calculations","charge","electrons","workers","quadrature","cutoff",
    "images","kmesh","k_shift","band_path","divisions",
    "dos_sigma","dos_points","output","outputs","weighted","k",
    "bonds","bond_analysis","ase_format","frame","fragments","kpoints",
    "parameters","plots"
}
ALIASES = {
    "calculate":"calculations","calculation":"calculations","tarefas":"calculations",
    "charge":"charge","electrons":"electrons","n_electrons":"electrons",
    "workers":"workers","parallel":"workers","workers":"workers",
    "cutoff":"cutoff","kmesh":"kmesh","kgrid":"kmesh",
    "band":"band_path","bandpath":"band_path",
    "band_path":"band_path","band_divisions":"divisions",
    "kpoints":"kpoints","k_points":"kpoints","fragments":"fragments",
    "smearing":"dos_sigma","dos_sigma":"dos_sigma",
    "output":"output","output_dir":"output","output_format":"outputs",
    "ase_format":"ase_format", "bond_analysis":"bond_analysis",
    "parameters":"parameters","parms":"parameters","parameters_atomicos":"parameters",
    "plots":"plots","plot":"plots","plots_png":"plots"
}
TASK_ALIASES = {
    "energy":"energy","energy":"energy","e":"energy",
    "charges":"charges","charge":"charges","charges":"charges",
    "only_charges":"charges","only_charges":"charges","mulliken":"charges",
    "dos":"dos","densidade_estados":"dos","density_of_states":"dos",
    "bands":"bands","bands":"bands","band":"bands",
    "bonds":"bonds","bond":"bonds","coop":"bonds","cohp":"bonds",
    "bonds":"bonds","fmo":"fmo","fragments":"fmo",
    "orbitals":"orbitals","levels":"orbitals","levels":"orbitals",
    "matrices":"matrices","matrices":"matrices",
    "all":"all","all":"all","all":"all",
    "basic":"basic","basic":"basic"
}


def _normalize(s):
    return "".join(c for c in unicodedata.normalize("NFD",str(s))
                   if unicodedata.category(c)!="Mn").lower().replace("-","_")


def _path(base,value):
    raw=Path(value.strip().strip(chr(34)+chr(39))).expanduser()
    return raw if raw.is_absolute() else base/raw


def _numbers(s,type_func=float):
    return tuple(type_func(x) for x in re.split(r"[\s,]+",s.strip()) if x)


def _yesno(v):
    val=_normalize(v.strip())
    if val in ("yes","s","true","yes","1","on"):return True
    if val in ("no","n","false","no","0","off"):return False
    raise ValueError("Expected yes/no: "+v)


def _k_path(value):
    result=[]
    for entry in value.split(";"):
        entry=entry.strip()
        if not entry:
            continue
        if ":" not in entry:
            raise ValueError("band_path: use 'G:0,0,0; X:0.5,0,0'")
        label,coords=entry.split(":",1)
        p=_numbers(coords)
        if len(p)!=3 or not label.strip():
            raise ValueError("Each band point needs a name and three coordinates")
        result.append((label.strip(),np.array(p,dtype=float)))
    if len(result)<2:
        raise ValueError("band_path requires at least two special points")
    return result


def _bonds(value,atoms,tasks):
    """Sintaxe: bonds = 1-2; 2-3@1,0,0; 1-1@0,1,0."""
    wanted=set(tasks)
    items=[]
    for j,entry in enumerate(value.split(";"),1):
        entry=entry.strip()
        if not entry:
            continue
        mat=re.fullmatch(r"(\d+)\s*-\s*(\d+)(?:\s*@\s*([-+\d,\s]+))?",entry)
        if not mat:
            raise ValueError("bonds: use '1-2; 2-3@1,0,0' (atoms are numbered starting at 1)")
        a,b=int(mat.group(1)),int(mat.group(2))
        if not (1<=a<=len(atoms) and 1<=b<=len(atoms)):
            raise ValueError("Atom index in bonds is outside the structure")
        translation=_numbers(mat.group(3) or "0 0 0",int)
        if len(translation)!=3:
            raise ValueError("Bond translation requires three integers")
        if a==b and translation==(0,0,0):
            raise ValueError("A bond between an atom and itself in the same cell is undefined")
        # COOP and COHP are evaluated together for each pair, although only
        # the selected results are exported to the filtered JSON.
        kind="cohp" if "cohp" in wanted and "coop" not in wanted else "coop"
        items.append((kind,j,a-1,b-1,translation))
    if not items:
        raise ValueError("bonds is empty")
    return items


@dataclass
class EasyJob:
    config_path: Path
    source_type: str
    source_path: Path
    options: dict
    tasks: frozenset[str]
    output_dir: Path
    outputs: frozenset[str]
    workers: int
    quadrature: tuple[int,int,int]
    cutoff: float
    dos_sigma: float
    dos_points: int
    structure: object
    bands_kpoints: object = None
    bands_labels: list | None = None
    notes: tuple[str,...] = ()


def is_easy_input(path):
    """Safe autodetection without disrupting traditional YAeHMOP inputs."""
    path=Path(path)
    try:
        lines=path.read_text(encoding="utf-8-sig").splitlines()
    except UnicodeDecodeError:
        return False
    first=next((x.strip().lower() for x in lines if x.strip() and not x.lstrip().startswith(("#",";"))),"")
    if first.startswith("format"):
        first=re.sub(r"^format\s*[:=]\s*","",first)
    return first in KINDS


def read_easy_job(path):
    """Read and validate the entire file; paths are relative to the configuration folder."""
    config=Path(path).expanduser().resolve()
    lines=[]
    for lineno,raw in enumerate(config.read_text(encoding="utf-8-sig").splitlines(),1):
        value=raw.strip()
        if value and not value.startswith(("#",";")):
            lines.append((lineno,value))
    if len(lines)<2:
        raise ValueError("Friendly input: first line is format; second line is file path")
    format_line=lines[0][1]
    if "=" in format_line or ":" in format_line:
        parts=re.split(r"\s*[:=]\s*",format_line,maxsplit=1)
        if _normalize(parts[0])!="format":
            raise ValueError("The first line must be yaehmop, xyz, cif or ase")
        format_line=parts[1]
    kind=_normalize(format_line)
    if kind not in KINDS:
        raise ValueError("Unknown format: "+kind)
    source_line=lines[1][1]
    if re.match(r"^(file|file)\s*[:=]",source_line,re.I):
        source_line=re.split(r"\s*[:=]\s*",source_line,maxsplit=1)[1]
    source=_path(config.parent,source_line)
    if not source.is_file():
        raise FileNotFoundError(f"Structure file not found: {source}")
    opts={}
    for lineno,raw in lines[2:]:
        if "=" not in raw and ":" not in raw:
            raise ValueError(f"Line {lineno}: expected 'parameter = value': {raw!r}")
        key,value=re.split(r"\s*[:=]\s*",raw,maxsplit=1)
        key=ALIASES.get(_normalize(key.strip()),_normalize(key.strip()))
        if key not in SUPPORTED:
            raise ValueError(f"Line {lineno}: unknown option: {key!r}. "
                             f"Options: {', '.join(sorted(SUPPORTED))}")
        if key in opts:
            raise ValueError(f"Line {lineno}: repeated option: {key}")
        opts[key]=value.strip()
    if "charge" in opts and "electrons" in opts:
        raise ValueError("Specify either charge OR electrons, not both")
    tasks_raw=opts.get("calculations","all" if kind=="yaehmop" else "basic")
    task_values=[_normalize(x.strip()) for x in tasks_raw.split(",") if x.strip()]
    if not task_values:
        raise ValueError("calculations cannot be empty")
    for val in task_values:
        if val not in TASK_ALIASES:
            raise ValueError(f"Unknown calculation: {val}. Escolha {', '.join(TASKS)}")
    selected={TASK_ALIASES[x] for x in task_values}
    if "all" in selected:
        if len(selected)!=1:
            raise ValueError("Do not mix 'all' with individual calculations")
        # Task expansion depends on whether the geometry is molecular or periodic.
    elif "basic" in selected:
        selected.discard("basic")
        selected.update(("energy","charges"))

    workers=int(opts.get("workers",1))
    if workers<1 or workers>128:
        raise ValueError("workers must be an integer between 1 and 128")
    quadrature=_numbers(opts.get("quadrature","24 20 12"),int)
    if len(quadrature)!=3 or any(x<4 for x in quadrature):
        raise ValueError("quadrature: supply three integers >=4")
    cutoff=float(opts.get("cutoff",11.0))
    if cutoff<=0:
        raise ValueError("cutoff must be positive")
    sigma=float(opts.get("dos_sigma",0.15))
    points=int(opts.get("dos_points",600))
    if sigma<=0 or not 10<=points<=100000:
        raise ValueError("dos_sigma must be >0 and dos_points must be between 10 and 100000")
    images_value=opts.get("images","auto")
    images="auto" if _normalize(images_value)=="auto" else _numbers(images_value,int)
    from .parameter_files import load_parameter_selection
    parameter_table=load_parameter_selection(opts.get("parameters","default"),config.parent)
    if kind=="yaehmop":
        structure=read_input(source)
        if images!="auto" and structure.vectors is not None:
            if len(images)!=len(structure.vectors):
                raise ValueError("images dimensionality is incompatible with the lattice")
            structure.bounds=images
    elif kind=="xyz":
        structure=read_xyz(source,parameter_table=parameter_table)
        if images!="auto":
            raise ValueError("images is not applicable to XYZ molecules")
    else:
        structure=read_ase_structure(
            source,kind,cutoff,images,ase_format=opts.get("ase_format"),
            index=int(opts.get("frame",0)),parameter_table=parameter_table)
    return _finish_job(config,kind,source,opts,selected,workers,
                       quadrature,cutoff,sigma,points,structure)



def _special_kpoints(value,dim):
    points=[]
    for entry in value.split(";"):
        entry=entry.strip()
        if not entry:
            continue
        coordinates=_numbers(entry)
        if len(coordinates) not in (3,4):
            raise ValueError("kpoints: use '0 0 0; 0.5 0 0.5' ou 'kx ky kz weight'")
        k=np.array(coordinates[:3],dtype=float)
        if dim<3 and np.any(k[dim:]):
            raise ValueError("kpoints: nonperiodic direction must have zero coordinate")
        weight=coordinates[3] if len(coordinates)==4 else 1.0
        if weight<0:
            raise ValueError("kpoint weights cannot be negative")
        points.append((k,float(weight)))
    if not points or sum(w for _,w in points)<=0:
        raise ValueError("kpoints is empty or total weight is zero")
    return points


def _parse_fragments(value):
    from .fragments import parse_integer_ranges
    result=[]
    for entry in value.split(";"):
        entry=entry.strip()
        if not entry:
            continue
        if ":" not in entry:
            raise ValueError("fragments: use '1:4; 2-5:32' (atoms:electrons)")
        atoms,ne=entry.rsplit(":",1)
        result.append((parse_integer_ranges(atoms),float(ne)))
    if not result:
        raise ValueError("fragments is empty")
    return result


def _finish_job(config,kind,source,opts,selected,workers,quad,cutoff,sigma,
                points,structure):
    notes=[]
    if "parameters" in opts:
        elements=_parameter_overrides(opts["parameters"],config.parent,
                                      structure,kind)
        if elements:
            notes.append("Custom EHT parameters: "+", ".join(elements))
    periodic=structure.vectors is not None
    tasks=set(selected)
    if "all" in tasks:
        tasks=set(TASKS)
        if not periodic:
            tasks.discard("bands")
        if not opts.get("bonds") and not opts.get("bond_analysis") and not structure.bond_requests:
            tasks.discard("bonds")
        if not opts.get("fragments") and not structure.fmo_fragments:
            tasks.discard("fmo")
        if periodic:
            tasks.discard("orbitals")  # Bloch orbitals are represented through the bands
            tasks.discard("matrices")  # Per-k matrices are not exported by the core engine
    if "bands" in tasks and not periodic:
        raise ValueError("Electronic bands require a periodic structure (CIF/ASE/YAeHMOP)")
    if periodic and "matrices" in tasks:
        raise NotImplementedError("Per-k matrices are not exported yet; use calculations=bands")
    if periodic and "orbitals" in tasks:
        raise NotImplementedError("Per-k Bloch orbitals are not exported yet; use calculations=bands")
    if "fmo" in tasks and periodic:
        raise NotImplementedError("Periodic FMO requires FCO, which is not available")
    if "charge" in opts:
        structure.electrons=sum(x.valence for x in structure.atoms)-float(opts["charge"])
    if "electrons" in opts:
        structure.electrons=float(opts["electrons"])
    if not 0 <= structure.electrons:
        raise ValueError("Negative electron count after specifying charge")
    if "weighted" in opts:
        structure.weighted=_yesno(opts["weighted"])
    if "k" in opts:
        structure.k_constant=float(opts["k"])

    if opts.get("fragments"):
        structure.fmo_fragments=_parse_fragments(opts["fragments"])
        atomids=[i for frag,el in structure.fmo_fragments for i in frag]
        if sorted(atomids)!=list(range(1,len(structure.atoms)+1)):
            raise ValueError("fragments: each atom must appear exactly once")
    if "fmo" in tasks and not structure.fmo_fragments:
        raise ValueError("FMO requested without fragments; define 'fragments = 1:...; 2-3:...'")
    if "fmo" not in tasks:
        structure.fmo_fragments=[]

    if "bonds" in opts or "bond_analysis" in opts:
        structure.bond_requests=_bonds(
            opts.get("bonds",opts.get("bond_analysis","")),structure.atoms,tasks)
    if "bonds" in tasks and not structure.bond_requests:
        raise ValueError("Bond analysis requested without pairs. Define 'bonds = 1-2; 1-2@1,0,0'")
    if "bonds" not in tasks:
        structure.bond_requests=[]

    overrides=None
    labels=None
    if periodic:
        dim=len(structure.vectors)
        if "kpoints" in opts and "kmesh" in opts:
            raise ValueError("Choose kpoints OR kmesh (not both)")
        average_needed=bool(tasks & {"energy","charges","dos","bonds"})
        if average_needed:
            if "kpoints" in opts:
                structure.kpoints=_special_kpoints(opts["kpoints"],dim)
            elif "kmesh" in opts:
                structure.kpoints=build_kmesh(_numbers(opts["kmesh"],int),dim,
                    _numbers(opts["k_shift"]) if "k_shift" in opts else None)
            elif not structure.kpoints:
                mesh=(2,)*dim
                structure.kpoints=build_kmesh(mesh,dim)
                notes.append("Default k mesh 2x2x2 (or 2 per periodic direction): "
                             "preliminary result; check convergence with kmesh")
        else:
            structure.kpoints=[]
        if "bands" in tasks:
            if "band_path" in opts:
                spec=opts["band_path"].strip()
                if _normalize(spec)=="auto":
                    overrides,labels=bandpath_automatic(
                        structure.vectors,int(opts.get("divisions",80)))
                    notes.append("Band path suggested automatically by ASE")
                    structure.band_path=[]
                else:
                    structure.band_path=_k_path(spec)
                    structure.band_divisions=int(opts.get("divisions",20))
            elif not structure.band_path:
                overrides,labels=bandpath_automatic(structure.vectors,80)
                notes.append("Band path suggested automatically by ASE")
            elif "divisions" in opts:
                structure.band_divisions=int(opts["divisions"])
        else:
            structure.band_path=[]
    elif any(x in opts for x in ("kmesh","kpoints","k_shift",
                                 "band_path","divisions","images")):
        raise ValueError("kmesh, kpoints, band_path e images requerem cristal periodic")

    outputs_raw=opts.get("outputs","json,csv")
    outputs=set(_normalize(x.strip()) for x in outputs_raw.split(",") if x.strip())
    if "plots" in opts and _yesno(opts["plots"]):
        outputs.add("png")
    if not outputs or not outputs <= {"json","csv","png"}:
        raise ValueError("outputs: choose json, csv, png, or a combination")
    outdir=_path(config.parent,opts.get("output",f"results/{config.stem}")).resolve()
    if outdir.suffix.lower()==".json":
        raise ValueError("output must be a DIRECTORY; JSON will be named summary.json")
    return EasyJob(config,kind,source,opts,frozenset(tasks),outdir,
                   frozenset(outputs),workers,tuple(quad),cutoff,sigma,
                   points,structure,overrides,labels,tuple(notes))


def calculate_easy_job(job):
    """Executa a tarefa selecionada e devolve result filtrado + metadados."""
    from .engine import calculate
    from .spectra import gaussian_dos
    s=job.structure
    selected=set(job.tasks)
    report=calculate(
        s,quadrature=job.quadrature,workers=job.workers,cutoff=job.cutoff,
        kpoints=job.bands_kpoints,bands="bands" in selected,
        dos_sigma=job.dos_sigma,dos_points=job.dos_points,
        compute_dos="dos" in selected)
    if job.bands_labels is not None:
        report["special_kpoints"]=job.bands_labels
    data={
        "general_input":{
            "configuration_file":str(job.config_path),
            "source_format":job.source_type,
            "structure_file":str(job.source_path),
            "requested_calculations":sorted(selected),
            "workers":job.workers,
            "cutoff_angstrom":job.cutoff,
            "quadrature":job.quadrature,
            "integration_kpoint_count":len(s.kpoints),
            "image_bounds":s.bounds},
        "type":report["type"],"title":report["title"],
        "atoms":report["atoms"],"orbitals":report["orbitals"],
        "electrons":report["electrons"],
        "species":[atom.symbol for atom in s.atoms],
        "warnings":list(report.get("warnings",[]))+list(job.notes)
    }
    if "energy" in selected and s.vectors is None:
        for key in ("electronic_energy_eV","homo_eV","lumo_eV","gap_homo_lumo_eV"):
            if key in report:
                data[key]=report[key]
    if "charges" in selected and s.vectors is None:
        for key in ("mulliken_charges","orbital_populations"):
            if key in report:
                data[key]=report[key]
    if "orbitals" in selected and s.vectors is None:
        for key in ("orbital_energies_eV","occupations","orbital_names"):
            if key in report:
                data[key]=report[key]
    if "matrices" in selected and s.vectors is None:
        for key in ("overlap_matrix","hamiltonian_matrix","coefficients"):
            data[key]=report[key]
    if "dos" in selected and s.vectors is None:
        en=np.asarray(report["orbital_energies_eV"]).reshape(1,-1)
        grid,dos=gaussian_dos(en,np.ones(1),points=job.dos_points,sigma=job.dos_sigma)
        data["dos"]={"energy_eV":grid,"states_per_eV":dos,
                     "sigma_eV":job.dos_sigma}
    if "bands" in selected and s.vectors is not None:
        for key in ("kpoints","band_energies_eV","special_kpoints",
                    "band_path_gap_eV","minimum_direct_gap_eV",
                    "vbm_eV","cbm_eV"):
            if key in report:
                data[key]=report[key]
    if s.vectors is not None:
        props=report.get("sampled_properties")
        if props is not None:
            part={}
            if "energy" in selected:
                part.update({key:props[key] for key in
                    ("fermi_eV","average_electronic_energy_eV")})
            if "charges" in selected:
                part.update({key:props[key] for key in
                    ("mulliken_charges","orbital_populations","charge_sum")})
            if "dos" in selected:
                part["dos"]=props["dos"]
            if part:
                part["total_kpoints"]=props["total_kpoints"]
                data["sampled_properties"]=part
    if "bonds" in selected:
        for key in ("bond_analyses","bond_groups"):
            if key in report:
                data[key]=report[key]
    if "fmo" in selected and "fmo" in report:
        data["fmo"]=report["fmo"]
    return data



def _parameter_overrides(value,config_dir,structure,source_type):
    """Apply built-in defaults, classic YAeHMOP or eht_parms.dat records."""
    from .parameter_files import apply_parameter_selection
    return apply_parameter_selection(value,config_dir,structure,source_type)
