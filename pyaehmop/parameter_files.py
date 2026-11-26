"""EHT parameter overrides from classic YAeHMOP and eht_parms.dat files.

This module only changes data loading; the Hamiltonian, overlap,
occupation, band, Stoner and SOC algorithms remain untouched.
"""
from pathlib import Path
import math
import re

from .basis import Shell, shells_from_custom


def _number(text):
    value=float(text.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError(f"Non-finite EHT parameter: {text}")
    return value


def _table_record(parts, line_number):
    """YAeHMOP eht_parms.dat record:
    label atomic_number valence n_zeta n shell energy zeta1 zeta2 c1 c2.
    """
    if len(parts)!=11 or parts[5].lower() not in {"s","p","d","f"}:
        raise ValueError(f"Line {line_number}: invalid eht_parms.dat record")
    symbol=parts[0].upper()
    atomic_number=int(parts[1])
    valence=_number(parts[2])
    n_zeta=int(parts[3])
    n=int(parts[4])
    l="spdf".index(parts[5].lower())
    energy=_number(parts[6])
    zeta1,zeta2=_number(parts[7]),_number(parts[8])
    c1,c2=_number(parts[9]),_number(parts[10])
    if not 1<=atomic_number<=118 or not 0<=valence<=atomic_number:
        raise ValueError(f"Line {line_number}: invalid atomic number/valence")
    if n_zeta not in (1,2) or n<1 or n<=l:
        raise ValueError(f"Line {line_number}: invalid principal quantum number or zeta count")
    # YAeHMOP tables sometimes include disabled virtual shells (zeta1=0).
    if zeta1<=0 or c1==0:
        return symbol,atomic_number,valence,None
    if n_zeta==2:
        if zeta2<=0 or c2==0:
            raise ValueError(f"Line {line_number}: invalid double-zeta parameters")
        shell=Shell(n,l,energy,(zeta1,zeta2),(c1,c2))
    else:
        shell=Shell(n,l,energy,(zeta1,),(1.,))
    return symbol,atomic_number,valence,shell


def parse_parameter_file(path):
    """Read classic PARAMETERS lines or native YAeHMOP eht_parms.dat rows.

    Classic: C 6 4 2 1.71 -21.4 2 1.625 -11.4
    Table:   C 6 4 1 2 s -21.4 1.625 0.0 1.0 0.0
    Comments starting with #, ; or !, and optional PARAMETERS/END
    section markers, are supported. Only species actually present in
    the system are applied. Duplicate shell rows follow YAeHMOP's
    last-occurrence-wins rule.
    """
    path=Path(path)
    records={}
    formats={}
    for line_number,line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(),1):
        clean=re.split(r"[;#!]",line,1)[0].strip()
        if not clean or clean.lower() in ("parameters","end","end parameters","endparameters"):
            continue
        parts=clean.split()
        if not re.fullmatch(r"[A-Za-z]{1,3}",parts[0]):
            raise ValueError(f"{path.name}:{line_number}: expected an element symbol")
        # Original tables contain a historical NO (nobelium) row
        # with its valence column missing; like the bundled reader,
        # keep the valence from earlier rows of the same element.
        if len(parts)==10 and parts[4].lower() in ("s","p","d","f"):
            earlier=records.get(parts[0].upper())
            if earlier is None:
                raise ValueError(f"{path.name}:{line_number}: missing valence column")
            parts.insert(2,str(earlier["valence"]))
        table=len(parts)==11 and parts[5].lower() in ("s","p","d","f")
        try:
            if table:
                symbol,z,valence,shell=_table_record(parts,line_number)
                new_format="eht_parms.dat"
            else:
                symbol,entry=shells_from_custom(clean)
                z,valence,shell=entry["z"],entry["valence"],None
                if not entry["shells"]:
                    raise ValueError("no active orbital shells")
                new_format="classic"
            if symbol in formats and formats[symbol]!=new_format:
                raise ValueError("mixed table and classic records for one element")
            formats[symbol]=new_format
            if symbol not in records:
                records[symbol]={"z":z,"valence":valence,"shells":[]}
            target=records[symbol]
            if target["z"]!=z:
                raise ValueError("inconsistent atomic number across shell lines")
            # Native YAeHMOP tables contain historical rows with differing
            # Nvalen for one element; match the bundled reader by retaining
            # the valence declared on its first encountered shell row.
            if new_format=="classic":
                target["shells"]=entry["shells"]
            elif shell is not None:
                target["shells"]=[s for s in target["shells"] if s.l!=shell.l]
                target["shells"].append(shell)
        except (ValueError,NotImplementedError,IndexError) as exc:
            raise ValueError(f"{path.name}:{line_number}: {exc}") from exc
    if not records or not any(v["shells"] for v in records.values()):
        raise ValueError(f"No active EHT shell parameters in {path}")
    return records


def load_parameter_selection(value,config_dir):
    """Return selected species table without changing geometry or orbitals."""
    value=str(value).strip().strip("\"'")
    if not value or value.lower()=="default":
        return {}
    if re.match(r"^[A-Za-z]{1,2}\s+\d+\s+",value):
        symbol,entry=shells_from_custom(value)
        return {symbol:entry}
    path=Path(value).expanduser()
    if not path.is_absolute():
        path=Path(config_dir)/path
    if not path.is_file():
        raise FileNotFoundError(f"EHT parameters file not found: {path}")
    return parse_parameter_file(path)


def apply_parameter_selection(value,config_dir,structure,source_type):
    """Override shell data, preserving the existing scientific algorithms."""
    overrides=load_parameter_selection(value,config_dir)
    if not overrides:
        return []
    species={atom.symbol.upper() for atom in structure.atoms}
    matching=species.intersection(overrides)
    if not matching:
        raise ValueError("Parameter file does not contain any species in the structure: "
                         +", ".join(sorted(species)))
    for symbol in matching:
        if not overrides[symbol]["shells"]:
            raise ValueError(f"Element {symbol} has no active orbital shells")
    original_valence=sum(atom.valence for atom in structure.atoms)
    for atom in structure.atoms:
        symbol=atom.symbol.upper()
        if symbol in overrides:
            entry=overrides[symbol]
            atom.valence=entry["valence"]
            atom.shells=tuple(entry["shells"])
    if source_type!="yaehmop":
        structure.electrons+=sum(atom.valence for atom in structure.atoms)-original_valence
    return sorted(matching)
