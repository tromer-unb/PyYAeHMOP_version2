"""Editable native YAeHMOP EHT table model for the desktop GUI.

The editor only creates parameter data files; it never modifies the
original bundled table, EHT formulas, SCF equations or orbital integrals.
"""
from pathlib import Path
import re
from .basis import read_parameter_table
from .parameter_files import parse_parameter_file, _table_record

COLUMNS=("symbol","atomic_number","valence","n_zeta","n","shell",
         "energy_eV","zeta1","zeta2","coeff1","coeff2")
SHELL_ORDER="spdf"


def _records(table,species=None):
    if species is not None:
        wanted={str(s).upper() for s in species}
    else:
        wanted=None
    result=[]
    for element,info in sorted(table.items()):
        if wanted is not None and element not in wanted:
            continue
        for shell in sorted(info["shells"],key=lambda s:s.l):
            nz=len(shell.zetas)
            result.append({
                "symbol":element,
                "atomic_number":str(info["z"]),
                "valence":str(int(info["valence"]) if info["valence"].is_integer()
                              else info["valence"]),
                "n_zeta":str(nz),
                "n":str(shell.n),
                "shell":SHELL_ORDER[shell.l],
                "energy_eV":str(shell.energy),
                "zeta1":str(shell.zetas[0]),
                "zeta2":str(shell.zetas[1] if nz>1 else 0.0),
                "coeff1":str(shell.coeffs[0] if nz>1 else 1.0),
                "coeff2":str(shell.coeffs[1] if nz>1 else 0.0)})
    return result


def default_rows(species=None):
    """Use only active shells in the bundled YAeHMOP table."""
    return _records(read_parameter_table(),species)


def file_rows(path,species=None):
    """Load a classic or native table, preserving its physical shell data."""
    return _records(parse_parameter_file(path),species)


def validate_row(row):
    """Validate an editor row using the same native table parser."""
    values=[str(row.get(key,"")).strip() for key in COLUMNS]
    symbol=values[0].upper()
    if not re.fullmatch(r"[A-Za-z]{1,3}",symbol):
        raise ValueError("Element label must contain one to three letters")
    values[0]=symbol
    _,_,_,shell=_table_record(values,0)
    if shell is None:
        raise ValueError("A shell with zero zeta or zero coefficient is disabled; "
                         "delete it instead of saving an inactive shell")
    return dict(zip(COLUMNS,values))


def native_table_text(rows):
    """Serialize rows as native 11-column YAeHMOP eht_parms.dat format."""
    if not rows:
        raise ValueError("Add at least one active shell")
    validated=[validate_row(row) for row in rows]
    unique=set()
    per_species={}
    for row in validated:
        key=(row["symbol"],row["shell"])
        if key in unique:
            raise ValueError(f"Duplicate orbital shell: {key}")
        unique.add(key)
        symbol=row["symbol"]
        spec=(int(row["atomic_number"]),float(row["valence"]))
        if symbol in per_species and per_species[symbol]!=spec:
            raise ValueError(f"Inconsistent atomic number/valence for {symbol}")
        per_species[symbol]=spec
    lines=[
        "; Created by PyYAeHMOP Studio atomic parameter editor.",
        "; Label AtNo Nval Nzeta Nquant Shell Hii zeta1 zeta2 coeff1 coeff2"
    ]
    for row in sorted(validated,
                      key=lambda r:(int(r["atomic_number"]),SHELL_ORDER.index(r["shell"]))):
        lines.append(" ".join(row[key] for key in COLUMNS))
    return "\n".join(lines)+"\n"


def save_native_table(rows,path):
    """Write native parameter file and validate it by reading it back."""
    target=Path(path).expanduser()
    text=native_table_text(rows)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(text,encoding="utf-8")
    parsed=parse_parameter_file(target)
    if not parsed:
        raise RuntimeError("Saved file contains no active shells")
    return target.resolve()
