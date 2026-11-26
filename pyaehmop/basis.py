"""Slater orbital basis: s, p, d and f; Hückel parameters."""
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from math import factorial, sqrt, pi
import numpy as np

BOHR_PER_ANGSTROM = 1.889644746
ORDER = {"s": ("s",), "p": ("px", "py", "pz"),
         "d": ("dx2-y2", "dz2", "dxy", "dxz", "dyz"),
         "f": ("fz3", "fxz2", "fyz2", "fxyz", "fz(x2-y2)",
               "fx(x2-3y2)", "fy(3x2-y2)")}


@dataclass(frozen=True)
class Shell:
    n: int
    l: int
    energy: float
    zetas: tuple[float, ...]
    coeffs: tuple[float, ...]


@dataclass(frozen=True)
class Orbital:
    atom: int
    label: str
    n: int
    l: int
    energy: float
    zetas: tuple[float, ...]
    coeffs: tuple[float, ...]


@dataclass
class Atom:
    symbol: str
    position: np.ndarray
    valence: float
    shells: tuple[Shell, ...]


def radial_norm(n, zeta):
    """N=(2*zeta)^(n+1/2)/sqrt((2*n)!), for normalized Y_lm."""
    return (2*zeta)**(n+0.5) / sqrt(factorial(2*n))


def double_zeta_normalization(shell):
    """Normalize the combination of Slater functions with the same n,l at one center."""
    c = np.array(shell.coeffs)
    z = np.array(shell.zetas)
    ss = (2*np.sqrt(z[:,None]*z[None,:]) /
          (z[:,None]+z[None,:]))**(2*shell.n+1)
    return tuple((c/np.sqrt(c@ss@c)).tolist())


@lru_cache(maxsize=1)
def read_parameter_table():
    """Includes a locally redistributed eht_parms.dat file; never loads YAeHMOP C."""
    result = {}
    for line in (Path(__file__).parent / "eht_parms.dat").read_text().splitlines():
        if not line.strip() or line.lstrip().startswith(";"):
            continue
        vals = line.split()
        if len(vals) < 11:
            continue
        symbol, z, ne, nz, n, letter, en, z1, z2, c1, c2 = vals[:11]
        if letter.lower() not in "spdf":
            continue
        if float(z1) <= 0 or float(c1) == 0:
            continue
        zetas = (float(z1),float(z2)) if int(nz)==2 and float(z2)>0 else (float(z1),)
        coeffs = (float(c1),float(c2)) if len(zetas)==2 else (1.,)
        shell = Shell(int(n),"spdf".index(letter.lower()),float(en),zetas,coeffs)
        e = result.setdefault(symbol.upper(), {"valence":float(ne),"z":int(z),"shells":[]})
        # YAeHMOP has a set of s/p/d/f shells per species. Certain historical
        # rows occur repeatedly in the table; the last occurrence takes precedence.
        previous=[ix for ix,s in enumerate(e["shells"]) if s.l==shell.l]
        if previous:
            e["shells"][previous[-1]]=shell
        else:
            e["shells"].append(shell)
    return result


def shells_from_custom(fields):
    """Classical YAeHMOP format: X Z Nval Ns zs Es Np zp Ep [Nd zd Ed]."""
    v = fields.split()
    if len(v) < 6:
        raise ValueError("Invalid PARAMETERS line: "+fields)
    symbol = v[0].upper()
    atomic_no = int(v[1]); valence = float(v[2])
    groups = []
    rest = v[3:]
    for l in range(min(len(rest)//3, 4)):
        n, z, energy = int(float(rest[3*l])),float(rest[3*l+1]),float(rest[3*l+2])
        if n and z > 0:
            groups.append(Shell(n,l,energy,(z,),(1.,)))
    if len(rest)%3:
        raise NotImplementedError("Custom extended-format parameters are not yet supported")
    return symbol, {"z":atomic_no,"valence":valence,"shells":groups}


def make_basis(atoms):
    orbitals = []
    for atom_index, atom in enumerate(atoms):
        for shell in atom.shells:
            if shell.l > 3:
                raise NotImplementedError("Only s, p, d and f orbitals are supported")
            coeffs = double_zeta_normalization(shell)
            for label in ORDER["spdf"[shell.l]]:
                orbitals.append(Orbital(atom_index,label,shell.n,shell.l,
                                        shell.energy,shell.zetas,coeffs))
    if not orbitals:
        raise ValueError("The geometry contains no orbitals with valid parameters")
    return tuple(orbitals)


def angular_functions(x, l):
    """Normalized real spherical harmonics, ordering s/px/py/pz/d."""
    r2 = np.einsum("...i,...i->...",x,x)
    r = np.sqrt(np.maximum(r2,1e-24))
    xx,yy,zz = np.moveaxis(x,-1,0)
    if l == 0:
        return [np.full_like(r,1/np.sqrt(4*pi))]
    if l == 1:
        return [np.sqrt(3/(4*pi))*component/r for component in (xx,yy,zz)]
    if l == 2:
        return [np.sqrt(15/(16*pi))*(xx**2-yy**2)/r2,
                np.sqrt(5/(16*pi))*(3*zz**2-r2)/r2,
                np.sqrt(15/(4*pi))*xx*yy/r2,
                np.sqrt(15/(4*pi))*xx*zz/r2,
                np.sqrt(15/(4*pi))*yy*zz/r2]
    if l == 3:
        r3=r2*r
        return [np.sqrt(7/(16*pi))*zz*(5*zz*zz-3*r2)/r3,
                np.sqrt(21/(32*pi))*xx*(5*zz*zz-r2)/r3,
                np.sqrt(21/(32*pi))*yy*(5*zz*zz-r2)/r3,
                np.sqrt(105/(4*pi))*xx*yy*zz/r3,
                np.sqrt(105/(16*pi))*zz*(xx*xx-yy*yy)/r3,
                np.sqrt(35/(32*pi))*xx*(xx*xx-3*yy*yy)/r3,
                np.sqrt(35/(32*pi))*yy*(3*xx*xx-yy*yy)/r3]
    raise NotImplementedError("l>3 not yet supported")
