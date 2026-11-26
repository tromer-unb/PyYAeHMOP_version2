"""Convert internal coordinates (YAeHMOP Z-matrix) to Cartesian coordinates (Å).

Convention: first atom at origin, second on +z, third in xz plane;
dihedrals use the right-hand convention from the YAeHMOP manual.
Dummy atoms (&) remain in the geometry as reference sites.
"""
import math
import numpy as np


def _unit(vec):
    v=np.asarray(vec,dtype=float)
    norm=np.linalg.norm(v)
    if norm<1e-12:
        raise ValueError("Z-matrix references are coincident or collinear")
    return v/norm


def zmat_to_cartesian(records):
    """records: lines '[index] [symbol] [ref1] [distance] [ref2] [angle] [ref3] [dihedral]'.

    References are file atom IDs (not matrix indices) and must point
    to previously defined atoms. Returns [(symbol, np.array([x,y,z]))].
    """
    result=[]
    by_id={}
    for number,raw in enumerate(records):
        v=raw.replace(","," ").split()
        if len(v)<2:
            raise ValueError("Invalid Z-matrix line: "+raw)
        atom_id=int(v[0])
        symbol=v[1].upper()
        if atom_id in by_id:
            raise ValueError(f"Duplicate Z-matrix atom ID: {atom_id}")
        def reference(pos):
            if len(v)<=pos:
                raise ValueError("Missing geometry reference in "+raw)
            idx=int(v[pos])
            if idx not in by_id:
                raise ValueError(f"Z-matrix reference is missing or points forward: {idx}")
            return by_id[idx]
        if number==0:
            if len(v)>2:
                raise ValueError("First Z-matrix atom must not specify a bond")
            xyz=np.zeros(3)
        elif number==1:
            if len(v)<4:
                raise ValueError("Second atom requires a reference and bond length")
            reference(2)
            bond=float(v[3])
            if bond<=0:
                raise ValueError("Z-matrix bond length must be positive")
            xyz=np.array([0.,0.,bond])
        elif number==2:
            if len(v)<6:
                raise ValueError("Third atom requires a reference, distance and bond angle")
            a=reference(2);b=reference(4)
            bond=float(v[3]);theta=math.radians(float(v[5]))
            if bond<=0:
                raise ValueError("Z-matrix bond length must be positive")
            e=_unit(b-a)
            # Positive x for the third atom even if ref1 is the second atom.
            xyz=a+bond*(math.cos(theta)*e+math.sin(theta)*np.array([1.,0.,0.]))
        else:
            if len(v)<8:
                raise ValueError("Subsequent atoms require three references and a dihedral")
            a=reference(2);b=reference(4);c=reference(6)
            bond=float(v[3]);theta=math.radians(float(v[5]))
            phi=math.radians(float(v[7]))
            if bond<=0:
                raise ValueError("Z-matrix bond length must be positive")
            e=_unit(b-a)
            transverse=c-a-np.dot(c-a,e)*e
            transverse=_unit(transverse)
            normal=_unit(np.cross(e,transverse))
            xyz=a+bond*(math.cos(theta)*e+math.sin(theta)*
                        (transverse*math.cos(phi)-normal*math.sin(phi)))
        by_id[atom_id]=xyz
        result.append((symbol,xyz))
    return result
