"""Two-center STO integrals via prolate spheroidal Gauss-Laguerre/Legendre quadrature.

Independent mathematical implementation; does not use or call C code.
Adjust --quadrature to check convergence; this is a numerical integral, not an exact closed form.
"""
from math import factorial
from functools import lru_cache
import numpy as np
from scipy.special import roots_laguerre
from .basis import BOHR_PER_ANGSTROM, angular_functions, radial_norm, ORDER

@lru_cache(maxsize=20)
def _nodes(nu, nv, nphi):
    u,w_u=roots_laguerre(nu)
    v,w_v=np.polynomial.legendre.leggauss(nv)
    phi=(np.arange(nphi)+0.5)*2*np.pi/nphi
    return u,w_u,v,w_v,phi


def _frame(direction):
    axis=np.asarray(direction,dtype=float)
    axis/=np.linalg.norm(axis)
    guess=np.array([0.,0.,1.]) if abs(axis[2])<0.8 else np.array([0.,1.,0.])
    ex=np.cross(guess,axis);ex/=np.linalg.norm(ex)
    ey=np.cross(axis,ex)
    return ex,ey,axis


def overlap_block(orbitals_a, orbitals_b, distance_angstrom,
                  quadrature=(24,20,12)):
    """S_ab=<chi_a(0)|chi_b(R)> for every pair of functions on atoms A and B.

    Angularly normalized orbitals with exponents in bohr^-1.
    For R=0, enforce exact angular orthogonality at coincident centers.
    """
    n1,n2=len(orbitals_a),len(orbitals_b)
    if not n1 or not n2:
        return np.zeros((n1,n2))
    if min(quadrature)<4:
        raise ValueError("Each quadrature dimension must contain at least 4 nodes")
    d=np.asarray(distance_angstrom,dtype=float)
    R=np.linalg.norm(d)*BOHR_PER_ANGSTROM
    if R<1e-9:
        # At the same center, normalized Y_lm functions are orthonormal.
        result=np.zeros((n1,n2))
        for i,a in enumerate(orbitals_a):
            for j,b in enumerate(orbitals_b):
                if a.label==b.label:
                    total=0.
                    for z1,c1 in zip(a.zetas,a.coeffs):
                        for z2,c2 in zip(b.zetas,b.coeffs):
                            total+=(c1*c2*radial_norm(a.n,z1)*radial_norm(b.n,z2)
                                    *factorial(a.n+b.n) /
                                    (z1+z2)**(a.n+b.n+1))
                    result[i,j]=total
        return result
    nu,nv,nphi=quadrature
    x,wu,v,wv,phi=_nodes(nu,nv,nphi)
    min_z1=min(z for o in orbitals_a for z in o.zetas)
    min_z2=min(z for o in orbitals_b for z in o.zetas)
    scale=0.5*R*(min_z1+min_z2)
    if scale<=0:
        raise ValueError("Expoentes STO devem ser positivos")
    u=(1+x/scale)[:,None,None]
    vv=v[None,:,None]
    ph=phi[None,None,:]
    hh=R/2
    rho=hh*np.sqrt(np.maximum(0,(u*u-1)*(1-vv*vv)))
    zz=hh*u*vv
    e1,e2,e3=_frame(d)
    cos,sin=np.cos(ph),np.sin(ph)
    transverse=rho[...,None]*(cos[...,None]*e1+sin[...,None]*e2)
    ra_vec=transverse+(zz+hh)[...,None]*e3
    rb_vec=transverse+(zz-hh)[...,None]*e3
    r1=(hh*(u+vv))+0*ph
    r2=(hh*(u-vv))+0*ph
    # u integration: exp(-x) du -> w_u dx / scale;
    # each STO is evaluated with exp(+x/2), cancelling Gauss-Laguerre weight.
    w=wu[:,None,None]*wv[None,:,None]*np.ones_like(ph)*(2*np.pi/nphi)
    weight=(w*(hh**3)*(u*u-vv*vv)/scale).ravel()
    angleA={}
    angleB={}
    for l in set(o.l for o in orbitals_a):
        angleA[l]=angular_functions(ra_vec,l)
    for l in set(o.l for o in orbitals_b):
        angleB[l]=angular_functions(rb_vec,l)

    exp_part=0.5*x[:,None,None]
    def values(orbs, distances, angles):
        value=np.empty((weight.size,len(orbs)))
        for i,orb in enumerate(orbs):
            radial=np.zeros_like(distances)
            for z,c in zip(orb.zetas,orb.coeffs):
                radial += c*radial_norm(orb.n,z)*distances**(orb.n-1) * np.exp(-z*distances+exp_part)
            angle=angles[orb.l][ORDER["spdf"[orb.l]].index(orb.label)]
            value[:,i]=(radial*angle).ravel()
        return value
    va=values(orbitals_a,r1,angleA)
    vb=values(orbitals_b,r2,angleB)
    return va.T@(weight[:,None]*vb)
