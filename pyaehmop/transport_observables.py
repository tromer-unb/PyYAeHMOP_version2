"""Landauer observables for a frozen, per-spin transmission spectrum."""
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import expit
G0=7.748091729e-5  # S, including two spin channels
KB=8.617333262145e-5 # eV/K

def _grid(energy,T):
    e=np.asarray(energy,float); t=np.asarray(T,float)
    if e.ndim!=1 or t.shape!=e.shape or len(e)<2 or not np.isfinite(e).all() or not np.isfinite(t).all() or np.any(np.diff(e)<=0) or np.any(t<0):
        raise ValueError("Expected increasing finite energies and nonnegative transmission")
    return e,t

def conductance(energy,T):
    """Zero-temperature conductance versus chemical potential, in siemens."""
    e,t=_grid(energy,T)
    return G0*t

def current_voltage(energy,T,bias,temperature=300.,mu=0.):
    """Current in A; symmetric bias, frozen equilibrium T(E), spin degeneracy 2."""
    e,t=_grid(energy,T); v=np.atleast_1d(np.asarray(bias,float))
    if not np.isfinite(v).all() or not np.isfinite(temperature) or temperature<0 or not np.isfinite(mu):
        raise ValueError("Invalid temperature, bias or chemical potential")
    margin=20*KB*temperature
    if mu-max(abs(v))/2-margin<e[0] or mu+max(abs(v))/2+margin>e[-1]:
        raise ValueError("Energy grid must cover bias window and thermal tails (20 kBT)")
    out=[]
    for b in v:
        lo,hi=mu-abs(b)/2,mu+abs(b)/2
        if b==0: out.append(0.); continue
        if temperature==0:
            x=np.r_[lo,e[(e>lo)&(e<hi)],hi]
            val=trapezoid(np.interp(x,e,t),x)
        else:
            # Refine the quadrature if input energy spacing exceeds kBT.
            step=KB*temperature/20
            left=max(e[0],lo-margin); right=min(e[-1],hi+margin)
            x=np.linspace(left,right,int(np.ceil((right-left)/step))+1)
            val=trapezoid(np.interp(x,e,t)*(expit((hi-x)/(KB*temperature))-expit((lo-x)/(KB*temperature))),x)
        out.append(np.sign(b)*G0*val)
    return np.array(out)
