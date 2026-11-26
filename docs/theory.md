# Scientific methodology

[← Overview](../README.md) · [Architecture](architecture.md) · [Validation](validation.md)

PyYAeHMOP implements **Extended Hückel Theory** (EHT) in Python, including a numerical Slater-type orbital basis, molecular and periodic electronic structure and selected electronic analyses. A separate research module adds Stoner spin polarization and onsite atomic spin–orbit coupling.

This document describes the **model that is actually implemented** and distinguishes demonstrated functionality from physical validation.

## 1. Atomic-orbital basis

The electronic wave functions are expanded in normalized Slater-type atomic orbitals (STOs), including real spherical harmonics for `s`, `p`, `d` and `f` shells. Each function is centered on an atom and carries the element-specific exponent(s), principal quantum number and onsite energy from the chosen EHT parameter set.

Double-zeta shells are represented as a normalized linear combination of two STO radial factors. The basis is generally **nonorthogonal**, even for atoms with the same angular momentum on different sites.

Two-center overlap elements

```math
S_{\mu\nu}=\langle\chi_\mu\mid\chi_\nu\rangle
```

are evaluated numerically in prolate spheroidal coordinates with configurable quadrature orders. Integration precision is a numerical consideration, not an exact analytic assertion.

## 2. Extended Hückel Hamiltonian

Onsite diagonal elements are the semiempirical valence orbital energies:

```math
H_{\mu\mu}=\varepsilon_\mu.
```

The off-diagonal elements are constructed from overlap `S_{μν}`, onsite energies, and a modified Wolfsberg–Helmholtz-type coupling. The implemented weighted form is:

```math
H_{\mu\nu}=\frac{\varepsilon_\mu+\varepsilon_\nu}{2}
\left[K+d^2+(1-K)d^4\right]S_{\mu\nu},\quad
d=\frac{\varepsilon_\mu-\varepsilon_\nu}{\varepsilon_\mu+\varepsilon_\nu}.
```

The default `K` is **1.75** unless a supported input explicitly changes it. Different choices of EHT parameters and `K` can materially change the spectrum.

The resulting generalized eigenproblem is

```math
HC=SC\varepsilon.
```

PyYAeHMOP uses a Hermitian generalized eigensolver from SciPy/LAPACK. It does **not** replace the physical overlap matrix by the identity.

## 3. Electron occupations, charges and energies

For molecular systems, orbital energies are sorted and occupied using the supplied electron count. The electronic EHT energy is the sum of occupied one-electron energies:

```math
E_{\mathrm{el}}=\sum_n f_n\varepsilon_n.
```

Mulliken populations are computed using the orbital coefficients, density matrix and overlap, then summed by atom. The resulting net charges depend on the basis and the partitioning convention.

This electronic energy is **not a density-functional total energy or cohesive energy**. It does not incorporate all electrostatic, ion–ion, double-counting and correlation terms needed for direct comparison with experimental formation energies.

## 4. Periodic systems and band structures

For periodic cells, real-space overlap and Hamiltonian blocks are summed with Bloch phases. For a reciprocal-coordinate k point:

```math
S(\mathbf k)=\sum_{\mathbf R}S(\mathbf R)
  e^{2\pi i\mathbf k\cdot\mathbf R}.
```

The same construction applies to `H(k)`, followed by the generalized eigensolver at each k point.

Two different k-point selections serve different purposes:

- **Integration grid** (`kmesh` or `kpoints`): occupations, Fermi level and Brillouin-zone properties.
- **Band path** (`band_path`, `divisions`): graphical band dispersion along selected reciprocal-space directions.

A band gap extracted from a path is a *path-sampled gap*, which need not equal the fundamental gap over the whole Brillouin zone. Convergence depends on the k mesh, lattice images, overlap cutoff and quadrature.

## 5. Density of states

For a finite set of band eigenvalues, the DOS is approximated with Gaussian broadening:

```math
D(E)=\sum_{n,\mathbf k}w_{\mathbf k}\,
      g_\sigma(E-\varepsilon_{n\mathbf k}).
```

The width is specified by `dos_sigma` in eV, and `dos_points` sets energy-grid resolution. The DOS is broadened numerical data, not an analytic continuum limit. Choosing `sigma` and converging the k grid are essential for interpretation.

## 6. Bonding and fragment analysis

**COOP** uses overlap-weighted orbital contributions; **COHP** uses Hamiltonian-weighted contributions. Both require a specified pair of atomic sites and, for bonds across cells, a lattice translation. These analyses depend on the conventions of the implemented basis and energy sign.

**Fragment molecular orbital (FMO)** analysis diagonalizes molecular fragments and transforms orbital coefficients to a fragment basis. This implementation is molecular; periodic FCO is **not** claimed to be supported.

## 7. Spin polarization and spin–orbit coupling

The experimental research module at `pyaehmop/SPIN_TEST/` introduces a site/shell-resolved Stoner exchange potential and a spinor basis.

### Collinear self-consistent Stoner exchange

For a group of atomic orbitals sharing an atom and angular momentum, the local moment is iterated as

```math
m_{A,l}=n_{\uparrow,A,l}-n_{\downarrow,A,l}.
```

A site/shell diagonal Stoner matrix has elements `D_μ = I_(A,l) m_(A,l) / 2`. To preserve Hermiticity in the nonorthogonal AO metric, the implementation uses a symmetrized potential:

```math
H_\uparrow=H_0-\tfrac12(SD+DS),\qquad
H_\downarrow=H_0+\tfrac12(SD+DS).
```

The Fermi level is common to both spin channels, occupations use Fermi–Dirac smearing and moments are iterated with linear mixing. This is **not DFT exchange–correlation self-consistency**.

### Onsite atomic spin–orbit coupling

A local `L·S` interaction is added in a real-orbital, two-component spinor basis:

```math
H_{SOC}=\xi\,\mathbf L\cdot\mathbf S
       =\frac{\xi}{2}\begin{pmatrix}
       L_z&L_x-iL_y\\
       L_x+iL_y&-L_z
       \end{pmatrix}.
```

This is an **onsite** approximation. It does not include fully relativistic interatomic SOC integrals.

### Combined spin+SOC

With `spin=yes` and `soc=yes`, the spinor Hamiltonian includes both Stoner exchange and SOC. The Stoner iteration updates only its **z component**; a fully noncollinear vector-magnetization SCF is not implemented.

With SOC, spin is generally **not** a conserved label: a band cannot necessarily be called pure ↑ or ↓. The spin-resolved visualization uses the eigenstate expectation `⟨σ_z⟩`, while DOS uses projections `(1±⟨σ_z⟩)/2`.

The Fe-wire demonstration is geometrically motivated by Autès *et al.* (2006), but its semiempirical parameters were **not fitted to reproduce that publication quantitatively**.

## 8. Parallel computation

The independent tasks are distributed to Python worker processes using `concurrent.futures.ProcessPoolExecutor`. Parallelizable tasks include selected interatomic overlap blocks and independent k-point eigenproblems. This can improve throughput on multicore CPUs, but overhead and memory limit speedup, especially for small molecules. See [Architecture](architecture.md).

## 9. Supported vs unsupported physics

| Implemented | Not claimed |
|---|---|
| Independent EHT Hamiltonian and nonorthogonal generalized eigensolver | First-principles Kohn–Sham DFT |
| Periodic Bloch band paths and k-grid integration | Automatic self-consistent charge redistribution |
| Semiempirical Mulliken, DOS, COOP/COHP | Experimental accuracy without fitting/convergence |
| Collinear Stoner exchange SCF | General noncollinear spin-density SCF |
| Onsite `L·S` in spinor basis | Fully relativistic interatomic SOC |
| FMO for molecules | Every original YAeHMOP feature or periodic FCO |

## References

1. R. Hoffmann, “An Extended Hückel Theory. I. Hydrocarbons,” *Journal of Chemical Physics* **39**, 1397 (1963), [DOI: 10.1063/1.1734456](https://doi.org/10.1063/1.1734456).
2. Original [YAeHMOP source code and scientific parameter file](https://github.com/greglandrum/yaehmop).
3. G. Autès *et al.*, “Magnetism of iron: from the bulk to the monatomic wire,” *Journal of Physics: Condensed Matter* **18**, 6785 (2006), [DOI: 10.1088/0953-8984/18/29/018](https://doi.org/10.1088/0953-8984/18/29/018).

The formulas and their limitations must be interpreted together with the numerical implementation and chosen EHT parameters.
