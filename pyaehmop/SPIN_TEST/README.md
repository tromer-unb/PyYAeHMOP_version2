# SPIN_TEST — self-consistent Stoner spin polarization and spin-orbit coupling

**English experimental edition.** This is an isolated research prototype.
Its input, output, comments, figures, and diagnostic messages are in English.
The physical model, numerical parameters, algorithms and quadratures are
**unchanged** from the preceding SPIN_TEST implementation.

Location:

```text
pyaehmop/SPIN_TEST/
```

This module is developed separately from the conventional EHT engine,
so its experimental spin/SOC limitations remain clearly documented.

## Three mutually selectable physical modes

| Input | Mode | Meaning of up and down |
|---|---|---|
| `spin = yes`, `soc = no` | Self-consistent collinear Stoner | Exact spin-up and spin-down bands |
| `spin = no`, `soc = yes` | Nonmagnetic onsite atomic SOC | Spinor eigenstates; Kramers pairs where required |
| `spin = yes`, `soc = yes` | Self-consistent collinear-z Stoner with SOC | Spin is mixed; display projections ⟨σz⟩ |

**Important:** with SOC, an individual eigenstate generally does not
belong to a pure spin-up or spin-down channel. The colored bands represent
the expectation value `<sigma_z>`, ranging from −1 to +1.

## Reproduce the Fe-wire demonstration

The monatomic Fe wire has lattice constant
`a = 2.2013 angstroms ≈ 4.16 bohr`. The geometry is based on the
example analyzed by Autès and coauthors (2006), but the curves are
produced independently by **EHT + Stoner + SOC** in this repository.
The figure is **not** a quantitative reproduction of the paper's bands.

```bash
cd PyYAeHMOP
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -m pyaehmop.SPIN_TEST --demo
```

Run only one mode with:

```bash
python3 -m pyaehmop.SPIN_TEST \\
    pyaehmop/SPIN_TEST/examples/Fe_wire_spin.spin

python3 -m pyaehmop.SPIN_TEST \\
    pyaehmop/SPIN_TEST/examples/Fe_wire_soc.spin

python3 -m pyaehmop.SPIN_TEST \\
    pyaehmop/SPIN_TEST/examples/Fe_wire_both.spin
```

## Friendly spin configuration

Example `examples/Fe_wire_both.spin`:

```ini
ase
Fe_wire.extxyz

spin = yes
soc = yes

stoner_fe_d = 0.88
soc_fe_d = 0.06

initial_moment = 2.5
mixing = 0.3
temperature_ev = 0.025
tolerance = 0.00001
max_iter = 140

kmesh = 24
band_path = G:0,0,0; X:0.5,0,0
divisions = 75

images = 3
cutoff = 5.6
quadrature = 20 18 12
workers = 2

plots = yes
output = ../results/Fe_wire_both
```

This uses the same two-line structure format as the friendly main engine:
`yaehmop`, `xyz`, `cif` or `ase`, then the structure path.
All subsequent scientific options are in English.

Other available examples:

- `Fe_wire_spin.spin`: spin-polarized Stoner SCF without SOC.
- `Fe_wire_soc.spin`: SOC without magnetic Stoner exchange.
- `Fe_wire_both.spin`: simultaneous Stoner SCF and SOC.
- `Fe_wire_soc_exaggerated.spin`: deliberately enhanced SOC strength,
  **for visualization only**; not physical Fe calibration.
- `C_atom_both.spin`: a single carbon atom using a molecular XYZ input.

### Main parameters

| Key | Example | Meaning |
|---|---|---|
| `spin` | `yes` | Enable self-consistent collinear exchange |
| `soc` | `yes` | Include atomic `ξ L·S` |
| `stoner_fe_d` | `0.88` | Stoner I parameter for Fe d shells, eV |
| `soc_fe_d` | `0.06` | Onsite Fe d-shell SOC strength, eV |
| `initial_moment` | `2.5` | Starting moment in μB per Stoner group |
| `mixing` | `0.3` | SCF linear-mixing fraction |
| `temperature_ev` | `0.025` | Fermi–Dirac broadening kBT, eV |
| `tolerance` | `0.00001` | SCF convergence threshold |
| `max_iter` | `140` | Maximum magnetic SCF steps |
| `kmesh` | `24` | K-point grid in each periodic direction |
| `band_path` | `G:0,0,0; X:0.5,0,0` | Band structure path |
| `divisions` | `75` | Number of divisions per segment |
| `images` | `3` | Periodic image range |
| `cutoff` | `5.6` | STO overlap cutoff (angstroms) |
| `quadrature` | `20 18 12` | Numerical integration orders |
| `workers` | `2` | CPU processes |
| `plots` | `yes` | Save band and spin-DOS figures |

The physical constants in the demonstration are examples taken from
common tight-binding Fe models, not a complete combined
EHT/Stoner/SOC fitting procedure. Do not use them for quantitative
claims without a new parameterization and independent validation.

## Underlying unchanged equations

For each k point the nonmagnetic engine supplies the extended Hückel
Hamiltonian `H0(k)` and overlap `S(k)`.

For each atom/angular-momentum Stoner group:

```text
m_(A,l) = sum_mu [(P_up S)_(mu,mu) - (P_down S)_(mu,mu)]

D_mu = I_(A,l) * m_(A,l) / 2

H_up   = H0 - (S D + D S)/2
H_down = H0 + (S D + D S)/2
```

The potential is symmetrized in the nonorthogonal overlap metric.
The SCF uses one Fermi level for both spin channels, Fermi–Dirac
occupations, and linear mixing of the spin moments.

The onsite atomic spin-orbit term is:

```text
H_SOC = xi * (Lx*Sx + Ly*Sy + Lz*Sz)
      = xi/2 * [      Lz     Lx-i*Ly ]
               [   Lx+i*Ly     -Lz   ]
```

The Hermitian spinor generalized eigenproblem is

```text
H_spinor C = S_spinor C E
S_spinor = diag(S,S)
```

Spin-projected DOS uses weights
`w_up=(1+<sigma_z>)/2` and `w_down=(1-<sigma_z>)/2`.

## Produced files

Run `--demo` to generate:

- `results/Fe_wire_literature_comparison.png`: overview of spin-only,
  SOC-only, spin+SOC and energy-sorted eigenvalue differences.
- `results/Fe_wire_SOC_spin_mixing_detail.png`: zoomed avoided crossing
  and change of spin character.
- `results/Fe_wire_comparison_summary.json`: numerical mode comparison.
- `results/convergence_study.json`: k-mesh, cutoff and quadrature study.
- `results/Fe_wire_spin/`, `Fe_wire_soc/`, `Fe_wire_both/`:
  `spin_soc_result.json`, `bands_spin.csv`, `spin_dos.csv`,
  `historico_scf.csv` and `bands_spin_soc.png`.

Each JSON stores the exchange/SOC input, SCF convergence, total magnetic
moment, occupation numbers, k-point band energies, and spin projections.

## Scientific test results

The Fe-wire example produces:

| Mode | Iterations | Total z-moment, μB/Fe |
|---|---:|---:|
| Spin only | 34 | 2.000056 |
| SOC only | 1 (not a magnetic SCF) | approximately zero |
| Spin + SOC | 34 | 1.996872 |

In spin+SOC mode 22 sampled band states had
`abs(<sigma_z>) < 0.8`, demonstrating spin mixing.
The maximum energy-sorted shift along the displayed path was
approximately 48.7 meV. These changes are **not a state-tracked
adiabatic band comparison** at crossings.

The SCF and orbital operator tests validate:

- Hermitian angular momenta and `[Lx,Ly] = i Lz` for s/p/d/f;
- atomic p-shell SOC multiplets, including correct degeneracies;
- electron conservation and nonorthogonal overlap treatment;
- convergence to Fe-wire ferromagnetic and nonmagnetic solutions;
- spin-reversal symmetry without SOC;
- the zero-SOC limit;
- serial and parallel numerical agreement;
- Kramers degeneracy of this inversion-symmetric nonmagnetic wire;
- valid PNG figures and a molecular carbon test.

The original 17 SPIN_TEST tests and the 34 engine tests have been
rerun in the English copy.

## Convergence limitations

A separate study on Fe-wire spin+SOC gave:

| Calculation | z-moment μB/Fe | EHT+Stoner model energy (eV/cell) |
|---|---:|---:|
| kmesh 12, cutoff 5.6 Å | 1.996781 | −103.171294 |
| kmesh 24, cutoff 5.6 Å | 1.996872 | −103.155644 |
| kmesh 48, cutoff 5.6 Å | 1.996827 | −103.158984 |
| kmesh 24, cutoff 7.0 Å | 1.996873 | −103.155610 |

The energy varies by several meV with k-point sampling. The model is
**not sufficiently converged for quantitative magnetic anisotropy
energy predictions**, even if its moment appears numerically stable.

Crucially, the current SCF updates the Stoner exchange only along z
when SOC is enabled. It is not a fully noncollinear self-consistent
spin-density model. It omits SOC between different atomic centers,
charge self-consistency, electron correlation beyond the chosen
semiempirical exchange, and other relativistic effects.

## References

- Autès, Barreteau, Spanjaard and Desjonquères,
  *Magnetism of iron: from the bulk to the monatomic wire*,
  Journal of Physics: Condensed Matter **18**, 6785–6813 (2006).
  DOI: https://doi.org/10.1088/0953-8984/18/29/018
- Li, Smogunov, Barreteau, Ducastelle and Spanjaard,
  *Magnetocrystalline anisotropy energy of Fe(001) and Fe(110) slabs
  and nanoclusters*, Physical Review B **88**, 214413 (2013).
  DOI: https://doi.org/10.1103/PhysRevB.88.214413

The figures generated here are **new simulations**, not copies or
quantitative reconstructions of published band structures.

## Run tests

```bash
cd PyYAeHMOP
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \\
  python3 -m unittest discover -s pyaehmop/SPIN_TEST/tests -p 'test_*.py' -v
```


## Selecting external YAeHMOP EHT parameters

Set `parameters = default` or `parameters = ../../eht_parms.dat`
in any `.spin` input. The file can be a classic PARAMETERS block or
a complete native YAeHMOP EHT parameter table, including double-zeta
coefficients. This overrides the **nonmagnetic orbital EHT basis**;
independent `stoner_*` and `soc_*` parameters must still be supplied
explicitly, because they describe different physical interactions.

Example: `examples/Fe_wire_external_parameters.spin`.
A graphical input alternative is available via
`python3 -m pyaehmop.gui` from the English edition directory.
