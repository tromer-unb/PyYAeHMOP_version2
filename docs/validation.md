# Verification, regression tests and scientific limitations

[← Overview](../README.md) · [Examples](examples.md) · [Theory](theory.md)

## Verification status

The English **0.4.0** source was evaluated locally on Linux before the GitHub publication. For this source snapshot:

| Validation activity | Observed result |
|---|---:|
| Core engine, parameter files and GUI tests | 57 passed |
| Experimental spin/SOC tests | 17 passed |
| **Total automated tests** | **74 passed** |
| Numerical values matched against previous Python implementation | **10,314** |
| Maximum absolute change for those sampled values | **0** |
| Native external EHT table vs default Fe-wire spin+SOC | **0 difference** in sampled bands, spin projections, model energy and moment |
| pip wheel build and isolated installation | Successful in local environment |
| Linux desktop GUI construction and end-to-end calculations | Successful locally |

**No test runs have been initiated on GitHub Actions for this publication.** A manual workflow is provided so maintainers can trigger checks later if desired.

### Numerical regression coverage

The sampled comparisons included molecular systems (XeF₄, H–Fe and water), diamond periodic calculations and the three Fe-wire spin/SOC modes. Quantities included orbital/overlap/Hamiltonian matrices, electronic bands, Mulliken charge vectors, SCF magnetic moments and spin-projection arrays.

The equality refers only to those fixtures and numerical settings, not all theoretically possible inputs.

### What the tests actually exercise

- Integral normalization, overlap quadrature, generalized eigensolver and Hermiticity.
- Molecular occupation and Mulliken analysis.
- Brillouin-zone k-point bands and integration.
- COOP/COHP and FMO branches where implemented.
- Classic and native YAeHMOP orbital parameter file parsing.
- Double-zeta shell data preservation and custom-element files.
- XYZ/CIF/ASE friendly configuration and desktop-GUI job generation.
- Local GUI window construction and EHT/spin–SOC subprocess execution.
- Angular momentum operator algebra, selected analytic SOC multiplets and Kramers-related properties in relevant test fixtures.
- Spin reversal without SOC, nonconvergence signalling and parallel/serial consistency.

The tested conditions are defined by the source files in `pyaehmop/tests/` and `pyaehmop/SPIN_TEST/tests/`.

## Run tests locally

```bash
# From the repository root:
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

python -m unittest discover -s pyaehmop/tests -p "test_*.py" -v
python -m unittest discover -s pyaehmop/SPIN_TEST/tests -p "test_*.py" -v
```

A Tk desktop or suitable display server is required for real GUI window tests. In headless environments, GUI window tests may be skipped as intended.

## Convergence must be studied case-by-case

A local Fe-wire Stoner+SOC convergence check gave:

| K mesh | Cutoff (Å) | Moment (μB/Fe) | Model energy (eV/cell) |
|---|---:|---:|---:|
| 12 | 5.6 | 1.996781 | −103.171294 |
| 24 | 5.6 | 1.996872 | −103.155644 |
| 48 | 5.6 | 1.996827 | −103.158984 |
| 24 | 7.0 | 1.996873 | −103.155610 |

The moment is rather stable in this example, but the model energy shifts by **several meV**, indicating that this run is **not converged enough for quantitative magnetocrystalline anisotropy calculations**. SOC-driven energy differences can be much smaller than ordinary k-point errors.

Recommended practical procedure: vary k mesh, periodic image bounds/cutoff, numerical quadrature and Fermi broadening independently, then check SCF stability under changes of initial moment and mixing.

## Known model limitations

1. **No DFT:** Extended Hückel is a semiempirical one-electron method. Its electronic energy is not a complete first-principles total energy.
2. **Parameter dependence:** Onsite orbital energies, STO exponents, double-zeta coefficients, Stoner `I` and SOC `ξ` require material-specific validation.
3. **Collinear-z magnetism:** With combined spin+SOC, the exchange potential remains constrained to z; a general noncollinear spin-density SCF is not implemented.
4. **Onsite SOC:** The experimental operator omits full interatomic relativistic SOC matrix elements.
5. **K-path vs actual gap:** Finite-band-path gaps do not necessarily equal the fundamental electronic gap.
6. **Convergence:** Quadrature, neighbor cutoffs, lattice images, k-point integration and Fermi smearing can change numerical results.
7. **Finite input coverage:** The parser supports a documented subset of native YAeHMOP; not all Walsh, charge-iteration or periodic FCO options are implemented.
8. **Parallel scaling:** Multiple processes may help but are not guaranteed to be faster on tiny systems.

## Reporting a scientific issue

Use the GitHub [Issues](https://github.com/tromer-unb/PyYAeHMOP/issues) page and include a **minimal input**, reference parameter file, platform, package versions, settings, observed spectrum/charge and the expected physical or numerical property. Avoid sharing proprietary structures or unpublished sensitive data.

The goal is a transparent reproducible implementation whose approximations and caveats are documented alongside its results.
