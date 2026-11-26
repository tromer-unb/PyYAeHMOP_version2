# PyYAeHMOP — English edition 0.4.0

**Scientific methodology unchanged.** This directory contains the English
source package for PyYAeHMOP. The formulas, Slater basis functions,
overlap quadratures, EHT parameter table, generalized eigensolver, numerical
cutoffs, integration grids, and parallel algorithms have not been changed.
The associated spin-polarization and spin-orbit pilot is in `SPIN_TEST/`.

The source package lives in `pyaehmop/` within this repository.
See the root `README.md` for installation and licensing details.

This is an independent **extended Hückel theory (EHT)** implementation.
It does **not** call the YAeHMOP C executable or invoke its libraries.
It uses Python, NumPy and SciPy (whose BLAS/LAPACK backends can use native code).
It is not density-functional theory (DFT).

## Run a calculation

From the English edition's parent directory:

```bash
cd PyYAeHMOP
python3 -m pip install -r pyaehmop/requirements.txt

# Original YAeHMOP-format input
OPENBLAS_NUM_THREADS=1 python3 -m pyaehmop pyaehmop/inputs/XeF4

# Friendly English input: XYZ, selected calculations, CSV/JSON/PNG
OPENBLAS_NUM_THREADS=1 python3 -m pyaehmop \
    pyaehmop/examples_friendly/water_full.job

# Periodic calculation from CIF
OPENBLAS_NUM_THREADS=1 python3 -m pyaehmop \
    pyaehmop/examples_friendly/diamond_cif.job

# Spin and SOC demonstration in the same English edition
OPENBLAS_NUM_THREADS=1 python3 -m pyaehmop.SPIN_TEST --demo
```

The friendly-input guide is **[EASY_INPUT.md](EASY_INPUT.md)**. The
spin/SOC guide is **[SPIN_TEST/README.md](SPIN_TEST/README.md)**.
Input filenames inside job files are relative to the job's directory.

Use `python3 -m pyaehmop --help` and
`python3 -m pyaehmop.SPIN_TEST --help` for command-line options.

## Scientific model

The Slater-type orbitals are normalized and include s, p, d, and f
functions, including double-zeta shells where the source parameters
supply them. Two-center overlaps are integrated numerically using
prolate spheroidal coordinates and quadrature.

For the diagonal elements of the EHT Hamiltonian,

\[
H_{ii} = \varepsilon_i .
\]

For different orbitals, the weighted EHT expression is

\[
H_{ij}=\frac{\varepsilon_i+\varepsilon_j}{2}
\left[K+d^2+(1-K)d^4\right]S_{ij},
\qquad d=\frac{\varepsilon_i-\varepsilon_j}
{\varepsilon_i+\varepsilon_j},
\]

with default `K=1.75`. The generalized eigenproblem is
`H C = S C E`, solved by `scipy.linalg.eigh(H,S)`.
The overlap matrix **is not replaced by an identity matrix**.

For periodic solids, real-space overlap blocks are Fourier summed:

\[
S(\mathbf k)=\sum_{\mathbf R}S(\mathbf R)
e^{2\pi i\mathbf k\cdot\mathbf R}.
\]

Both real-space overlap integrals and independent k-point eigenproblems
can be distributed through `concurrent.futures.ProcessPoolExecutor`.

### Features

- Cartesian molecular, Z-matrix, crystal, CIF, and ASE structure input;
- EHT electronic levels, HOMO/LUMO, band paths and density of states;
- Mulliken atomic charges;
- COOP and COHP for molecular and periodic bonds, including neighboring cells;
- Fragment molecular orbital (FMO) analysis of molecules;
- Custom atomic EHT parameters;
- Parallel CPU processing and JSON, CSV and optional PNG output;
- Experimental self-consistent Stoner spin polarization and on-site
  spin-orbit coupling through the **separate SPIN_TEST package**.

### English field names

For consistency, the translated edition uses English keys throughout
the friendly input, API, JSON and CSV outputs. Examples:

| Earlier Portuguese key | Current English key |
|---|---|
| `calculos` | `calculations` |
| `cargas` | `charges` |
| `energia_eletronica_eV` | `electronic_energy_eV` |
| `energias_bandas_eV` | `band_energies_eV` |
| `cargas_mulliken` | `mulliken_charges` |
| `amostragem_propriedades` | `sampled_properties` |
| `processos` | `workers` |
| `malha_k` | `kmesh` |
| `caminho_bandas` | `band_path` |
| `saida` | `output` |

Native YAeHMOP input syntax is intentionally retained as defined by
that external file format, so existing YAeHMOP files continue to work.

## Test and numeric verification

```bash
cd PyYAeHMOP
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

python3 -m unittest discover -s pyaehmop/tests -p 'test_*.py' -v
python3 -m unittest discover -s pyaehmop/SPIN_TEST/tests -p 'test_*.py' -v
python3 compare_numerics.py
```

The numeric regression probe independently computed XeF4, FeH,
diamond, water through the friendly XYZ interface, and all three Fe-wire
spin/SOC cases with both editions. **10,314 scalar results and matrix
elements matched exactly**, with a maximum absolute difference of zero.
This validates that translating the interface did not alter the sampled
numerical results; it is not a guarantee of predictive accuracy for
materials outside the validated examples.

## Important scientific limitations

- EHT electronic energies are not total DFT or experimental energies.
- Periodic energies, DOS, and moments require k-mesh convergence.
- A band gap along a selected path is not necessarily the fundamental gap.
- Semiempirical EHT, Stoner and SOC parameters require validation per material.
- Not all original YAeHMOP input features are supported: Walsh diagrams,
  charge iteration, electrostatics and periodic FCO remain unavailable.
- SOC in SPIN_TEST is an onsite atomic L·S model with collinear-z Stoner SCF,
  not fully noncollinear spin-density self-consistency.

## License and provenance

The original YAeHMOP project and scientific parameter data are credited in
the bundled `LICENSE.txt` file.
YAeHMOP reference implementation: https://github.com/greglandrum/yaehmop

The English edition modifies language and naming only; it does not
reparametrize the scientific implementation.


## Desktop graphical interface (new in 0.4.0)

Run `python3 -m pyaehmop.gui` or, after installation, `pyaehmop-gui`.
Use the buttons to select a structure, default or external YAeHMOP
parameter file or an interactive shell editor, EHT/spin/SOC mode, CPU workers, quadrature, bands, DOS,
charges and output formats. Click **Preview structure**, then **Run
calculation**. The seven tabs display 3D geometry, band structure,
density of states, Mulliken charges, magnetic SCF, summaries and logs.
The GUI transparently saves the generated `.job` or `.spin` input and
runs the **existing Python numerical engine in a subprocess**. It does
not introduce a separate EHT or SOC solver.

Complete instructions: [GUI_GUIDE.md](GUI_GUIDE.md).

### External YAeHMOP parameter files

The friendly input now supports either:

```ini
parameters = default
```

or:

```ini
parameters = yaehmop.txt
```

The external file may contain classic YAeHMOP `PARAMETERS` records
(one row per species), or the full native `eht_parms.dat` table with
s/p/d/f and double-zeta coefficients. The code applies the specified
atomic shells without changing mathematical methodology; for
XYZ/CIF/ASE it can supply shells absent from the built-in table.
The same `parameters` setting works in `.spin` files.

Examples:
`examples_friendly/water_default_parameters.job`,
`examples_friendly/water_custom_parameters.job`, and
`SPIN_TEST/examples/Fe_wire_external_parameters.spin`.
