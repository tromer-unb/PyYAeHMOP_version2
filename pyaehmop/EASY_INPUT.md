# EASY_INPUT — English friendly configuration format

This input format requires **no Python programming**. It chooses which
properties to calculate and how to read a structure (YAeHMOP, XYZ, CIF,
or other formats through ASE).

## Basic format

The first line is the structure format; the second is the structure
filename or path. The remaining lines are optional `key = value`
settings.

Example: `my_calculation.job`

```ini
xyz
water.xyz

calculations = energy, charges, dos
charge = 0
workers = 4
dos_sigma = 0.15
plots = yes
outputs = json, csv
output = results/water
```

The first two lines may instead be written:

```ini
format = xyz
file = water.xyz
```

Paths are interpreted relative to the **.job file's directory**,
not necessarily the terminal's working directory. Lines starting with
`#` or `;` are comments.

Run from the English edition's parent directory:

```bash
cd PyYAeHMOP
python3 -m pyaehmop pyaehmop/examples_friendly/water_full.job
```

The engine auto-detects the friendly input. Traditional YAeHMOP input
files still work without conversion.

## Supported calculation types

The `calculations` option accepts a comma-separated list.

| Value | Output |
|---|---|
| `energy` | Electronic energy; molecular HOMO, LUMO and orbital gap |
| `charges` | Mulliken charges and orbital populations |
| `orbitals` | Molecular orbital energy levels and occupations |
| `dos` | Gaussian-broadened density of states |
| `bands` | Electronic bands along a reciprocal-space path |
| `bonds` | COOP/COHP analyses of specified atom pairs |
| `fmo` | Molecular fragment orbital analysis |
| `matrices` | Molecular overlap, Hamiltonian and coefficient matrices |
| `basic` | Equivalent to `energy, charges` |
| `all` | All relevant implemented properties for the given input |

Requesting only `charges` suppresses unnecessary band and DOS outputs.
The underlying Hamiltonian must still be diagonalized to compute electron
occupations and charges.

## Configuration options

| Key | Example | Notes |
|---|---|---|
| `workers` | `4` | Number of CPU processes |
| `charge` | `0` | Total molecular/cell charge |
| `electrons` | `36` | Alternative to specifying charge |
| `weighted` | `yes` | Weighted EHT coupling |
| `k` | `1.75` | EHT coupling constant |
| `quadrature` | `24 20 12` | Quadrature orders for STO integrals |
| `cutoff` | `11` | Periodic overlap cutoff in angstroms |
| `images` | `auto` or `5 5 5` | Lattice images in each periodic direction |
| `kmesh` | `4 4 4` | Uniform k-point integration grid |
| `kpoints` | `0 0 0; 0.5 0 0.5` | Explicit integration k points; optional fourth field = weight |
| `band_path` | `G:0,0,0; X:0.5,0,0.5` | Band path in reduced coordinates, or `auto` |
| `divisions` | `40` | Points per path segment (or total for ASE auto path) |
| `dos_sigma` | `0.15` | DOS broadening width in eV |
| `dos_points` | `600` | Number of DOS energy samples |
| `bonds` | `1-2; 2-3@1,0,0` | Pair indices and optional lattice translation |
| `fragments` | `1:4; 2-5:32` | Atom groups and electron counts for FMO |
| `parameters` | `C 6 4 2 1.71 -21.4 2 1.625 -11.4` | Custom YAeHMOP-style EHT shell record |
| `outputs` | `json,csv,png` | Output formats |
| `plots` | `yes` | Generate available PNG plots |
| `output` | `results/example` | Output directory |
| `ase_format` | `vasp` | Optional explicit ASE input format |
| `frame` | `0` | Select a single ASE trajectory frame |

Boolean settings use **yes** or **no**.
Custom parameter data can also be provided as a path to a text file
with one atom-species record per line.

## Four structure formats

### 1. YAeHMOP (legacy structure with English controls)

```ini
yaehmop
../inputs/XeF4

calculations = energy, charges, orbitals, fmo, bonds
workers = 2
outputs = json, csv
output = ../results/friendly_xef4
```

The standard YAeHMOP file is read by the Python parser and processed by
the independent Python engine. Existing scientific syntax such as
`GEOMETRY`, `LATTICE`, `PARAMETERS`, `BAND` and `COOP` is
intentionally not translated; it is part of the external file standard.

### 2. XYZ molecules

Source `water.xyz`:

```text
3
Water molecule
O  0.0000000000 0.0000000000 0.0000000000
H  0.0000000000 0.0000000000 0.9572000000
H  0.9266272065 0.0000000000 -0.2399872084
```

Job file:

```ini
xyz
water.xyz
calculations = energy, charges, dos
dos_sigma = 0.20
plots = yes
outputs = json, csv
```

Standard XYZ has no lattice vectors, so this is treated as a molecule.
For periodic extended-XYZ, use `ase`.

### 3. CIF crystals

```ini
cif
diamond_primitive.cif

calculations = energy, charges, bands, dos
kmesh = 4 4 4
band_path = auto
divisions = 80
cutoff = 11
images = auto
workers = 4
plots = yes
outputs = json, csv
```

CIF structures are read through ASE, which expands space-group operations
when applicable. Disordered or partially occupied sites are rejected
when detected. Review whether the chosen primitive cell and band path
are appropriate for the intended comparison.

### 4. ASE-supported structures

```ini
ase
diamond_primitive.POSCAR

ase_format = vasp
calculations = energy, charges
kmesh = 2 2 2
workers = 4
```

Supported formats depend on the installed ASE version. Common examples
include POSCAR, CONTCAR, extended XYZ, CIF, and trajectory files.

## Reproducing the official diamond example

The original diamond example uses custom carbon EHT parameters that
are not identical to the generic carbon parameter table.

The included `diamond_cif.job` explicitly defines:

```ini
parameters = C 6 4 2 1.71 -21.4 2 1.625 -11.4
```

The translated CIF and translated YAeHMOP reference examples produce
the same 81 band energies up to about 1.3×10⁻⁷ eV.

## Output files

Depending on selected calculations:

| File | Content |
|---|---|
| `summary.json` | Selected properties with English field names |
| `molecular_levels.csv` | Orbital levels and occupations |
| `energies.csv` | Key energy properties |
| `mulliken.csv` | Atomic Mulliken charges |
| `bands.csv` | K points and band energies |
| `dos.csv` | Gaussian DOS |
| `bonds.csv` | COOP/COHP data |
| `fmo_levels.csv` | Fragment orbital levels |
| `bands.png` | Band diagram |
| `dos.png` | Density-of-states plot |
| `charges.png` | Mulliken charge diagram |

Only requested properties are exported.

## Dependencies

The basic engine needs NumPy and SciPy; CIF/ASE need ASE; PNG plots need
Matplotlib:

```bash
python3 -m pip install './pyaehmop[all]'
```

Or install NumPy/SciPy directly with the provided `requirements.txt`.

## Reproducibility and limitations

- Default k meshes for new periodic inputs are preliminary, not converged.
- EHT parameters must be appropriate for each chemical species/material.
- The model is not DFT and does not provide reliable total energies without
  the qualifications in the main README.
- Features outside the supported YAeHMOP subset fail explicitly or are
  ignored with warnings. Check warnings in each calculation.
- For experimental spin and spin-orbit calculations, use the isolated
  English `SPIN_TEST/` interface with `spin=yes` and/or `soc=yes`.


## Built-in or external YAeHMOP orbital parameters (0.4.0)

In either a `.job` input or a `.spin` input, specify exactly one of:

```ini
parameters = default
```

```ini
parameters = water_parameters.txt
```

The parameter file path is relative to the input file or may be
absolute. Supported parameter files are (1) classic YAeHMOP PARAMETERS
records, with or without `PARAMETERS`/`END` section markers, and
(2) native YAeHMOP `eht_parms.dat` rows including double-zeta shells.

Example classic parameter file:

```text
H 1 1 1 1.3000 -13.6000
O 8 6 2 2.2750 -32.3000 2 2.2750 -14.8000
```

The packaged native table is `pyaehmop/eht_parms.dat`; use its path
as the value of `parameters` to load the entire table. Multiple species
are allowed; entries for elements absent from your structure are
ignored. The supplied species' shells are replaced by those in the
selected file, while other species retain the built-in parameters.

A third, backward-compatible form is an inline classic YAeHMOP record:

```ini
parameters = C 6 4 2 1.71 -21.4 2 1.625 -11.4
```

When the geometry is a native YAeHMOP file, `default` retains the
parameters already defined by that source file, rather than overriding
them unexpectedly. For a GUI without manual editing, launch
`python3 -m pyaehmop.gui` (see `GUI_GUIDE.md`).
