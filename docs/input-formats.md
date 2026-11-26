# Input formats and calculation recipes

[← Overview](../README.md) · [Installation](installation.md) · [Examples](examples.md)

Use a friendly `.job` for standard EHT or a `.spin` for the experimental spin/SOC solver. You may also pass a supported **native YAeHMOP input file** directly to the EHT CLI.

The friendly format begins with two nonempty lines: **structure format** and **structure file path**. Remaining lines are `key = value`. Comment lines start with `#` or `;`.

For examples copied from this page, put the corresponding structure files in the **same directory as the .job/.spin** unless an absolute path is provided. Relative file paths and output paths are resolved from the configuration file directory.

## 1. Molecule from XYZ: energy, DOS and charges

`water.job`:

```ini
xyz
water.xyz

parameters = default
calculations = energy, charges, dos
charge = 0
dos_sigma = 0.20
dos_points = 500
workers = 2
outputs = json, csv, png
plots = yes
output = results/water
```

`water.xyz`:

```text
3
Water molecule
O 0.0000000000 0.0000000000 0.0000000000
H 0.0000000000 0.0000000000 0.9572000000
H 0.9266272065 0.0000000000 -0.2399872084
```

Run:

```bash
python -m pyaehmop water.job
```

A bundled, validated variant is at [`pyaehmop/examples_friendly/water_full.job`](../pyaehmop/examples_friendly/water_full.job).

## 2. Only Mulliken charges

```ini
xyz
water.xyz
parameters = default
calculations = charges
workers = 2
outputs = json, csv
output = results/water_charges
```

The Hamiltonian still needs to be solved to calculate occupied populations, but unnecessary bands/DOS reports are not requested.

## 3. Periodic CIF: diamond bands and DOS

```ini
cif
diamond_primitive.cif

parameters = C 6 4 2 1.71 -21.4 2 1.625 -11.4
calculations = energy, charges, bands, dos
kmesh = 4 4 4
band_path = auto
divisions = 80
cutoff = 11
images = auto
workers = 4
outputs = json, csv, png
plots = yes
output = results/diamond
```

The carbon parameter line reproduces the custom EHT onsite parameters used by the bundled diamond comparison example; it is not equivalent to all alternative carbon tables.

**Important:** a band-path gap from a finite path is not a fully converged fundamental electronic gap.

Built-in example: [`diamond_cif.job`](../pyaehmop/examples_friendly/diamond_cif.job).

## 4. ASE structure: POSCAR / extended XYZ

```ini
ase
diamond_primitive.POSCAR

ase_format = vasp
calculations = energy, charges
kmesh = 2 2 2
workers = 2
outputs = json, csv
```

ASE must be installed; supported file formats vary with the installed ASE version. Structures with unsupported partial occupancies should be rejected rather than interpreted as ordinary integer-occupancy atoms.

## 5. Existing YAeHMOP input

To run a native YAeHMOP-style input through the independent Python engine:

```bash
python -m pyaehmop pyaehmop/inputs/XeF4
```

Or make a friendly wrapper:

```ini
yaehmop
../inputs/XeF4

parameters = default
calculations = energy, charges, orbitals
workers = 2
outputs = json, csv
output = results/xef4
```

The native keywords such as `GEOMETRY`, `LATTICE`, `PARAMETERS`, `BAND` and `COOP` belong to the original input format and are not translated. Native YAeHMOP input compatibility covers the **supported documented subset**, not every historical feature.

## 6. EHT parameter selection

### Built-in parameters

```ini
parameters = default
```

This uses the bundled parameter table for XYZ/CIF/ASE and preserves source-defined parameters when an original YAeHMOP input already specifies them.

### Classic YAeHMOP parameters in a separate file

`custom_yaehmop.txt`:

```text
; Classic YAeHMOP PARAMETERS lines
PARAMETERS
H 1 1 1 1.3000 -13.6000
O 8 6 2 2.2750 -32.3000 2 2.2750 -14.8000
END
```

Job:

```ini
xyz
water.xyz
parameters = custom_yaehmop.txt
calculations = energy, charges
```

The optional `PARAMETERS` and `END` markers and comments are recognized. Multi-species files are supported. Unspecified atomic species retain their available default shells.

### Native `eht_parms.dat` tables (double-zeta)

The bundled example [`pyaehmop/eht_parms.dat`](../pyaehmop/eht_parms.dat) contains standard YAeHMOP parameter rows. Each complete table row consists of:

```text
Element  Z  Valence  Nzeta  n  shell  Hii  zeta1  zeta2  coeff1  coeff2
```

For example, a user can include a locally available parameter table:

```ini
parameters = yaehmop_parameters.txt
```

The parser supports s/p/d/f shells and double-zeta coefficients. Repeated native shell rows follow the last-occurrence-wins rule. A file with no suitable parameters for the structure causes a clear error.

### Legacy inline record

```ini
parameters = C 6 4 2 1.71 -21.4 2 1.625 -11.4
```

This existing syntax is still accepted. A human-friendly alternative is the **GUI atomic-shell editor**: select an element, edit onsite energies and exponents with controls, and save directly in native YAeHMOP format.

### Physical caution

External EHT files specify **orbital energies and radial functions**, *not* Stoner exchange constants or SOC strengths. In a `.spin` input, `stoner_* = ...` and `soc_* = ...` must be configured independently.

## 7. Spin polarization, SOC, or both

`Fe_wire_both.spin`:

```ini
ase
Fe_wire.extxyz
parameters = default

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
output = results/Fe_wire_both
```

Run:

```bash
python -m pyaehmop.SPIN_TEST Fe_wire_both.spin
```

Choose one or both:

| Settings | Magnetic calculation |
|---|---|
| `spin=yes; soc=no` | Collinear Stoner spin-up/spin-down |
| `spin=no; soc=yes` | Nonmagnetic onsite SOC |
| `spin=yes; soc=yes` | Combined Stoner+SOC with `⟨σz⟩` projection |

The parameter strengths above are **illustrative** and not a universal Fe calibration. See [theory](theory.md#7-spin-polarization-and-spinorbit-coupling) for the collinear-z limitation.

[Validated spin examples](../pyaehmop/SPIN_TEST/examples).

## 8. Main friendly-input options

| Keyword | Example | Description |
|---|---|---|
| `calculations` | `energy,charges,dos` | Requested properties for EHT |
| `parameters` | `default` / `param.txt` | Atomic shell parameter source |
| `charge` | `0` | System charge |
| `electrons` | `36` | Explicit electron count (instead of charge) |
| `workers` | `4` | CPU worker processes |
| `quadrature` | `24 20 12` | STO integration orders |
| `cutoff` | `11` | Overlap cutoff, Å |
| `images` | `auto` / `3 3 3` | Periodic image extent |
| `kmesh` | `4 4 4` | K-grid integration |
| `kpoints` | `0 0 0; 0.5 0 0.5` | Explicit k points |
| `band_path` | `auto` | ASE high-symmetry path or explicit coordinates |
| `divisions` | `80` | Band path sampling |
| `dos_sigma` | `0.15` | Gaussian broadening, eV |
| `dos_points` | `600` | DOS mesh resolution |
| `bonds` | `1-2; 2-3@1,0,0` | COOP/COHP atom pairs and translations |
| `fragments` | `1:4; 2-5:32` | Fragment groups and electron counts for FMO |
| `outputs` | `json,csv,png` | Export formats |
| `plots` | `yes` | Save available PNG figures |
| `output` | `results/case` | Results directory |

For additional settings and validation details consult the [complete friendly-input manual](../pyaehmop/EASY_INPUT.md). The desktop GUI exposes the same parameters through fields, selectors and checkboxes.
