# Extended Hückel parameter sources and formats

[PDF tutorial](PyYAeHMOP_v2_tutorial.pdf) · [BandFit reproduction](../examples/tutorials/bandfit/README.md)

PyYAeHMOP computes Slater-orbital overlaps and the extended Hückel
Hamiltonian directly. It does not invoke the YAeHMOP executable.
The Hamiltonian uses onsite energies Hii and Hückel off-diagonal elements
proportional to Sij, with K and optional weighted Hij. The following
parameter sources are alternative choices, not sequential requirements.

## 1. Bundled YAeHMOP table

    parameters = default
    K = 1.75
    weighted_hij = true

The bundled table is pyaehmop/eht_parms.dat. K=1.75 is a typical EHT
choice, not a universal optimum. Default parameters are a starting point,
not a guarantee of transferability for a particular material.

## 2. Your own YAeHMOP parameter file

    parameters = my_parameters.dat

Relative paths are resolved against the job directory. The native
eht_parms.dat row has 11 fields:

    element Z valence n_zeta n shell Hii_eV zeta1 zeta2 c1 c2

Example (illustrative conventional C values; not the fitted graphene table):

    C 6 4 1 2 s -21.4 1.625 0 1 0
    C 6 4 1 2 p -11.4 1.625 0 1 0

For a contracted double-zeta p shell, use n_zeta=2 and supply both positive
exponents and both nonzero contraction coefficients. The solver normalizes
the contraction. Exponents use inverse Bohr and energies eV. The classic
YAeHMOP PARAMETERS format is also accepted, e.g.:

    C 6 4 2 1.71 -21.4 2 1.625 -11.4

That record is element, atomic number, valence, Ns, zeta_s, Hss,
Np, zeta_p, Hpp. Use the 11-column table for double-zeta shells.
External files replace the shell set for each included species; omitted
species retain bundled defaults. Include all required shells for overridden
species. For a C/Mo/S device, a carbon-only file therefore does not
calibrate the Mo or S parameters. The device uses one common K and weighted
formula; do not silently mix separate fitted K values at an interface.

## 3. YAeHMOP BandFit parameters

BandFit is a separate tool for fitting EHT bands to an electronic-structure
reference. Its outputs must be converted to the native table (or use the
already exported .dat examples). Use the matching K in the job; K is not
read from a native parameter row. For the graphene p-DZ7 four-band example:

    parameters = ../transport/carbon_pdz7_4B.dat
    K = 2.0068527138311567
    weighted_hij = true
    energy_reference_eV = -18.800896640571224

This is the exact reference in the original 4B job; do not reuse it for a
new material automatically. For all six sets and their matching absolute
energy references see [BandFit jobs](../examples/tutorials/bandfit/README.md).
The example fits 2, 4, or 6 selected frontier bands, not 2, 4 or 6 orbitals.
Double-zeta describes radial contraction, not two independent angular basis
functions. The 4B calibration kept p onsite energy and all radial data
fixed and optimized s onsite, K and a band-alignment offset. No transport
data were fitted. The alignment shift used in fitting is distinct from
the absolute energy reference supplied for transport.

To reproduce the fit, the repository includes the six frontier DFT
energies at 338 k points and the original initial parameters. Clone
BandFit separately and run refit.py (see linked tutorial). TranSIESTA and
SIESTA are not required to run any shipped PyYAeHMOP calculation.

For new materials, validate band ordering, dispersion, gap, wavefunction
character and transferability. For heterojunctions, also validate absolute
band alignment, cross-species couplings and charge redistribution.
