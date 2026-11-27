# Reproducible examples: PyYAeHMOP Version 2

This is the entry point to all runnable examples. No computation is triggered
by GitHub. Each linked job is a real input file, not pseudocode. Numerical
reference outputs were generated on a PBS-enabled research cluster; no access to that cluster is required.

## Installation and common commands

From a Linux terminal:

    git clone https://github.com/tromer-unb/PyYAeHMOP_version2.git
    cd PyYAeHMOP_version2
    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -e ".[all]"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

Run a single input from the repository root:

    python -m pyaehmop examples/tutorials/graphene_x/run.device.job --workers 4

For PBS, submit from the same directory. Set PYTHON to your environment's
absolute Python executable; use your site's queue/account settings if required:

    qsub -v SUITE=devices,PYTHON="$(command -v python)" examples/tutorials/run_suite.pbs
    qsub -v SUITE=electronic,PYTHON="$(command -v python)" examples/tutorials/run_suite.pbs
    qsub -v SUITE=bandfit,PYTHON="$(command -v python)" examples/tutorials/run_suite.pbs
    qsub -v SUITE=full_bandfit,PYTHON="$(command -v python)" examples/tutorials/run_suite.pbs

The last command is substantially heavier (301 transverse k points, 401
energies, six parameterizations). Do not run on a cluster login node.
Results are written to the directories specified in each job (ignored by Git).
PBS status files are written to this directory's results/. The archive script
copies compact spectra and summary files to reference/:

    qsub -v PYTHON="$(command -v python)" examples/tutorials/archive.pbs

For detailed electrode placement see [geometry guide](../../docs/device-geometry.md); for parameter selection see [parameter guide](../../docs/parameters.md). Download the [illustrated PDF tutorial](../../docs/PyYAeHMOP_v2_tutorial.pdf).

## 1. Molecular and periodic electronic structure

All inputs and structure files are in [examples_friendly](../../pyaehmop/examples_friendly).
The command is python -m pyaehmop PATH --workers 4.

| Feature | Executable input |
|---|---|
| Molecular energy, levels, orbitals, Mulliken charges, DOS | [water_full.job](../../pyaehmop/examples_friendly/water_full.job) |
| Charges only | [water_charges.job](../../pyaehmop/examples_friendly/water_charges.job) |
| Default parameter table | [water_default_parameters.job](../../pyaehmop/examples_friendly/water_default_parameters.job) |
| External parameter table | [water_custom_parameters.job](../../pyaehmop/examples_friendly/water_custom_parameters.job) |
| Crystal bands, DOS and charges; CIF | [diamond_cif.job](../../pyaehmop/examples_friendly/diamond_cif.job) |
| ASE extended XYZ | [diamond_ase.job](../../pyaehmop/examples_friendly/diamond_ase.job) |
| POSCAR | [diamond_poscar.job](../../pyaehmop/examples_friendly/diamond_poscar.job) |
| Native YAeHMOP crystal input | [diamond_reference.job](../../pyaehmop/examples_friendly/diamond_reference.job) |
| Native molecular input, FMO, bond analysis | [xef4_original.job](../../pyaehmop/examples_friendly/xef4_original.job) |

For periodic COOP/COHP see the existing [native inputs](../../pyaehmop/inputs)
and [extension tests](../../pyaehmop/tests/test_extensions.py).
The GUI can load the same ordinary molecular/crystal inputs, but has no
device-building tab. Start it on a graphical Linux desktop with
python -m pyaehmop.gui; see [GUI guide](../../docs/gui.md).

## 2. Self-consistent periodic spin and atomic SOC

Command: python -m pyaehmop.SPIN_TEST PATH.

| Feature | Input |
|---|---|
| Collinear Stoner SCF in Fe wire | [Fe_wire_spin.spin](../../pyaehmop/SPIN_TEST/examples/Fe_wire_spin.spin) |
| Atomic SOC only | [Fe_wire_soc.spin](../../pyaehmop/SPIN_TEST/examples/Fe_wire_soc.spin) |
| Spin + SOC | [Fe_wire_both.spin](../../pyaehmop/SPIN_TEST/examples/Fe_wire_both.spin) |
| External EHT parameters | [Fe_wire_external_parameters.spin](../../pyaehmop/SPIN_TEST/examples/Fe_wire_external_parameters.spin) |
| Exaggerated SOC demonstration | [Fe_wire_soc_exaggerated.spin](../../pyaehmop/SPIN_TEST/examples/Fe_wire_soc_exaggerated.spin) |
| Isolated C atom | [C_atom_both.spin](../../pyaehmop/SPIN_TEST/examples/C_atom_both.spin) |

These use the original constrained collinear-z Stoner solver, not open-system
NEGF spin self-consistency. [Scalar reference outputs](reference/electronic_summary.json)
and [PBS status](reference/electronic_status.txt) include all 15 electronic examples.

## 3. Uniform graphene and BandFit

See [bandfit tutorial](bandfit/README.md) for all six parameter files, individual
full-resolution jobs and quick jobs, and the reproducible partial fit to the
338-point SIESTA DFT reference. The original TranSIESTA comparison remains
[here](../transport/README.md).

## 4. Independent electrodes and central region

Each device uses a .device.job input and three XYZ files in the SAME coordinate
frame. Lead repeat vectors point from left to right. The common transverse
periodicity can be absent (1D) or nonzero (2D/3D). The solver rejects omitted
long-range lead couplings; increase the thickness of the contact layers when
raising cutoff. See [the device manual](../devices/README.md) for every key.

| Implementation | Input | Archived data |
|---|---|---|
| Pristine graphene / transport Z | [pristine](../devices/pristine.device.job) | [T/G/I](../devices/benchmarks) |
| Independent central vacancy | [vacancy](../devices/vacancy.device.job) | [T/G/I](../devices/benchmarks) |
| CO only in the device | [adsorbate](../devices/adsorbate.device.job) | [T/G/I](../devices/benchmarks) |
| Dissimilar leads, C/Mo/S chain | [asymmetric](../devices/asymmetric.device.job) | [T/G/I](../devices/benchmarks) |
| Lateral graphene–MoS2 | [graphene_mos2](../devices/graphene_mos2.device.job) | [T/G/I](../devices/benchmarks) |
| Frozen spin-polarized transport | [spin](../devices/spin.device.job) | [T/G/I](../devices/benchmarks) |
| Spinor SOC transport | [soc](../devices/soc.device.job) | [T/G/I](../devices/benchmarks) |
| Frozen spin + SOC | [spin_soc](../devices/spin_soc.device.job) | [T/G/I](../devices/benchmarks) |
| Graphene transport X | [graphene_x](graphene_x/run.device.job) | [T/G/I](graphene_x/reference/) |
| Finite-width nanoribbon (1D) | [nanoribbon](nanoribbon/run.device.job) | [T/G/I](nanoribbon/reference/) |
| One vacancy in 4x wider transverse cell | [wide_vacancy](wide_vacancy/run.device.job) | [T/G/I](wide_vacancy/reference/) |
| Central frozen 0.5 eV potential | [potential_barrier](potential_barrier/run.device.job) | [T/G/I](potential_barrier/reference/) |
| Parallel magnetic contacts | [spin_parallel](spin_parallel/run.device.job) | [T/G/I](spin_parallel/reference/) |
| Antiparallel magnetic contacts | [spin_antiparallel](spin_antiparallel/run.device.job) | [T/G/I](spin_antiparallel/reference/) |
| Graphene–MoS2–graphene channel | [graphene_mos2_graphene](graphene_mos2_graphene/run.device.job) | [T/G/I](graphene_mos2_graphene/reference/) |

The seven new examples have individual directories with run.device.job,
left.xyz, device.xyz, right.xyz and magnetic/potential JSON when needed.
Each reference directory contains transmission.dat, conductance.dat,
current_voltage.dat and spectra.svg (three panels). Collinear examples also
provide current_spin.dat. The [summary](reference/devices_summary.json) records
atom counts, transmission and current. [PBS status](reference/devices_status.txt).

The original eight jobs have shared geometry and parameter files in
../devices, and individual curves in ../devices/benchmarks; do not copy just
the input file without its dependencies. All jobs reference the native
carbon p-DZ7 4-band parameter file in ../transport, or default parameters
where specified. Relative paths resolve against the job directory.

## 5. Interpretation and limitations

These new examples are engineering demonstrations using illustrative,
unrelaxed structures and reduced numerical grids. They are not independent
DFT-NEGF predictions. The nanoribbon is unpassivated and nonmagnetic;
vacancy/adsorbate examples remain periodically repeated transversely.
The wider vacancy reduces the concentration but is not an isolated-defect
convergence study. The graphene X and Z cells have different transverse
periods: compare transmission/current per transverse width, not raw
currents. Heterojunctions use strained commensurate cells, default Mo/S
parameters, and an explicit reference energy without validated band
alignment or charge transfer. Magnetic contact examples have prescribed
exchange potentials, not spontaneous graphene magnetism. SOC reports total
spinor transmission, not spin-flip-resolved channels. I(V) is a frozen
equilibrium Landauer approximation, not self-consistent finite-bias NEGF.

All normal EHT and original BandFit/TranSIESTA examples remain unchanged.
No new engine functionality has been introduced in this tutorial update.

## PBS execution record (2026-11-26)

| Suite | PBS job | Outcome |
|---|---|---|
| 15 independent device examples | 62855 | 15/15 passed, 53 s |
| Six quick BandFit transport examples | 62856 | 6/6 passed, 17 s |
| Nine electronic + six periodic/molecular spin examples | 62857 | 15/15 passed, 25 s |
| Six full-resolution transport re-runs | 62858 | 6/6 passed, 9 min 29 s |
| Six BandFit partial optimizations | 62859 | 6/6 reproduced exactly, 31 s |
| Compact tutorial reference archive | 62860 | Passed, 3 s |

The full-resolution rerun reproduced all six archived curves within 5e-11 absolute T at the 41 shared sample energies. See [full spectra](bandfit/reference/full_transmission.csv) and [reproduction metrics](bandfit/reference/reproduction_metrics.json). The original comparison with TranSIESTA remains in ../transport.

## Rebuild the PDF

The [PDF manual](../../docs/PyYAeHMOP_v2_tutorial.pdf) is generated from shipped inputs and archived numerical tables. With NumPy and Matplotlib installed, run python docs/build_tutorial.py from the checkout. No electronic or transport calculation is required for this documentation step.

## CIF alternatives for transport

[CIF and XYZ transport guide](../../docs/transport-cif.md) explains the paired uniform inputs and how to convert the three-region CIF sources for all explicit device examples. The original XYZ examples remain unchanged.
