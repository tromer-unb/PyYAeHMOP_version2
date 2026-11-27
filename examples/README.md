# All runnable examples — start here

Clone the repository, install with pip or Conda ([installation](../docs/installation.md)), then run the following commands **from the repository root**. No PBS or external transport software is required. Each transport row has two real executable job files, not just geometry files.

## Uniform periodic transport (XYZ and CIF)

Run either linked job with `python -m pyaehmop PATH --workers 4`. The CIF reads a 3D cell with 30 Å vacuum; the XYZ uses two explicit periodic vectors. Same material, parameters and energy reference. See [format comparison](../docs/transport-cif.md).

| Model | XYZ full | CIF full | XYZ quick | CIF quick |
|---|---|---|---|---|
| pdz7_2B | [pdz7_2B.transport.job](tutorials/bandfit/pdz7_2B.transport.job) | [pdz7_2B_cif.transport.job](tutorials/bandfit/pdz7_2B_cif.transport.job) | [pdz7_2B_smoke.transport.job](tutorials/bandfit/pdz7_2B_smoke.transport.job) | [pdz7_2B_smoke_cif.transport.job](tutorials/bandfit/pdz7_2B_smoke_cif.transport.job) |
| pdz7_4B | [pdz7_4B.transport.job](tutorials/bandfit/pdz7_4B.transport.job) | [pdz7_4B_cif.transport.job](tutorials/bandfit/pdz7_4B_cif.transport.job) | [pdz7_4B_smoke.transport.job](tutorials/bandfit/pdz7_4B_smoke.transport.job) | [pdz7_4B_smoke_cif.transport.job](tutorials/bandfit/pdz7_4B_smoke_cif.transport.job) |
| pdz7_6B | [pdz7_6B.transport.job](tutorials/bandfit/pdz7_6B.transport.job) | [pdz7_6B_cif.transport.job](tutorials/bandfit/pdz7_6B_cif.transport.job) | [pdz7_6B_smoke.transport.job](tutorials/bandfit/pdz7_6B_smoke.transport.job) | [pdz7_6B_smoke_cif.transport.job](tutorials/bandfit/pdz7_6B_smoke_cif.transport.job) |
| shared_2B | [shared_2B.transport.job](tutorials/bandfit/shared_2B.transport.job) | [shared_2B_cif.transport.job](tutorials/bandfit/shared_2B_cif.transport.job) | [shared_2B_smoke.transport.job](tutorials/bandfit/shared_2B_smoke.transport.job) | [shared_2B_smoke_cif.transport.job](tutorials/bandfit/shared_2B_smoke_cif.transport.job) |
| shared_4B | [shared_4B.transport.job](tutorials/bandfit/shared_4B.transport.job) | [shared_4B_cif.transport.job](tutorials/bandfit/shared_4B_cif.transport.job) | [shared_4B_smoke.transport.job](tutorials/bandfit/shared_4B_smoke.transport.job) | [shared_4B_smoke_cif.transport.job](tutorials/bandfit/shared_4B_smoke_cif.transport.job) |
| shared_6B | [shared_6B.transport.job](tutorials/bandfit/shared_6B.transport.job) | [shared_6B_cif.transport.job](tutorials/bandfit/shared_6B_cif.transport.job) | [shared_6B_smoke.transport.job](tutorials/bandfit/shared_6B_smoke.transport.job) | [shared_6B_smoke_cif.transport.job](tutorials/bandfit/shared_6B_smoke_cif.transport.job) |

Shared geometries: [graphene.xyz](transport/graphene.xyz) · [graphene.cif](transport/graphene.cif). Full p-DZ7 4B shortcut: [graphene.transport.job](transport/graphene.transport.job) · [graphene_cif.transport.job](transport/graphene_cif.transport.job).

## Independent electrodes and scattering region (15 paired cases)

Run **either** job directly with `python -m pyaehmop PATH --workers 4`. No manual conversion is required. The CIF jobs use ordered P1 region files with a shared padded box and origin.json, not three independent bulk CIF cells. The physical repeats remain in the job. The XYZ jobs use the original Cartesian geometries.

| Case | XYZ input | CIF input | CIF region files | Explanation / results |
|---|---|---|---|---|
| pristine | [pristine.device.job](devices/pristine.device.job) | [pristine_cif.device.job](devices/pristine_cif.device.job) | [L/D/R + origin](devices/cif_sources/pristine) | [Guide](devices/README.md) |
| vacancy | [vacancy.device.job](devices/vacancy.device.job) | [vacancy_cif.device.job](devices/vacancy_cif.device.job) | [L/D/R + origin](devices/cif_sources/vacancy) | [Guide](devices/README.md) |
| adsorbate | [adsorbate.device.job](devices/adsorbate.device.job) | [adsorbate_cif.device.job](devices/adsorbate_cif.device.job) | [L/D/R + origin](devices/cif_sources/adsorbate) | [Guide](devices/README.md) |
| asymmetric | [asymmetric.device.job](devices/asymmetric.device.job) | [asymmetric_cif.device.job](devices/asymmetric_cif.device.job) | [L/D/R + origin](devices/cif_sources/asymmetric) | [Guide](devices/README.md) |
| graphene_mos2 | [graphene_mos2.device.job](devices/graphene_mos2.device.job) | [graphene_mos2_cif.device.job](devices/graphene_mos2_cif.device.job) | [L/D/R + origin](devices/cif_sources/graphene_mos2) | [Guide](devices/README.md) |
| spin | [spin.device.job](devices/spin.device.job) | [spin_cif.device.job](devices/spin_cif.device.job) | [L/D/R + origin](devices/cif_sources/spin) | [Guide](devices/README.md) |
| soc | [soc.device.job](devices/soc.device.job) | [soc_cif.device.job](devices/soc_cif.device.job) | [L/D/R + origin](devices/cif_sources/soc) | [Guide](devices/README.md) |
| spin_soc | [spin_soc.device.job](devices/spin_soc.device.job) | [spin_soc_cif.device.job](devices/spin_soc_cif.device.job) | [L/D/R + origin](devices/cif_sources/spin_soc) | [Guide](devices/README.md) |
| graphene_x | [run.device.job](tutorials/graphene_x/run.device.job) | [run_cif.device.job](tutorials/graphene_x/run_cif.device.job) | [L/D/R + origin](tutorials/graphene_x/cif_sources/run) | [Guide](tutorials/graphene_x/README.md) |
| nanoribbon | [run.device.job](tutorials/nanoribbon/run.device.job) | [run_cif.device.job](tutorials/nanoribbon/run_cif.device.job) | [L/D/R + origin](tutorials/nanoribbon/cif_sources/run) | [Guide](tutorials/nanoribbon/README.md) |
| wide_vacancy | [run.device.job](tutorials/wide_vacancy/run.device.job) | [run_cif.device.job](tutorials/wide_vacancy/run_cif.device.job) | [L/D/R + origin](tutorials/wide_vacancy/cif_sources/run) | [Guide](tutorials/wide_vacancy/README.md) |
| potential_barrier | [run.device.job](tutorials/potential_barrier/run.device.job) | [run_cif.device.job](tutorials/potential_barrier/run_cif.device.job) | [L/D/R + origin](tutorials/potential_barrier/cif_sources/run) | [Guide](tutorials/potential_barrier/README.md) |
| spin_parallel | [run.device.job](tutorials/spin_parallel/run.device.job) | [run_cif.device.job](tutorials/spin_parallel/run_cif.device.job) | [L/D/R + origin](tutorials/spin_parallel/cif_sources/run) | [Guide](tutorials/spin_parallel/README.md) |
| spin_antiparallel | [run.device.job](tutorials/spin_antiparallel/run.device.job) | [run_cif.device.job](tutorials/spin_antiparallel/run_cif.device.job) | [L/D/R + origin](tutorials/spin_antiparallel/cif_sources/run) | [Guide](tutorials/spin_antiparallel/README.md) |
| graphene_mos2_graphene | [run.device.job](tutorials/graphene_mos2_graphene/run.device.job) | [run_cif.device.job](tutorials/graphene_mos2_graphene/run_cif.device.job) | [L/D/R + origin](tutorials/graphene_mos2_graphene/cif_sources/run) | [Guide](tutorials/graphene_mos2_graphene/README.md) |

## Electronic structure and spin/SOC

These are the original molecular/crystalline inputs, not transport jobs. Standard XYZ is Cartesian and nonperiodic for molecules, CIF provides crystal cell and fractional sites, ASE/POSCAR retain their lattice. All remain in the full tutorial.

| Case | Executable input |
|---|---|
| diamond_ase | [diamond_ase.job](../pyaehmop/examples_friendly/diamond_ase.job) |
| diamond_cif | [diamond_cif.job](../pyaehmop/examples_friendly/diamond_cif.job) |
| diamond_poscar | [diamond_poscar.job](../pyaehmop/examples_friendly/diamond_poscar.job) |
| diamond_reference | [diamond_reference.job](../pyaehmop/examples_friendly/diamond_reference.job) |
| water_charges | [water_charges.job](../pyaehmop/examples_friendly/water_charges.job) |
| water_custom_parameters | [water_custom_parameters.job](../pyaehmop/examples_friendly/water_custom_parameters.job) |
| water_default_parameters | [water_default_parameters.job](../pyaehmop/examples_friendly/water_default_parameters.job) |
| water_full | [water_full.job](../pyaehmop/examples_friendly/water_full.job) |
| xef4_original | [xef4_original.job](../pyaehmop/examples_friendly/xef4_original.job) |
| C_atom_both | [C_atom_both.spin](../pyaehmop/SPIN_TEST/examples/C_atom_both.spin) |
| Fe_wire_both | [Fe_wire_both.spin](../pyaehmop/SPIN_TEST/examples/Fe_wire_both.spin) |
| Fe_wire_external_parameters | [Fe_wire_external_parameters.spin](../pyaehmop/SPIN_TEST/examples/Fe_wire_external_parameters.spin) |
| Fe_wire_soc | [Fe_wire_soc.spin](../pyaehmop/SPIN_TEST/examples/Fe_wire_soc.spin) |
| Fe_wire_soc_exaggerated | [Fe_wire_soc_exaggerated.spin](../pyaehmop/SPIN_TEST/examples/Fe_wire_soc_exaggerated.spin) |
| Fe_wire_spin | [Fe_wire_spin.spin](../pyaehmop/SPIN_TEST/examples/Fe_wire_spin.spin) |

Run ordinary jobs with `python -m pyaehmop PATH`; run .spin jobs with `python -m pyaehmop.SPIN_TEST PATH`.

## Full illustrated manual and validation

- [Complete PDF with all input files and device plots](../docs/PyYAeHMOP_v2_tutorial.pdf)
- [Original tutorial index and numerical reference results](tutorials/README.md)
- [Direct CIF device validation](transport/cif_reference/direct_device_validation.json)
- [Uniform XYZ/CIF validation](transport/cif_reference/validation.json)
- [Parameters](../docs/parameters.md) and [electrode construction](../docs/device-geometry.md)

All original examples are retained. New devices are illustrative frozen-Hamiltonian models; I(V) is not self-consistent finite-bias NEGF.
