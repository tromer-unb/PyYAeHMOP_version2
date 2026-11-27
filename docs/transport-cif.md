# CIF and XYZ transport examples

[PDF tutorial](PyYAeHMOP_v2_tutorial.pdf) · [Device geometry guide](device-geometry.md)

## 1. Uniform graphene: two equivalent inputs

Clone and install Version 2 as described in [installation](installation.md).
Both commands use the same independent EHT-NEGF engine, identical carbon
parameters, energy reference, cutoff and k mesh. No external transport code
is needed.

From the repository root:

    python -m pyaehmop examples/transport/graphene.transport.job --workers 4
    python -m pyaehmop examples/transport/graphene_cif.transport.job --workers 4

These are full-resolution calculations (301 transverse k points, 401
energies). Start with the faster pair:

    python -m pyaehmop examples/tutorials/bandfit/pdz7_4B_smoke.transport.job --workers 4
    python -m pyaehmop examples/tutorials/bandfit/pdz7_4B_smoke_cif.transport.job --workers 4

| Setting | XYZ input | CIF input |
|---|---|---|
| File | graphene.xyz | graphene.cif |
| Periodic lattice | x, z | a=x, b=y (vacuum), c=z |
| Cell lengths (Å) | 2.459512, 4.26 | 2.459512, 30, 4.26 |
| lattice_vectors | 2.459512 0 0; 0 0 4.26 | omit (read from CIF) |
| images | 5 4 | 5 1 4 |
| transport_axis | 1 | 2 |
| transverse_kpoints | 301 | 301 1 |
| cutoff (Å) | 11 | 11 |

XYZ stores Cartesian coordinates; CIF stores fractional positions and cell
information. ASE reads the CIF as a 3D periodic crystal. The 30 Å vacuum
along y is longer than the 11 Å interaction cutoff, so there is no coupling
between the periodic sheets in this example. One k point is sufficient in
that noninteracting direction. Do not use a single transverse k point in
a physically dispersive direction. The transport_axis index is zero-based
and selects a *lattice vector*, not a Cartesian component. With an arbitrary
CIF, check the actual cell and vector order. For transport along x, use
transport_axis=0 and transverse_kpoints=1 301 in this CIF (y then z); XYZ
uses transport_axis=0 and transverse_kpoints=301 (z).

Do not copy the CIF's images/axis settings into a 2D XYZ job, or specify
lattice_vectors together with CIF. The CIF example is P1, fully occupied
and contains the same four carbon sites. For other materials use a
fully ordered structure with physically appropriate parameterization.
Validate cutoff, mesh, vacuum and energy reference.

## 2. All six BandFit models

In [bandfit](../examples/tutorials/bandfit/) every XYZ full and quick job
has a matching _cif.transport.job with a distinct output directory.
The carbon .dat files and K values are identical for each pair.
Use the corresponding model's archived absolute energy reference.
No refitting is needed to switch file format.

## 3. Explicit devices: CIF source files, XYZ solver input

The .device.job parser still requires three Cartesian XYZ files.
A CIF does not by itself define the left electrode, device and right
electrode or their absolute relative translations.

Each of the 15 shipped explicit-device examples now has a matching
cif_sources/<job-name>/ directory containing left.cif, device.cif,
right.cif and origin.json. The seven tutorials using run.device.job
use cif_sources/run/. These P1 CIFs share a padded, artificial coordinate
box and preserve all sites. They are **coordinate interchange files,
not physical bulk electrode CIF cells**. The true repeat vectors and
transverse periodicity are specified in the original .device.job.
origin.json stores the common Cartesian offset so the converter can
restore the three absolute positions and atom ordering. This also
preserves zero-based magnetic atom indices.

Example: pristine graphene

    python examples/transport/device_from_cif.py examples/devices/pristine.device.job
    python -m pyaehmop examples/devices/results/cif_roundtrip/pristine/converted.device.job --workers 4

Example: different electrodes, graphene–MoS2

    python examples/transport/device_from_cif.py examples/devices/graphene_mos2.device.job
    python -m pyaehmop examples/devices/results/cif_roundtrip/graphene_mos2/converted.device.job --workers 4

The converted XYZ files and job are generated locally under results/,
leaving the original inputs untouched. The calculations are written to
the generated job's calculated/ directory. Use the same converter for
vacancy, adsorbate, nanoribbon, spin, SOC and other shipped device jobs.
For example, for nanoribbon pass examples/tutorials/nanoribbon/run.device.job
and run examples/tutorials/nanoribbon/results/cif_roundtrip/run/converted.device.job.

For a new device, start with physically meaningful bulk CIFs, read
them using ASE, construct sufficiently thick principal layers, choose a
common transverse cell, and translate the contacts and device into one
global coordinate frame before exporting XYZ. The shipped converter
is only for the included CIF sets with their matching origin.json; it
does not automatically create an interface or align arbitrary CIFs.
See [device geometry](device-geometry.md).

## 4. Reproduce the validation

The [validation script](../examples/transport/validate_cif.py) runs all
six quick CIF jobs and 15 converted device jobs and compares the spectra
with the corresponding XYZ results (generating the XYZ results if absent).
From the installed checkout:

    python -m examples.transport.validate_cif

On a PBS cluster (optional):

    qsub -v PYTHON="$(command -v python)" examples/transport/validate_cif.pbs
    qsub -v PYTHON="$(command -v python)" examples/transport/validate_cif_full.pbs

[Archived numerical checks](../examples/transport/cif_reference/validation.json)
report maximum absolute differences in T, G (S) and I (A), and coordinate
round-trip error in Å. This is numerical equivalence of file formats,
not a new physical validation of the illustrative devices.

The full p-DZ7 4B comparison (301 k points and 401 energies) completed
with max |T_XYZ - T_CIF| = 8.60423e-15. See
[full spectra](../examples/transport/cif_reference/graphene_xyz_cif_full.dat)
and [full validation](../examples/transport/cif_reference/full_validation.json).
All six parameter sets have been tested at the quick resolution; the
other five full-resolution CIF jobs are supplied but have not been rerun.

All 15 explicit device round trips passed with zero recorded coordinate,
transmission, conductance and current differences using the 12-decimal
XYZ exporter. Final PBS job 62870 completed with exit 0. The validation
runner uses a separate Python process for each job to avoid accumulated
native-library state. An earlier single-process rerun was interrupted
by a segmentation fault; the fresh-process validation completed.
