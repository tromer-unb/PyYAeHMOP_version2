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

## 3. Explicit devices: direct XYZ and CIF jobs

Each of the 15 devices now has two runnable jobs in the same directory.
See the [complete catalog](../examples/README.md) for direct links.
Run either job with the same command, without manual conversion:

    python -m pyaehmop examples/devices/pristine.device.job --workers 4
    python -m pyaehmop examples/devices/pristine_cif.device.job --workers 4

Different electrodes:

    python -m pyaehmop examples/devices/graphene_mos2.device.job --workers 4
    python -m pyaehmop examples/devices/graphene_mos2_cif.device.job --workers 4

The CIF job differs only in geometry selection and output path. Example:

    structure_format = cif
    cif_origin = cif_sources/pristine/origin.json
    electrode_left = cif_sources/pristine/left.cif
    scattering_region = cif_sources/pristine/device.cif
    electrode_right = cif_sources/pristine/right.cif
    transverse_vectors = 2.459512 0 0
    left_repeat = 0 0 8.52
    right_repeat = 0 0 8.52

In the XYZ job, structure_format defaults to xyz, cif_origin is absent,
and the three paths end in .xyz. All remaining physical and numerical
parameters are identical. The CIF reader restores global Cartesian
positions internally (including the common offset), preserves atom
order, and ignores the artificial box periodicity. Only the explicit
repeat/transverse vectors define the physical periodicity of devices.

### Why the common box matters

All shipped three-region CIFs are fully occupied P1 **coordinate
interchange boxes**, not physical bulk electrode unit cells. They share
one orthogonal padded box and origin.json. The box keeps all atom sites
away from periodic wrapping; its lengths must not be substituted for
left_repeat or right_repeat. The parser rejects non-P1, partially occupied
sites and mismatched cells. For magnetic jobs preserve the site ordering:
exchange_eV indexes atoms starting at zero.

The seven individual tutorials have run.device.job and run_cif.device.job,
with geometries in cif_sources/run/. The eight shared examples have
<name>.device.job and <name>_cif.device.job, with geometries in
cif_sources/<name>/. Results of CIF jobs are written under results/cif_direct/.

For a new device, start with physically meaningful bulk CIFs, expand
symmetry and build sufficiently thick principal layers, choose a common
transverse cell, place all three regions in one Cartesian frame, and
export either XYZ or three ordered P1 CIFs in the same padded box with
matching origin metadata. A bulk CIF alone does not build a heterojunction.
See [device geometry](device-geometry.md). The optional previous converter
remains available for comparing exported XYZ files:

    python examples/transport/device_from_cif.py examples/devices/pristine.device.job

### Verification

The direct input validation compares the CIF and XYZ coordinates, all
three spectra and collinear current channels where present. See
[direct validation](../examples/transport/cif_reference/direct_device_validation.json).
Run from the repository root:

    python -m examples.transport.validate_direct_cif

Or use optional PBS:

    qsub -v PYTHON="$(command -v python)" examples/transport/validate_direct_cif.pbs

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

Native CIF device validation completed for all 15 cases (PBS 62874, exit 0). Maximum recorded coordinate, T, G and I differences were zero; collinear current channels also matched. Only the geometry reader was added; the NEGF Hamiltonian and solver are unchanged.
