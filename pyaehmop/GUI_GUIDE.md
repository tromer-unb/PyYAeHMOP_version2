# PyYAeHMOP Studio — Desktop graphical interface

The 0.4.0 edition adds a **real Tkinter desktop application**. It uses
the original numerical EHT engine and the unchanged Stoner/SOC scientific
methods. It is not a browser prototype or a static interface mockup.

## Start the application

On the authorized Linux desktop:

```bash
cd PyYAeHMOP
python3 -m pyaehmop.gui
```

If the package is installed through pip, the following command is also
available:

```bash
pyaehmop-gui
```

Dependencies: Python, NumPy, SciPy, Matplotlib and ASE. Tkinter must be
available from your operating system. On Debian/Ubuntu, if necessary:
`sudo apt install python3-tk`. To install Python dependencies:

```bash
python3 -m pip install './pyaehmop[gui]'
```

A graphical desktop session (an X11 or compatible display) is required.

## GUI workflow: no manual input file needed

1. Click **Open structure** and select a YAeHMOP input, XYZ, CIF, POSCAR,
   extended XYZ or any ASE-readable structure. Select the format in the
   **Structure & parameters** group.
2. Set **EHT basis parameters** to **default**, or select **file** and
   click **…** to browse a classic YAeHMOP parameter file or native
   `eht_parms.dat` table. Alternatively click **Edit atomic shell table…**
   to edit every Slater-shell parameter (Hii, n, zeta1/2 and c1/2) in
   a table and save it as a native YAeHMOP file. No hand-written
   parameter file is required.
3. Select **Standard EHT**, **Spin SCF**, **Spin-orbit only**, or
   **Spin SCF + SOC**.
4. For Standard EHT, check **Energy**, **Mulliken charges**, **Bands**,
   **DOS**, **COOP/COHP bonds**, **FMO**, **Orbitals** or **Matrices**.
   For spin/SOC, supply the Stoner/SOC shell strengths; its existing
   model calculates spin DOS automatically and bands if a path is given.
5. Select k mesh, band path, cutoffs, quadrature, CPU workers and
   DOS broadening as appropriate. Set output directory and export options.
6. Click **Preview structure** for a 3D view of atoms and the cell.
   Click **Run calculation** to compute. A spinner indicates progress.
7. Switch between **Structure**, **Bands**, **DOS**, **Charges**,
   **SCF**, **Summary** and **Console** tabs. Every plot has toolbar
   controls for zooming, panning and saving a PNG.

The **Edit atomic shell table…** window loads default shells for the
selected structure, permits importing a YAeHMOP table, and includes
**Add new shell**, **Remove selected shell**, **Apply changes** and
**Save YAeHMOP table and use** buttons. The saved custom table is
selected automatically in the main GUI. The original parameter table
is never overwritten unless you explicitly choose that same path in
an operating-system save dialog.

The GUI transparently creates `gui_calculation.job` or
`gui_calculation.spin` in your selected results directory. You
**do not need to edit or write** these files. Use **Save input** if you
want an editable, portable record of your settings; use **Open .job /
.spin** to reload one.

The calculation runs in a separate Python process. Tkinter stays
responsive and the **Stop** button terminates the child process on
request. The GUI does not implement a second Hamiltonian solver.

## Native YAeHMOP parameters in the friendly input

The two available cases are:

```ini
xyz
water.xyz
parameters = default
calculations = energy, charges, dos
```

or

```ini
xyz
water.xyz
parameters = water_parameters.txt
calculations = energy, charges, dos
```

The parameter file may use either of the following YAeHMOP conventions.
For compatibility with older jobs, the GUI also offers an **inline**
parameter source, where you can paste a classic `C 6 4 ...` record
without choosing a file.

### Classic YAeHMOP PARAMETERS block

The `PARAMETERS` and `END` markers are optional. One line may
contain s, p, d and f triplets (n, zeta, onsite energy):

```text
; Comments are allowed
PARAMETERS
H 1 1 1 1.3000 -13.6000
O 8 6 2 2.2750 -32.3000 2 2.2750 -14.8000
END
```

The columns are:

```text
Symbol AtomicNumber Valence  ns zeta_s Hss  np zeta_p Hpp  [nd zeta_d Hdd ...]
```

This is the same classic format used by the original YAeHMOP input
section. Classic records use one STO exponent per angular-momentum shell.

### Native YAeHMOP eht_parms.dat table (supports double-zeta)

```text
; Label AtNo Nval Nzeta Nquant Shell Hii zeta1 zeta2 coeff1 coeff2
FE 26 8 1 4 s  -9.100 1.9000 0.0000 1.0000 0.0000
FE 26 8 1 4 p  -5.320 1.9000 0.0000 1.0000 0.0000
FE 26 8 2 3 d -12.600 5.3500 2.0000 0.5505 0.6260
```

The complete bundled native file is available at:

`pyaehmop/eht_parms.dat`

The parser recognizes shells `s,p,d,f`, preserves double-zeta
exponents and coefficients, ignores comments beginning with
`#`, `;` or `!`, and permits entries for many species; only
species found in your structure are applied. Repeated native shell rows
use the final row, matching the source implementation's precedent.

If `parameters = default`, the bundled EHT table is used for
XYZ/CIF/ASE files. When the geometry is an original YAeHMOP input,
its existing native parameter assignments are preserved. Thus
`default` never unexpectedly overrides parameters contained in
that original input.

If your structure contains an element without usable built-in orbitals,
supplying an external parameter file can enable that species for XYZ,
CIF and ASE geometries, provided enough orbitals are defined.

**Changing the onsite energies or radial exponents changes the
calculated results.** The parameter selection feature is designed to
make this choice explicit, not to silently apply or optimize a fit.

## Integration with spin and SOC

The same file choice works in a `.spin` job:

```ini
ase
Fe_wire.extxyz
parameters = ../../eht_parms.dat
spin = yes
soc = yes
stoner_fe_d = 0.88
soc_fe_d = 0.06
initial_moment = 2.5
kmesh = 24
band_path = G:0,0,0; X:0.5,0,0
workers = 2
```

The EHT parameter file specifies the nonmagnetic orbital basis and
onsite energies. The `stoner_*` and `soc_*` values are **separate
physical parameters**, entered under **Spin / spin-orbit settings**.
They are not derived automatically from an EHT table.

For SOC, bands are colored by projected spin `<sigma_z>`; they are
not rigorously separable into pure spin-up/down states. The underlying
SCF remains collinear along z and is not a fully noncollinear solver.

## Output files and figures

The GUI produces the same scientific results as the CLI. The results
folder may contain `summary.json`, `mulliken.csv`, `bands.csv`,
`dos.csv`, `bands.png`, `dos.png` and `charges.png`.

For spin/SOC the engine produces `spin_soc_result.json`,
`bands_spin_soc.png`, the spin-resolved DOS and SCF history CSVs,
according to the selected model.

The **Open results JSON** button can display saved numerical results
without re-running them. Atom colors and simple bonds in the 3D
structure preview are visualization aids; they are **not** a bond-order
calculation. Electron-level COOP and COHP are separate numerical tasks.

## Verify the implementation

From `PyYAeHMOP`:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

python3 -m unittest discover -s pyaehmop/tests -p 'test_*.py' -v
python3 -m unittest discover -s pyaehmop/SPIN_TEST/tests -p 'test_*.py' -v
python3 compare_numerics.py
```

The parameter/GUI suite checks classic YAeHMOP records, full native
tables, molecular and spin-job parameter overrides, new elements,
entry generation/round trips, actual GUI-runner calculations, and
Matplotlib plotting. The existing core and spin/SOC numerical regression
tests remain part of the test process.

For scientific accuracy, always check k mesh, numerical quadrature,
overlap cutoff and, where relevant, the suitability of semiempirical
EHT/Stoner/SOC parameters for your material.
