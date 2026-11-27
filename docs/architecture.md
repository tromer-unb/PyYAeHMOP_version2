# Architecture and parallel execution

[← Overview](../README.md) · [Theory](theory.md)

## Design goals

1. Reimplement the supported **EHT numerical method in Python**, without invoking the original YAeHMOP C executable.
2. Separate **physical engine** from **input formats** and **user interfaces**.
3. Allow CPU parallelism where tasks are independent.
4. Keep semiempirical spin/SOC research extensions isolated from standard EHT.
5. Make regression tests and example inputs inspectable and reproducible.

## Main package modules

| Module | Responsibility |
|---|---|
| `basis.py` | Atomic parameters, normalized Slater orbitals and real harmonics |
| `overlap.py` | Numerical two-center overlap integrals |
| `engine.py` | EHT Hamiltonian construction, occupations, molecular and periodic eigenproblems |
| `input.py` | Supported original YAeHMOP syntax |
| `structure_io.py` | XYZ, CIF and ASE structures |
| `friendly.py` | Human-readable `.job` input and task selection |
| `parameter_files.py` | Native YAeHMOP and classic external EHT parameters |
| `spectra.py` | DOS, k-point properties and COOP/COHP |
| `fragments.py` | Molecular fragment orbital analysis |
| `plots.py` and `export.py` | PNG and CSV presentation |
| `parameter_editor_model.py` | Typed shell table and native export |
| `SPIN_TEST/solver.py` | Magnetic SCF, spinor Hamiltonian and SOC |
| `SPIN_TEST/operators.py` | s/p/d/f angular momentum and onsite L·S |
| `SPIN_TEST/config.py` | Spin/SOC friendly input |
| `SPIN_TEST/plotting.py` | Spin- and SOC-colored band diagrams |

## Data flow

```text
YAeHMOP input ─┐
XYZ / CIF ─────┼──▶ structure + element parameters ─▶ STO overlaps
ASE / POSCAR ──┘                    │                    │
        .job options       ────────┘                    ▼
                                            EHT Hamiltonian H, overlap S
                                                        │
                                             generalized eigensolver
                                                        │
                              ┌─────────────────────────┴─────────────┐
                              ▼                                       ▼
                      molecular occupations                    periodic k bands
                              │                                       │
                       charges / FMO                          DOS / COOP / COHP
                              └─────────────────────────┬─────────────┘
                                                        ▼
                                             JSON / CSV / PNG

Optional research extension:
 EHT H,S + Stoner I / SOC ξ ─▶ spinor H,S ─▶ iterative magnetic SCF
                                          └▶ spin-projected bands / DOS
```

## Parallel execution

`concurrent.futures.ProcessPoolExecutor` partitions appropriate independent tasks among Python worker processes. Examples include:

- Pairwise, two-center atomic overlap calculations.
- Distinct reciprocal-space k-point generalized eigenproblems.

Specify `workers = 4`. Independent per-k calculations are suitable for process-level concurrency, but memory footprint, Python interprocess communication, BLAS threading and small problem sizes affect performance.

A conservative option is:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m pyaehmop pyaehmop/examples_friendly/diamond_cif.job
```

Use serial `workers=1` as a reference and compare results before benchmarking. Parallelism should not change the Hamiltonian formula.

## Scope and contribution

The EHT package and experimental magnetic module are independently testable. Proposed changes to EHT physical formulas, weighting conventions or SOC models should arrive in separate feature branches with tests and scientific references; documentation and interface updates should not silently change model behavior.

See [CONTRIBUTING.md](../CONTRIBUTING.md).
