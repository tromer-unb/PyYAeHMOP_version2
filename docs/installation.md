# Installation

[← Overview](../README.md) · [Inputs](input-formats.md)

PyYAeHMOP Version 2 (0.5.0.dev0) uses **Python 3.10 or newer**, NumPy and SciPy. Structure import from CIF/POSCAR uses ASE; plotting uses Matplotlib.

> **Distribution status:** This repository is the authoritative source. The instructions below use **pip installation from GitHub** or **Conda to install dependencies and then the local project**. We do **not** claim that a separate PyPI or conda-forge release exists.

## Quick install (Linux, macOS or Windows command line)

On systems with a supported Python installation and Git:

```bash
git clone https://github.com/tromer-unb/PyYAeHMOP_version2.git
cd PyYAeHMOP_version2
python3 -m pip install -e .
```

The basic installation provides EHT calculations and standard YAeHMOP/XYZ inputs. For optional ASE structures and Matplotlib plotting:

```bash
python3 -m pip install -e ".[all]"
```

Use `python` instead of `python3` if that is your environment's Python command.

**Install without cloning** (core engine):

```bash
python3 -m pip install "git+https://github.com/tromer-unb/PyYAeHMOP_version2.git"
```

For a deterministic code revision, replace the Git URL with `git+https://github.com/tromer-unb/PyYAeHMOP_version2.git@<commit-sha>`; record the SHA in your publication.

## Isolated pip environment (recommended)

```bash
git clone https://github.com/tromer-unb/PyYAeHMOP_version2.git
cd PyYAeHMOP_version2
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[all]"
```

## Conda installation (conda-forge)

The root [`environment.yml`](../environment.yml) creates a reproducible dependency environment and installs PyYAeHMOP from the current checkout with pip:

```bash
git clone https://github.com/tromer-unb/PyYAeHMOP_version2.git
cd PyYAeHMOP_version2
conda env create -f environment.yml
conda activate pyaehmop
```

To update the environment after changing the source:

```bash
conda env update -f environment.yml --prune
```

Alternatively, create an environment manually:

```bash
conda create -n pyaehmop -c conda-forge python=3.11 numpy scipy ase matplotlib pip
conda activate pyaehmop
python -m pip install -e .
```

This is **not** `conda install pyaehmop`: there is no claim of a published conda-forge recipe.

## Check the installed package

```bash
python -c "import pyaehmop; print(pyaehmop.__version__)"
python -m pyaehmop --help
python -m pyaehmop.SPIN_TEST --help
```

### Command-line smoke examples

```bash
python -m pyaehmop pyaehmop/inputs/XeF4
python -m pyaehmop pyaehmop/examples_friendly/water_full.job
python -m pyaehmop pyaehmop/examples_friendly/diamond_cif.job
python -m pyaehmop.SPIN_TEST --demo
```

Run all example commands from the **repository root**, so the relative paths above resolve correctly. In an installed package outside a source checkout, copy the bundled sample inputs to your own working directory.

## Scientific numerical performance

For multi-process calculations, avoid *nested* BLAS threading consuming more CPU cores than expected:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python -m pyaehmop pyaehmop/examples_friendly/diamond_cif.job
```

Then set `workers = 2`, `workers = 4`, etc. in the job. Speedup depends on the structure, workload, number of k points, memory and CPU architecture; the software does not promise linear speedup.

## Troubleshooting

| Symptom | Suggested action |
|---|---|
| `ModuleNotFoundError: numpy/scipy` | Install from `pip install -e .` in the selected environment |
| `ImportError: ase` | Install `pip install -e ".[formats]"` or Conda ASE |
| `ImportError: matplotlib` | Install `pip install -e ".[plots]"` |
| Missing atomic parameters | Use `parameters = custom.txt` |
| Non-positive overlap matrix | Check basis, geometry, integration orders and overlap cutoff |
| Slow periodic calculations | Reduce diagnostic mesh, then perform a deliberate convergence study |

## Validation and releases

You can run the local test suite:

```bash
python -m unittest discover -s pyaehmop/tests -p 'test_*.py' -v
python -m unittest discover -s pyaehmop/SPIN_TEST/tests -p 'test_*.py' -v
```

Version 2 has no GitHub Actions workflows. Run validation only in your own environment. No PyPI/conda-forge package or GitHub Release asset is implied by the source installation instructions.

## Transport examples

Clone the source checkout (recommended, since the example geometry/data folders are not all installed as package data). Start with [the tutorial index](../examples/tutorials/README.md) and [PDF manual](PyYAeHMOP_v2_tutorial.pdf). A workstation is sufficient for the quick cases; PBS is entirely optional. On Windows use `.venv\Scripts\Activate.ps1` in PowerShell instead of the Unix activation command.
