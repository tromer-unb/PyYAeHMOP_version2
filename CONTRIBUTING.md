# Contributing to PyYAeHMOP

Thank you for helping improve this independent research software.

[Overview](README.md) · [Theory](docs/theory.md) · [Tests](docs/validation.md) · [Security](SECURITY.md)

## Contribution principles

- Scientific correctness and explicit approximations take priority over feature count.
- Numerical changes must have regression tests and a short explanation of the physical model.
- Additions to EHT, Stoner exchange or SOC must cite their methodological basis and define their validity limits.
- New input keywords should be documented in [input formats](docs/input-formats.md), with a runnable fixture.
- Never silently introduce breaking input or JSON-output changes; describe them in the changelog.
- Do not commit private experimental structures, login credentials, temporary outputs or machine-specific directories.

## Branching model

| Branch | Purpose |
|---|---|
| `main` | Reviewed, stable publication snapshots |
| `develop` | Integration of proposed changes before a release |
| `feature/<topic>` | Scientific or user-interface enhancements |
| `fix/<topic>` | Reproducible bug fixes |
| `docs/<topic>` | Documentation-only changes |
| `release/<version>` | Release preparation when needed |

Proposed changes should normally start from `develop`. Open a pull request targeting `develop`, describe scope and validation, and integrate into `main` when ready to publish. Do not rewrite public branch history.

## Development setup

```bash
git clone https://github.com/tromer-unb/PyYAeHMOP.git
cd PyYAeHMOP
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[gui]"
```

## Running the local scientific tests

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python -m unittest discover -s pyaehmop/tests -p "test_*.py" -v
python -m unittest discover -s pyaehmop/SPIN_TEST/tests -p "test_*.py" -v
```

The Tk window test requires a desktop display; on a headless runner it may be skipped. The optional GitHub workflow runs **only when manually dispatched**. No GitHub Action is required or triggered automatically by a pull request or push.

## Pull request checklist

- [ ] Describe the motivation and what changed.
- [ ] Confirm no third-party or confidential file is being redistributed without permission.
- [ ] Include input cases and exact parameters to reproduce the behavior.
- [ ] Explain whether the **scientific equations** changed.
- [ ] Report local test status, including skipped tests.
- [ ] Update README, docs, changelog and parameter examples as applicable.
- [ ] Verify a user can install from a clean checkout.

## Scientific issue reports

A helpful issue includes: a minimal input, the parameter table, expected versus observed numbers, Python and package versions, relevant machine/OS details, and convergence controls (`kmesh`, `cutoff`, `quadrature`, `temperature_ev`, and so forth).

If results differ physically from DFT or experiment, first ensure both models, bases and electronic populations are comparable. An EHT-vs-DFT discrepancy alone does not prove an implementation bug.

## Attribution

Retain citations for Hoffmann's EHT methodology and the original YAeHMOP implementation. Preserve the original license notice for `eht_parms.dat`. This project is unofficial; avoid names or descriptions that imply endorsement from its original maintainers.
