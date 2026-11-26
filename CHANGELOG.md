# Changelog

This project uses a [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) style and semantic source-version labels. Version numbers do not imply a PyPI, conda-forge or GitHub Release publication.

## [0.4.0] — Source snapshot

### Added

- Human-friendly input for original YAeHMOP, XYZ, CIF and ASE structures.
- Visual desktop GUI for Linux, with 3D structure, band, DOS, charge and SCF views.
- Native YAeHMOP parameter-table import and interactive orbital-shell editing.
- Dedicated input modes for EHT, self-consistent Stoner spin, onsite SOC and their combined spinor model.
- External parameter selection through `parameters = default` or `parameters = <filename>`.
- Local validation of 74 automated tests, including spin/SOC and GUI/parameter regression tests.
- Repository documentation, screenshots/figures, installation paths, attribution and contributor guidance.

### Important limitations

- Spin+SOC is collinear-z Stoner SCF with onsite SOC, not full noncollinear spin-density SCF.
- Extended Hückel electronic energies are not full DFT total energies.
- Existing native YAeHMOP inputs are supported only within a documented subset.

## [0.3.1] — English interface snapshot

- English naming, user-facing diagnostics, friendly parameters, documentation and plots.
- Compared 10,314 sampled numerical values with the preceding Python version, with no observed differences.

## Earlier development snapshots

- Independent EHT engine in Python and CPU-parallel overlap/k-point tasks.
- Molecular levels, Mulliken charges, crystal bands, DOS, COOP/COHP and molecular FMO.
- Original-format YAeHMOP compatibility and tested examples.

Prior snapshots are described here historically; no earlier public Git tag is implied.
