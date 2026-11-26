## Summary

Describe the purpose and scope of this change.

## Type of change

- [ ] Scientific model / new method
- [ ] Bug fix
- [ ] Input format / CLI
- [ ] Desktop GUI
- [ ] Documentation or packaging
- [ ] Tests / reproducibility

## Scientific impact

Does this change any EHT, overlap, generalized eigenproblem, Stoner, SOC, DOS, COOP/COHP or FMO formula? If yes, explain the method, approximations and references. If no, state that the scientific method is unchanged.

## Reproducibility

Provide a minimal input file, parameter source, k mesh, cutoff and expected numerical outcome.

## Local checks

- [ ] Core tests passed (`python -m unittest discover -s pyaehmop/tests -p 'test_*.py'`)
- [ ] Spin tests passed if affected
- [ ] Example files are reproducible
- [ ] Package installs from a clean checkout
- [ ] No generated results, credentials or proprietary data committed
- [ ] Documentation and changelog updated if relevant

**Test environment / results:**

## Notes

Issues, possible numerical tolerances and known limitations:
