# Third-party notices and project independence

PyYAeHMOP is an **independent, unofficial** Python implementation that follows Extended Hückel Theory and is inspired by the capabilities of YAeHMOP. It is not the official YAeHMOP distribution, is not maintained by its original authors and does not invoke the original YAeHMOP C executable as a runtime dependency.

## Original YAeHMOP project and scientific data

**Original project:** YAeHMOP (Yet Another Extended Hückel Molecular Orbital Package)  
**Original author / project attribution:** Greg Landrum and original YAeHMOP contributors  
**Source:** https://github.com/greglandrum/yaehmop

The repository includes a copy of the YAeHMOP EHT parameter table `pyaehmop/eht_parms.dat` and certain compatibility/testing examples or reference structures associated with the original code or its documented conventions.

The original project carries a **BSD-style, two-condition redistribution license**. Its license is preserved in full at [`pyaehmop/LICENSE.txt`](pyaehmop/LICENSE.txt), including copyright statements, conditions and warranty disclaimer. Please preserve those notices when redistributing this material.

## New Python project

The new Python source, UI code and project-specific documentation are offered under the repository [MIT license](LICENSE), subject to the above retained third-party notices.

An MIT license on new files does **not** alter the licensing requirements of redistributed YAeHMOP materials or any additional content obtained from other authors. Users are responsible for reviewing the applicable notices when bundling data into derivative works.

## Scientific context

Scientific ideas, methods and experimental parameter choices are described and referenced in [the theory documentation](docs/theory.md) and [the example validation notes](docs/validation.md). Those citations distinguish methods from the new implementation and avoid suggesting endorsement by the cited researchers.

The Fe-wire examples are independent computations based on a geometry discussed in the literature, and their figures are **not reproduced figures from a copyrighted publication**.
