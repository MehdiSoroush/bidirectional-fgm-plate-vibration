# Research provenance and verification records

## Timeline

- **2018–2020:** Research on free and nonlinear vibrations of bidirectionally functionally graded plates, using analytical mechanics, Galerkin reduction, perturbation methods, and Mathematica calculations.
- **2026:** Modernized Python Ritz implementation, verification of linear finite-element benchmarks, and preparation of a reproducible research software package.

The Python implementation reproduces selected nonlinear coefficients and reduced-order equations from the 2018–2020 notebooks. Complete Mathematica re-execution remains a separate verification task.

## Provenance categories

1. **Linear computational models:** Ritz eigenfrequencies, material gradation, symmetry/thickness checks, nodal MAC comparisons, and Abaqus input-deck generation.
2. **Recorded finite-element results:** Abaqus 2021 frequencies, mesh refinement, thickness series, graded cases and derived nodal MAC values. The numerical tables label their origin explicitly.
3. **Analytical nonlinear model:** The 2018–2020 Mathematica calculations and the Python reconstruction of six P/Q coefficient cases, cubic time-domain responses, and conservative period/energy checks.
4. **Spatial nonlinear extension:** A 2026 variational membrane-energy reduction using spatial material gradation, immovable in-plane edges, and numerical convergence and symmetry checks. Geometrically nonlinear finite-element verification remains planned.

The original solver outputs used for the recorded frequencies were reported to 0.01 Hz. Higher precision should not be inferred.

## Key technical limitations

- An idealized effective isotropic rule-of-mixtures material law is assumed.
- The modern Python model is classical Kirchhoff–Love bending; Abaqus uses S4R shell kinematics.
- The shell/classical frequency difference is measured, but the individual physical contributions were not isolated.
- GitHub Actions evaluates the Python benchmarks; independent Abaqus solver execution requires a licensed Abaqus installation.
- The repository includes reference tables, Abaqus U3 nodal datasets, model-generation scripts, scientific figures, and an independent technical research report.
- The nonlinear single-transverse-mode approximations require further amplitude-range and matched-boundary-condition finite-element investigation.
- The repository contains analytical and computational research and no experimental dataset.

## Licensing

Code and original repository documentation are released under the **MIT License**. Publications and cited third-party literature retain their own terms. Project author: **Mehdi Soroush**.
