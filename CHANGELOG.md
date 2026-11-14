# Changelog

## Unreleased — nonlinear reconstruction (2026-11-14)

- Reconstructed analytical single-mode P/Q expressions from the 2018–2020 Mathematica source.
- Added six source-coefficient regression cases and an undamped cubic time-response solver.
- Implemented conservative energy and exact-period verification.
- Added an initial eight-notebook source audit and reproducible figure generation.
- Derived an independent spatially graded variational von Kármán single-mode model with immovable in-plane edges and interior membrane relaxation.
- Added reproducible 10×10 Ritz / 8×8 membrane benchmark data and regression tests.
- Documented source-equation rotary-inertia and forcing-normalization consistency checks.
- Analytical derivation checks, Mathematica re-execution and matched-restraint nonlinear finite-element comparison remain planned research tasks.

## 0.1.0 — 2026-11-13

Reproducible linear vibration benchmark associated with research conducted in 2018–2020.

- Python Kirchhoff–Love double-sine Rayleigh–Ritz solver.
- Material gradation model and regression tests.
- Recorded Abaqus frequencies, thickness and mesh comparisons.
- Abaqus 2021 / Python 2.7-compatible input generator and U3 exporter.
- Nodal MAC implementation and recorded FEM–Ritz mode-shape correlation benchmarks.
- GitHub Actions checks and MIT licensing.
- Documentation of the model assumptions, numerical verification scope, and nonlinear development priorities.

Supplementary data and documentation were added in commit `5bbbe3d`: recorded Abaqus U3 datasets, scientific PNG figures, and the independent technical research report. The versioned datasets are used in modal regression tests.
