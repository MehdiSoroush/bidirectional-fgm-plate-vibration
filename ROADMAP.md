# Research software roadmap

## v0.1 — Linear computational benchmark
- [x] Classical Python Ritz stiffness/mass assembly and generalized eigenvalue solver
- [x] Homogeneous exact frequency and plate-thickness verification
- [x] Five graded frequency benchmark cases
- [x] Abaqus 2021 S4R generator and U3 CSV export code
- [x] Recorded Abaqus mesh and gradation benchmark tables
- [x] Modal assurance criterion implementation and tests
- [x] Continuous integration for the self-contained Python benchmarks
- [x] Publish retained nodal Abaqus U3 datasets and report PDF as research assets
- [x] Enable the nodal-data MAC regression tests in CI

## v0.2 — Audited nonlinear analytical research
- [ ] Complete symbolic review of the 2018–2020 Mathematica derivations
- [x] Inventory the eight notebooks and reconstruct active P/Q coefficient expressions
- [ ] Re-execute and verify original nonlinear response calculations
- [x] Reconstruct six P/Q coefficient sets in Python
- [x] Integrate the source-form cubic ODE and verify energy and period
- [x] Derive an independent spatially graded single-mode von Kármán membrane-energy reduction
- [x] Implement and test immovable-edge in-plane relaxation, nonlinear stiffness and pressure gain
- [ ] Independently compare trajectories with rerun Mathematica source
- [ ] Perform matched-restraint nonlinear Abaqus validation of the new variational ROM

## v1.0 — Extended nonlinear reduced-order framework
- [ ] Evaluate movable and mixed in-plane supports in addition to immovable-edge conditions
- [ ] Evaluate amplitude-dependent frequencies and mode coupling
- [ ] Validate against nonlinear Abaqus solutions across mesh and amplitude

The roadmap lists computational development and validation priorities.
