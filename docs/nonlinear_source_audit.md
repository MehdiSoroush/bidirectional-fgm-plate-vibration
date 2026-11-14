# Nonlinear analytical coefficients and computational reconstruction

## Scope and chronology

The analytical investigation of bidirectionally graded plate vibration dates to **2018–2020**. The modern Python implementation was developed in **2026**. Eight Mathematica 9 notebooks record the analytical and numerical development. The reconstructed computational quantities are specified below, while complete independent symbolic verification remains a separate stage.

## Notebook inventory

| Notebook | Analytical or computational content | Verification status |
|---|---|---|
| relations.nb | Volume fractions and mixture relations | Examined |
| finah_eq.nb | Nonlinear transverse plate operator and Galerkin projection | Inspected; symbolic check pending |
| final_coeff.nb | Linear, cubic and inertial projected terms | Inspected; symbolic check pending |
| initil_def_test.nb | Initial-deflection parameters and response plots | Examined |
| main_code_naderifar.nb | Material, geometry, analytical coefficient functions P/Q | **Active P and Q reconstructed** |
| time_response_numerical.nb | NDSolve cubic time responses | **Reduced ODE reconstructed** |
| validate_linear.nb | Alternative validation material and parameters | Kept separate |
| validate_nonlinear.nb | Alternative nonlinear validation parameters | Kept separate; independent check pending |

## Source equation

The time-domain notebook solves

\[
\ddot q+Pq+Qq^3=F\sin(\Omega t),\quad
q(0)=0.02,\quad \dot q(0)=0.
\]

A recorded comparison uses F = 5 × 10⁶ and Ω = 1500 rad/s, over 0–0.03 s. There is **no damping term**. This F is inserted directly into the mass-normalized acceleration equation. Its conversion from an applied pressure has **not yet been verified**; a separate load-projection function appears in the analytical notebook.

## Model mapping

The analytical definitions use a = δh and b = ra, where δ = a/h and r = b/a. They average the ceramic volume fraction over the plate area:

\[
\overline V_c=1/[(1+\alpha)(1+\beta)],
\quad\overline E=E_m+(E_c-E_m)\overline V_c,
\quad\overline\rho=\rho_m+(\rho_c-\rho_m)\overline V_c.
\]

Area-integrated constitutive coefficients, bending stiffness, inertia and Airy coupling terms enter P and Q. Their numerical port is in the Python module **src/fgm_plate/analytical_single_mode.py**. The area-averaged **single-mode** reduction is not identical to the spatially resolved **10×10 Ritz** model. At gradation (1,1), their respective linear estimates are approximately **417.52 Hz** and **402.57 Hz**. Neither value validates the other model's assumptions.

## Reproduced coefficients

Six source coefficient pairs are in **data/reference/nonlinear_source_coefficients.csv**. The notebook reports approximately six significant digits; Python reproduces them within relative tolerance 3 × 10⁻⁶.

| Case | a/h | b/a | α | β | P (s⁻²) | Q (m⁻² s⁻²) |
|---|---:|---:|---:|---:|---:|---:|
| s00 | 20 | 1 | 0 | 0 | 1.38808e7 | 1.14191e10 |
| s10 | 20 | 1 | 1 | 0 | 9.60979e6 | 7.90555e9 |
| s11 | 20 | 1 | 1 | 1 | 6.88208e6 | 5.66158e9 |
| s201 | 20 | 1 | 0.5 | 0.5 | 9.04429e6 | 7.44034e9 |
| s151 | 15 | 1 | 0.5 | 0.5 | 2.84936e7 | 2.34405e10 |
| s2015 | 20 | 1.5 | 0.5 | 0.5 | 4.72292e6 | 4.26359e9 |

The coefficient agreement establishes numerical consistency with the recorded analytical formulation. The underlying continuum mechanics requires independent derivation checks.

## Structural feature of the averaged reduction

For fixed geometry, the reconstructed source formulation yields a constant Q/P across the material-gradation exponents. For the baseline square plate, **Q/P = 822.65625 m⁻²**. The free Duffing frequency normalized by its own linear frequency therefore has the same amplitude dependence for (0,0), (1,0) and (1,1), even though their absolute frequencies differ. This conclusion applies **only to the audited source reduction**, not to arbitrary spatially varying nonlinear plates.

For the unforced reduced ODE, independent energy-integral verification uses

\[
T(A)=4\int_0^{\pi/2}
\left[P+\frac12QA^2(1+\sin^2\theta)\right]^{-1/2}d\theta,\quad
\omega(A)=2\pi/T(A).
\]

The exact reduced-ODE period is checked against a numerically integrated quarter-cycle; the free-motion energy is also checked. These checks evaluate the reduced oscillator; nonlinear finite-element comparison is a separate verification task.

## Reproduction

Install dependencies and run tests, then run **python scripts/reproduce_nonlinear_transients.py**. That script generates reconstructed time-response and free-frequency PNGs in the results/figures directory.

## Open verification questions

1. Re-execute the 2018–2020 Mathematica notebooks and compare the time-history arrays numerically.
2. Independently verify the von Kármán/Airy reduction and physical forcing projection.
3. Develop spatially resolved nonlinear coefficients and membrane-consistent boundary conditions.
4. Assess nonlinear modal truncation and amplitude validity.
5. Build geometrically nonlinear Abaqus benchmarks and compare them independently.

The numerical verification covers reduced-order coefficients and trajectories. Geometrically nonlinear finite-element comparison remains a planned research stage.
