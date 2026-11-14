# Comparison of analytical notebook formulations

The 2018–2020 Mathematica notebooks contain alternative rotary-inertia and generalized-force expressions. The following comparisons identify their algebraic differences and numerical implications.

## Rotary inertia

The transverse operators in finah_eq.nb and final_coeff.nb explicitly include the term

\[
-\rho h^3 w_{,xxtt}/12,
\]

without an accompanying \(-\rho h^3w_{,yytt}/12\). The one-sine Galerkin projection of these particular operators consequently has an inertia factor

\[
1+h^2k_x^2/12.
\]

By contrast, the later P coefficient in main_code_naderifar.nb contains

\[
1+h^2(k_x^2+k_y^2)/12.
\]

For a=b=0.8 m, h=0.04 m, and mode (1,1), the x-only formulation gives P **0.2052% larger** than the symmetric-rotary formulation when all other parameters are held fixed. Its linear frequency is **0.1025% larger** (about **0.608 Hz**). The distinction is small in this geometry but exact algebraic equivalence cannot be claimed.

The separate spatially graded Kirchhoff–Love implementation intentionally **omits rotary inertia**; its frequencies are not expected to match either rotary-corrected coefficient exactly.

## External excitation

The notebook main_code_naderifar.nb contains a generalized load expression Ft, whereas time_response_numerical.nb inserts F directly into the mass-normalized ODE. The conversion between a physical pressure and the numerical F remains unresolved.

For the new unit-center-normalized homogeneous sine mode under uniformly distributed transverse pressure p(t), the generalized acceleration gain is

\[
G=\frac{\int\phi\,d\Omega}{\int\rho h\phi^2\,d\Omega}
=\frac{16}{\pi^2\rho h}.
\]

For the square fundamental mode, the literal source Ft/f expression reduces to \(4/[\pi^2\rho h(1+h^2(k_x^2+k_y^2)/12)]\). Assuming a uniformly distributed physical pressure, the expressions differ by a factor of four in the thin-plate limit. Their comparison requires matching forcing distributions and generalized-coordinate normalization. The numerical F=5×10⁶ is retained as a mass-normalized equation parameter.

## Generalized coordinate

The symbolic Galerkin notebooks use \(w=h\,r(t)\sin(k_xx)\sin(k_yy)\) with dimensionless r(t). The numerical P/Q model uses a generalized amplitude q with dimensions of length. Comparing cubic coefficients requires the explicit transformation \(q=hr\), including the h² scaling in the cubic term.

A complete symbolic and dimensional reconciliation of these expressions is listed among the remaining analytical verification tasks.
