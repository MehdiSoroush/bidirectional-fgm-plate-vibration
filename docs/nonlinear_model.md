# Nonlinear time-domain computational model

The analytical single-mode model is governed by

\[
\ddot q+Pq+Qq^3=F\sin(\Omega t).
\]

The state-space form is \(\dot y=[\dot q,-Pq-Qq^3+F\sin(\Omega t)]\). The implementation uses an adaptive DOP853 integrator. The numerical example uses q(0) = 0.02 m, dq/dt(0) = 0, F = 5 × 10⁶ in the reduced ODE and Ω = 1500 rad/s. The coefficient F enters the mass-normalized equation; its conversion to distributed pressure requires the modal forcing projection.

For zero forcing, the conserved energy is

\[
H=\frac12\dot q^2+\frac12 Pq^2+\frac14 Qq^4.
\]

An independent conservative-period quadrature gives the exact free-oscillator frequency. The first-order weak-nonlinearity approximation is

\[
\omega(A)\simeq\sqrt P\left(1+\frac{3QA^2}{8P}\right).
\]

The weak-nonlinearity approximation is applicable when the cubic correction is small. The period and energy checks verify the reduced oscillator implementation. Its coefficients employ area-averaged properties, whereas the spatial Ritz model retains coordinate-dependent properties.

The [analytical coefficient study](nonlinear_source_audit.md) gives the coefficient regression cases, notebook correspondence, and further verification tasks.
