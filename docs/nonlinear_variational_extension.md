# Spatially graded von Kármán plate: variational single-mode extension

**Computational study: November 2026.** The spatially resolved model extends the 2018–2020 nonlinear analytical investigation. It differs from the area-averaged model in constitutive treatment, mode-shape representation, and in-plane restraints.

## Physical system and restraints

The plate has \(\Omega=[0,a]\times[0,b]\), thickness \(h\), \(V_c(x,y)=(x/a)^\alpha(y/b)^\beta\), spatial Young's modulus \(E=E_m+(E_c-E_m)V_c\), and density \(\rho=\rho_m+(\rho_c-\rho_m)V_c\). The plane-stress law is locally isotropic with constant Poisson ratio.

The transverse Ritz displacement is \(w(x,y,t)=q(t)\phi(x,y)\), with \(\phi(a/2,b/2)=1\). The essential transverse restraint is **w=0 at all four edges**; natural bending conditions are satisfied in the weak Ritz sense.

For this nonlinear extension, **u=v=0 on every plate edge**. The in-plane support is *immovable*, whereas the linear Abaqus benchmark uses minimal pins to suppress rigid-body motion. In-plane displacement is allowed in the interior and is condensed quasi-statically.

Assumptions: Kirchhoff–Love bending, von Kármán membrane strains, small in-plane strain, no damping, no initial membrane stress, no imperfections, no in-plane inertia, no rotary inertia, and one prescribed transverse mode shape. These assumptions exclude a full finite-strain or multi-transverse-mode nonlinear FE interpretation.

## Variational derivation

The mid-surface engineering strain vector is

\[
\boldsymbol\epsilon_m=
[u_{,x}+w_{,x}^2/2,\ 
 v_{,y}+w_{,y}^2/2,\ 
 u_{,y}+v_{,x}+w_{,x}w_{,y}]^T.
\]

The spatial plane-stress membrane stiffness is

\[
\mathbf A(x,y)=\frac{E(x,y)h}{1-\nu^2}
\begin{bmatrix}1&\nu&0\\ \nu&1&0\\ 0&0&(1-\nu)/2\end{bmatrix}.
\]

Let \(u,v=q^2\,\mathbf N(x,y)\mathbf z\), using independent double-sine bases for each in-plane component; this satisfies u=v=0 on the boundary. Then

\[
\boldsymbol\epsilon_m=q^2(\mathbf B\mathbf z+\mathbf g),\qquad
\mathbf g=[\phi_{,x}^2/2,\phi_{,y}^2/2,\phi_{,x}\phi_{,y}]^T.
\]

Define

\[
\mathbf K_m=\int_\Omega\mathbf B^T\mathbf A\mathbf B\,d\Omega,\qquad
\mathbf f_m=\int_\Omega\mathbf B^T\mathbf A\mathbf g\,d\Omega,\qquad
C_0=\int_\Omega\mathbf g^T\mathbf A\mathbf g\,d\Omega.
\]

Minimizing the membrane potential yields

\[
\mathbf z_*=-\mathbf K_m^{-1}\mathbf f_m,\qquad
C=C_0-\mathbf f_m^T\mathbf K_m^{-1}\mathbf f_m,\qquad
U_m=\tfrac12 Cq^4.
\]

The bending stiffness and transverse mass for the chosen Ritz shape are

\[
K_b=\int_\Omega D(x,y)[\phi_{,xx}^2+\phi_{,yy}^2+
2\nu\phi_{,xx}\phi_{,yy}+2(1-\nu)\phi_{,xy}^2]\,d\Omega,
\]

\[
M=\int_\Omega \rho(x,y)h\phi^2\,d\Omega,\quad
D(x,y)=E(x,y)h^3/[12(1-\nu^2)].
\]

Consequently,

\[
M\ddot q+K_bq+2Cq^3=0,\qquad P=K_b/M,\quad Q=2C/M.
\]

For a *spatially uniform physical pressure* \(p(t)\), the new model's forcing projection is

\[
\ddot q+Pq+Qq^3=p(t)\,
\frac{\int_\Omega\phi\,d\Omega}{M}.
\]

The uniform-pressure projection belongs to the spatial formulation. The numerical forcing coefficient in the 2018–2020 single-mode model uses a separate normalization whose pressure mapping remains to be established.

## Reproducible benchmark

Default numerical setup: a=b=0.8 m, h=0.04 m, **10×10 transverse Ritz modes**, **8×8 in-plane modes per component**, and **48×48 Gauss quadrature**. The table compares the spatially resolved and area-averaged reductions under their respective assumptions.

| Gradation | Linear f (Hz) | Spatial Q (10⁹ m⁻²s⁻²) | Spatial Q/P (m⁻²) | Averaged Q/P (m⁻²) | Relaxation |
|---|---:|---:|---:|---:|---:|
| (0,0) | 594.18074 | 11.52876 | 827.153 | 822.656 | 29.42% |
| (1,0) | 477.31102 | 7.64476 | 849.965 | 822.656 | 28.90% |
| (1,1) | 402.57036 | 5.45674 | 852.884 | 822.656 | 28.78% |
| (1,5) | 334.65096 | 3.50175 | 792.030 | 822.656 | 32.48% |
| (5,1) | 334.65096 | 3.50175 | 792.030 | 822.656 | 32.48% |
| (5,5) | 316.53602 | 3.09578 | 782.644 | 822.656 | 29.96% |

The area-averaged model has gradation-invariant Q/P for fixed geometry, while the spatial model yields gradation-dependent values. Since the formulations also differ in mode enrichment and boundary conditions, the effect of each individual assumption cannot be isolated from this comparison alone.

## Numerical checks

The test suite verifies homogeneous plate frequencies, coordinate-swap symmetry, positivity of membrane energy, monotonic non-increase of minimized membrane energy under enrichment of the in-plane trial space, the linear Ritz eigenfrequency recovered from the enriched transverse mode, thickness scaling \(P\propto h^2,\ Q\propto h^0\) for the declared restraints, homogeneous uniform-pressure gain \(16/(\pi^2\rho h)\), and conservative reduced-ODE energy. These checks establish numerical consistency of the reduced-order model. Matched-restraint nonlinear finite-element comparison is a separate stage.

For gradation (1,1), the 10×10 transverse mode gives Q=7.35968×10⁹ at one in-plane term per direction, decreasing to Q=5.45674×10⁹ with eight terms per direction. Higher-order in-plane refinement and alternate boundary conditions remain relevant.

To regenerate the numerical CSV and figures, execute:

~~~bash
python scripts/sweep_spatial_nonlinear.py
pytest -q
~~~

Source: src/fgm_plate/nonlinear_membrane.py. Data: data/reference/spatial_nonlinear_immovable.csv.

## Outstanding verification

1. Perform a matched-restraint geometrically nonlinear Abaqus comparison, with consistent displacement amplitude normalization.
2. Extend in-plane and transverse convergence studies across graded cases and amplitudes.
3. Investigate movable and mixed in-plane supports; stiffness Q is boundary-sensitive.
4. Include additional transverse modes and nonlinear modal interaction as needed.
5. Independently verify the 2018–2020 Airy-stress derivations and forcing conventions.

The current results comprise analytical derivations, reduced-order simulations, and numerical regression checks. Validation against geometrically nonlinear Abaqus solutions remains outstanding.

## Background

M. Amabili, “Nonlinear vibrations of rectangular plates with different boundary conditions: theory and experiments,” *Computers & Structures* 82 (2004), 2587–2605. DOI: 10.1016/j.compstruc.2004.03.077. The present computational reduction and numerical results were independently developed.
