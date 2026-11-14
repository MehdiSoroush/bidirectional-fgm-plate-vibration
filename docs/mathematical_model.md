# Mathematical model: linear plate benchmark

The linear free-vibration benchmark is formulated using Kirchhoff–Love plate theory and Rayleigh–Ritz discretization. Nonlinear reductions are documented in the [single-mode oscillator](nonlinear_model.md) and [variational plate extension](nonlinear_variational_extension.md).

## Domain and gradation

The square plate occupies $0\le x\le a$, $0\le y\le b$ with $a=b=0.8\,\mathrm m$ and a baseline thickness of $h=0.04\,\mathrm m$.

$$
V_c(x,y)=\left(\frac{x}{a}\right)^\alpha\left(\frac{y}{b}\right)^\beta,\qquad \alpha,\beta\ge0.
$$

$$
E(x,y)=E_m+(E_c-E_m)V_c(x,y),\qquad
\rho(x,y)=\rho_m+(\rho_c-\rho_m)V_c(x,y).
$$

The properties are $E_c=380$ GPa, $E_m=70$ GPa, $\rho_c=3800$ kg/m3, $\rho_m=2700$ kg/m3, and $\nu=0.30$. The $(0,0)$ case is homogeneous ceramic. These constitutive assumptions are idealized, not experimentally calibrated.

## Kirchhoff–Love Ritz approximation

The transverse field is approximated by a double-sine basis:

$$
w(x,y,t)=\sum_{m=1}^{N}\sum_{n=1}^{N}q_{mn}(t)
\sin(m\pi x/a)\sin(n\pi y/b).
$$

With $D(x,y)=E(x,y)h^3/[12(1-\nu^2)]$, the mass and stiffness matrices are

$$
M_{ij}=\iint_\Omega \rho h\phi_i\phi_j\,d\Omega,
$$

$$
K_{ij}=\iint_\Omega D\left[\phi_{i,xx}\phi_{j,xx}
+\phi_{i,yy}\phi_{j,yy}
+\nu(\phi_{i,xx}\phi_{j,yy}+\phi_{i,yy}\phi_{j,xx})
+2(1-\nu)\phi_{i,xy}\phi_{j,xy}\right]d\Omega.
$$

Tensor-product Gauss–Legendre quadrature evaluates the matrices, followed by the symmetric generalized eigenproblem

$$
K\mathbf q=\omega^2M\mathbf q,\qquad f=\omega/(2\pi).
$$

The homogeneous simply supported reference is

$$
f_{11}=\frac{\pi}{2}\left(\frac{1}{a^2}+\frac{1}{b^2}\right)
\sqrt{\frac{D}{\rho_c h}}.
$$

## Abaqus implementation and limitations

S4R shell elements represent the plate; element-centroid properties approximate the continuous gradation. Transverse displacement U3 is set to zero along all edges, and minimal in-plane constraints eliminate rigid-body motions. The linear modal restraint uses minimal in-plane pins; the variational nonlinear model instead adopts immovable in-plane edges.

The first *transverse-bending* mode is identified by examining its U3 field, not by assuming Abaqus Mode 1 is always the bending mode. A model-form difference from the classical Ritz solution is expected because the two kinematic formulations are not identical.

## Modal assurance criterion

For corresponding nodal displacement samples $u$ and $v$,

$$
\operatorname{MAC}(u,v)=\frac{(u^Tv)^2}{(u^Tu)(v^Tv)}.
$$

For best-fit scale $c=(u^Tv)/(v^Tv)$, the relative L2 residual is $\|u-cv\|_2/\|u\|_2$. A MAC near unity indicates strong correlation of the sampled linear displacement fields. Nonlinear stiffness and response require independent checks.
