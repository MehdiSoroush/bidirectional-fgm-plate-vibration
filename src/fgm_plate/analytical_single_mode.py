"""Coefficients from the 2018–2020 analytical single-mode plate model.

The P and Q coefficients follow the definitions in main_code_naderifar.nb.
The model employs area-averaged constitutive properties, in contrast to the
spatially resolved stiffness and inertia of the multi-term Ritz formulation.
"""
from dataclasses import dataclass
from math import pi


@dataclass(frozen=True)
class AnalyticalMaterial:
    e_ceramic: float = 380e9
    e_metal: float = 70e9
    rho_ceramic: float = 3800.0
    rho_metal: float = 2700.0
    nu: float = 0.30
    thickness: float = 0.04

    def __post_init__(self):
        if min(self.e_ceramic, self.e_metal, self.rho_ceramic,
               self.rho_metal, self.thickness) <= 0:
            raise ValueError("Moduli, densities and thickness must be positive")
        if not -1.0 < self.nu < 0.5:
            raise ValueError("Invalid Poisson ratio for stable isotropic material")


@dataclass(frozen=True)
class AnalyticalCoefficients:
    p: float  # linear coefficient, s^-2
    q: float  # cubic coefficient, m^-2 s^-2


def source_coefficients(alpha: float, beta: float, *, slenderness: float = 20.0,
                        aspect_ratio: float = 1.0, m: int = 1, n: int = 1,
                        material: AnalyticalMaterial = AnalyticalMaterial()
                        ) -> AnalyticalCoefficients:
    """Evaluate the analytical single-mode P and Q coefficients.

    slenderness = a/h; aspect_ratio = b/a. The formulation uses integrated,
    area-averaged constitutive and inertia coefficients.
    """
    if alpha < 0 or beta < 0:
        raise ValueError("Gradation exponents must be nonnegative")
    if slenderness <= 0 or aspect_ratio <= 0:
        raise ValueError("Plate dimensions must be positive")
    if any(isinstance(k, bool) or int(k) != k or k < 1 for k in (m, n)):
        raise ValueError("m and n must be positive integers")

    h, nu = material.thickness, material.nu
    a = slenderness * h
    b = aspect_ratio * a
    area = a * b
    vc_mean = 1.0 / ((1.0 + alpha) * (1.0 + beta))
    e_mean = material.e_metal + (material.e_ceramic - material.e_metal) * vc_mean
    rho_mean = material.rho_metal + (material.rho_ceramic - material.rho_metal) * vc_mean
    c11 = area * e_mean / (1.0 - nu**2)
    c22 = c11
    c12 = nu * c11
    c66 = area * e_mean / (2.0 * (1.0 + nu))
    d11, d22 = c11 * h**3 / 12, c22 * h**3 / 12
    d12, d66 = c12 * h**3 / 12, c66 * h**3 / 12
    i0 = area * h * rho_mean
    s11 = c22 / (h * (c11*c22 - c12*c12))
    s22 = s11
    s12 = -c12 / (h * (c11*c22 - c12*c12))
    am, bn = m*pi/a, n*pi/b
    psi3 = (s11*bn**2 - s12*am**2) / (16*(s11*s22 - s12**2))
    psi4 = (s22*am**2 - s12*bn**2) / (16*(s11*s22 - s12**2))

    p = 12*pi**4 * (
        2*m**2*n**2*a**2*b**2*(d12 + 2*d66)
        + d22*n**4*a**4 + d11*b**4*m**4
    ) / (i0 * (12*a**2*b**2 + h**2*pi**2*(b**2*m**2 + n**2*a**2))
         * (a**2*b**2))

    q = (2*am**2*psi4 + 2*bn**2*psi3
         + (am**4/s11 + bn**4/s22)/16)
    q /= i0 * (1 + h**2*(am**2 + bn**2)/12)
    return AnalyticalCoefficients(float(p), float(q))
