"""Material and geometry definitions for the in-plane graded plate.

All dimensions and properties use SI units. The material law is an
assumed rule-of-mixtures model, not a measured microstructure model.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Plate:
    a: float = 0.8
    b: float = 0.8
    h: float = 0.04
    nu: float = 0.30
    E_c: float = 380e9
    E_m: float = 70e9
    rho_c: float = 3800.0
    rho_m: float = 2700.0

    def __post_init__(self):
        if min(self.a, self.b, self.h, self.E_c, self.E_m, self.rho_c, self.rho_m) <= 0:
            raise ValueError('Lengths and material properties must be positive')
        if not (-1.0 < self.nu < 0.5):
            raise ValueError('Invalid isotropic Poisson ratio')


def ceramic_fraction(x, y, alpha, beta, plate=Plate()):
    """V_c = (x/a)^alpha (y/b)^beta; alpha=beta=0 is all ceramic."""
    if alpha < 0 or beta < 0:
        raise ValueError('Exponents must be nonnegative')
    x, y = np.asarray(x), np.asarray(y)
    if np.any((x < -1e-12) | (x > plate.a+1e-12)) or np.any((y < -1e-12) | (y > plate.b+1e-12)):
        raise ValueError('Coordinates outside the plate')
    return (x / plate.a)**alpha * (y / plate.b)**beta


def properties(x, y, alpha, beta, plate=Plate()):
    vc = ceramic_fraction(x, y, alpha, beta, plate)
    E = plate.E_m + (plate.E_c - plate.E_m) * vc
    rho = plate.rho_m + (plate.rho_c - plate.rho_m) * vc
    return E, rho
