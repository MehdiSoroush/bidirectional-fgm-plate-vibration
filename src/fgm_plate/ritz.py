"""Classical Kirchhoff--Love Rayleigh--Ritz solution, linear free vibration.

The sine basis enforces w=0 on the edges; it does not impose in-plane
membrane constraints. Do not interpret this as a nonlinear plate solver.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss
from .material import Plate, properties


@dataclass
class RitzResult:
    alpha: float
    beta: float
    frequency_hz: float
    modes: list
    coefficients: np.ndarray
    K: np.ndarray
    M: np.ndarray

    @property
    def relative_coefficients(self):
        first = self.coefficients[self.modes.index((1, 1))]
        if abs(first) < 1e-14:
            raise ValueError('First-mode coefficient is essentially zero')
        return self.coefficients / first


def matrices(alpha, beta, nmax=10, quadrature=48, plate=Plate()):
    """Integrate symmetric plate stiffness/mass using tensor Gauss quadrature."""
    if nmax < 1 or quadrature < 4:
        raise ValueError('Invalid basis size or quadrature order')
    pts, wt = leggauss(quadrature)
    x, wx = (pts + 1) / 2, wt / 2
    y, wy = (pts + 1) / 2, wt / 2
    X, Y = np.meshgrid(x, y, indexing='ij')
    W = np.outer(wx, wy) * plate.a * plate.b
    E, rho = properties(X * plate.a, Y * plate.b, alpha, beta, plate)
    D = E * plate.h**3 / (12 * (1 - plate.nu**2))
    modes = [(m, n) for m in range(1, nmax+1) for n in range(1, nmax+1)]
    phi, xx, yy, xy = [], [], [], []
    for m, n in modes:
        kx, ky = np.pi*m/plate.a, np.pi*n/plate.b
        shape = np.sin(m*np.pi*X) * np.sin(n*np.pi*Y)
        phi.append(shape)
        xx.append(-kx*kx*shape)
        yy.append(-ky*ky*shape)
        xy.append(kx*ky*np.cos(m*np.pi*X)*np.cos(n*np.pi*Y))
    shape = (len(modes), -1)
    P = np.asarray(phi).reshape(shape)
    XX = np.asarray(xx).reshape(shape)
    YY = np.asarray(yy).reshape(shape)
    XY = np.asarray(xy).reshape(shape)
    d = (D*W).ravel()
    mass = (rho*plate.h*W).ravel()
    K = ((XX*d) @ XX.T + (YY*d) @ YY.T
         + plate.nu*((XX*d) @ YY.T + (YY*d) @ XX.T)
         + 2*(1-plate.nu)*(XY*d) @ XY.T)
    M = (P*mass) @ P.T
    return (K+K.T)/2, (M+M.T)/2, modes


def solve(alpha, beta, nmax=10, quadrature=48, plate=Plate()):
    K, M, modes = matrices(alpha, beta, nmax=nmax, quadrature=quadrature, plate=plate)
    eigenvalues, vec = eigh(K, M, subset_by_index=[0, 0], check_finite=True)
    if eigenvalues[0] <= 0:
        raise ValueError('Nonpositive eigenvalue')
    return RitzResult(alpha, beta, np.sqrt(eigenvalues[0])/(2*np.pi),
                      modes, vec[:, 0].copy(), K, M)


def evaluate_shape(result, x, y, plate=Plate()):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    shape = np.zeros(np.broadcast_shapes(x.shape, y.shape), dtype=float)
    for (m, n), c in zip(result.modes, result.coefficients):
        shape += c*np.sin(m*np.pi*x/plate.a)*np.sin(n*np.pi*y/plate.b)
    return shape


def exact_homogeneous_frequency(plate=Plate()):
    """Closed-form (1,1) frequency for a homogeneous ceramic Kirchhoff plate."""
    D = plate.E_c*plate.h**3/(12*(1-plate.nu**2))
    omega = np.pi**2 * np.sqrt(D/(plate.rho_c*plate.h))*(1/plate.a**2+1/plate.b**2)
    return omega/(2*np.pi)
