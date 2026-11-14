"""Undamped cubic oscillator from the analytical time-response notebook.

Equation: q'' + P*q + Q*q^3 = F*sin(Omega*t).
The forcing amplitude is defined in the mass-normalized ODE; conversion
from a physical distributed pressure depends on the modal load projection.
"""
from dataclasses import dataclass
from math import sin, sqrt
from typing import Tuple
import numpy as np
from scipy.integrate import quad, solve_ivp


@dataclass(frozen=True)
class DuffingModel:
    p: float
    q: float
    forcing_amplitude: float = 0.0
    forcing_frequency: float = 0.0

    def __post_init__(self):
        if self.p <= 0 or self.q < 0 or self.forcing_frequency < 0:
            raise ValueError("Require P>0, Q>=0, and Omega>=0")

    @property
    def linear_angular_frequency(self) -> float:
        return sqrt(self.p)

    def rhs(self, t: float, state: np.ndarray) -> np.ndarray:
        displacement, velocity = state
        acceleration = (-self.p*displacement - self.q*displacement**3
                        + self.forcing_amplitude*sin(self.forcing_frequency*t))
        return np.array([velocity, acceleration])

    def mechanical_energy(self, displacement: np.ndarray,
                          velocity: np.ndarray) -> np.ndarray:
        """Conserved in the absence of forcing."""
        return (0.5*np.asarray(velocity)**2
                + 0.5*self.p*np.asarray(displacement)**2
                + 0.25*self.q*np.asarray(displacement)**4)

    def integrate(self, *, t_end: float = 0.03, q0: float = 0.02,
                  v0: float = 0.0, num: int = 1501,
                  rtol: float = 1e-10, atol: float = 1e-12
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        if t_end <= 0 or num < 2:
            raise ValueError("Positive t_end and at least two samples required")
        t = np.linspace(0, t_end, num)
        solution = solve_ivp(self.rhs, (0, t_end), (q0, v0),
                             method="DOP853", t_eval=t, rtol=rtol, atol=atol,
                             max_step=min(t_end/20, 1/10000))
        if not solution.success:
            raise RuntimeError("Numerical integration failed: " + solution.message)
        return t, solution.y[0], solution.y[1]


def leading_order_free_frequency(p: float, q: float, amplitude: float) -> float:
    """Weak-nonlinearity approximation, not exact at finite amplitude."""
    if p <= 0 or q < 0 or amplitude < 0:
        raise ValueError("Require P>0, Q>=0, amplitude>=0")
    return sqrt(p)*(1 + 3*q*amplitude**2/(8*p))


def exact_free_angular_frequency(p: float, q: float, amplitude: float) -> float:
    """Conservative oscillator period from the energy integral."""
    if p <= 0 or q < 0 or amplitude < 0:
        raise ValueError("Require P>0, Q>=0, amplitude>=0")
    if amplitude == 0 or q == 0:
        return sqrt(p)
    integrand = lambda theta: (
        p + 0.5*q*amplitude**2*(1+np.sin(theta)**2)
    )**(-0.5)
    quarter_period = quad(integrand, 0, np.pi/2,
                          epsabs=1e-12, epsrel=1e-12)[0]
    return 2*np.pi/(4*quarter_period)
