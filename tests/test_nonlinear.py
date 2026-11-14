"""ODE verification: linear limit, conserved energy and independent period."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from fgm_plate.analytical_single_mode import source_coefficients
from fgm_plate.nonlinear import (
    DuffingModel, leading_order_free_frequency, exact_free_angular_frequency,
)


def test_linear_limit_matches_closed_form():
    m = DuffingModel(p=4.0, q=0.0)
    t,x,v = m.integrate(t_end=8, q0=.02, v0=-.03, num=400)
    np.testing.assert_allclose(x,.02*np.cos(2*t)-.015*np.sin(2*t),atol=3e-10)
    np.testing.assert_allclose(v,-.04*np.sin(2*t)-.03*np.cos(2*t),atol=8e-10)


def test_energy_conservation_unforced():
    c=source_coefficients(1,1)
    model=DuffingModel(c.p,c.q)
    t,x,v=model.integrate(t_end=.04,q0=.02,num=1801)
    energy=model.mechanical_energy(x,v)
    assert np.max(np.abs(energy-energy[0]))/energy[0] < 2e-8


def test_forced_time_response_is_finite_and_initial_acceleration():
    model=DuffingModel(1.38808e7,1.14191e10,5e6,1500)
    assert model.rhs(0.0,[.02,0.0])[1] == pytest.approx(
        -model.p*.02-model.q*.02**3)
    _,x,v=model.integrate(t_end=.006,q0=.02,num=300)
    assert np.isfinite(x).all() and np.isfinite(v).all()


def test_free_frequency_small_amplitude():
    p,q=9e6,1e9
    assert leading_order_free_frequency(p,q,0)==pytest.approx(np.sqrt(p))
    assert leading_order_free_frequency(p,q,.001)>np.sqrt(p)


def test_exact_energy_period_against_time_integration():
    c=source_coefficients(1,1)
    amplitude=.008
    model=DuffingModel(c.p,c.q)
    expected=exact_free_angular_frequency(c.p,c.q,amplitude)
    def descending_zero(t,state):
        return state[0]
    descending_zero.direction=-1
    descending_zero.terminal=True
    solution=solve_ivp(model.rhs,(0,.01),(amplitude,0),
        events=descending_zero,rtol=1e-11,atol=1e-12,method="DOP853")
    assert solution.status==1
    measured=2*np.pi/(4*solution.t_events[0][0])
    assert measured == pytest.approx(expected,rel=5e-8)


def test_exact_free_frequency_limits():
    c=source_coefficients(0,0)
    omega0=np.sqrt(c.p)
    assert exact_free_angular_frequency(c.p,c.q,0)==pytest.approx(omega0)
    assert exact_free_angular_frequency(c.p,0,.03)==pytest.approx(omega0)
    a=.0002
    exact=exact_free_angular_frequency(c.p,c.q,a)
    approx=leading_order_free_frequency(c.p,c.q,a)
    assert exact == pytest.approx(approx,rel=1e-7)
    assert exact>omega0


def test_invalid_ode_inputs():
    with pytest.raises(ValueError): DuffingModel(-1,2)
    with pytest.raises(ValueError): DuffingModel(1,-2)
    with pytest.raises(ValueError): DuffingModel(1,2).integrate(t_end=0)
