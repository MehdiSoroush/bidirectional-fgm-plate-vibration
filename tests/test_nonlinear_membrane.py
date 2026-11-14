"""Independent variational membrane reduction checks."""
import numpy as np
import pytest
from fgm_plate.material import Plate
from fgm_plate.ritz import solve
from fgm_plate.nonlinear_membrane import spatial_nonlinear_coefficients as calc
from fgm_plate.nonlinear import DuffingModel


def test_homogeneous_sine_frequency_analytical():
    r=calc(0,0,inplane_terms=3,quadrature=24)
    exact=solve(0,0,nmax=1,quadrature=24).frequency_hz
    assert r.linear_frequency_hz==pytest.approx(exact,rel=3e-13)
    assert r.cubic_coefficient>0
    assert r.uniform_pressure_gain==pytest.approx(16/(np.pi**2*3800*.04),rel=1e-13)
    assert 0<r.membrane_relaxed<r.membrane_unrelaxed


def test_square_gradation_swap_symmetry():
    x=calc(1,5,inplane_terms=3,quadrature=28)
    y=calc(5,1,inplane_terms=3,quadrature=28)
    assert x.linear_coefficient==pytest.approx(y.linear_coefficient,rel=3e-12)
    assert x.cubic_coefficient==pytest.approx(y.cubic_coefficient,rel=3e-12)


def test_condensed_membrane_energy_nonincreasing():
    r=[calc(1,1,inplane_terms=n,quadrature=28) for n in (1,2,3,4)]
    assert all(v.cubic_coefficient>0 for v in r)
    assert all(r[i+1].membrane_relaxed<=r[i].membrane_relaxed*(1+1e-10)
               for i in range(3))
    assert r[-1].membrane_relaxed<r[0].membrane_unrelaxed


def test_spatial_gradation_breaks_average_invariant():
    rs=[calc(a,b,inplane_terms=3,quadrature=28) for a,b in [(0,0),(1,0),(1,1),(5,5)]]
    ratios=[r.ratio_cubic_to_linear for r in rs]
    assert max(ratios)-min(ratios)>0.02*ratios[0]


def test_enriched_shape_linear_frequency_matches_ritz():
    r=calc(1,1,mode_shape='ritz',transverse_terms=3,
           inplane_terms=3,quadrature=28)
    ref=solve(1,1,nmax=3,quadrature=28)
    assert r.linear_frequency_hz==pytest.approx(ref.frequency_hz,rel=2e-13)


def test_thickness_scaling_for_immovable_membrane_model():
    thick=calc(1,1,plate=Plate(h=.04),inplane_terms=2,quadrature=24)
    thin=calc(1,1,plate=Plate(h=.02),inplane_terms=2,quadrature=24)
    assert thin.linear_coefficient==pytest.approx(thick.linear_coefficient/4,rel=1e-12)
    assert thin.cubic_coefficient==pytest.approx(thick.cubic_coefficient,rel=1e-12)


def test_energy_model_is_valid_duffing_ode():
    result=calc(1,1,inplane_terms=2,quadrature=20)
    model=DuffingModel(result.linear_coefficient,result.cubic_coefficient)
    t,x,v=model.integrate(t_end=.002,q0=.004,num=100)
    en=model.mechanical_energy(x,v)
    assert np.max(np.abs(en-en[0]))/en[0]<1e-8


def test_invalid_parameters():
    with pytest.raises(ValueError): calc(-1,1)
    with pytest.raises(ValueError): calc(1,1,inplane_terms=0)
    with pytest.raises(ValueError): calc(1,1,mode_shape='unphysical')
