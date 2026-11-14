import pytest
from fgm_plate import Plate, solve, exact_homogeneous_frequency


def test_homogeneous_closed_form():
    exact = exact_homogeneous_frequency()
    result = solve(0,0,nmax=3,quadrature=20)
    assert result.frequency_hz == pytest.approx(exact, rel=1e-10)
    assert result.frequency_hz == pytest.approx(594.1807407444,abs=0.00001)


@pytest.mark.parametrize('case,frequency',[
    ((1,1),402.57036),((5,5),316.5360),((1,5),334.65096),((5,1),334.65096)])
def test_ten_by_ten_benchmarks(case, frequency):
    result=solve(*case,nmax=10,quadrature=48)
    assert result.frequency_hz == pytest.approx(frequency,abs=0.002)


def test_axis_exchange_invariance():
    a=solve(1,5,nmax=4,quadrature=24)
    b=solve(5,1,nmax=4,quadrature=24)
    assert a.frequency_hz==pytest.approx(b.frequency_hz,abs=1e-8)


def test_thickness_scaling_in_classical_model():
    a=solve(0,0,nmax=1,quadrature=10,plate=Plate(h=0.04))
    b=solve(0,0,nmax=1,quadrature=10,plate=Plate(h=0.02))
    assert a.frequency_hz/b.frequency_hz==pytest.approx(2.0,rel=1e-10)
