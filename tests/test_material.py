import numpy as np
import pytest
from fgm_plate import Plate, ceramic_fraction, properties


def test_ceramic_limit_and_edge_values():
    p = Plate()
    assert ceramic_fraction(0, 0, 0, 0, p) == pytest.approx(1.0)
    assert ceramic_fraction(0.4, 0.4, 1, 1, p) == pytest.approx(0.25)
    assert ceramic_fraction(0, 0.6, 1, 5, p) == pytest.approx(0.0)
    E, rho = properties(0.4, 0.4, 1, 1, p)
    assert E == pytest.approx(147.5e9)
    assert rho == pytest.approx(2975.0)


def test_out_of_domain_and_negative_gradation():
    with pytest.raises(ValueError):
        ceramic_fraction(-0.2, 0.4, 1, 1)
    with pytest.raises(ValueError):
        ceramic_fraction(0.2, 0.4, -1, 1)


def test_symmetry_under_coordinate_swap():
    x = np.array([0.1, 0.2, 0.7])
    y = np.array([0.5, 0.6, 0.3])
    assert np.allclose(ceramic_fraction(x,y,1,5),ceramic_fraction(y,x,5,1))
