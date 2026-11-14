"""Regression against coefficient literals in the numerical Mathematica notebook."""
import csv
from pathlib import Path
import pytest
from fgm_plate.analytical_single_mode import source_coefficients

ROOT = Path(__file__).resolve().parents[1]


def load_references():
    with (ROOT/"data/reference/nonlinear_source_coefficients.csv").open(
        newline=""
    ) as f:
        return list(csv.DictReader(f))


@pytest.mark.parametrize("row", load_references(), ids=lambda r: r["case"])
def test_source_coefficients_match_notebook(row):
    c = source_coefficients(float(row["alpha"]), float(row["beta"]),
                            slenderness=float(row["slenderness"]),
                            aspect_ratio=float(row["aspect_ratio"]))
    # Original notebook values were printed to approximately six digits.
    assert c.p == pytest.approx(float(row["source_P"]), rel=3e-6)
    assert c.q == pytest.approx(float(row["source_Q"]), rel=3e-6)


def test_coordinate_swap_in_square_plate():
    for a,b in [(1,5),(0.5,2),(0,1)]:
        x,y=source_coefficients(a,b),source_coefficients(b,a)
        assert x.p == pytest.approx(y.p, rel=1e-14)
        assert x.q == pytest.approx(y.q, rel=1e-14)


def test_parameters():
    with pytest.raises(ValueError): source_coefficients(-0.5,1)
    with pytest.raises(ValueError): source_coefficients(1,1,aspect_ratio=0)
    with pytest.raises(ValueError): source_coefficients(1,1,m=0)


def test_normalized_nonlinearity_independent_of_gradation_in_source_model():
    ref = source_coefficients(0,0)
    ratio = ref.q/ref.p
    for a,b in [(1,0),(1,1),(0.5,0.5),(5,5),(1,5)]:
        c=source_coefficients(a,b)
        assert c.q/c.p == pytest.approx(ratio, rel=1e-14)
    assert ratio == pytest.approx(822.65625, rel=1e-13)
