from pathlib import Path
import pytest
from fgm_plate import mac, solve, compare_to_ritz, load_u3

DATA=Path(__file__).resolve().parents[1]/'data'/'abaqus_u3'


def test_mac_sign_and_shape():
    assert mac([1,2,3],[-1,-2,-3])==pytest.approx(1.0)
    assert mac([1,0],[0,1])==pytest.approx(0.0)


@pytest.mark.parametrize('a,b,filename,sine,full',[
    (0,0,'mode_00_n40.csv',0.999987165,0.999987165),
    (1,1,'mode_11_n40.csv',0.986894443,0.999986420),
    (1,5,'mode_15_n40.csv',0.993918368,0.999973790),
])
def test_recorded_mode_shapes(a,b,filename,sine,full):
    path=DATA/filename
    if not path.is_file():
        pytest.skip('Nodal Abaqus U3 data are not part of the current hosted snapshot')
    result=solve(a,b,nmax=10,quadrature=48)
    stats=compare_to_ritz(path,result)
    assert stats['nodes']==1681
    assert stats['single_sine_mac']==pytest.approx(sine,abs=3e-7)
    assert stats['full_ritz_mac']==pytest.approx(full,abs=3e-7)
    assert stats['full_ritz_mac']>=stats['single_sine_mac']-1e-8


def test_export_node_count_80_by_80():
    path=DATA/'mode_11_n80.csv'
    if not path.is_file():
        pytest.skip('Nodal Abaqus U3 data are not part of the current hosted snapshot')
    x,y,u=load_u3(path)
    assert len(x)==6561
