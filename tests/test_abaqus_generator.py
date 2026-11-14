from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def test_generates_minimal_s4r_deck(tmp_path):
    out=tmp_path/'inp'
    subprocess.run([sys.executable,str(ROOT/'abaqus'/'generate_abaqus2021.py'),
                    '--n','3','--alpha','1','--beta','5','--outdir',str(out)],check=True)
    paths=list(out.glob('*.inp'))
    assert len(paths)==1
    text=paths[0].read_text()
    assert text.count('*Elset, elset=E')==9
    assert text.count('*Shell Section,')==9
    assert '*Element, type=S4R' in text
    assert 'EDGE, 3, 3' in text
    assert '*Frequency, eigensolver=Lanczos' in text
    assert '1, 1, 2' in text
