"""Calculate unweighted MAC with Abaqus 2021 U3 CSV exports.

Uses versioned Abaqus displacement datasets in data/abaqus_u3/.
"""
from pathlib import Path
from fgm_plate import solve, compare_to_ritz

DATA = Path(__file__).resolve().parents[1]/'data'/'abaqus_u3'
CASES = [(0,0,'mode_00_n40.csv'),(1,1,'mode_11_n40.csv'),(1,5,'mode_15_n40.csv')]


def main():
    for alpha, beta, filename in CASES:
        path = DATA/filename
        if not path.exists():
            raise SystemExit('Required recorded U3 data missing: %s' % path)
        res = solve(alpha, beta, nmax=10)
        stats = compare_to_ritz(path, res)
        print('(alpha,beta)=(%g,%g), n=%d; Ritz=%.5f Hz; sine MAC=%.9f; full MAC=%.9f; residual=%.3f%%'
              % (alpha,beta,stats['nodes'],res.frequency_hz,stats['single_sine_mac'],
                 stats['full_ritz_mac'],100*stats['full_ritz_residual']))


if __name__ == '__main__':
    main()
