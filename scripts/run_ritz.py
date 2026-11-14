"""Run the documented 10x10 Rayleigh--Ritz reference cases."""
import argparse
import csv
from pathlib import Path
from fgm_plate import solve


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nmax', type=int, default=10)
    parser.add_argument('--quadrature', type=int, default=48)
    parser.add_argument('--output', type=Path, default=Path('results/ritz_frequencies.csv'))
    args = parser.parse_args()
    rows = []
    for alpha, beta in [(0,0),(1,1),(5,5),(1,5),(5,1)]:
        r = solve(alpha, beta, nmax=args.nmax, quadrature=args.quadrature)
        rows.append({'alpha':alpha,'beta':beta,'nmax':args.nmax,'quadrature':args.quadrature,
                     'frequency_hz':'%.10f' % r.frequency_hz})
        print('(%g,%g): %.6f Hz' % (alpha, beta, r.frequency_hz))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
