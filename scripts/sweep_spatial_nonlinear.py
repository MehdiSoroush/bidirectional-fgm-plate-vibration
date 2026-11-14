"""Reproduce the spatially graded immovable-edge von Karman ROM benchmarks.

Writes the benchmark CSV and two PNG/SVG scientific figures.
The model uses immovable in-plane edges; linear Abaqus employs minimal pins.
"""
import argparse
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fgm_plate.analytical_single_mode import source_coefficients
from fgm_plate.nonlinear_membrane import spatial_nonlinear_coefficients

CASES=[(0,0),(1,0),(1,1),(1,5),(5,1),(5,5)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--transverse-terms',type=int,default=10)
    parser.add_argument('--inplane-terms',type=int,default=8)
    parser.add_argument('--quadrature',type=int,default=48)
    args=parser.parse_args()
    data_path=ROOT/'data/reference/spatial_nonlinear_immovable.csv'
    fig_dir=ROOT/'results/figures'
    fig_dir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for a,b in CASES:
        r=spatial_nonlinear_coefficients(a,b,mode_shape='ritz',
            transverse_terms=args.transverse_terms,
            inplane_terms=args.inplane_terms,quadrature=args.quadrature)
        s=source_coefficients(a,b)
        rows.append(dict(alpha=a,beta=b,transverse_terms=args.transverse_terms,
            inplane_terms=args.inplane_terms,quadrature=args.quadrature,
            frequency_hz=r.linear_frequency_hz,P_s_minus2=r.linear_coefficient,
            Q_m_minus2_s_minus2=r.cubic_coefficient,
            Q_over_P_m_minus2=r.ratio_cubic_to_linear,
            relaxed_membrane_energy_fraction=r.energy_relaxation_fraction,
            pressure_gain_m_per_s2_per_Pa=r.uniform_pressure_gain,
            archived_average_Q_over_P_m_minus2=s.q/s.p))
    with data_path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    x=np.arange(len(rows))
    labels=[f'({r["alpha"]},{r["beta"]})' for r in rows]
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    ax.plot(x,[r['Q_over_P_m_minus2'] for r in rows],
            'o-',label='Spatial Ritz + immovable-edge membrane')
    ax.plot(x,[r['archived_average_Q_over_P_m_minus2'] for r in rows],
            's--',label='Archived averaged single-mode reduction')
    ax.set_xticks(x,labels)
    ax.set(xlabel='Gradation (alpha, beta)',
           ylabel=r'Normalized cubic stiffness Q/P (m$^{-2}$)',
           title='Comparison of two distinct nonlinear plate reductions')
    ax.legend(frameon=False)
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(fig_dir/'nonlinear_spatial_gradation.png',dpi=200)
    fig.savefig(fig_dir/'nonlinear_spatial_gradation.svg')
    plt.close(fig)
    nvals=list(range(1,args.inplane_terms+1))
    qr=[]
    for n in nvals:
        r=spatial_nonlinear_coefficients(1,1,mode_shape='ritz',
            transverse_terms=args.transverse_terms,
            inplane_terms=n,quadrature=args.quadrature)
        qr.append(r.cubic_coefficient/1e9)
    fig,ax=plt.subplots(figsize=(8.0,4.6))
    ax.plot(nvals,qr,'o-')
    ax.set(xlabel='In-plane sine modes per direction',
           ylabel=r'Cubic stiffness Q ($10^9$ m$^{-2}$ s$^{-2}$)',
           title='In-plane Ritz relaxation convergence: (1,1) gradation')
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(fig_dir/'nonlinear_inplane_convergence.png',dpi=200)
    fig.savefig(fig_dir/'nonlinear_inplane_convergence.svg')
    plt.close(fig)
    print('Wrote:',data_path)
    print('Figures in:',fig_dir)
    print('Spatial variational model benchmarks and figures generated.')


if __name__=='__main__':
    main()
