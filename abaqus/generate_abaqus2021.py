# -*- coding: ascii -*-
"""Write Abaqus 2021 S4R frequency input decks (Python 2.7 compatible).

NOTE: This script only GENERATES input decks. Running an Abaqus job
requires a licensed installation, and the model has not been executed here.
SI units; one material/section per element, based on centroid properties.
"""
from __future__ import print_function
import os
import argparse

A = 0.8
B = 0.8
NU = 0.30
EC = 380.0e9
EA = 70.0e9
RC = 3800.0
RA = 2700.0


def node(i, j, n):
    return j*(n+1)+i+1


def make(alpha, beta, n, h, outpath):
    if min(n, h) <= 0 or n < 2 or alpha < 0 or beta < 0:
        raise ValueError('Mesh, thickness and exponents must be valid')
    lines = ['*Heading',
             '** Al2O3/Al in-plane graded plate: alpha=%s beta=%s; mesh=%dx%d; thickness=%g m' % (alpha,beta,n,n,h),
             '** Element-centroid mixing; SI units; S4R linear eigenmodes',
             '*Node']
    for j in range(n+1):
        for i in range(n+1):
            lines.append('%d, %.12g, %.12g, 0.' % (node(i,j,n),A*i/float(n),B*j/float(n)))
    lines.append('*Element, type=S4R')
    for j in range(n):
        for i in range(n):
            eid = j*n+i+1
            lines.append('%d, %d, %d, %d, %d' % (eid,node(i,j,n),node(i+1,j,n),node(i+1,j+1,n),node(i,j+1,n)))
    edge = sorted(set([node(i,0,n) for i in range(n+1)] +
                      [node(i,n,n) for i in range(n+1)] +
                      [node(0,j,n) for j in range(n+1)] +
                      [node(n,j,n) for j in range(n+1)]))
    lines.append('*Nset, nset=EDGE')
    for k in range(0,len(edge),16):
        lines.append(', '.join(str(v) for v in edge[k:k+16]))
    for j in range(n):
        for i in range(n):
            eid = j*n+i+1
            vc = ((i+0.5)/float(n))**alpha*((j+0.5)/float(n))**beta
            E = EA+(EC-EA)*vc
            rho = RA+(RC-RA)*vc
            lines.extend(['*Elset, elset=E%d' % eid,str(eid),
                          '*Material, name=MAT%d' % eid,
                          '*Elastic', '%.12g, %.12g' % (E,NU),
                          '*Density','%.12g' % rho,
                          '*Shell Section, elset=E%d, material=MAT%d' % (eid,eid),
                          '%.12g' % h])
    # w=0 along edges; minimal in-plane pins to eliminate RB modes only.
    # Linear modal restraints suppress rigid-body motion with minimal in-plane pins.
    lines.extend(['*Boundary','EDGE, 3, 3',
                  '%d, 1, 2' % node(0,0,n),
                  '%d, 2, 2' % node(n,0,n),
                  '*Step, name=FREQUENCY, perturbation',
                  '*Frequency, eigensolver=Lanczos','12,',
                  '*Output, field, frequency=1','*Node Output','U',
                  '*End Step'])
    f = open(outpath, 'w')
    try:
        f.write('\n'.join(lines)+'\n')
    finally:
        f.close()
    return outpath


def main():
    parser = argparse.ArgumentParser(description='Generate Abaqus 2021 S4R inputs')
    parser.add_argument('--n',type=int,default=40)
    parser.add_argument('--h',type=float,default=0.04)
    parser.add_argument('--alpha',type=float,default=1.0)
    parser.add_argument('--beta',type=float,default=1.0)
    parser.add_argument('--all',action='store_true',help='generate five baseline cases')
    parser.add_argument('--outdir',default='inputs')
    args = parser.parse_args()
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    cases = [(0,0),(1,1),(5,5),(1,5),(5,1)] if args.all else [(args.alpha,args.beta)]
    for alpha,beta in cases:
        atext = '%g' % alpha
        btext = '%g' % beta
        path = os.path.join(args.outdir,'fgm_a%s_b%s_n%d_h%g.inp' % (atext,btext,args.n,args.h))
        print(make(alpha,beta,args.n,args.h,path))


if __name__ == '__main__':
    main()
