"""Variational single-transverse-mode von Karman ROM with spatial E(x,y).

Transverse w=q(t)*phi(x,y) has zero edge deflection. In-plane u,v vanish
at every edge and relax quasi-statically through a sine--sine Ritz basis.
The model employs a spatially varying constitutive law and a variational,
quasi-static condensation of in-plane membrane degrees of freedom.
"""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve
from .material import Plate, properties
from .ritz import solve as solve_linear_ritz


@dataclass(frozen=True)
class MembraneReduction:
    alpha: float
    beta: float
    mode_shape: str
    transverse_terms: int
    inplane_terms: int
    transverse_mass: float
    bending_stiffness: float
    membrane_unrelaxed: float
    membrane_relaxed: float
    linear_coefficient: float
    cubic_coefficient: float
    maximum_inplane_coeff: float
    uniform_pressure_gain: float

    @property
    def linear_frequency_hz(self):
        return np.sqrt(self.linear_coefficient)/(2*np.pi)

    @property
    def ratio_cubic_to_linear(self):
        return self.cubic_coefficient/self.linear_coefficient

    @property
    def energy_relaxation_fraction(self):
        return 1-self.membrane_relaxed/self.membrane_unrelaxed


def _transverse_coefficients(alpha,beta,plate,shape,nmax,quadrature):
    if shape == 'sine':
        return [(1,1)],np.array([1.0])
    if shape != 'ritz':
        raise ValueError("mode_shape must be 'sine' or 'ritz'")
    r=solve_linear_ritz(alpha,beta,nmax=nmax,quadrature=quadrature,plate=plate)
    center=sum(c*np.sin(m*np.pi/2)*np.sin(n*np.pi/2)
               for (m,n),c in zip(r.modes,r.coefficients))
    if abs(center)<1e-14:
        raise ValueError('Ritz bending mode has zero center displacement')
    return r.modes,r.coefficients/center


def spatial_nonlinear_coefficients(alpha,beta,*,plate=Plate(),
        mode_shape='sine',transverse_terms=4,inplane_terms=4,quadrature=36):
    """Return P=K_b/M and Q=2*C/M for the immovable-edge plate.

    q has units of metres; phi(a/2,b/2)=1. The strain vector is
      [u_x + w_x^2/2, v_y + w_y^2/2, u_y+v_x+w_x*w_y].
    Membrane energy is minimized over an in-plane sine--sine Ritz space.
    """
    if alpha<0 or beta<0 or inplane_terms<1 or transverse_terms<1 or quadrature<8:
        raise ValueError('Invalid material, basis or quadrature parameter')
    if not isinstance(inplane_terms,int) or not isinstance(transverse_terms,int):
        raise ValueError('Basis sizes must be integers')
    modes,coeff=_transverse_coefficients(alpha,beta,plate,mode_shape,
                                         transverse_terms,quadrature)
    roots,weights=leggauss(quadrature)
    x,wx=(roots+1)/2,weights/2
    y,wy=(roots+1)/2,weights/2
    X,Y=np.meshgrid(x,y,indexing='ij')
    wt=np.outer(wx,wy).reshape(-1)*plate.a*plate.b
    E,rho=properties(X*plate.a,Y*plate.b,alpha,beta,plate)
    E,rho=E.ravel(),rho.ravel()
    X,Y=X.ravel(),Y.ravel()
    phi=np.zeros_like(X);px=np.zeros_like(X);py=np.zeros_like(X)
    pxx=np.zeros_like(X);pyy=np.zeros_like(X);pxy=np.zeros_like(X)
    for (m,n),v in zip(modes,coeff):
        sx,cx=np.sin(m*np.pi*X),np.cos(m*np.pi*X)
        sy,cy=np.sin(n*np.pi*Y),np.cos(n*np.pi*Y)
        kx=m*np.pi/plate.a;ky=n*np.pi/plate.b
        phi+=v*sx*sy
        px+=v*kx*cx*sy;py+=v*ky*sx*cy
        pxx-=v*kx*kx*sx*sy;pyy-=v*ky*ky*sx*sy
        pxy+=v*kx*ky*cx*cy
    nu=plate.nu
    mass=np.sum(rho*plate.h*wt*phi*phi)
    pressure_gain=np.sum(wt*phi)/mass
    D=E*plate.h**3/(12*(1-nu**2))
    bending=np.sum(D*wt*(pxx*pxx+pyy*pyy+2*nu*pxx*pyy+2*(1-nu)*pxy*pxy))
    if mass<=0 or bending<=0:
        raise ArithmeticError('Nonpositive mass or bending stiffness')
    g0=.5*px*px;g1=.5*py*py;g2=px*py
    A=E*plate.h/(1-nu**2)
    shear=(1-nu)/2
    C0=np.sum(wt*A*(g0*g0+g1*g1+2*nu*g0*g1+shear*g2*g2))
    N=inplane_terms**2;dof=2*N
    Bx=np.zeros((dof,len(wt)));By=np.zeros_like(Bx);Bxy=np.zeros_like(Bx)
    for i,(m,n) in enumerate((m,n) for m in range(1,inplane_terms+1)
                                   for n in range(1,inplane_terms+1)):
        sx,cx=np.sin(m*np.pi*X),np.cos(m*np.pi*X)
        sy,cy=np.sin(n*np.pi*Y),np.cos(n*np.pi*Y)
        phix=m*np.pi/plate.a*cx*sy
        phiy=n*np.pi/plate.b*sx*cy
        Bx[i]=phix;Bxy[i]=phiy
        By[N+i]=phiy;Bxy[N+i]=phix
    AW=wt*A;KW=AW*shear
    K=((Bx*AW)@Bx.T+(By*AW)@By.T
       +nu*((Bx*AW)@By.T+(By*AW)@Bx.T)+(Bxy*KW)@Bxy.T)
    K=(K+K.T)/2
    force=(Bx@(AW*(g0+nu*g1))+By@(AW*(g1+nu*g0))
           +Bxy@(KW*g2))
    inplane=solve(K,force,assume_a='pos',check_finite=True)
    C=C0-float(force@inplane)
    if not 0<=C<=C0*(1+1e-10):
        raise ArithmeticError('Membrane relaxation violates energy bounds')
    if C<=0:
        raise ArithmeticError('Degenerate membrane model')
    return MembraneReduction(alpha,beta,mode_shape,transverse_terms,
                             inplane_terms,float(mass),float(bending),
                             float(C0),float(C),float(bending/mass),
                             float(2*C/mass),float(np.max(np.abs(inplane))),
                             float(pressure_gain))
