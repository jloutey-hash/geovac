"""Validation gate: two-center geminal (Slater-TC w-kernel) ERI vs an INDEPENDENT
direct r12 double-grid Cartesian quadrature. Backs the LiH molecular xTC sparsity
run (xtc_lih_molecular.py). The w-ERI uses the same multipole assembly as the
validated TwoCenterLM Coulomb ERI with the radial kernel replaced by the Legendre
multipole of w(r12) = 1/r + D (Slater geminal, Kato cusp)."""
import sys; sys.path.insert(0, 'debug')
import numpy as np
from math import factorial
from scipy.special import lpmv
from two_center_grid_lm import TwoCenterLM
import xtc_lih_molecular as X

def Ylm_real(l, m, u, phi):
    mm = abs(m)
    if m == 0: return np.sqrt((2*l+1)/(4*np.pi))*lpmv(0, l, u)
    norm = np.sqrt((2*l+1)/(2*np.pi)*factorial(l-mm)/factorial(l+mm))
    ang = np.cos(mm*phi) if m > 0 else np.sin(mm*phi)
    return norm*lpmv(mm, l, u)*ang
def sto_norm(z, l): return 1.0/np.sqrt(factorial(2*l+2)/(2*z)**(2*l+3))
def chi_at(orb, XX, YY, ZZ, R):
    z, l, m, c = orb; zc = ZZ if c == 'A' else ZZ-R
    rc = np.sqrt(XX*XX+YY*YY+zc*zc)+1e-300; uc = np.clip(zc/rc, -1, 1); ph = np.arctan2(YY, XX)
    return sto_norm(z, l)*rc**l*np.exp(-z*rc)*Ylm_real(l, m, uc, ph)
def quad_grid(nr=30, nu=14, nphi=16, rmax=40.0):
    t = np.linspace(0, 1, nr); r = rmax*t**2+1e-6
    wr = np.zeros_like(r); wr[1:-1] = (r[2:]-r[:-2])/2; wr[0] = (r[1]-r[0])/2; wr[-1] = (r[-1]-r[-2])/2
    u, wu = np.polynomial.legendre.leggauss(nu); phi = (np.arange(nphi)+0.5)*2*np.pi/nphi
    R3, U3, P3 = np.meshgrid(r, u, phi, indexing='ij'); s = np.sqrt(np.clip(1-U3**2, 0, 1))
    XX = (R3*s*np.cos(P3)).ravel(); YY = (R3*s*np.sin(P3)).ravel(); ZZ = (R3*U3).ravel()
    WR, WU, WP = np.meshgrid(wr, wu, np.full(nphi, 2*np.pi/nphi), indexing='ij')
    return XX, YY, ZZ, (R3**2*WR*WU*WP).ravel()
def direct_eri(a, b, c, d, R, grid, g):
    XX, YY, ZZ, vol = grid
    rho1 = chi_at(a, XX, YY, ZZ, R)*chi_at(b, XX, YY, ZZ, R)*vol
    rho2 = chi_at(c, XX, YY, ZZ, R)*chi_at(d, XX, YY, ZZ, R)*vol
    tot = 0.0
    for i0 in range(0, len(XX), 400):
        i1 = min(i0+400, len(XX))
        dx = XX[i0:i1, None]-XX[None, :]; dy = YY[i0:i1, None]-YY[None, :]; dz = ZZ[i0:i1, None]-ZZ[None, :]
        tot += float(rho1[i0:i1] @ (X.w_kernel(np.sqrt(dx*dx+dy*dy+dz*dz), g) @ rho2))
    return tot

if __name__ == '__main__':
    R, g = 3.0, 1.2
    eng = TwoCenterLM(R, nr=1000, nu=28, nphi=28, rmax=45.0, Lmax=10, real=True)
    FLw = X.multipole_kernels(eng.r, lambda rr: X.w_kernel(rr, g), 10)
    FLc = X.multipole_kernels(eng.r, lambda rr: 1.0/np.maximum(rr, 1e-30), 10)
    grid = quad_grid(30, 14, 16, 40.0)
    cases = [((2.7,0,0,'A'),(2.7,0,0,'A'),(1.0,0,0,'B'),(1.0,0,0,'B')),
             ((2.7,0,0,'A'),(1.0,0,0,'B'),(2.7,0,0,'A'),(1.0,0,0,'B')),
             ((0.65,0,0,'A'),(1.0,0,0,'B'),(0.65,0,0,'A'),(1.0,0,0,'B')),
             ((0.65,1,0,'A'),(0.65,1,0,'A'),(1.0,0,0,'B'),(1.0,0,0,'B'))]
    print("w-ERI: multipole vs DIRECT r12 quadrature (g=%.1f) | Coulomb-limit multipole vs engine-native"%g)
    for c in cases:
        wmp = X.eri_kernel_chemist(eng, *c, FLw, 10); wdq = direct_eri(*c, R, grid, g)
        cmp = X.eri_kernel_chemist(eng, *c, FLc, 10); cnat = eng.eri(*c)
        print("  w: mp=%+.6f direct=%+.6f d=%+.1e || coul: mp=%+.6f native=%+.6f d=%+.1e"
              % (wmp, wdq, wmp-wdq, cmp, cnat, cmp-cnat))
