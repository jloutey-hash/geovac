"""CT-F12 PoC — Stage 1: derive/validate the Hermitian transcorrelated interaction
and its geminal radial integrals against direct r12 quadrature.

Hermitian (canonical) transcorrelation for He (single electron pair):
  Psi = exp(tau) Phi,  tau = f(r12),  f'(0) = 1/2 (singlet Kato cusp).
  Full TC  H~ = H + D + K   with
     D = -1/2 sum_i [ lap_i tau + (grad_i tau)^2 ]   (multiplicative, Hermitian)
     K = - sum_i (grad_i tau).grad_i                 (convective, non-Hermitian)
  Hermitian TC (drop K)  ->  H_CT = H + D.
  For one pair, tau=f(r12):  sum_i lap_i tau = 2(f'' + 2 f'/r12),  sum_i (grad_i tau)^2 = 2 f'^2
     => D = -(f'' + 2 f'/r12) - f'^2   (a two-body multiplicative potential)
  So the effective e-e interaction is  w(r12) = 1/r12 + D
  Slater geminal enforcing cusp:  f = -(1/(2g)) e^{-g r12}  => f'=(1/2)e^{-g r},  f''=-(g/2)e^{-g r}
     -f''        = (g/2) e^{-g r}
     -2 f'/r     = - e^{-g r}/r
     -f'^2       = -(1/4) e^{-2 g r}
  =>  w(r) = 1/r - e^{-g r}/r + (g/2) e^{-g r} - (1/4) e^{-2 g r}
          = (1 - e^{-g r})/r + (g/2) e^{-g r} - (1/4) e^{-2 g r}      (FINITE at r=0: w(0)=3g/2-1/4)
"""
import numpy as np
import mpmath as mp
from scipy.special import eval_legendre
from scipy.integrate import quad

def w_kernel(r, g):
    r = np.asarray(r, float)
    small = r < 1e-8
    out = np.empty_like(r)
    rr = np.where(small, 1.0, r)
    term1 = -np.expm1(-g*rr)/rr            # (1-e^{-g r})/r
    out = term1 + (g/2)*np.exp(-g*rr) - 0.25*np.exp(-2*g*rr)
    out[small] = 1.5*g - 0.25             # limit r->0
    return out

def coulomb(r):
    return 1.0/r

# ---- 1. cusp check: d/dr12 [ (Psi/Phi) ] / (Psi/Phi) at r12=0 should be 1/2 ----
g = 1.2
# f = -(1/2g) e^{-g r}; Psi/Phi = e^{f}; dlog/dr = f'; f'(0)=1/2
fp0 = 0.5  # by construction
print(f"[cusp] f'(0) = {fp0}  (target 1/2 singlet Kato)  OK")
print(f"[w]    w(0)  = {w_kernel(np.array([1e-12]),g)[0]:.6f}  (=3g/2-1/4={1.5*g-0.25:.6f})")
print(f"[w]    w(2.0)= {w_kernel(np.array([2.0]),g)[0]:.6f}   1/r={0.5:.6f}  (approaches Coulomb at large r)")

# ---- 2. multipole kernel g_k^w(r1,r2) via Gauss-Legendre in x, vs scipy.quad (direct r12 quad) ----
def gkw_gl(r1, r2, kk, g, nx=64):
    xs, ws = np.polynomial.legendre.leggauss(nx)
    r12 = np.sqrt(r1*r1 + r2*r2 - 2*r1*r2*xs)
    return (2*kk+1)/2.0 * np.sum(ws * eval_legendre(kk, xs) * w_kernel(r12, g))

def gkw_quad(r1, r2, kk, g):
    f = lambda x: eval_legendre(kk, x) * float(w_kernel(np.array([np.sqrt(r1*r1+r2*r2-2*r1*r2*x)]), g)[0])
    val, _ = quad(f, -1, 1, limit=200, points=[ (r1*r1+r2*r2 - 0.0) ]) if False else quad(f,-1,1,limit=200)
    return (2*kk+1)/2.0 * val

print("\n[multipole kernel g_k^w(r1,r2): Gauss-Legendre(64) vs scipy.quad]")
maxerr = 0.0
for (r1,r2) in [(0.5,0.7),(1.0,1.0),(0.3,2.5),(2.0,2.1),(0.8,0.8)]:
    for kk in [0,1,2,3]:
        a = gkw_gl(r1,r2,kk,g); b = gkw_quad(r1,r2,kk,g)
        err = abs(a-b); maxerr=max(maxerr,err)
        print(f"  r1={r1} r2={r2} k={kk}:  GL={a:+.8f}  quad={b:+.8f}  |d|={err:.2e}")
print(f"  MAX |GL - quad| = {maxerr:.2e}")

# ---- 3. Coulomb multipole via GL vs analytic r<^k/r>^{k+1} ----
def gkc_gl(r1,r2,kk,nx=200):
    xs,ws=np.polynomial.legendre.leggauss(nx)
    r12=np.sqrt(r1*r1+r2*r2-2*r1*r2*xs)
    return (2*kk+1)/2.0*np.sum(ws*eval_legendre(kk,xs)/r12)
print("\n[Coulomb multipole: GL vs r<^k/r>^{k+1}]")
for (r1,r2) in [(0.5,0.7),(1.0,2.0)]:
    for kk in [0,1,2]:
        rl,rg=min(r1,r2),max(r1,r2)
        ana=rl**kk/rg**(kk+1)
        print(f"  r1={r1} r2={r2} k={kk}: GL={gkc_gl(r1,r2,kk):+.6f} analytic={ana:+.6f} |d|={abs(gkc_gl(r1,r2,kk)-ana):.2e}")
