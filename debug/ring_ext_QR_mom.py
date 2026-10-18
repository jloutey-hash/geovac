"""QR step 2b: independent momentum-space route for Coulomb rings.
Kernel momenta k_e on the n Coulomb edges (radial measure d^3k/k^2 = dk dOmega), densities on the
vertices between consecutive edges, rho~_v(|k - k'|).  The partial waves of the fingerprint
transform are geometric:  (1 + c(k^2 + k'^2) - 2 c k k' x)^(-3/2)
   = (h/(c k k'))^{3/2} (1 - h^2)^{-1} sum_l (2l+1) h^l P_l(x),   z = (c(k^2+k'^2)+1)/(2ckk'),
h = z - sqrt(z^2 - 1); the ring's angular integrals give (4pi)^n (2l+1), and sum (2l+1) H^l
= (1+H)/(1-H)^2.  Result:
   Ring_n = (4pi)^{2n} (2pi)^{-3n} int_{R+^n} d^n k  prod_v w_v  (1 + H)/(1 - H)^2,
   w_v = (1/(2 c_v k k')) [1/Delta_v(k - k') - 1/Delta_v(k + k')],
   h_v = (Delta_v(k + k') - Delta_v(k - k'))^2 / (4 c_v k k'),  H = prod h_v.
Tensor Gauss-Legendre with k = s u/(1-u)."""
import numpy as np
import json, sys
from pathlib import Path
from numpy.polynomial.legendre import leggauss

def ring_mom(cs, N, s=1.0):
    n = len(cs)
    x, w = leggauss(N); u = (x + 1)/2; wu = w/2
    k = s*u/(1 - u); wk = wu*s/(1 - u)**2
    shape = [N]*n
    K = [k.reshape([N if i == j else 1 for i in range(n)]) for j in range(n)]
    W = [wk.reshape([N if i == j else 1 for i in range(n)]) for j in range(n)]
    prodw = 1.0; H = 1.0
    for v in range(n):
        a, b, c = K[v-1], K[v], cs[v]
        Dm = np.sqrt(c*(a - b)**2 + 1); Dp = np.sqrt(c*(a + b)**2 + 1)
        prodw = prodw*(1/Dm - 1/Dp)/(2*c*a*b)
        H = H*(Dp - Dm)**2/(4*c*a*b)
    f = prodw*(1 + H)/(1 - H)**2
    for v in range(n):
        f = f*W[v]
    return (4*np.pi)**(2*n)/(2*np.pi)**(3*n)*np.sum(f)

if __name__ == "__main__":
    rec = {}
    # triangle control (n = 3): scales (1,1,t) -> c = (1, 1, 1/t^2)
    for t in [0.7, 1.6]:
        vals = [ring_mom([1.0, 1.0, 1/t**2], N) for N in (80, 120, 160)]
        print("triangle t", t, vals, flush=True)
        rec[f"tri_{t}"] = vals
    for t in [0.7, 1.6]:
        vals = [ring_mom([1.0, 1.0, 1.0, 1/t**2], N) for N in (40, 60, 80, 96)]
        print("square t", t, vals, flush=True)
        rec[f"sq_{t}"] = vals
    Path("debug/data/ring_ext_QR_mom.json").write_text(json.dumps(rec, indent=1))
