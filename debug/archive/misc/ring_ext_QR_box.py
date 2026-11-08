"""QR step 0: the 3-D massless one-loop box in closed form.
Melrose reduction (N = D + 1, the dimension-shift coefficient N - D - 1 vanishes):
Box(P0..P3) = int d^3l prod_i 1/|l - P_i|^2 = sum_i b_i Tri_i, b = S^{-1} 1, S_ij = |P_i-P_j|^2,
Tri = pi^3/(d12 d23 d31) (star-triangle).  sympy simplifies the sum to
   Box = pi^3 [sum_i 1/prod_{j != i} d_ij] / (d01 d23 + d02 d13 + d03 d12).
Check: direct Feynman-parameter evaluation, primary sectors with x_j = y_j^2 (removes the
integrable vertex singularity), at random configurations."""
import numpy as np
import json
from pathlib import Path
from scipy import integrate
from scipy.special import gamma

def box_closed(P):
    d = np.array([[np.linalg.norm(P[i]-P[j]) for j in range(4)] for i in range(4)])
    num = sum(1/np.prod([d[i, j] for j in range(4) if j != i]) for i in range(4))
    return np.pi**3*num/(d[0,1]*d[2,3]+d[0,2]*d[1,3]+d[0,3]*d[1,2])

def box_feyn(P):
    s = np.array([[np.sum((P[i]-P[j])**2) for j in range(4)] for i in range(4)])
    tot = 0.0
    for i in range(4):
        o = [j for j in range(4) if j != i]
        def f(y3, y2, y1):
            x = [y1*y1, y2*y2, y3*y3]
            Q = sum(x[a]*s[i, o[a]] for a in range(3)) + sum(x[a]*x[b]*s[o[a], o[b]] for a in range(3) for b in range(a+1, 3))
            return (1+sum(x))*Q**-2.5*8*y1*y2*y3
        tot += integrate.tplquad(f, 0, 1, 0, 1, 0, 1, epsabs=0, epsrel=1e-11)[0]
    return gamma(2.5)*np.pi**1.5*tot

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    rec = []
    for k in range(4):
        P = rng.normal(size=(4, 3))
        a, b = box_feyn(P), box_closed(P)
        print(k, a, b, abs(a/b-1), flush=True)
        rec.append({"feyn": a, "closed": b, "rel": abs(a/b-1)})
    Path("debug/data/ring_ext_QR_box.json").write_text(json.dumps(rec, indent=1))
