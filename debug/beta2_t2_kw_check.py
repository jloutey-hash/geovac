"""Float64 validation of the NEW exact (k,w) factorization of the Paper 59 collinear T2.

Identity (derived this sprint):
   j0(k b) = int_0^1 cos(k b w) dw ,  b = s+t
   => T2 = (8/pi) int_0^1 int_0^1 J(s,t) ds dt
         = (8/pi) int_0^inf dk int_0^1 dw  Re[ Q(k,w)^2 ],  Q = int_0^1 e^{i k w s} P(s,k) ds
   and since P(s,k)=P(1-s,k),  Q = e^{i k w/2} R  with R real:
   ** T2 = (8/pi) int_0^inf dk int_0^1 dw  cos(k w) R(k,w)^2,
      R(k,w) = int_0^1 cos(k w (s-1/2)) P(s,k) ds **
The outer (s,t) DOUBLE integral becomes the SQUARE of a 1D integral => O(N) not O(N^2) per k,
and the (s,t) complex off-axis singularity (b = +- i(sqrt(c_s)+sqrt(c_t))) never appears.
"""
import numpy as np

def P(s, k):
    c = s*(1.0-s); D = np.sqrt(c*k*k+1.0)
    return c*np.exp(-D)*(1/D**3 + 3/D**4 + 3/D**5)

def gl(N, a, b):
    x, w = np.polynomial.legendre.leggauss(N)
    return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w

def F_of_k(k, Ns, Nw):
    xs, ws = gl(Ns, 0.0, 1.0)
    Pv = P(xs, k)
    xw, ww = gl(Nw, 0.0, 1.0)
    # R(k,w) for all w at once
    ph = np.outer(k*xw, xs-0.5)                 # (Nw, Ns)
    R = np.cos(ph) @ (ws*Pv)
    return float(np.sum(ww*np.cos(k*xw)*R*R))

def T2_kw(K, Nk_per, Ns, Nw, panels=None):
    tot = 0.0
    if panels is None:
        # panels of width ~pi/2 to resolve k-oscillation
        edges = [0.0]
        while edges[-1] < K: edges.append(min(K, edges[-1]+1.5))
    else:
        edges = panels
    for a, b in zip(edges[:-1], edges[1:]):
        xk, wk = gl(Nk_per, a, b)
        for kk, wkk in zip(xk, wk):
            nw = max(20, int(2.2*kk)+20)
            tot += wkk*F_of_k(kk, Ns, min(nw, Nw))
    return (8/np.pi)*tot

if __name__ == '__main__':
    ref = 0.3953557659017139641
    for K in (20, 40, 60, 80):
        v = T2_kw(K, 12, 300, 900)
        print(f"K={K:3d}: {v:.15f}   diff vs anchor {v-ref:+.3e}")
