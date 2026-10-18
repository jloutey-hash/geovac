"""QR step 1: the K_v kernel K(x, y) = int rho(z) / (|x - z| |z - y|) d^3z for the fingerprint density
rho(z) = a^3 K0(a|z|) / (2 pi^2)  (transform (1 + k^2/a^2)^(-3/2), unit charge).
1-D form (Fourier + Feynman parameter u, then the q-integral in closed form):
   K(x,y) = (2/pi) int_0^{pi/2} dth  a Im F1(a (b - i X)) / X,   u = sin^2 th,
   X = |(1-u) x + u y|,  b = |y - x| sqrt(u (1-u)),
   F1(z) = int_0^inf e^{-zq} (1+q^2)^{-3/2} dq = -z + (pi z/2)(H_1(z) - Y_1(z))  (Struve H, Bessel Y).
Check against direct 3-D quadrature in prolate spheroidal coordinates with foci x, y
(d^3z/(|z-x||z-y|) = (R/2) dxi deta dphi), five random (x, y), scales a = 1 and a = 0.6."""
import json
import numpy as np
import mpmath as mp
from pathlib import Path
from scipy import integrate
from scipy.special import k0

def K_1d(x, y, a, dps=30):
    with mp.workdps(dps):
        x = [mp.mpf(v) for v in x]; y = [mp.mpf(v) for v in y]; a = mp.mpf(a)
        R = mp.sqrt(sum((p - q)**2 for p, q in zip(x, y)))
        def f(th):
            u = mp.sin(th)**2
            X = mp.sqrt(sum(((1 - u)*p + u*q)**2 for p, q in zip(x, y)))
            b = R*mp.sqrt(u*(1 - u))
            z = a*(b - 1j*X)
            F1 = -z + mp.pi*z/2*(mp.struveh(1, z) - mp.bessely(1, z))
            return a*mp.im(F1)/X
        return 2/mp.pi*mp.quad(f, [0, mp.pi/4, mp.pi/2])

def K_3d(x, y, a):
    x = np.array(x, float); y = np.array(y, float)
    R = np.linalg.norm(y - x); e = (y - x)/R; m = (x + y)/2
    mpar = m @ e; mperp = np.linalg.norm(m - mpar*e)
    rho = lambda r: a**3*k0(a*r)/(2*np.pi**2)
    def zabs(xi, eta, ph):
        s = R/2*xi*eta; rp = R/2*np.sqrt(max((xi*xi - 1)*(1 - eta*eta), 0.0))
        z2 = (mpar + s)**2 + mperp**2 + rp**2 - 2*rp*mperp*np.cos(ph)
        return np.sqrt(max(z2, 1e-300))
    # zero of |z|: phi = 0 (with the minus sign above), on the axis point where m + s e + rp e1 = 0
    nx, ny = np.linalg.norm(x), np.linalg.norm(y)
    xi0, eta0 = (nx + ny)/R, (nx - ny)/R
    def f_eta(eta, xi):
        return 2*integrate.quad(lambda ph: rho(zabs(xi, eta, ph)), 0, np.pi, epsabs=0, epsrel=1e-12, limit=200)[0]
    def f_xi(xi):
        pts = [eta0] if -1 < eta0 < 1 else None
        return integrate.quad(f_eta, -1, 1, args=(xi,), points=pts, epsabs=0, epsrel=1e-12, limit=200)[0]
    I = integrate.quad(f_xi, 1, xi0, epsabs=0, epsrel=1e-12, limit=200)[0] + \
        integrate.quad(f_xi, xi0, xi0 + 60/(a*R), epsabs=0, epsrel=1e-12, limit=200)[0] + \
        integrate.quad(f_xi, xi0 + 60/(a*R), np.inf, epsabs=0, epsrel=1e-12, limit=200)[0]
    return R/2*I

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    rec = []
    for k in range(5):
        x = rng.normal(size=3)*0.9; y = rng.normal(size=3)*1.3
        a = 1.0 if k < 3 else 0.6
        A = K_1d(x, y, a); B = K_3d(x, y, a)
        rel = abs(float(A)/B - 1)
        print(k, "a", a, "|x|", np.linalg.norm(x), "|y|", np.linalg.norm(y), "K_1d", mp.nstr(A, 20), "K_3d", B, "rel", rel, flush=True)
        rec.append({"x": list(x), "y": list(y), "a": a, "K_1d": mp.nstr(A, 25), "K_3d": B, "rel": rel})
    Path("debug/data/ring_ext_QR_kernel.json").write_text(json.dumps(rec, indent=1))
