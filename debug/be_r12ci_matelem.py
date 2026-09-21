"""Be R12-CI, Part 2 (robust): the 2x2 {Phi_0, F Phi_0} generalized eigenproblem by
importance-sampled MC of SMOOTH matrix-element integrands (gradient-form kinetic, analytic
gradients).  Avoids the VMC local-energy node divergence entirely.

Basis:  phi_0 = Phi_0 = D_up(r1,r2) D_dn(r3,r4)   (unnormalized; gen-eig handles the norm)
        phi_1 = F Phi_0 ,  F = sum_{i<j} f(r_ij) , f(r)=1-e^{-r}
Matrix elements (V = -Z sum 1/r_i + sum_{i<j} 1/r_ij):
   S_ij = INT phi_i phi_j
   H_ij = 1/2 INT grad phi_i . grad phi_j  +  INT phi_i phi_j V
Solve H c = E S c ; lowest E is the variational energy of the linear ansatz.

The genuinely-4-body Coulomb term lives in H_11's V_ee piece INT (F Phi_0)^2 sum 1/r_kl
(the <f_ab (1/r_kl) f_cd> with 4 distinct electrons); Part 3 validates the reduced RI-free
form against this.  Here all terms are done by the same brute MC (reference).

GATE: H_00 / S_00 must equal the analytic reference E0 = -14.539 Ha.
"""
import numpy as np

d = np.load("debug/data/be_r12ci_ref.npz")
Z1, Z2, EREF = float(d["z1"]), float(d["z2"]), float(d["Eref"])
Zn = 4.0
N1 = 2.0 * Z1 ** 1.5
OV = N1 * 6.0 / (Z1 + Z2) ** 4
from numpy.polynomial.legendre import leggauss
_xg, _wg = leggauss(600); _rr = 12.5 * (_xg + 1.0); _wr = 12.5 * _wg
_u = _rr * np.exp(-Z2 * _rr) - OV * N1 * np.exp(-Z1 * _rr)
N2 = 1.0 / np.sqrt(np.sum(_u * _u * _rr ** 2 * _wr))


# --- orbitals and radial derivatives ---------------------------------------- #
def s1_val(r):
    return N1 * np.exp(-Z1 * r)


def s1_der(r):
    return -Z1 * N1 * np.exp(-Z1 * r)


def s2_val(r):
    return N2 * (r * np.exp(-Z2 * r) - OV * N1 * np.exp(-Z1 * r))


def s2_der(r):
    return N2 * (np.exp(-Z2 * r) - Z2 * r * np.exp(-Z2 * r) + OV * N1 * Z1 * np.exp(-Z1 * r))


GEM_G = 2.0  # short-range geminal f = r e^{-g r} : cusp f'(0)=1, ->0 at large r (well-conditioned)


def f_gem(r):
    return r * np.exp(-GEM_G * r)


def f_der(r):
    return (1.0 - GEM_G * r) * np.exp(-GEM_G * r)


# --- Phi_0 and its gradient (nw,4,3) ---------------------------------------- #
def phi0_and_grad(R):
    r = np.linalg.norm(R, axis=2)                       # (nw,4)
    rhat = R / r[:, :, None]
    a = s1_val(r); b = s2_val(r); ap = s1_der(r); bp = s2_der(r)
    Dup = a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1]
    Ddn = a[:, 2] * b[:, 3] - b[:, 2] * a[:, 3]
    phi = Dup * Ddn
    g = np.zeros_like(R)
    # d Dup / d r1 , r2  (times Ddn) ; d Ddn / d r3, r4 (times Dup)
    g[:, 0] = (ap[:, 0] * b[:, 1] - bp[:, 0] * a[:, 1])[:, None] * rhat[:, 0] * Ddn[:, None]
    g[:, 1] = (a[:, 0] * bp[:, 1] - b[:, 0] * ap[:, 1])[:, None] * rhat[:, 1] * Ddn[:, None]
    g[:, 2] = (ap[:, 2] * b[:, 3] - bp[:, 2] * a[:, 3])[:, None] * rhat[:, 2] * Dup[:, None]
    g[:, 3] = (a[:, 2] * bp[:, 3] - b[:, 2] * ap[:, 3])[:, None] * rhat[:, 3] * Dup[:, None]
    return phi, g


def F_and_grad(R):
    nw = R.shape[0]
    F = np.zeros(nw); gF = np.zeros_like(R)
    for i in range(4):
        for j in range(i + 1, 4):
            dij = R[:, i] - R[:, j]; rij = np.linalg.norm(dij, axis=1)
            F += f_gem(rij)
            gv = (f_der(rij) / rij)[:, None] * dij
            gF[:, i] += gv; gF[:, j] -= gv
    return F, gF


def potential(R):
    r = np.linalg.norm(R, axis=2)
    V = -Zn * np.sum(1.0 / r, axis=1)
    for i in range(4):
        for j in range(i + 1, 4):
            V += 1.0 / np.linalg.norm(R[:, i] - R[:, j], axis=1)
    return V


# --- importance sampling: each electron ~ mixture of two Gammas (1s & 2s scales) --- #
TH1, TH2 = 1.0 / (2 * Z1), 1.0 / (2 * Z2)


def g3d(r):
    """product-density per-electron pdf in 3D at radius r (isotropic)."""
    from scipy.special import gamma as G
    pdf1 = r ** 2 * np.exp(-r / TH1) / (TH1 ** 3 * G(3))
    pdf2 = r ** 2 * np.exp(-r / TH2) / (TH2 ** 3 * G(3))
    grad = 0.5 * pdf1 + 0.5 * pdf2                       # radial pdf of r
    return grad / (4 * np.pi * r ** 2)                  # 3D density


def sample(nw, rng):
    pick = rng.uniform(size=(nw, 4)) < 0.5
    r = np.where(pick, rng.gamma(3, TH1, size=(nw, 4)), rng.gamma(3, TH2, size=(nw, 4)))
    u = rng.normal(size=(nw, 4, 3)); u /= np.linalg.norm(u, axis=2, keepdims=True)
    R = r[:, :, None] * u
    w = 1.0 / np.prod(g3d(r), axis=1)                   # 1/Pi g3d(r_i)
    return R, w


def matrix_elements(nw=2_000_000, batch=500_000, seed=1):
    rng = np.random.default_rng(seed)
    acc = np.zeros(6)                                    # S00,S01,S11, plus H pieces below
    S = np.zeros((2, 2)); H = np.zeros((2, 2)); n = 0
    while n < nw:
        m = min(batch, nw - n)
        R, w = sample(m, rng)
        phi, gph = phi0_and_grad(R)
        F, gF = F_and_grad(R)
        V = potential(R)
        p0 = phi; p1 = F * phi
        gp0 = gph
        gp1 = gF * phi[:, None, None] + F[:, None, None] * gph
        def E(x):
            return np.sum(x * w) / nw * (nw)             # accumulate; normalize by count later
        # accumulate sums (weighted)
        S[0, 0] += np.sum(w * p0 * p0); S[0, 1] += np.sum(w * p0 * p1); S[1, 1] += np.sum(w * p1 * p1)
        gdot00 = np.sum(gp0 * gp0, axis=(1, 2))
        gdot01 = np.sum(gp0 * gp1, axis=(1, 2))
        gdot11 = np.sum(gp1 * gp1, axis=(1, 2))
        H[0, 0] += np.sum(w * (0.5 * gdot00 + p0 * p0 * V))
        H[0, 1] += np.sum(w * (0.5 * gdot01 + p0 * p1 * V))
        H[1, 1] += np.sum(w * (0.5 * gdot11 + p1 * p1 * V))
        n += m
    S /= nw; H /= nw
    S[1, 0] = S[0, 1]; H[1, 0] = H[0, 1]
    return S, H


if __name__ == "__main__":
    S, H = matrix_elements()
    print("=" * 68)
    print("Be R12-CI 2x2  {Phi_0, F Phi_0}  (matrix-element MC, gradient kinetic)")
    print("=" * 68)
    print("S =", S.ravel())
    print("H =", H.ravel())
    E0_mc = H[0, 0] / S[0, 0]
    print(f"\nGATE  H_00/S_00 = {E0_mc:.5f} Ha   (analytic E0 = {EREF:.5f} ; diff {E0_mc-EREF:+.5f})")
    from scipy.linalg import eigh
    w, v = eigh(H, S)
    print(f"\n2x2 lowest eigenvalue E_R12 = {w[0]:.5f} Ha")
    print(f"   correlation captured vs Phi_0 : {w[0]-EREF:+.5f} Ha")
    print(f"   Be exact (nonrel) = -14.6674 ; variational (E>=exact): {w[0] >= -14.6674}")
