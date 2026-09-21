"""Be R12-CI, Part 2: variational energy of  Psi = (1 + c F) Phi_0  by Metropolis VMC.

Phi_0 = D_up(r1,r2) * D_dn(r3,r4)  (closed-shell 1s^2 2s^2; electrons 1,2 = up / 3,4 = dn),
  D_up = 1s(r1)2s(r2) - 2s(r1)1s(r2),  D_dn similarly.
F = sum_{i<j} f(r_ij),  f(r) = 1 - e^{-r}  (bounded, cusp f'(0)=1).
Psi = (1 + c F) Phi_0.  E(c) = <Psi|H|Psi>/<Psi|Psi> by VMC (local energy, FD Laplacian).

The many-body correlation terms -- INCLUDING the genuinely-4-body ones -- are ALL contained
in |Psi|^2 and E_local automatically (VMC needs no integral decomposition).  This gives the
reference correlated Be energy for the ansatz; Part 3 validates that the reduced RI-free
form reproduces the 4-body piece with exchange.

VALIDATION GATE: at c=0, VMC(Phi_0) must reproduce the analytic reference E0 = -14.539 Ha.
"""
import numpy as np

d = np.load("debug/data/be_r12ci_ref.npz")
Z1, Z2, EREF = float(d["z1"]), float(d["z2"]), float(d["Eref"])
Zn = 4.0

# analytic orbital constants (match Part 1's numerically-orthogonalized 2s)
N1 = 2.0 * Z1 ** 1.5                                  # 1s norm
OV = N1 * 6.0 / (Z1 + Z2) ** 4                        # <1s|chi>, chi = r e^{-Z2 r}
# 2s = N2 (chi - OV*1s) ; get N2 by the same radial quadrature as Part 1
from numpy.polynomial.legendre import leggauss
_xg, _wg = leggauss(600); _rr = 12.5 * (_xg + 1.0); _wr = 12.5 * _wg
_chi = _rr * np.exp(-Z2 * _rr); _s1 = N1 * np.exp(-Z1 * _rr)
_u = _chi - OV * _s1
N2 = 1.0 / np.sqrt(np.sum(_u * _u * _rr ** 2 * _wr))


def orb_1s(r):
    return N1 * np.exp(-Z1 * r)


def orb_2s(r):
    return N2 * (r * np.exp(-Z2 * r) - OV * N1 * np.exp(-Z1 * r))


def f_gem(r):
    return 1.0 - np.exp(-r)


def radii(R):
    """R: (nw,4,3) -> per-electron |r_i| (nw,4)."""
    return np.linalg.norm(R, axis=2)


def Fsum(R):
    """F = sum_{i<j} f(r_ij) : (nw,)."""
    nw = R.shape[0]
    out = np.zeros(nw)
    for i in range(4):
        for j in range(i + 1, 4):
            out += f_gem(np.linalg.norm(R[:, i] - R[:, j], axis=1))
    return out


def Phi0(R):
    """closed-shell determinant D_up(r1,r2) D_dn(r3,r4) : (nw,)."""
    r = radii(R)
    a = orb_1s(r); b = orb_2s(r)                      # (nw,4)
    Dup = a[:, 0] * b[:, 1] - b[:, 0] * a[:, 1]
    Ddn = a[:, 2] * b[:, 3] - b[:, 2] * a[:, 3]
    return Dup * Ddn


def Psi(R, c):
    return (1.0 + c * Fsum(R)) * Phi0(R)


def potential(R):
    """V = -Zn sum 1/r_i + sum_{i<j} 1/r_ij : (nw,)."""
    r = radii(R)
    V = -Zn * np.sum(1.0 / np.maximum(r, 1e-12), axis=1)
    for i in range(4):
        for j in range(i + 1, 4):
            V += 1.0 / np.maximum(np.linalg.norm(R[:, i] - R[:, j], axis=1), 1e-12)
    return V


def local_energy(R, c, h=2e-3):
    """E_L = -1/2 sum_i lap_i Psi / Psi + V , FD Laplacian."""
    psi0 = Psi(R, c)
    lap = np.zeros(R.shape[0])
    for i in range(4):
        for k in range(3):
            Rp = R.copy(); Rp[:, i, k] += h
            Rm = R.copy(); Rm[:, i, k] -= h
            lap += (Psi(Rp, c) + Psi(Rm, c) - 2.0 * psi0) / h ** 2
    return -0.5 * lap / psi0 + potential(R)


def vmc(c, nw=3000, nsteps=2600, nequil=700, step=0.45, seed=1, measure_every=10):
    rng = np.random.default_rng(seed)
    # init walkers near the atom
    R = rng.normal(scale=0.8, size=(nw, 4, 3))
    p = Psi(R, c) ** 2
    Es = []
    for it in range(nsteps):
        Rn = R + rng.normal(scale=step, size=R.shape)
        pn = Psi(Rn, c) ** 2
        acc = rng.uniform(size=nw) < (pn / np.maximum(p, 1e-300))
        R[acc] = Rn[acc]; p[acc] = pn[acc]
        if it >= nequil and (it - nequil) % measure_every == 0:
            Es.append(local_energy(R, c))
    Es = np.concatenate(Es)
    # drop non-finite (node-straddling FD outliers) and clip extreme tails
    Es = Es[np.isfinite(Es)]
    lo, hi = np.percentile(Es, [0.5, 99.5])
    Ecl = Es[(Es >= lo) & (Es <= hi)]
    return Ecl.mean(), Ecl.std() / np.sqrt(Ecl.size), Es.size


if __name__ == "__main__":
    print("=" * 70)
    print("Be VMC :  Psi = (1 + c F) Phi_0 ,  f(r)=1-e^{-r}")
    print(f"  analytic reference E0(Phi_0) = {EREF:.5f} Ha   (validation target at c=0)")
    print("=" * 70)
    print(f"{'c':>7}{'E(c) [Ha]':>14}{'+/-':>10}{'vs E0':>11}", flush=True)
    m0, e0, n0 = vmc(0.0)
    print(f"{0.0:>7.3f}{m0:>14.5f}{e0:>10.5f}{m0-EREF:>+11.5f}   <- c=0 gate", flush=True)
    best = (m0, 0.0)
    for c in (0.05, 0.10, 0.15, 0.20, 0.30):
        m, e, n = vmc(c)
        flag = ""
        if m < best[0]:
            best = (m, c); flag = "  *best"
        print(f"{c:>7.3f}{m:>14.5f}{e:>10.5f}{m-EREF:>+11.5f}{flag}", flush=True)
    print("-" * 70)
    print(f"best: E={best[0]:.5f} Ha at c={best[1]:.3f}  "
          f"(correlation captured vs Phi_0: {best[0]-EREF:+.5f} Ha)")
    print(f"Be exact (nonrel) = -14.6674 ; variational bound: E >= exact = {best[0]>=-14.6674}")
