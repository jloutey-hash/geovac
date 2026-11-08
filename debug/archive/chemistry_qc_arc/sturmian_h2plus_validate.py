"""
Paper 60 -- VALIDATION for the molecular resource estimate.

(1) Solve the ONE-ELECTRON molecular isoenergetic Shibuya-Wulfman secular equation
    [ W(kR) - k S(kR) ] C = 0,  E = -k^2/2
    for H2+ (two protons, Z=1, separation R along z), s-orbital Coulomb-Sturmian basis,
    by the isoenergetic self-consistency scan (build at scale k, find eigenvalue k'=k).
    Compare the sigma_g ground root to the known H2+ electronic energy at R=2 bohr
    (E_elec = -1.1026 Ha  =>  k = 1.485).  This grounds the resource estimate in a
    working algorithm and pins lambda_eff = k_max (subnormalization of the whitened matrix).

(2) Robustness of the Gaussian metric conditioning across even-tempered ratios
    (the resource head-to-head must not hinge on one ratio).

W and S are built on the SAME 2D cylindrical grid used in sturmian_sw_conditioning.py.
Diagnostic only.
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from numpy.linalg import cond
from scipy.linalg import eigh
from scipy.special import genlaguerre

# ----- 2D cylindrical grid -----------------------------------------------------------
rho = np.linspace(1e-4, 40.0, 900)
z = np.linspace(-30.0, 40.0, 1500)
RHO, ZZ = np.meshgrid(rho, z, indexing="ij")
drho = rho[1] - rho[0]
dz = z[1] - z[0]


def integ(g):
    return 2 * np.pi * np.sum(g * RHO) * drho * dz


def cs_at(n, zc, kk):
    """CS s-orbital n at scale kk centered at z=zc. L2-normalized."""
    rr = np.sqrt(RHO ** 2 + (ZZ - zc) ** 2)
    f = np.exp(-kk * rr) * genlaguerre(n - 1, 1)(2 * kk * rr)
    nrm = np.sqrt(integ(f * f))
    return f / nrm


def build_WS(nmax, R, kk):
    """Wulfman potential W and Shibuya-Wulfman metric S at basis scale kk, separation R.
    Centers at z=0 and z=R.  v = -1/r_A - 1/r_B (Z=1).  W = (1/kk) * <i|(1/r_A+1/r_B)|j>."""
    zA, zB = 0.0, R
    rA = np.sqrt(RHO ** 2 + (ZZ - zA) ** 2)
    rB = np.sqrt(RHO ** 2 + (ZZ - zB) ** 2)
    invr = 1.0 / rA + 1.0 / rB
    basis = [(n, zA) for n in range(1, nmax + 1)] + [(n, zB) for n in range(1, nmax + 1)]
    fs = [cs_at(n, zc, kk) for (n, zc) in basis]
    N = len(fs)
    W = np.zeros((N, N))
    S = np.zeros((N, N))
    grads = [np.gradient(f, rho, z) for f in fs]
    for i in range(N):
        for j in range(N):
            W[i, j] = (1.0 / kk) * integ(fs[i] * invr * fs[j])
            gir, giz = grads[i]
            gjr, gjz = grads[j]
            S[i, j] = (1.0 / (2 * kk ** 2)) * integ(gir * gjr + giz * gjz) + 0.5 * integ(fs[i] * fs[j])
    return W, S


def iso_roots(nmax, R, kk):
    """Generalized eigenvalues k' of [W - k' S]C=0 built at scale kk (all real, S>0)."""
    W, S = build_WS(nmax, R, kk)
    S = 0.5 * (S + S.T)
    W = 0.5 * (W + W.T)
    ev = eigh(W, S, eigvals_only=True)
    return np.sort(ev)[::-1], cond(S)   # descending; deepest bound = largest k'


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)
    R = 2.0
    nmax = 4

    print("=" * 78)
    print(f"(1) H2+ isoenergetic self-consistency scan  (R={R} bohr, nmax={nmax}, N={2*nmax})")
    print("=" * 78)
    print("  build basis at scale k -> solve [W-k'S]C=0 -> physical roots where k'(k)=k")
    print(f"  {'k(basis)':>9} | {'largest k-roots k_i(k)':>40} | {'cond(S)':>8}")
    print("  " + "-" * 66)
    ks = np.linspace(0.8, 2.2, 15)
    track = []
    kmax_seen = 0.0
    for kk in ks:
        roots, cS = iso_roots(nmax, R, kk)
        top = roots[:4]
        kmax_seen = max(kmax_seen, roots.max())
        track.append((kk, roots))
        print(f"  {kk:>9.3f} | {np.array2string(top, precision=3):>40} | {cS:>8.2f}")

    # self-consistent sigma_g: find k where the largest root k'_0(k) == k
    def largest_root(kk):
        return iso_roots(nmax, R, kk)[0][0]
    # bisection on g(k)=largest_root(k)-k
    lo, hi = 1.2, 1.8
    glo = largest_root(lo) - lo
    ghi = largest_root(hi) - hi
    root_k = None
    if glo * ghi < 0:
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            gm = largest_root(mid) - mid
            if glo * gm <= 0:
                hi = mid
                ghi = gm
            else:
                lo = mid
                glo = gm
        root_k = 0.5 * (lo + hi)
    print()
    if root_k:
        E = -0.5 * root_k ** 2
        print(f"  self-consistent sigma_g:  k = {root_k:.4f}  ->  E_elec = -k^2/2 = {E:.4f} Ha")
        print(f"  reference H2+ at R=2:      k = 1.4850      E_elec = -1.1026 Ha "
              f"(total -0.6026 + 1/R)")
        print(f"  s-only basis error: {abs(E-(-1.1026)):.4f} Ha ({100*abs(E+1.1026)/1.1026:.1f}%) "
              f"-- s-only underbinds (no p-polarization), expected")
    else:
        print("  (no sign change bracketed in [1.2,1.8]; see scan above)")
    print(f"\n  lambda_eff = k_max over all roots/scales in scan = {kmax_seen:.3f}")
    print("  (subnormalization of the whitened matrix M~ = S^-1/2 W S^-1/2; sets Q_QPE)")

    print("\n" + "=" * 78)
    print("(2) Gaussian metric conditioning robustness across even-tempered ratios (R=2)")
    print("=" * 78)

    def gauss_overlap(a, b, dAB):
        return (4 * a * b / (a + b) ** 2) ** 0.75 * np.exp(-a * b / (a + b) * dAB ** 2)

    def gaussian_metric(nm, RR, a0, ratio):
        exps = a0 * ratio ** np.arange(nm)
        centers = [0.0] * nm + [RR] * nm
        allexp = list(exps) + list(exps)
        NN = 2 * nm
        Sg = np.zeros((NN, NN))
        for i in range(NN):
            for j in range(NN):
                Sg[i, j] = gauss_overlap(allexp[i], allexp[j], abs(centers[i] - centers[j]))
        return Sg

    Ns = np.array([2 * nm for nm in range(2, 11)], float)
    print(f"  {'ratio':>6} | {'cond growth fit':>26} | cond @ N=8 , N=16")
    print("  " + "-" * 60)
    for ratio in (1.6, 2.0, 2.5, 3.0):
        c = np.array([cond(gaussian_metric(nm, R, 0.10, ratio)) for nm in range(2, 11)])
        pp = np.polyfit(np.log(Ns), np.log(c), 1)
        pe = np.polyfit(Ns, np.log(c), 1)
        rp = 1 - np.sum((np.log(c) - np.polyval(pp, np.log(Ns))) ** 2) / np.sum((np.log(c) - np.log(c).mean()) ** 2)
        re = 1 - np.sum((np.log(c) - np.polyval(pe, Ns)) ** 2) / np.sum((np.log(c) - np.log(c).mean()) ** 2)
        form = f"N^{pp[0]:.2f} (R2={rp:.3f})" if rp >= re else f"e^({pe[0]:.3f}N) (R2={re:.3f})"
        print(f"  {ratio:>6.1f} | {form:>26} | {c[3]:.2e} , {c[7]:.2e}")
    print("\n  (SW metric for the same problem: cond ~ N^1.81, cond@N=8=25.8, N=16=91.8)")
    print("\nDONE.")
