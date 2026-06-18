"""DE-RISK PROTOTYPE (group1 Bite-A, Paper 40 rank>=2 genuine test).

Reconstruct gamma_Lambda(G) for SU(2) (calibration) and SU(3) (rank-2 test)
via numerical Weyl integration, checking whether 4/pi emerges (Reading A) over
the Weyl-formula decoy 16/pi^2 (Reading B).

Key identity (avoids 0/0 on Weyl walls): with chi_lambda = A_{lambda+rho}/A_rho
and |A_rho|^2 = |Delta|^2,
    K_Lambda |Delta|^2 = |D|^2 |Delta|^2 / Z = |sum sqrt(dim) A_{lambda+rho}|^2 / Z,
so integrate S(theta) = sum_{Cas<=Lam^2} sqrt(dim) A_{lambda+rho}(theta), where
A_mu(theta) = sum_{w in W} sgn(w) exp(i <w.mu, theta>) is the alternating sum
(evaluated as a permutation sum of x_i^{l_j}, x_i = e^{i theta_i}).

    gamma = (1/(Z |W|)) int_T |S|^2 d_G(e,theta) dmu_T,
    d_G(e,theta) = sqrt(h^vee) |theta_principal|.
"""
import numpy as np
from itertools import permutations

_S3 = list(permutations(range(3)))
def _sgn(perm):
    s = 1
    for i in range(3):
        for j in range(i+1, 3):
            if perm[i] > perm[j]:
                s = -s
    return s


def principal(t):
    return (t + np.pi) % (2*np.pi) - np.pi


# ---------------------------------------------------------------------------
# SU(2) calibration (rank 1) via the numerator identity
# ---------------------------------------------------------------------------
# Torus diag(e^{i phi}, e^{-i phi}); A_{j+rho} ~ sin((2j+1) phi) (drop common 2i),
# |Delta|^2 = 4 sin^2 phi.  |W| = 2.  d = sqrt(2)|theta|, theta=(phi,-phi),
# = 2 |phi_principal|.  Cutoff j <= j_max => n_max = 2 j_max + 1.

def su2_gamma_grid(n_max, Nphi=40000):
    phi = (np.arange(Nphi) + 0.5) * (2*np.pi) / Nphi
    js = [k/2 for k in range(0, n_max)]
    S = np.zeros_like(phi)                 # sum sqrt(2j+1) sin((2j+1) phi)
    Z = 0.0
    for j in js:
        n = 2*j + 1
        S += np.sqrt(n) * np.sin(n*phi)
        Z += n
    # |D|^2 |Delta|^2 = |2i S|^2 = 4 S^2 ; weight = |Delta|^2/|W|; K weight = 4 S^2/(Z*2)
    Kw = 4*S**2 / (Z*2.0)
    dmu = (2*np.pi)/Nphi / (2*np.pi)
    dist = 2*np.abs(principal(phi))
    norm = np.sum(Kw*dmu)
    gamma = np.sum(Kw*dist*dmu)
    return gamma, norm


def su2_calibrate():
    from geovac.central_fejer_su2 import gamma_n_via_sum_rule
    print("=== SU(2) calibration: grid gamma vs exact sum-rule ===")
    ok = True
    for n in (8, 16, 32, 64):
        g_grid, norm = su2_gamma_grid(n)
        g_exact = float(gamma_n_via_sum_rule(n))
        rel = abs(g_grid-g_exact)/g_exact
        print(f"n={n:3d}  grid={g_grid:.6f}  sumrule={g_exact:.6f}  rel={rel:.2e}  "
              f"norm={norm:.6f}  n*g/log n={n*g_grid/np.log(n):.4f}")
        if rel > 1e-3:
            ok = False
    print("  asymptote 4/pi =", 4/np.pi, " | calibration", "PASS" if ok else "FAIL")
    return ok


# ---------------------------------------------------------------------------
# SU(3) rate (rank 2)
# ---------------------------------------------------------------------------
def su3_dim(p, q):
    return (p+1)*(q+1)*(p+q+2)//2

def su3_casimir(p, q):
    l = np.array([p+q+2, q+1, 0], dtype=float); l -= l.mean()
    r = np.array([2, 1, 0], dtype=float);       r -= r.mean()
    return float(l@l - r@r)

def _alt(x, powers):
    tot = np.zeros_like(x[0], dtype=complex)
    for perm in _S3:
        term = np.ones_like(x[0], dtype=complex)
        for i in range(3):
            term = term * x[i]**powers[perm[i]]
        tot = tot + _sgn(perm)*term
    return tot

def su3_gamma_grid(Lam2, Ngrid=256):
    a = (np.arange(Ngrid) + 0.5) * (2*np.pi) / Ngrid
    A, B = np.meshgrid(a, a, indexing='ij')
    t1, t2, t3 = A, B, -A-B
    x = (np.exp(1j*t1), np.exp(1j*t2), np.exp(1j*t3))
    # collect irreps with Cas <= Lam2
    irreps = []
    pmax = int(np.ceil(np.sqrt(Lam2))) + 6
    for p in range(pmax):
        for q in range(pmax):
            if su3_casimir(p, q) <= Lam2:
                irreps.append((p, q))
    S = np.zeros_like(A, dtype=complex)        # sum sqrt(dim) A_{lambda+rho}
    Z = 0.0
    for (p, q) in irreps:
        l = [p+q+2, q+1, 0]
        S += np.sqrt(su3_dim(p, q)) * _alt(x, l)
        Z += su3_dim(p, q)
    Kw = np.abs(S)**2 / (Z*6.0)                # |D|^2|Delta|^2/(Z|W|), |W|=6
    dmu = ((2*np.pi)/Ngrid)**2 / (2*np.pi)**2
    p1, p2, p3 = principal(t1), principal(t2), principal(t3)
    dist = np.sqrt(3.0) * np.sqrt(p1**2 + p2**2 + p3**2)
    norm = np.sum(Kw*dmu)
    gamma = np.sum(Kw*dist*dmu)
    Lam = np.sqrt(Lam2)
    return gamma, norm, len(irreps), Lam, Z


def su3_run():
    print("\n=== SU(3) rate via Weyl integration ===")
    print(" Lam2 n_irr   Z       norm     gamma    Lam    Lam*g/logLam")
    rows = []
    for Lam2 in (16, 25, 36, 49, 64, 100, 144, 196, 256):
        g, norm, ni, Lam, Z = su3_gamma_grid(Lam2, Ngrid=256)
        est = Lam*g/np.log(Lam)
        rows.append((Lam, g, est, norm))
        print(f" {Lam2:4d} {ni:4d}  {Z:8.1f}  {norm:.4f}  {g:.5f}  {Lam:.3f}  {est:.4f}")
    print("  Reading A 4/pi =", round(4/np.pi, 4), " Reading B 16/pi^2 =", round(16/np.pi**2, 4))
    return rows


def _fit_c(Lams, gammas):
    """4-param subleading fit Lam*gamma = c*logLam + b + c1/Lam + c2 logLam/Lam.
    Returns c (leading log coefficient)."""
    Lams = np.asarray(Lams, float); g = np.asarray(gammas, float)
    y = Lams*g
    X = np.column_stack([np.log(Lams), np.ones_like(Lams),
                         1.0/Lams, np.log(Lams)/Lams])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef[0]


def ratio_test():
    """Convention-independent cross-group test: c(SU3)/c(SU2) -> 1 under the
    SAME dual-Coxeter normalisation Cas(adjoint)=h^vee for both groups."""
    from geovac.central_fejer_su2 import gamma_n_via_sum_rule
    print("\n=== Cross-group ratio test (Cas(ad)=h^vee for both) ===")
    # SU(2): Cas_dC(j)=j(j+1) (adjoint j=1 -> 2 = h^vee). Lam=sqrt(j_max(j_max+1)).
    su2 = []
    for n_max in range(6, 90, 2):
        jmax = (n_max-1)/2
        Lam = np.sqrt(jmax*(jmax+1))
        g = float(gamma_n_via_sum_rule(n_max))
        su2.append((Lam, g))
    # SU(3): Cas_dC = Cas_basic/2 (adjoint (1,1): Cas_basic=6 -> 3 = h^vee).
    su3 = []
    for Lam2 in (16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196):
        # retain Cas_basic <= 2*Lam2
        g, norm, ni, _, Z = su3_gamma_grid_dc(2*Lam2, Ngrid=420)
        su3.append((np.sqrt(Lam2), g))
    # fit on the asymptotic tail (drop the few smallest Lam)
    c2 = _fit_c([L for L, _ in su2 if L > 6], [g for L, g in su2 if L > 6])
    c3 = _fit_c([L for L, _ in su3 if L > 4], [g for L, g in su3 if L > 4])
    print(f"  c(SU2) = {c2:.4f}")
    print(f"  c(SU3) = {c3:.4f}")
    print(f"  ratio c(SU3)/c(SU2) = {c3/c2:.4f}   (Reading A -> 1.0)")
    print(f"  4/pi = {4/np.pi:.4f}")


def su3_gamma_grid_dc(cas_basic_cut, Ngrid=420):
    """SU(3) gamma with cutoff on the BASIC Casimir <= cas_basic_cut."""
    a = (np.arange(Ngrid) + 0.5) * (2*np.pi) / Ngrid
    A, B = np.meshgrid(a, a, indexing='ij')
    t1, t2, t3 = A, B, -A-B
    x = (np.exp(1j*t1), np.exp(1j*t2), np.exp(1j*t3))
    irreps = []
    pmax = int(np.ceil(np.sqrt(cas_basic_cut))) + 8
    for p in range(pmax):
        for q in range(pmax):
            if su3_casimir(p, q) <= cas_basic_cut:
                irreps.append((p, q))
    S = np.zeros_like(A, dtype=complex)
    Z = 0.0
    for (p, q) in irreps:
        S += np.sqrt(su3_dim(p, q)) * _alt(x, [p+q+2, q+1, 0])
        Z += su3_dim(p, q)
    Kw = np.abs(S)**2 / (Z*6.0)
    dmu = ((2*np.pi)/Ngrid)**2 / (2*np.pi)**2
    p1, p2, p3 = principal(t1), principal(t2), principal(t3)
    dist = np.sqrt(3.0) * np.sqrt(p1**2 + p2**2 + p3**2)
    norm = np.sum(Kw*dmu)
    gamma = np.sum(Kw*dist*dmu)
    return gamma, norm, len(irreps), None, Z


if __name__ == "__main__":
    su2_calibrate()
    su3_run()
    ratio_test()
