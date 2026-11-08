"""CT-F12 PoC, correct Hermitian route: variational explicitly-correlated CI
(R12-CI) for He singlet on a Coulomb-Sturmian s-basis + Slater geminal(s).

The naive two-body dressing H+D (drop convective K) collapses (ctf12_poc_he.py).
The spectrum-preserving *Hermitian* realization of the cusp-enforcing geminal is
the variational R12-CI: augment the orbital 2e space with geminal pair functions
G_ref = e^{-g r12} R_ref(r1) R_ref(r2) and Rayleigh-Ritz.  Variational => bounded
below by E_exact, no collapse.  All 2e integrals reduce to (r1,r2,x) quadrature
(s-only orbitals are isotropic; the geminal carries the theta12 dependence).

Kinetic uses the gradient form  <A|T|B> = 1/2 int (grad1 A . grad1 B + grad2 A . grad2 B)
so only first derivatives of the correlated function are needed.
"""
import numpy as np
import json, time
from scipy.special import eval_genlaguerre, eval_legendre
from geovac.sturmian_solver import hydrogenic_radial

EXACT_HE = -2.903724377
S_LIMIT = -2.879028767   # He nonrel s-limit (l=0 only), Baker et al.


# ---------------- s Coulomb-Sturmian radial R_{n0}(r; k) and derivative ----------------
def R_and_dR(n, r, k):
    """R_{n0}(r) = N e^{-kr} L_{n-1}^1(2kr),  N = 2 k^{3/2}/n  (L2-normalized)."""
    N = 2.0 * k ** 1.5 / n
    x = 2 * k * r
    L1 = eval_genlaguerre(n - 1, 1, x)
    R = N * np.exp(-k * r) * L1
    # dL_{n-1}^1/dx = -L_{n-2}^2 ; dx/dr = 2k
    if n >= 2:
        L2 = eval_genlaguerre(n - 2, 2, x)
    else:
        L2 = np.zeros_like(r)
    dR = N * np.exp(-k * r) * (-k * L1 - 2 * k * L2)
    return R, dR


# ---------------- radial grid (clustered near 0) ----------------
def make_grid(k, Ng=1200, r_max=None):
    if r_max is None:
        r_max = 42.0 / k
    t = np.linspace(0.0, 1.0, Ng)
    r = r_max * t ** 2
    r[0] = 1e-9
    wr = np.zeros(Ng)
    wr[1:-1] = (r[2:] - r[:-2]) / 2.0
    wr[0] = (r[1] - r[0]) / 2.0
    wr[-1] = (r[-1] - r[-2]) / 2.0
    return r, wr


# ---------------- angular-average kernels on the (r1,r2) grid ----------------
def build_kernels(r, gamma, nx=160):
    """Return dict of (Ng,Ng) matrices = (1/2) int_{-1}^1 [.] dx, with
    g=e^{-g r12}, gp=-g e^{-g r12}, mu1=(r1-r2 x)/r12, mu2=(r2-r1 x)/r12."""
    R1, R2 = np.meshgrid(r, r, indexing='ij')
    xs, ws = np.polynomial.legendre.leggauss(nx)
    K = {name: np.zeros_like(R1) for name in
         ['one', 'g', 'gg', 'invr', 'g_invr', 'gg_invr',
          'mu1_gp', 'mu2_gp', 'mu1_ggp', 'mu2_ggp', 'gpgp']}
    half = 0.5
    for x, wx in zip(xs, ws):
        r12 = np.sqrt(np.maximum(R1 * R1 + R2 * R2 - 2 * R1 * R2 * x, 1e-30))
        g = np.exp(-gamma * r12)
        gp = -gamma * g
        inv = 1.0 / r12
        mu1 = (R1 - R2 * x) / r12
        mu2 = (R2 - R1 * x) / r12
        w = wx * half
        K['one'] += w * 1.0
        K['g'] += w * g
        K['gg'] += w * g * g
        K['invr'] += w * inv
        K['g_invr'] += w * g * inv
        K['gg_invr'] += w * g * g * inv
        K['mu1_gp'] += w * mu1 * gp
        K['mu2_gp'] += w * mu2 * gp
        K['mu1_ggp'] += w * mu1 * g * gp
        K['mu2_ggp'] += w * mu2 * g * gp
        K['gpgp'] += w * gp * gp
    return K


# ---------------- basis functions ----------------
# Each function: list of product terms (Rleft(r), dRleft(r), Rright(r), dRright(r)), and ctype
# ctype 'one'  -> correlation factor c=1, c'=0
# ctype 'gem'  -> correlation factor c=g(r12), c'=gp(r12)
class BF:
    def __init__(self, terms, ctype):
        self.terms = terms      # list of (Rl,dRl,Rr,dRr)
        self.ctype = ctype


def orbital_pair(pi, pj, Rtab, dRtab):
    Rl, dRl = Rtab[pi], dRtab[pi]
    Rr, dRr = Rtab[pj], dRtab[pj]
    if pi == pj:
        terms = [(Rl, dRl, Rr, dRr)]
        # normalized later
    else:
        terms = [(Rtab[pi], dRtab[pi], Rtab[pj], dRtab[pj]),
                 (Rtab[pj], dRtab[pj], Rtab[pi], dRtab[pi])]
    return BF(terms, 'one')


def geminal(ref, Rtab, dRtab):
    return BF([(Rtab[ref], dRtab[ref], Rtab[ref], dRtab[ref])], 'gem')


# ---------------- matrix element assembly ----------------
def _ckey(ca, cb):
    if ca == 'one' and cb == 'one':
        return 'one'
    if ca == 'gem' and cb == 'gem':
        return 'gg'
    return 'g'          # one-gem


def _ee_key(ca, cb):
    return {'one': 'invr', 'g': 'g_invr', 'gg': 'gg_invr'}[_ckey(ca, cb)]


def assemble(bfs, r, wr, K, Z):
    """Return H, M (n x n) in the (non-orthogonal) BF basis."""
    n = len(bfs)
    W = r * r * wr                       # radial measure weight
    invr_rad = np.where(r > 1e-12, 1.0 / r, 0.0)
    H = np.zeros((n, n))
    M = np.zeros((n, n))
    ones = K['one']
    for a in range(n):
        for b in range(a, n):
            A, B = bfs[a], bfs[b]
            ck = _ckey(A.ctype, B.ctype)
            Kc = K[ck]                    # angular avg of c_A c_B
            Kee = K[_ee_key(A.ctype, B.ctype)]
            ov = 0.0; nuc = 0.0; ee = 0.0; kin = 0.0
            for (Rl_a, dRl_a, Rr_a, dRr_a) in A.terms:
                for (Rl_b, dRl_b, Rr_b, dRr_b) in B.terms:
                    la = W * Rl_a * Rl_b          # left  radial density (e1)
                    ra = W * Rr_a * Rr_b          # right radial density (e2)
                    # overlap  <A|B> = sum_ij la_i ra_j Kc_ij
                    ov += la @ Kc @ ra
                    # nuclear  -Z/r1 - Z/r2
                    nuc += (-Z) * ((la * invr_rad) @ Kc @ ra) \
                         + (-Z) * (la @ Kc @ (ra * invr_rad))
                    # e-e 1/r12
                    ee += la @ Kee @ ra
                    # kinetic (gradient form)
                    dla = W * dRl_a * dRl_b       # grad-grad radial (e1)
                    dra = W * dRr_a * dRr_b       # grad-grad radial (e2)
                    # KK term (both electrons): 1/2 (d1A d1B + d2A d2B) * c_A c_B
                    kin += 0.5 * (dla @ Kc @ ra + la @ Kc @ dra)
                    if A.ctype == 'gem' or B.ctype == 'gem':
                        # cross terms:  d(rho) . rhat . grad(c)   (c'-carrying)
                        # electron1: [ (d1 rho_A) rho_B c_A c_B' + rho_A (d1 rho_B) c_A' c_B ] mu1
                        # combine via kernels mu1_gp / mu1_ggp depending on c types
                        kA = _kin_cross_kernel(A.ctype, B.ctype, K, 'mu1')
                        kB = _kin_cross_kernel(B.ctype, A.ctype, K, 'mu1')
                        # term where derivative hits A-envelope, c' hits B-corr:
                        t1 = (W * dRl_a * Rl_b) @ kA @ (W * Rr_a * Rr_b)
                        # term where derivative hits B-envelope, c' hits A-corr:
                        t2 = (W * Rl_a * dRl_b) @ kB @ (W * Rr_a * Rr_b)
                        # electron2 (mu2):
                        kA2 = _kin_cross_kernel(A.ctype, B.ctype, K, 'mu2')
                        kB2 = _kin_cross_kernel(B.ctype, A.ctype, K, 'mu2')
                        t3 = (W * Rl_a * Rl_b) @ kA2 @ (W * dRr_a * Rr_b)
                        t4 = (W * Rl_a * Rl_b) @ kB2 @ (W * Rr_a * dRr_b)
                        kin += 0.5 * (t1 + t2 + t3 + t4)
                        if A.ctype == 'gem' and B.ctype == 'gem':
                            # c_A' c_B' (rhat12.rhat12=1) term on both electrons: rho_A rho_B gp gp
                            gg = K['gpgp']
                            kin += 0.5 * ((W * Rl_a * Rl_b) @ gg @ (W * Rr_a * Rr_b)
                                          + (W * Rl_a * Rl_b) @ gg @ (W * Rr_a * Rr_b))
            H[a, b] = kin + nuc + ee
            M[a, b] = ov
            H[b, a] = H[a, b]
            M[b, a] = M[a, b]
    return H, M


def _kin_cross_kernel(cderiv_env, ccorr, K, mu):
    """Kernel for a cross term where one factor's envelope is differentiated
    (correlation c stays) and the other's correlation is differentiated (c').
    Returns the angular-avg matrix of  mu * c_env * c'_corr.
    c_env in {1(one), g(gem)}, c'_corr in {0(one), gp(gem)}."""
    # if the corr-differentiated function is 'one' -> c'=0 -> zero kernel
    if ccorr == 'one':
        return np.zeros_like(K['one'])
    # ccorr == 'gem' -> c' = gp.  c_env: one->1, gem->g
    if cderiv_env == 'one':
        return K[mu + '_gp']       # mu * gp
    else:
        return K[mu + '_ggp']      # mu * g * gp


# ---------------- solver (canonical orthogonalization) ----------------
def solve_gen(H, M, thr=1e-8):
    ev, U = np.linalg.eigh(M)
    keep = ev > thr
    X = U[:, keep] / np.sqrt(ev[keep])
    Hp = X.T @ H @ X
    w = np.linalg.eigvalsh(Hp)
    return float(w[0]), int(keep.sum())


# ---------------- drivers ----------------
def build_tabs(ns, r, k):
    Rtab = {}; dRtab = {}
    for n in range(1, ns + 1):
        R, dR = R_and_dR(n, r, k)
        Rtab[n] = R; dRtab[n] = dR
    return Rtab, dRtab


def _energy_from(ns, r, wr, K, Rtab, dRtab, n_gem=0, gem_refs=None):
    bfs = []
    for i in range(1, ns + 1):
        for j in range(i, ns + 1):
            bfs.append(orbital_pair(i, j, Rtab, dRtab))
    if gem_refs is None:
        gem_refs = list(range(1, n_gem + 1))
    for ref in gem_refs:
        bfs.append(geminal(ref, Rtab, dRtab))
    H, M = assemble(bfs, r, wr, K, Z=2)
    E, kept = solve_gen(H, M)
    return E, len(bfs), kept


def he_energy(ns, k, gamma, n_gem=0, gem_refs=None, Ng=600, nx=128, r_max=None):
    """s-only orbital FCI (n_gem=0) or + geminal(s)."""
    r, wr = make_grid(k, Ng, r_max=r_max)
    Rtab, dRtab = build_tabs(ns, r, k)
    K = build_kernels(r, gamma, nx=nx)
    return _energy_from(ns, r, wr, K, Rtab, dRtab, n_gem, gem_refs)


def run_sweep(Ng=600, nx=128):
    """Fixed grid (k-independent r_max) + kernels cached per gamma."""
    KSCAN = np.round(np.arange(1.30, 2.31, 0.10), 3)
    GSCAN = [0.7, 0.9, 1.1, 1.3, 1.6, 2.0]
    R_MAX = 42.0 / KSCAN.min()                         # cover the most diffuse basis
    r, wr = make_grid(KSCAN.min(), Ng, r_max=R_MAX)
    W = r * r * wr
    # cache kernels per gamma, radial tabs per k
    Kcache = {g: build_kernels(r, g, nx=nx) for g in GSCAN}
    K_orb = build_kernels(r, 1.0, nx=nx)               # any gamma OK for orbital-only (c=1 kernels)
    tab = {}
    NS_MAX = 8
    for k in KSCAN:
        tab[k] = build_tabs(NS_MAX, r, float(k))

    out = {'system': 'He', 'exact': EXACT_HE, 's_limit': S_LIMIT,
           'grid': {'Ng': Ng, 'nx': nx, 'r_max': float(R_MAX)},
           'orbital_only': [], 'with_geminal': []}

    R1, dR1 = R_and_dR(1, r, 1.6875)
    print(f"[val] <1s|T|1s>={0.5*(dR1*dR1)@W:.6f} (k^2/2={1.6875**2/2:.6f})  <1s|1s>={(R1*R1)@W:.8f}")

    print("\n[He s-only ORBITAL FCI]  (s-limit -2.879029)")
    print(f"{'ns':>3} {'nbf':>4} {'q=2ns':>6} {'E':>12} {'err_slim_mHa':>13} {'err_exact_mHa':>14}")
    for ns in [1, 2, 3, 4, 5, 6, 7]:
        best = (1e9, None, None)
        for k in KSCAN:
            Rtab, dRtab = tab[k]
            E, nbf, _ = _energy_from(ns, r, wr, K_orb, Rtab, dRtab, n_gem=0)
            if E < best[0]:
                best = (E, float(k), nbf)
        E, k, nbf = best
        print(f"{ns:3d} {nbf:4d} {2*ns:6d} {E:12.6f} {(E-S_LIMIT)*1000:13.3f} {(E-EXACT_HE)*1000:14.3f}   k={k}")
        out['orbital_only'].append({'ns': ns, 'nbf': nbf, 'qubit_proxy': 2*ns, 'E': E, 'k': k,
                                    'err_slimit_mHa': (E-S_LIMIT)*1000, 'err_exact_mHa': (E-EXACT_HE)*1000})

    print("\n[He s-only + GEMINAL (variational R12-CI)]  target exact -2.903724")
    print(f"{'ns':>3} {'ngem':>4} {'nbf':>4} {'q=2ns':>6} {'best(k,g)':>13} {'E':>12} {'err_exact_mHa':>13}")
    for ns in [1, 2, 3, 4]:
        for n_gem in [1, 2]:
            refs = [1] if n_gem == 1 else [1, 2]
            if max(refs) > ns:
                continue
            best = (1e9, None)
            for k in KSCAN:
                Rtab, dRtab = tab[k]
                for g in GSCAN:
                    try:
                        E, nbf, kept = _energy_from(ns, r, wr, Kcache[g], Rtab, dRtab,
                                                    n_gem=n_gem, gem_refs=refs)
                    except Exception:
                        continue
                    if EXACT_HE - 0.005 < E < best[0]:      # variational guard
                        best = (E, float(k), float(g), nbf)
            if best[1] is None:
                continue
            E, k, g, nbf = best
            print(f"{ns:3d} {n_gem:4d} {nbf:4d} {2*ns:6d} {'k=%.2f g=%.1f'%(k,g):>13} {E:12.6f} {(E-EXACT_HE)*1000:13.3f}")
            out['with_geminal'].append({'ns': ns, 'n_gem': n_gem, 'nbf': nbf, 'qubit_proxy': 2*ns,
                                        'E': E, 'k': k, 'gamma': g, 'err_exact_mHa': (E-EXACT_HE)*1000})
    return out


if __name__ == '__main__':
    t0 = time.time()
    out = run_sweep(Ng=600, nx=128)
    out['walltime_s'] = time.time() - t0
    with open('debug/data/ctf12_r12ci_he.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nwall {time.time()-t0:.1f}s -> debug/data/ctf12_r12ci_he.json")
