# -*- coding: utf-8 -*-
"""B1: joint S^3 x S^1_T translation-seminorm verification (2026-06-10, v3.112.0).

Rebuilds the compact-temporal product-carrier convergence content of the descoped
Paper 45 in the sound (action-seminorm) architecture: spatial side = Peter-Weyl window
j <= 1 on SU(2) (dim 14, the P38 scalar-prototype window), temporal side = Toeplitz
Fourier window K modes (WH7 Step-1 algebra), joint seminorm = sup over unit directions
in the 4-dim product Lie algebra su(2) + u(1) of the commutator norm
L(T) = sup_{|X|=1} ||[G(X), T]||  (the s->0 form of the translation seminorm,
matching Paper 38's L_n definition).

Checks:
  T0  grid self-test: Peter-Weyl orthonormality of the 14-dim window (quadrature exact)
  T1  pure-factor exactness: L_joint(1 (x) S_q) = 2 pi q / T_t; L_joint(A (x) 1) = L_s(A)
  T2  Leibniz envelope: max(L_s(A)||B||, ||A|| L_t(B)) <= L_joint(A(x)B)
                         <= L_s(A)||B|| + ||A|| L_t(B)
  T3  joint kernel condition: L_joint(F) = 0 iff F constant (band-limited panel)
  T4  additive smoothing: ||Phi_s (x) Phi_t (F) - F|| <= (gamma_s + gamma_t) Lip(F)
      with gamma_s, gamma_t the first moments of the factor Fejer-type kernels
  T5  Riemannian-limit recovery: K_t = 0 (N_t = 1) reduces joint to spatial bit-exactly
      (the load-bearing P45 falsifier, preserved in the rebuilt architecture)

Companion: debug/wh7_toeplitz_temporal_probe.py (Step 1), Paper 38 (spatial side),
Paper 45 new product-carrier Proposition. Frozen falsifier: tests/test_wh7_b1_joint.py.
"""
import json
import numpy as np
from pathlib import Path
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parent
RNG = np.random.default_rng(20260611)
JS = [0.0, 0.5, 1.0]                       # Peter-Weyl window j <= 1, dim = 1+4+9 = 14
T_TIME = 2 * np.pi                          # temporal circumference


# ---------- SU(2) representation machinery ----------------------------------
def jmats(j):
    """Hermitian angular momentum matrices, basis m = +j..-j (descending)."""
    dim = int(round(2 * j + 1))
    m = j - np.arange(dim)
    Jz = np.diag(m)
    Jp = np.zeros((dim, dim))
    for a in range(1, dim):                 # J+ |m> = sqrt((j-m)(j+m+1)) |m+1>
        mm = m[a]
        Jp[a - 1, a] = np.sqrt((j - mm) * (j + mm + 1))
    Jx = (Jp + Jp.T) / 2
    Jy = (Jp - Jp.T) / 2j
    return Jx, Jy, Jz


def wigner_D(j, al, be, ga):
    """D^j(alpha, beta, gamma) = e^{-i al Jz} e^{-i be Jy} e^{-i ga Jz}."""
    Jx, Jy, Jz = jmats(j)
    return (expm(-1j * al * Jz) @ expm(-1j * be * Jy) @ expm(-1j * ga * Jz))


def su2_grid(n_al=8, n_be=16, n_ga=16):
    """Quadrature grid exact for the j <= 1 window products. Haar-normalized weights."""
    al = np.arange(n_al) * 2 * np.pi / n_al
    ga = np.arange(n_ga) * 4 * np.pi / n_ga
    x, wx = np.polynomial.legendre.leggauss(n_be)       # x = cos(beta)
    be = np.arccos(x)
    pts, wts = [], []
    for a in al:
        for b, wb in zip(be, wx):
            for g in ga:
                pts.append((a, b, g))
                wts.append(wb / (n_al * n_ga * 2))       # /2 normalizes int sin(b) db
    return pts, np.array(wts)


PTS, WTS = su2_grid()
NG = len(PTS)

# Peter-Weyl window: psi^j_{m'm} = sqrt(2j+1) D^j_{m'm}(g); columns = grid points
PSI, LABELS = [], []
DCACHE = {j: np.array([wigner_D(j, *p) for p in PTS]) for j in JS}   # (NG, dim, dim)
for j in JS:
    dim = int(round(2 * j + 1))
    for a in range(dim):
        for b in range(dim):
            PSI.append(np.sqrt(2 * j + 1) * DCACHE[j][:, a, b])
            LABELS.append((j, a, b))
PSI = np.array(PSI)                                      # (14, NG)
NW = len(PSI)

GRAM = (PSI * WTS) @ PSI.conj().T


def compress_spatial(fvals):
    """P M_f P on the window from grid values of f."""
    return (PSI.conj() * (WTS * fvals)) @ PSI.T


def spatial_generators():
    """Left-translation generators on the window: block diag J_a^{(j)} on m' index."""
    gens = []
    for axis in range(3):
        blocks = []
        for j in JS:
            Ja = jmats(j)[axis]
            dim = int(round(2 * j + 1))
            blocks.append(np.kron(Ja, np.eye(dim)))
        G = np.zeros((NW, NW), dtype=complex)
        off = 0
        for B in blocks:
            d = B.shape[0]
            G[off:off + d, off:off + d] = B
            off += d
        gens.append(G)
    return gens


A_GENS = spatial_generators()


def rep_window(al, be, ga):
    """R(u) = blockdiag D^j(u^{-1}) (x) I  (left translation on the window)."""
    blocks = []
    for j in JS:
        Dj = wigner_D(j, al, be, ga)
        dim = int(round(2 * j + 1))
        blocks.append(np.kron(np.conj(Dj.T), np.eye(dim)))   # D^j(u)^dagger = D^j(u^-1)
    R = np.zeros((NW, NW), dtype=complex)
    off = 0
    for B in blocks:
        d = B.shape[0]
        R[off:off + d, off:off + d] = B
        off += d
    return R


def geodesic_dist(al, be, ga):
    """d(e, u) = theta in [0, 2pi] from the defining rep: cos(theta/2) = Re tr(u)/2."""
    u = wigner_D(0.5, al, be, ga)
    c = np.clip(np.real(np.trace(u)) / 2, -1.0, 1.0)
    return 2 * np.arccos(c)


# ---------- temporal (Step-1 Toeplitz) machinery -----------------------------
def shift_matrix(K, q):
    N = 2 * K + 1
    S = np.zeros((N, N), dtype=complex)
    for k in range(-K, K + 1):
        if -K <= k + q <= K:
            S[k + q + K, k + K] = 1.0
    return S


def Dt(K, T=T_TIME):
    return np.diag(2 * np.pi * np.arange(-K, K + 1) / T)


# ---------- seminorms --------------------------------------------------------
def dirs4(n_rand=60):
    base = list(np.eye(4))
    extra = RNG.normal(size=(n_rand, 4))
    extra /= np.linalg.norm(extra, axis=1)[:, None]
    return base + list(extra)


DIRS4 = dirs4()


def L_spatial(A):
    return max(np.linalg.norm(G @ A - A @ G, 2) for G in A_GENS)


def L_temporal(B, K):
    D = Dt(K)
    return float(np.linalg.norm(D @ B - B @ D, 2))


def L_joint(Tm, K):
    n_t = 2 * K + 1
    gens = [np.kron(G, np.eye(n_t)) for G in A_GENS] + [np.kron(np.eye(NW), Dt(K))]
    best = 0.0
    for u in DIRS4:
        G = sum(c * g for c, g in zip(u, gens))
        best = max(best, np.linalg.norm(G @ Tm - Tm @ G, 2))
    return float(best)


# ---------- band functions and Lipschitz constants ---------------------------
def band_vals(j, a, b, real=True):
    v = DCACHE[j][:, a, b]
    return np.real(v) * np.sqrt(2) if real else v


def lip_spatial_band(j, a, b, real=True, n_dirs=40):
    """Lip of f = (sqrt2 Re) D^j_{ab} via exact directional derivative
    d/ds f(exp(s X_n)^{-1} g) |_0 = [i (n.J) D^j(g)]_{ab}."""
    Jx, Jy, Jz = jmats(j)
    dirs = [np.eye(3)[i] for i in range(3)]
    extra = RNG.normal(size=(n_dirs, 3))
    dirs += list(extra / np.linalg.norm(extra, axis=1)[:, None])
    best = 0.0
    for n in dirs:
        nJ = n[0] * Jx + n[1] * Jy + n[2] * Jz
        dv = 1j * np.einsum('ab,gbc->gac', nJ, DCACHE[j])[:, a, b]
        vals = np.real(dv) * np.sqrt(2) if real else dv
        best = max(best, float(np.max(np.abs(vals))))
    return best


# ---------- smoothing kernels -------------------------------------------------
def spatial_fejer():
    """K(g) = c |sum_{j<=1} (2j+1) chi_j(g)|^2, Haar-normalized; returns grid values
    plus first moment gamma_s."""
    dirich = sum((2 * j + 1) * np.einsum('gaa->g', DCACHE[j]) for j in JS)
    K = np.abs(dirich) ** 2
    K /= float(np.sum(WTS * K))
    dist = np.array([geodesic_dist(*p) for p in PTS])
    gamma = float(np.sum(WTS * K * dist))
    return K, gamma


def smooth_spatial(A, Kvals):
    """Phi_s(A) = int K(g) R(g) A R(g)^dagger dg (conjugation average)."""
    out = np.zeros_like(A)
    for p, w, k in zip(PTS, WTS, Kvals):
        if w * k < 1e-12:
            continue
        R = rep_window(*p)
        out += (w * k) * (R @ A @ R.conj().T)
    return out


def temporal_fejer(K, T=T_TIME, n_grid=2000):
    """Classical Fejer kernel with K+1 terms on the circle; symbols and 1st moment."""
    sym = {q: max(0.0, 1 - abs(q) / (K + 1)) for q in range(-K, K + 1)}
    s = np.linspace(0, T, n_grid, endpoint=False)
    th = 2 * np.pi * s / T
    with np.errstate(divide='ignore', invalid='ignore'):
        F = np.where(np.abs(np.sin(th / 2)) < 1e-12, (K + 1),
                     np.sin((K + 1) * th / 2) ** 2 / ((K + 1) * np.sin(th / 2) ** 2))
    w = (T / n_grid)
    F = F / np.sum(F * w)                        # normalize: int F ds = 1
    dist = np.minimum(s, T - s)
    gamma = float(np.sum(F * dist * w))
    return sym, gamma


def run():
    out = {}

    # T0 — grid self-test
    t0 = float(np.max(np.abs(GRAM - np.eye(NW))))
    out["T0_orthonormality"] = t0

    K = 2
    n_t = 2 * K + 1

    # T1 — pure-factor exactness
    Sq = shift_matrix(K, 1)
    L_pure_t = L_joint(np.kron(np.eye(NW), Sq), K)
    f_half = compress_spatial(band_vals(0.5, 0, 0))
    L_s_alone = L_spatial(f_half)
    L_pure_s = L_joint(np.kron(f_half, np.eye(n_t)), K)
    out["T1_pure_factor"] = {
        "temporal_L": L_pure_t, "temporal_target": 2 * np.pi * 1 / T_TIME,
        "temporal_err": abs(L_pure_t - 2 * np.pi / T_TIME),
        "spatial_L_joint": L_pure_s, "spatial_L_alone": L_s_alone,
        "spatial_err": abs(L_pure_s - L_s_alone)}

    # T2 — Leibniz envelope on a panel of products
    panel = []
    for (j, a, b) in [(0.5, 0, 0), (1.0, 1, 1), (1.0, 0, 2)]:
        A = compress_spatial(band_vals(j, a, b))
        for q in (1, 2):
            B = shift_matrix(K, q) + shift_matrix(K, -q)
            Tj = np.kron(A, B)
            lo = max(L_spatial(A) * np.linalg.norm(B, 2),
                     np.linalg.norm(A, 2) * L_temporal(B, K))
            hi = (L_spatial(A) * np.linalg.norm(B, 2)
                  + np.linalg.norm(A, 2) * L_temporal(B, K))
            Lj = L_joint(Tj, K)
            panel.append({"j": j, "q": q, "lower": lo, "L": Lj, "upper": hi,
                          "ok": bool(lo - 1e-9 <= Lj <= hi + 1e-9)})
    out["T2_leibniz"] = {"all_ok": all(r["ok"] for r in panel), "panel": panel}

    # T3 — joint kernel condition
    const_L = L_joint(np.kron(compress_spatial(np.full(NG, 2.2)), np.eye(n_t)), K)
    noncst = []
    for _ in range(10):
        j, a, b = [(0.5, 0, 1), (1.0, 0, 0), (1.0, 2, 1)][int(RNG.integers(3))]
        q = int(RNG.integers(1, K + 1))
        A = compress_spatial(band_vals(j, a, b))
        F = (np.kron(A, shift_matrix(K, q) + shift_matrix(K, -q))
             + RNG.normal() * np.kron(np.eye(NW), shift_matrix(K, 1) + shift_matrix(K, -1)))
        noncst.append(L_joint(F, K))
    out["T3_kernel"] = {"const_L": const_L, "min_nonconst_L": min(noncst)}

    # T4 — additive smoothing rate on product observables
    Ks_vals, gamma_s = spatial_fejer()
    sym_t, gamma_t = temporal_fejer(K)
    rows = []
    for (j, a, b) in [(0.5, 0, 0), (1.0, 1, 1)]:
        fv = band_vals(j, a, b)
        A = compress_spatial(fv)
        lipf = lip_spatial_band(j, a, b)
        for q in (1, 2):
            B = shift_matrix(K, q) + shift_matrix(K, -q)        # h = 2cos(q theta)
            liph = 2 * 2 * np.pi * q / T_TIME                   # max|h'|
            suph = 2.0
            supf = float(np.max(np.abs(fv)))
            lip_joint = np.sqrt((lipf * suph) ** 2 + (supf * liph) ** 2)  # upper bd
            F = np.kron(A, B)
            Fs = np.kron(smooth_spatial(A, Ks_vals), sym_t[q] * B)
            defect = float(np.linalg.norm(Fs - F, 2))
            bound = (gamma_s + gamma_t) * lip_joint
            rows.append({"j": j, "q": q, "defect": defect, "bound": bound,
                         "ok": bool(defect <= bound)})
    out["T4_smoothing"] = {"gamma_s": gamma_s, "gamma_t": gamma_t,
                           "all_ok": all(r["ok"] for r in rows), "rows": rows}

    # T5 — Riemannian-limit recovery at N_t = 1 (K = 0)
    errs = []
    for (j, a, b) in [(0.5, 0, 0), (1.0, 1, 1), (1.0, 0, 2)]:
        A = compress_spatial(band_vals(j, a, b))
        errs.append(abs(L_joint(np.kron(A, np.eye(1)), 0) - L_spatial(A)))
    out["T5_riemannian_limit"] = {"max_err": max(errs)}

    return out


if __name__ == "__main__":
    res = run()
    (ROOT / "data" / "wh7_b1_joint_product.json").write_text(
        json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(f"T0 PW orthonormality residual : {res['T0_orthonormality']:.2e}")
    t1 = res["T1_pure_factor"]
    print(f"T1 pure factors               : temporal err {t1['temporal_err']:.2e}, "
          f"spatial err {t1['spatial_err']:.2e}")
    print(f"T2 Leibniz envelope           : all_ok = {res['T2_leibniz']['all_ok']}")
    t3 = res["T3_kernel"]
    print(f"T3 kernel condition           : const {t3['const_L']:.2e}, "
          f"min nonconst {t3['min_nonconst_L']:.3f}")
    t4 = res["T4_smoothing"]
    print(f"T4 additive smoothing         : gamma_s={t4['gamma_s']:.3f}, "
          f"gamma_t={t4['gamma_t']:.3f}, all_ok = {t4['all_ok']}")
    print(f"T5 Riemannian limit (N_t=1)   : max err {res['T5_riemannian_limit']['max_err']:.2e}")
