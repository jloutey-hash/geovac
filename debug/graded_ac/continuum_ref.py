"""Continuum reference k_X^inf(t) on the smooth round S^3 (PI decision 2026-10-29).

(R1) Uncompressed spectral sum.  Matrix elements of the FULL functions a, b and
     phi = a*b (the true product) between Dirac eigenspaces, summed with the
     exact Duhamel / eps^2 kernels over all levels n <= N_ref, plus a rigorous
     Gaussian tail bound on the omitted part.
(R2) Closed-form anchor from the local Seeley-DeWitt a_4 density on round S^3.

HARD RULE: no S^3 WINDOW value of k_A, k_Phi or c is formed for n_max <= 7.
The only window values formed are c at t = 0.1 for n_max = 10, 12, 14 (allowed).

Construction.  Ambient L^2(SU(2)) (x) C^2 with basis sqrt(2j+1) D^j_{m m'} (x) e_s,
j <= J_amb; D = 2 sigma.J_R + 3/2.  The Dirac eigenbasis is the exact
Clebsch-Gordan coupling |j m; (j' = j, 1/2) J M>: J = j + 1/2 -> lambda = +(2j+3/2)
(level n = 2j), J = j - 1/2 -> lambda = -(2j+1/2) (level n = 2j - 1).  Multiplication
by D^{1/2}_{M M'} acts by CG on (m, m') (identity on spin).  All CG coefficients
here are j (x) 1/2 ones, used in closed form (checked against sympy).

Exactness of matrix elements.  A matrix element X_{ki} of a product of two
degree-1 multipliers (and D, which is diagonal) passes only through ambient
states with j_l in {j_i +- 1/2}; if j_i <= J_amb - 1 every such state exists, so
the element equals that of the untruncated operator.  The trusted set T is all
levels n <= N_ref = 2 J_amb - 3 (j <= J_amb - 1).

Tail bound (rigorous).  Write lambda_n = n + 3/2, mult(n) = 2(n+1)(n+2) (both
signs), Delta = maximal level jump of the band-limited operator (<= 3: one
multiplier changes j by 1/2, two by <= 1, and n in {2j, 2j-1}).  For an
eigenvector e_i at level n: sum_k |A_ik|^2 <= ||A||^2 and
sum_k |{D,A}_ik|^2 <= ((2 lambda_n + Delta) ||A||)^2, and the Duhamel weight
phi(x, y) <= exp(-min(x, y)) with min |lambda| >= lambda_n - Delta.  So the omitted
part of the eps^2 coefficient is bounded by
   B(N) = sum_{n > N} mult(n) [ t ||A||^2 e^{-t lambda_n^2}
                               + t^2 (2 lambda_n + Delta)^2 ||A||^2 e^{-t (lambda_n - Delta)^2} ],
and of the kinetic functional by t^2 sum_{n>N} mult(n) ||C||^2 e^{-t(lambda_n - Delta)^2}.
Operator norms, algebraically: |D^{1/2}_{M M'}| <= 1 pointwise; [D, f] =
2 sum_a sigma_a (x) (J_R^a f), J_R^a D^{1/2}_{M M'} = sum_{M''} D^{1/2}_{M M''} (S^a)_{M'' M'},
so |J_R^a f| <= ||S^a|| (sum_{M''} |D_{M M''}|^2)^{1/2} = 1/2 and ||[D, f]|| <= 2*3*(1/2) = 3.
Hence ||A|| <= 2 ||a|| ||[D,b]|| <= 6, ||C|| = ||[D, ab]|| <= ||a|| ||[D,b]|| + ||[D,a]|| ||b|| <= 6,
and for pure gauge ||[D,b] + h.c.|| <= 6.  The infinite tail sum is bounded by its
terms up to n = N + 4000 plus a geometric remainder (terms decreasing faster than 1/2 per step there).
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

LOG: List[str] = []


def log(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


# =============================================================================
# closed-form CG for j1 (x) 1/2  (all spins as doubled integers)
# =============================================================================


def cg_half(j1x2: np.ndarray, m1x2: np.ndarray, sx2: np.ndarray, Jx2: np.ndarray) -> np.ndarray:
    """<j1 m1; 1/2 s | J, m1+s> (Condon-Shortley), vectorized, doubled-integer inputs."""
    j1 = j1x2 / 2.0
    M = (m1x2 + sx2) / 2.0
    up = Jx2 == j1x2 + 1
    dn = Jx2 == j1x2 - 1
    with np.errstate(invalid="ignore", divide="ignore"):
        den = 2 * j1 + 1
        out = np.zeros(np.broadcast(j1x2, m1x2, sx2, Jx2).shape)
        plus = sx2 > 0
        v_up = np.where(plus, np.sqrt(np.maximum(j1 + M + 0.5, 0) / den),
                        np.sqrt(np.maximum(j1 - M + 0.5, 0) / den))
        v_dn = np.where(plus, -np.sqrt(np.maximum(j1 - M + 0.5, 0) / den),
                        np.sqrt(np.maximum(j1 + M + 0.5, 0) / den))
        out = np.where(up, v_up, np.where(dn, v_dn, 0.0))
        out = np.where(np.abs(M) <= Jx2 / 2.0 + 1e-12, out, 0.0)
    return out


def check_cg_against_sympy(jmax2: int = 7) -> float:
    from sympy import Rational
    from sympy.physics.wigner import clebsch_gordan
    worst = 0.0
    for j1x2 in range(0, jmax2 + 1):
        for m1x2 in range(-j1x2, j1x2 + 1, 2):
            for sx2 in (1, -1):
                for Jx2 in (j1x2 + 1, j1x2 - 1):
                    if Jx2 < 0 or abs(m1x2 + sx2) > Jx2:
                        continue
                    ref = float(clebsch_gordan(Rational(j1x2, 2), Rational(1, 2), Rational(Jx2, 2),
                                               Rational(m1x2, 2), Rational(sx2, 2), Rational(m1x2 + sx2, 2)))
                    val = float(cg_half(np.array(j1x2), np.array(m1x2), np.array(sx2), np.array(Jx2)))
                    worst = max(worst, abs(ref - val))
    return worst


# =============================================================================
# ambient + coupled bases, multipliers
# =============================================================================


class Sphere:
    """Ambient basis up to 2j = J2max and its exact coupled Dirac eigenbasis."""

    def __init__(self, J2max: int) -> None:
        self.J2max = J2max
        # ambient: (j2, m2, mp2, s) with s in {+1, -1} (doubled spin projection)
        offs = {}
        cnt = 0
        for j2 in range(J2max + 1):
            offs[j2] = cnt
            cnt += 2 * (j2 + 1) ** 2
        self.offs = offs
        self.n_amb = cnt
        # coupled states
        lam, lev, jx2, mx2, Jx2, MJx2 = [], [], [], [], [], []
        rows, cols, vals = [], [], []
        col = 0
        for j2 in range(J2max + 1):
            for Jc2 in (j2 + 1, j2 - 1):
                if Jc2 < 0:
                    continue
                for m2 in range(-j2, j2 + 1, 2):
                    for MJ2 in range(-Jc2, Jc2 + 1, 2):
                        for sx2 in (1, -1):
                            mp2 = MJ2 - sx2
                            if abs(mp2) > j2:
                                continue
                            c = float(cg_half(np.array(j2), np.array(mp2), np.array(sx2), np.array(Jc2)))
                            if c != 0.0:
                                rows.append(self.amb_index(j2, m2, mp2, sx2))
                                cols.append(col)
                                vals.append(c)
                        if Jc2 == j2 + 1:
                            lam.append(j2 + 1.5); lev.append(j2)
                        else:
                            lam.append(-(j2 + 0.5)); lev.append(j2 - 1)
                        jx2.append(j2); mx2.append(m2); Jx2.append(Jc2); MJx2.append(MJ2)
                        col += 1
        assert col == self.n_amb
        self.V = sp.csr_matrix((vals, (rows, cols)), shape=(self.n_amb, self.n_amb))
        self.lam = np.array(lam)
        self.level = np.array(lev)
        self.jx2 = np.array(jx2)

    def amb_index(self, j2: int, m2: int, mp2: int, sx2: int) -> int:
        a = (j2 - m2) // 2
        b = (j2 - mp2) // 2
        s = 0 if sx2 > 0 else 1
        return self.offs[j2] + (a * (j2 + 1) + b) * 2 + s

    def multiplier_ambient(self, M2: int, Mp2: int) -> sp.csr_matrix:
        """Multiplication by D^{1/2}_{M, M'} (doubled M2, Mp2 in {+1,-1})."""
        rows, cols, vals = [], [], []
        for j1 in range(self.J2max + 1):
            ms = np.arange(-j1, j1 + 1, 2)
            m1, mp1 = np.meshgrid(ms, ms, indexing="ij")
            m1 = m1.ravel(); mp1 = mp1.ravel()
            for j2 in (j1 + 1, j1 - 1):
                if j2 < 0 or j2 > self.J2max:
                    continue
                m2, mp2 = m1 + M2, mp1 + Mp2
                ok = (np.abs(m2) <= j2) & (np.abs(mp2) <= j2)
                c = (np.sqrt((j1 + 1) / (j2 + 1))
                     * cg_half(np.full_like(m1, j1), m1, np.full_like(m1, M2), np.full_like(m1, j2))
                     * cg_half(np.full_like(m1, j1), mp1, np.full_like(m1, Mp2), np.full_like(m1, j2)))
                ok &= c != 0
                for s in (1, -1):
                    for k in np.nonzero(ok)[0]:
                        rows.append(self.amb_index(j2, int(m2[k]), int(mp2[k]), s))
                        cols.append(self.amb_index(j1, int(m1[k]), int(mp1[k]), s))
                        vals.append(c[k])
        return sp.csr_matrix((vals, (rows, cols)), shape=(self.n_amb, self.n_amb))

    def coupled(self, X: sp.csr_matrix) -> sp.csr_matrix:
        return (self.V.T @ X @ self.V).tocsr()

    def D_amb(self) -> sp.csr_matrix:
        return (self.V @ sp.diags(self.lam) @ self.V.T).tocsr()


def commutator_D(lam: np.ndarray, X: sp.csr_matrix) -> sp.csr_matrix:
    L = sp.diags(lam)
    return (L @ X - X @ L).tocsr()


# =============================================================================
# eps^2 functionals for diagonal D (D = diag(lam))
# =============================================================================


def _phi(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    dlt = x - y
    small = np.abs(dlt) < 1e-12
    safe = np.where(small, 1.0, dlt)
    val = np.exp(-y) * (-np.expm1(-safe)) / safe
    return np.where(small, np.exp(-0.5 * (x + y)), val)


def gauge_eps2(lam: np.ndarray, A: sp.csr_matrix, t: float, T: np.ndarray) -> Tuple[float, float]:
    """(first, second) of the eps^2 heat coefficient of D + eps A, D = diag(lam):
    first = -t sum_{i in T} (A^2)_ii e^{-t lam_i^2} (full rows of A),
    second = (t^2/2) sum_{i,j in T} (lam_i + lam_j)^2 |A_ij|^2 phi."""
    Ac = A.tocoo()
    w = np.abs(Ac.data) ** 2
    inT = T[Ac.row]
    mu = lam ** 2
    first = -t * float(np.sum(w[inT] * np.exp(-t * mu[Ac.row[inT]])))
    both = inT & T[Ac.col]
    i, j = Ac.row[both], Ac.col[both]
    second = 0.5 * t * t * float(np.sum((lam[i] + lam[j]) ** 2 * w[both] * _phi(t * mu[i], t * mu[j])))
    return first, second


def kinetic(lam: np.ndarray, C: sp.csr_matrix, t: float, T: np.ndarray) -> float:
    Cc = C.tocoo()
    both = T[Cc.row] & T[Cc.col]
    i, j = Cc.row[both], Cc.col[both]
    mu = lam ** 2
    return 0.5 * t * t * float(np.sum(np.abs(Cc.data[both]) ** 2 * _phi(t * mu[i], t * mu[j])))


def tail_bound_gauge(N: int, t: float, normA: float, Delta: int) -> float:
    tot, n = 0.0, N + 1
    prev = None
    while True:
        lmb = n + 1.5
        m = 2 * (n + 1) * (n + 2)
        e2 = math.exp(-t * (lmb - Delta) ** 2) if lmb > Delta else 1.0
        term = m * (t * normA ** 2 * math.exp(-t * lmb ** 2) + t * t * (2 * lmb + Delta) ** 2 * normA ** 2 * e2)
        tot += term
        if term == 0.0:
            return tot * 2.0 + 1e-300   # later terms underflow (Gaussian decay)
        if prev is not None and term < 0.5 * prev and term < 1e-40 * max(tot, 1e-300):
            return tot * 2.0          # geometric remainder <= the last term sum doubling
        prev = term
        n += 1
        if n > N + 100000:
            return float("inf")


def tail_bound_kinetic(N: int, t: float, normC: float, Delta: int) -> float:
    tot, n, prev = 0.0, N + 1, None
    while True:
        lmb = n + 1.5
        e2 = math.exp(-t * (lmb - Delta) ** 2) if lmb > Delta else 1.0
        term = 2 * (n + 1) * (n + 2) * t * t * normC ** 2 * e2
        tot += term
        if term == 0.0:
            return tot * 2.0 + 1e-300   # later terms underflow (Gaussian decay)
        if prev is not None and term < 0.5 * prev and term < 1e-40 * max(tot, 1e-300):
            return tot * 2.0
        prev = term
        n += 1
        if n > N + 100000:
            return float("inf")


# =============================================================================
# (R2) closed-form a_4 anchor on the round unit S^3
# =============================================================================
#
# Hermitian gammas, D = -i gamma^a e_a + spin, A = gamma^a alpha_a with alpha real:
# D + A = -i gamma^a (nabla_a + i alpha_a), F = d alpha, D_A^2 = -(nabla_A^2 + E),
# E_F = i sum_{a<b} F_ab gamma^a gamma^b, Omega_F = i F.  Gilkey a_4 F-part:
# (1/360)(180 tr E^2 + 30 tr Omega^2) = (s/360)(90 - 30)|F|^2 = (s/6)|F|^2,
# |F|^2 = sum_{a,b} F_ab^2, s = tr_S 1 = 2.  No F-free alpha^2 term (gauge invariance).
#   k_A^inf(t) = (4 pi)^{-3/2} t^{1/2} (s/6) int |F|^2 + O(t^{3/2}).
# Kinetic: (t^2/2) int_0^1 Tr(C e C^dag e) ds = (t^2/2) Tr(C C^dag e^{-tD^2}) + ...
#   = (t^2/2)(4 pi t)^{-3/2} s int |d phi|^2 + O(t^{3/2})     (C = [D, phi] = -i c(d phi)).
# For the pair: (b, a) is a column of the unitary D^{1/2}(g), i.e. (z1, z2) up to an
# isometry, so alpha = 2 Im(a db) = 2 Im(z2 dz1) (up to isometry/sign), F = 2 Im(dz1 ^ dz2)
# up to sign: in Hopf coordinates |Im(dz1 ^ dz2)|^2_{a<b} = 1 identically, so
#   int |F|^2 = 4 * 2 * vol = 8 vol = 16 pi^2.
# phi = a b = z1 z2 (up to isometry), harmonic of degree 2 (Laplacian 8), int|z1 z2|^2 = vol/6:
#   int |d phi|^2 = 8 vol / 6 = (4/3) vol = 8 pi^2 / 3.
# Hence  k_A ~ (2/3) sqrt(pi) t^{1/2},  k_Phi ~ (1/3) sqrt(pi) t^{1/2},  R(0+) = 2.
VOL = 2 * math.pi ** 2
INT_F2 = Fraction(8) * 1          # times vol
INT_DPHI2 = Fraction(4, 3)        # times vol
S_SPIN = 2


def r2_leading() -> Tuple[float, float]:
    pref = (4 * math.pi) ** (-1.5)
    kA = pref * (S_SPIN / 6) * float(INT_F2) * VOL
    kP = pref * 0.5 * S_SPIN * float(INT_DPHI2) * VOL
    return kA, kP


def r2_symbolic_checks() -> Dict[str, object]:
    """Verify the Hopf-coordinate integrals symbolically (sympy)."""
    import sympy as s
    eta, x1, x2 = s.symbols("eta xi1 xi2", real=True)
    c, si = s.cos(eta), s.sin(eta)
    # orthonormal coframe e1 = d eta, e2 = cos eta d xi1, e3 = sin eta d xi2
    # w = dz1 ^ dz2 in that coframe (components e13, e21, e23), times e^{i(xi1+xi2)}
    ph = s.exp(s.I * (x1 + x2))
    w13, w21, w23 = -s.I * si * ph, s.I * c * ph, -ph
    im2 = sum(s.simplify(s.im(s.expand_complex(z)) ** 2) for z in (w13, w21, w23))
    im2 = s.simplify(im2)
    vol_el = c * si
    vol = s.integrate(s.integrate(s.integrate(vol_el, (eta, 0, s.pi / 2)), (x1, 0, 2 * s.pi)), (x2, 0, 2 * s.pi))
    z1z2sq = (c * si) ** 2
    m = s.integrate(s.integrate(s.integrate(z1z2sq * vol_el, (eta, 0, s.pi / 2)), (x1, 0, 2 * s.pi)), (x2, 0, 2 * s.pi))
    return {"|Im(dz1^dz2)|^2 (pointwise)": im2, "vol": vol, "int|z1 z2|^2 / vol": s.simplify(m / vol)}


# =============================================================================
# driver
# =============================================================================


def build(J2max: int):
    t0 = time.time()
    S = Sphere(J2max)
    a_amb = S.multiplier_ambient(-1, +1)     # a = D^{1/2}_{-1/2,+1/2}
    b_amb = S.multiplier_ambient(+1, +1)     # b = D^{1/2}_{+1/2,+1/2}
    a = S.coupled(a_amb)
    b = S.coupled(b_amb)
    lam = S.lam
    Db = commutator_D(lam, b)
    Da = commutator_D(lam, a)
    om = (a @ Db).tocsr()
    A = (om + om.getH()).tocsr()
    Apg = (Db + Db.getH()).tocsr()
    phi = (a @ b).tocsr()
    C = commutator_D(lam, phi)
    log(f"  built 2j <= {J2max}: dim {S.n_amb}, nnz(A) = {A.nnz}, nnz(C) = {C.nnz}, {time.time() - t0:.1f} s")
    return S, a, b, A, Apg, phi, C


def max_level_jump(S: Sphere, X: sp.csr_matrix, T: np.ndarray) -> int:
    Xc = X.tocoo()
    keep = T[Xc.col] & (np.abs(Xc.data) > 1e-13)
    return int(np.max(np.abs(S.level[Xc.row[keep]] - S.level[Xc.col[keep]])))


def main() -> None:
    t_start = time.time()
    log("continuum_ref.py -- smooth-S^3 reference k_X^inf(t)  (no S^3 window value at n_max <= 7)")
    log("=" * 78)
    log(f"CG closed form vs sympy, 2j <= 7: max |diff| = {check_cg_against_sympy():.1e}")

    # ---------------- consistency with the true-window prototype (norms only) ----
    from graded_product import true_window  # noqa: E402
    S2, a2, b2, *_ = build(5)
    for n in (2, 3):
        W = true_window(n)
        Wm = np.abs(S2.lam) <= n + 0.5
        idx = np.nonzero(Wm)[0]
        aw = a2[idx][:, idx].toarray(); bw = b2[idx][:, idx].toarray()
        lw = S2.lam[idx]
        ref_a = W.multiplier(Fraction(1, 2), Fraction(-1, 2), Fraction(1, 2))
        ref_b = W.multiplier(Fraction(1, 2), Fraction(1, 2), Fraction(1, 2))
        Dw = np.diag(lw)
        om_new = aw @ (Dw @ bw - bw @ Dw)
        Dr = W.D_float()
        om_ref = ref_a @ (Dr @ ref_b - ref_b @ Dr)
        sv1 = np.sort(np.linalg.svd(om_new, compute_uv=False))
        sv2 = np.sort(np.linalg.svd(om_ref, compute_uv=False))
        log(f"  window n_max={n}: singular values of P a P [D, P b P], new builder vs graded_product: "
            f"max diff {np.abs(sv1 - sv2).max():.1e}; ||.||_F = {np.linalg.norm(om_new):.6f} "
            f"(controls.log gauge 4.797118/2 at n_max=2)")

    # ---------------- reference build ------------------------------------------
    J_amb2 = 54                        # 2 J_amb
    S, a, b, A, Apg, phi, C = build(J_amb2)
    N_ref = J_amb2 - 3                 # levels n <= N_ref have j <= J_amb - 1
    T = S.level <= N_ref
    # the trusted rows are exact; check a, b unitary-column identity on the trusted block:
    # (a^dag a + b^dag b) = 1 as multiplication operators (|a|^2 + |b|^2 = 1)
    aa = (a.getH() @ a + b.getH() @ b).tocsr()
    idx = np.nonzero(S.level <= N_ref - 2)[0][:4000]
    sub = aa[idx][:, idx].toarray()
    log(f"  check |a|^2 + |b|^2 = 1 on trusted rows: max |.-1| = {np.abs(sub - np.eye(len(idx))).max():.1e}")
    Delta = max(max_level_jump(S, A, T), max_level_jump(S, C, T), max_level_jump(S, Apg, T))
    log(f"  N_ref = {N_ref} (levels), measured max level jump Delta = {Delta} (analytic bound 3); "
        "norm bounds ||A||, ||C||, ||A_pg|| <= 6 (algebraic, see docstring)")
    Delta_b = 3
    NORM = 6.0

    def ref_values(t: float):
        f, s2 = gauge_eps2(S.lam, A, t, T)
        kA = f + s2
        kP = kinetic(S.lam, C, t, T)
        fp, sp2 = gauge_eps2(S.lam, Apg, t, T)
        bA = tail_bound_gauge(N_ref, t, NORM, Delta_b)
        bP = tail_bound_kinetic(N_ref, t, NORM, Delta_b)
        return kA, kP, fp + sp2, (f, s2, fp, sp2), bA, bP

    kA0, kP0 = r2_leading()
    log("=" * 78)
    log("(R2) closed-form a_4 anchor on the unit round S^3")
    chk = r2_symbolic_checks()
    for k, v in chk.items():
        log(f"     sympy: {k} = {v}")
    log(f"     int|F|^2 = 8 vol = 16 pi^2;  int|d phi|^2 = (4/3) vol = 8 pi^2/3;  s = 2")
    log(f"     k_A^inf  ~ (4pi)^-3/2 (s/6) int|F|^2 t^1/2   = (2/3) sqrt(pi) t^1/2 = {kA0:.12f} t^1/2")
    log(f"     k_Phi^inf ~ (4pi)^-3/2 (s/2) int|dphi|^2 t^1/2 = (1/3) sqrt(pi) t^1/2 = {kP0:.12f} t^1/2")
    log(f"     ratio R(0+) = (1/3) int|F|^2 / int|dphi|^2 = 2 exactly")
    # flat-torus check of the same a_4 coefficients (independent continuum route of graded_product)
    from graded_product import torus_continuum_gauge, torus_continuum_higgs_kin  # noqa: E402
    log("     flat T^2 check of the a_4 constants (s/6 and s/2): predicted k_A = (4/3) pi t, k_Phi = 2 pi t")
    for tt in (0.1, 0.01, 0.001):
        ga = torus_continuum_gauge((1, 0), (0, 1), tt)
        hk = torus_continuum_higgs_kin((1, 0), (0, 1), tt)
        log(f"       t={tt}: k_A/(4pi t/3) = {ga / (4 * math.pi * tt / 3):.8f}, k_Phi/(2 pi t) = {hk / (2 * math.pi * tt):.8f}")

    log("=" * 78)
    log("(R1) uncompressed spectral sums; t -> 0 approach to (R2)")
    log("   t        k_A^inf            k_Phi^inf          R = kA/kPhi     (R-2)/t     kA/(R2 lead)  kPhi/(R2 lead)  tail bnd A (rel)  tail bnd Phi (rel)")
    for t in (0.64, 0.32, 0.16, 0.08, 0.04, 0.02):
        kA, kP, kpg, parts, bA, bP = ref_values(t)
        R = kA / kP
        log(f"   {t:<7} {kA:.14f}  {kP:.14f}  {R:.10f}  {(R - 2) / t:+.6f}  "
            f"{kA / (kA0 * math.sqrt(t)):.8f}    {kP / (kP0 * math.sqrt(t)):.8f}     {bA / abs(kA):.1e}           {bP / kP:.1e}")

    log("=" * 78)
    log("Pure gauge (a = 1): A_pg = [D,b] + h.c.;  continuum k_A^inf must be 0 to the tail bound")
    log("   t       k_pg^inf (sum)      first term          second term         tail bound (abs)")
    for t in (0.64, 0.16, 0.04, 0.02):
        kA, kP, kpg, parts, bA, bP = ref_values(t)
        fp, sp2 = parts[2], parts[3]
        log(f"   {t:<6}  {kpg:+.3e}          {fp:+.10e}   {sp2:+.10e}   {bA:.1e}")

    log("=" * 78)
    log("Denominators at the scored points t = tau/(n_max + 1/2)^2 (reference only, no window value)")
    log("   tau  n_max  t            k_A^inf            k_Phi^inf          tail rel A   tail rel Phi")
    for tau in (1.0, 4.0):
        for n in (2, 3, 4, 5):
            t = tau / (n + 0.5) ** 2
            kA, kP, kpg, parts, bA, bP = ref_values(t)
            log(f"   {tau:<4} {n:<5}  {t:.8f}   {kA:.14f}  {kP:.14f}  {bA / abs(kA):.1e}      {bP / kP:.1e}")

    log("=" * 78)
    log("S^3 window c at fixed t = 0.1, n_max = 10, 12, 14 ONLY (outside the scored range)")
    t = 0.1
    kA, kP, kpg, parts, bA, bP = ref_values(t)
    log(f"   reference t=0.1: k_A^inf = {kA:.14f}, k_Phi^inf = {kP:.14f}")
    for n in (10, 12, 14):
        W = np.abs(S.lam) <= n + 0.5
        idx = np.nonzero(W)[0]
        aw = a[idx][:, idx].tocsr(); bw = b[idx][:, idx].tocsr()
        lw = S.lam[idx]
        omw = (aw @ commutator_D(lw, bw)).tocsr()
        Aw = (omw + omw.getH()).tocsr()
        Cw = commutator_D(lw, (aw @ bw).tocsr())
        Tw = np.ones(len(idx), bool)
        f, s2 = gauge_eps2(lw, Aw, t, Tw)
        kAw = f + s2
        kPw = kinetic(lw, Cw, t, Tw)
        pgw = (commutator_D(lw, bw) + commutator_D(lw, bw).getH()).tocsr()
        fp, sp2 = gauge_eps2(lw, pgw, t, Tw)
        c = 4 * (kAw / kA) / (kPw / kP)
        log(f"   n_max={n:<3} dim {len(idx):<5}  kA/kA^inf = {kAw / kA:.12f}  kPhi/kPhi^inf = {kPw / kP:.12f}  "
            f"c = {c:.12f}  |c/4-1| = {abs(c / 4 - 1):.2e}  rho_pg = {abs((fp + sp2) / kAw):.2e}")
    log("=" * 78)
    log(f"wall time {time.time() - t_start:.1f} s")
    (HERE / "continuum_ref.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
