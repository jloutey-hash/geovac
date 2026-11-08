r"""Closed-form B-table seeds -- killing the `_seed_B` mpmath quadrature.

The B-moment  B_l^{m,s}(p,c) = int_1^inf xi^p (xi^2-1)^s [d^mQ_l/dxi^m] e^{-c xi} dxi
is the universal V_ee bottleneck (95% of vee_mp, 92% of the corr recurrence),
seeded by quadrature in `neumann_vee_general_m._seed_B`.  It has a closed form.

KEY (dissolves the term-divergence obstruction).  s >= m always holds physically
(s=(mu_i+mu_j+m)/2, m in {mu_i+mu_j,|mu_i-mu_j|} => mu_i+mu_j >= m => s >= m).  In
   d^mQ_l = sum_{a=0}^m C(m,a) d^aP_l d^{m-a}Q_0  -  d^mW_{l-1},
every pole term d^kQ_0 ~ (xi-1)^{-k} (k=m-a <= m <= s) is FULLY POLYNOMIALIZED by
the intact (xi^2-1)^s = (xi-1)^s (xi+1)^s weight:
   (xi^2-1)^s d^kQ_0 = coeff_k [ (xi-1)^{s-k}(xi+1)^s - (xi-1)^s(xi+1)^{s-k} ],
both terms polynomials.  So the ONLY transcendental piece is the a=m term (a log-
moment against Q_0); everything else is monomial moments A_n(c) (closed form).

   B_l^{m,s}(p,c) = <poly_A | Q_0>_c  +  sum_n poly_rest[n] A_n(c),
   poly_A = xi^p (xi^2-1)^s d^mP_l  ( = ngm._W_poly(l,m,s,p) ),
   poly_rest = xi^p (xi^2-1)^s [ sum_{a<m} C(m,a) d^aP_l d^{m-a}Q_0(polynomialized)
                                  - d^mW_{l-1} ].

The one transcendental primitive:
   L_n(c) = int_1^inf xi^n Q_0(xi) e^{-c xi} dxi = 0.5 (Lp_n - Lm_n),
   Lm_n = int xi^n ln(xi-1) e^{-c xi} = e^{-c} sum_{j=0}^n C(n,j) j! (H_j - g - ln c)/c^{j+1},
   Lp_n = int xi^n ln(xi+1) e^{-c xi} :  Lp_0 = (ln2 e^{-c})/c + (1/c) e^c E_1(2c),
          Lp_n = (e^{-c} ln2)/c + (n/c) Lp_{n-1} + (1/c) R_n,
          R_n  = sum_{i=0}^{n-1} (-1)^i A_{n-1-i}(c) + (-1)^n e^c E_1(2c).
(g = Euler-Mascheroni, H_j = harmonic number, A_k = _mono_moments.)  E_1 enters
once (the ln(xi+1) tail); no divergence.

This driver validates L_n and the closed-form B against mpmath quadrature.
"""
from __future__ import annotations
import time
import mpmath as mp
from geovac import neumann_vee_general_m as ngm

mp.mp.dps = 40


# ---------------------------------------------------------------------------
# the log-moment primitive  L_n(c) = int_1^inf xi^n Q_0 e^{-c xi} dxi
# ---------------------------------------------------------------------------
def L_moments(n_max: int, c) -> list:
    """L_n(c), n = 0..n_max, closed form (0.5 (Lp - Lm))."""
    c = mp.mpf(c)
    ec = mp.e ** (-c)          # e^{-c}
    ecp = mp.e ** c            # e^{+c}
    ln2 = mp.log(2)
    lnc = mp.log(c)
    g = mp.euler
    E1_2c = mp.e1(2 * c)
    A = ngm._mono_moments(c, n_max + 1)     # A_k(c), k=0..n_max

    # Lm_n
    Lm = []
    for n in range(n_max + 1):
        tot = mp.mpf(0)
        H = mp.mpf(0)
        for j in range(n + 1):
            if j >= 1:
                H += mp.mpf(1) / j
            tot += mp.binomial(n, j) * mp.factorial(j) * (H - g - lnc) / c ** (j + 1)
        Lm.append(ec * tot)

    # Lp_n by recurrence
    Lp = [(ln2 * ec) / c + (ecp * E1_2c) / c]    # Lp_0
    for n in range(1, n_max + 1):
        # R_n = sum_{i=0}^{n-1} (-1)^i A_{n-1-i} + (-1)^n e^c E1(2c)
        R = mp.mpf(0)
        for i in range(n):
            R += (-1) ** i * A[n - 1 - i]
        R += (-1) ** n * ecp * E1_2c
        Lp.append((ec * ln2) / c + (n * Lp[n - 1]) / c + R / c)

    return [mp.mpf('0.5') * (Lp[n] - Lm[n]) for n in range(n_max + 1)]


def L_quad(n, c):
    """L_n(c) by direct mpmath quadrature (ground truth)."""
    c = mp.mpf(c)

    def f(u):                                   # u = xi - 1
        xi = 1 + u
        Q0 = mp.mpf('0.5') * mp.log((xi + 1) / u)
        return xi ** n * Q0 * mp.e ** (-c * xi)
    return mp.quad(f, [0, mp.mpf('0.1'), mp.mpf('1'), mp.mpf(4), mp.inf])


# ---------------------------------------------------------------------------
# polynomial helpers for the polynomialized pole terms
# ---------------------------------------------------------------------------
def _xm1_pow(j: int) -> list:
    """(xi-1)^j coeffs (low->high)."""
    out = [mp.mpf(1)]
    base = [mp.mpf(-1), mp.mpf(1)]
    for _ in range(j):
        out = ngm._polymul(out, base)
    return out


def _xp1_pow(j: int) -> list:
    """(xi+1)^j coeffs."""
    out = [mp.mpf(1)]
    base = [mp.mpf(1), mp.mpf(1)]
    for _ in range(j):
        out = ngm._polymul(out, base)
    return out


# ---------------------------------------------------------------------------
# closed-form B_l^{m,s}(p,c)
# ---------------------------------------------------------------------------
def B_closed(l: int, m: int, s: int, p: int, c, Lmom=None, A=None) -> mp.mpf:
    """B_l^{m,s}(p,c) closed form.  Lmom/A optional precomputed tables."""
    c = mp.mpf(c)
    assert s >= m, f"closed form needs s>=m (got s={s}, m={m})"

    # --- log-moment of poly_A = xi^p (xi^2-1)^s d^mP_l  (=ngm._W_poly) ---
    poly_A = list(ngm._W_poly(l, m, s, p))
    if Lmom is None:
        Lmom = L_moments(len(poly_A), c)
    log_mom = sum((poly_A[n] * Lmom[n] for n in range(len(poly_A)) if poly_A[n] != 0),
                  mp.mpf(0))

    # --- polynomial remainder ---
    xi2m1s = list(ngm._xi2m1_poly(s))
    # -d^mW_{l-1}
    if l >= 1:
        rest = [-cc for cc in ngm._polyder(list(ngm._W_lm1_poly(l)), m)]
        rest = ngm._polymul(xi2m1s, rest)            # (xi^2-1)^s * (-d^mW)
    else:
        rest = [mp.mpf(0)]
    # sum_{a<m} C(m,a) d^aP_l (xi^2-1)^s d^{m-a}Q_0
    Pl = list(ngm._leg_coeffs(l))
    for a in range(m):
        k = m - a
        daP = ngm._polyder(Pl, a) if a > 0 else Pl
        coeff_k = -((-1) ** (k - 1)) * mp.factorial(k - 1) * mp.mpf('0.5')
        # (xi^2-1)^s d^kQ_0 = coeff_k [ (xi-1)^{s-k}(xi+1)^s - (xi-1)^s(xi+1)^{s-k} ]
        term1 = ngm._polymul(_xm1_pow(s - k), _xp1_pow(s))
        term2 = ngm._polymul(_xm1_pow(s), _xp1_pow(s - k))
        polyQ = [coeff_k * (t1 - (term2[i] if i < len(term2) else 0))
                 for i, t1 in enumerate(term1)]
        if len(term2) > len(term1):
            polyQ += [coeff_k * (-term2[i]) for i in range(len(term1), len(term2))]
        contrib = ngm._polymul(daP, polyQ)
        contrib = [mp.binomial(m, a) * cc for cc in contrib]
        rest = ngm._polyadd(rest, contrib)
    rest = ngm._shift(rest, p)                       # * xi^p

    if A is None:
        A = ngm._mono_moments(c, len(rest) + 1)
    mono_mom = sum((rest[nn] * A[nn] for nn in range(len(rest)) if rest[nn] != 0),
                   mp.mpf(0))
    return log_mom + mono_mom


if __name__ == "__main__":
    print("=== (1) log-moment primitive L_n(c) closed form vs quadrature ===")
    for c in [mp.mpf(2), mp.mpf('3.7'), mp.mpf(8)]:
        Lm = L_moments(6, c)
        worst = mp.mpf(0)
        for n in range(7):
            q = L_quad(n, c)
            worst = max(worst, abs(Lm[n] - q) / abs(q))
        print(f"  c={float(c):4.1f}: worst relerr(L_n, n=0..6) = {float(worst):.1e}")

    print("\n=== (2) closed-form B_l^{m,s}(p,c) vs ngm._seed_B quadrature ===")
    alpha = 1.0
    for (m, s) in [(0, 0), (1, 1), (2, 2), (1, 2), (2, 3)]:
        c = mp.mpf(2.0 * alpha)
        worst = mp.mpf(0)
        Lmom = L_moments(40, c)
        A = ngm._mono_moments(c, 60)
        for l in (m, m + 1):
            for p in range(0, 9):
                bc = B_closed(l, m, s, p, c, Lmom, A)
                bq = ngm._seed_B(l, m, s, p, c)
                if abs(bq) > 1e-40:
                    worst = max(worst, abs(bc - bq) / abs(bq))
        print(f"  (m,s)=({m},{s}): worst relerr(B, l={m}/{m+1}, p=0..8) = {float(worst):.1e}")
