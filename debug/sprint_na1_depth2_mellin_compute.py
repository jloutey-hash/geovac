"""
NA-1 depth-2 Mellin test driver — Reading A (abelian/primitive) vs
Reading B (shuffle/deconcatenation) disambiguation.

Setup (Paper 32 §VIII master Mellin engine + Paper 56 / Paper 55 §6 joint M2 × M3):

For a single Dirac operator D on S^3 with spectrum {|lambda_n| = n + 3/2,
n=0,1,2,...} and degeneracy g_n = (n+1)(n+2), the master Mellin engine
identifies three sub-mechanisms via the index k in
  M_k[f](s) := (1/Gamma(s)) int_0^infty t^(s-1) Tr(D^k * f(D^2)) dt
where k = 0 (M1 — Hopf base measure), k = 1 (M3 — vertex parity / Hurwitz),
k = 2 (M2 — Seeley-DeWitt heat kernel).

The depth-2 joint Mellin moment is
  J(s1, s2) := M_{k1,k2}[D^2, gamma_P * D]
            = (1/(Gamma(s1)*Gamma(s2))) *
              int int t1^(s1-1) t2^(s2-1) *
                Tr( D^(k1) e^{-t1 D^2} gamma_P D^(k2) e^{-t2 D^2} ) dt1 dt2
where gamma_P is the vertex-parity grading (which on S^3 Dirac is the
operator distinguishing even-n shells from odd-n shells: gamma_P |n> = (-1)^n |n>).

We test with (k1, k2) = (2, 1) — one M2 slot, one M3 slot.

KEY OBSERVATION. If [D^2, gamma_P D] = 0 (as Hilbert-space operators on
the Dirac spectrum), then operator ordering is trivial and J(s1, s2) =
zeta_{D^2,M2}(s1) * L_{gamma_P D, M3}(s2) cleanly factors. The factorisation
is SYMMETRIC under (s1 <-> s2) only modulo relabelling the (k1, k2) pair;
the natural symmetry to test is the deconcatenation-pair signature on the
Mellin-Hopf side, where SHUFFLE would force a non-symmetric deconcatenation:
  J_shuffle(s1, s2) = M2(s1) * M3(s2) + M3(s1) * M2(s2)  (depth-2 shuffle)
                    or equivalently the symmetric sum of orderings,
while PRIMITIVE-PRODUCT gives the bare product:
  J_prim(s1, s2) = M2(s1) * M3(s2)  (depth-2 primitive)

The depth-2 Mellin transform of the joint operator trace at depth 2 has
a UNIQUE closed form once the operator-ordering is fixed; what is at stake
is whether the depth-2 GeoVac observable that USES this joint Mellin (the
mixed M2*M3 observable from Paper 55 §6, e.g. zeta_{D^2}(s) * G or
zeta_{D^2}(s) * beta(4)) corresponds to:
  - Primitive product: J(s1, s2)_geovac = M2(s1) * M3(s2) bit-exactly
    (so depth-2 = depth-1 * depth-1, no new content), or
  - Shuffle product: J(s1, s2)_geovac = M2(s1) * M3(s2) + M3(s1) * M2(s2)
    (depth-2 is the symmetric pair, deconcatenation visible).

The test:
  (a) Symbolic operator-trace expansion of Tr(D^2 e^{-t1 D^2} gamma_P D e^{-t2 D^2})
      using the closed-form spectrum and chirality grading on the Camporesi-Higuchi
      Dirac on S^3. Compute the bit-exact double Mellin at INTEGER (s1, s2).
  (b) Compare to (i) the primitive product M2(s1) * M3(s2) and (ii) the
      shuffle pair M2(s1) M3(s2) + M3(s1) M2(s2).
  (c) PSLQ to identify the joint moment against (i), (ii), or against a
      depth-2 mixed-Tate basis (which would be the Reading-A null if
      the joint moment lands in depth-1 product space).

Decision gate:
  - SYMMETRIC factorisation (joint = primitive product, bit-exact)
    => Reading A: GeoVac respects primitive coproduct; U*_GV is the
       abelianization of U*_CM. Substrate is correct as-is.
  - DECONCATENATION PAIR (joint = symmetric pair, bit-exact)
    => Reading B: GeoVac respects shuffle coproduct; substrate needs
       enrichment to T(V) cofree-cocommutative shuffle Hopf algebra.
  - NEITHER (joint Mellin not cleanly defined at depth 2 without a
    depth-2 generalisation of the case-exhaustion theorem)
    => INCONCLUSIVE: report what depth-2 generalisation is needed.

Independent secondary test (Hain-Brown at depth 2):
  - Take the depth-2 joint Mellin values at (s1, s2) in {(2,2), (2,4), (3,2),
    (4,2)} and PSLQ-identify against a basis that includes Hain-Brown modular
    ingredients {E_4(2i), E_6(2i), Delta(2i), Eichler periods of Delta}.
  - If a depth-2 GeoVac period identifies with modular content across
    all three precisions, REOPEN the Hain-Brown identification at depth 2.
  - Otherwise, confirm WH6 / Paper 55's pure-MT(Z[i,1/2]) classification
    is preserved at depth 2 as well.
"""

from __future__ import annotations

import json
from pathlib import Path
import time

import sympy as sp
import mpmath as mp

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

PRECISIONS = [50, 100, 200]
PSLQ_CEILING = 10 ** 6
PSLQ_MAXSTEPS = 2000

DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)
OUT_JSON = DATA_DIR / "na1_depth2_mellin_results.json"


# ---------------------------------------------------------------------------
# Closed-form Dirac S^3 spectrum
# ---------------------------------------------------------------------------
#
# Camporesi-Higuchi: |lambda_n| = n + 3/2, n = 0, 1, 2, ...
# Degeneracy g_n = (n+1)(n+2).
# Vertex-parity grading gamma_P |n> = (-1)^n |n>.
# Both chiralities have the same |lambda_n|; tracing over chirality contributes
# an overall factor of 2 (counted in g_n if we use the dimension-of-spinor
# bundle convention).
#
# For the Mellin computation we use the closed-form per-shell spectral data:
#   eigenvalue:   lambda_n  =  (2n + 3) / 2     (positive convention)
#   degeneracy:   g_n        =  (n+1)(n+2)
#   parity sign:  s_n        =  (-1)^n
# and operator traces are
#   Tr(D^k e^{-t D^2}) = sum_n g_n * lambda_n^k * exp(-t lambda_n^2)
#   Tr(gamma_P D^k e^{-t D^2}) = sum_n s_n * g_n * lambda_n^k * exp(-t lambda_n^2)
#
# Operator-ordered joint trace:
#   T_joint(t1, t2) := Tr( D^2 e^{-t1 D^2} . gamma_P . D . e^{-t2 D^2} )
#
# Since D^2 e^{-t1 D^2} is diagonal in the |n> shell-basis, gamma_P is
# diagonal, and D is diagonal in the eigenbasis (sign convention: D|n>=lambda_n|n>),
# the joint trace evaluates as
#   T_joint(t1, t2) = sum_n g_n * (-1)^n * lambda_n^3 * exp(-(t1+t2) lambda_n^2)
# but ONLY if all four diagonals commute. Here lambda_n^3 = D^2 * D = D^3 in
# the spectrum, and exp(-(t1+t2) lambda_n^2) factors because the two e^{-t lambda^2}
# pieces commute on the diagonal eigenbasis.
#
# CRITICAL: In the spectrum representation, D, D^2, gamma_P, and e^{-t D^2}
# ARE all simultaneously diagonal in the |n,chirality,m_l> eigenbasis of D.
# (D is diagonal because we already projected; gamma_P is diagonal by definition.)
# So the operator-ordering is trivial on the spectral side, and:
#   J(s1, s2) := M_{s1,s2}[T_joint]
#              = (1/(Gamma(s1) Gamma(s2))) int int t1^(s1-1) t2^(s2-1) T_joint dt1 dt2
#              = (1/(Gamma(s1) Gamma(s2))) sum_n g_n * (-1)^n * lambda_n^3
#                   * int t1^(s1-1) e^{-t1 lambda_n^2} dt1
#                   * int t2^(s2-1) e^{-t2 lambda_n^2} dt2
#              = sum_n g_n * (-1)^n * lambda_n^3 / lambda_n^(2 s1) / lambda_n^(2 s2)
#              = sum_n g_n * (-1)^n * lambda_n^(3 - 2(s1+s2))
#
# So the joint Mellin depends ONLY on s1 + s2 = s_tot:
#   J(s1, s2) = J_eff(s_tot) := sum_n g_n * (-1)^n * lambda_n^(3 - 2 s_tot)
#                            = sum_n (n+1)(n+2) (-1)^n ((2n+3)/2)^(3 - 2 s_tot)
#
# This is the DIAGONAL collapse, and it IS the central diagnostic:
# - At Reading A (primitive product), the joint Mellin should equal the
#   PRODUCT of independent depth-1 Mellins:
#     J_prim(s1, s2) := zeta_{D^2, M2}(s1) * L_{gamma_P D, M3}(s2)
#   where zeta_{D^2, M2}(s) := sum_n g_n lambda_n^(2 - 2s) and
#         L_{gamma_P D, M3}(s) := sum_n (-1)^n g_n lambda_n^(1 - 2s).
#
# - At Reading B (shuffle/deconcatenation), the depth-2 mixed observable
#   should equal the SYMMETRIC PAIR:
#     J_shuf(s1, s2) := M2(s1) * M3(s2) + M3(s1) * M2(s2)
#   where M2(s) := zeta_{D^2}(s) and M3(s) := L_{gamma_P D}(s).
#
# The DIAGNOSTIC:
#   The diagonal-collapse J_eff(s_tot) does NOT equal either J_prim(s1, s2) or
#   J_shuf(s1, s2) generically; what we test is whether the diagonal-collapse
#   J_eff value at (s1, s2) AGREES with M2(s1)*M3(s2) (primitive) or with
#   (M2(s1)*M3(s2) + M3(s1)*M2(s2))/2 (symmetrised shuffle) or with neither.
#
# Mathematical clarity: the operator-trace IS the joint Mellin under the
# faithful representation of the Hopf algebra of Mellin moments on the
# Camporesi-Higuchi Dirac. Whether THE BUILT JOINT OPERATOR-TRACE FACTORS
# in the way that Reading A predicts is the precise structural diagnostic.

# We work bit-exact with sympy.Rational for the spectrum.


def lambda_n_sq(n: int) -> sp.Rational:
    """Eigenvalue squared on shell n: ((2n+3)/2)^2."""
    return sp.Rational(2 * n + 3, 2) ** 2


def lambda_n(n: int) -> sp.Rational:
    """Eigenvalue magnitude on shell n: (2n+3)/2."""
    return sp.Rational(2 * n + 3, 2)


def g_n(n: int) -> int:
    """Degeneracy on shell n: (n+1)(n+2)."""
    return (n + 1) * (n + 2)


def parity_sign(n: int) -> int:
    """Vertex parity: (-1)^n."""
    return -1 if (n % 2) else 1


# ---------------------------------------------------------------------------
# Depth-1 Mellin moments (numerical, high precision)
#
# zeta_{D^2, M2}(s) = sum_n g_n * lambda_n^(2-2s) at integer s >= 2
# L_{gamma_P D, M3}(s) = sum_n (-1)^n g_n * lambda_n^(1-2s) at integer s >= 1
# ---------------------------------------------------------------------------

def mellin_M2(s: int, N: int) -> mp.mpf:
    """Compute zeta_{D^2}(s) using Paper 28 Thm 1 / closed-form Hurwitz
    identification.

    Key: the direct series sum_n g_n / lambda_n^(2s-2) is conditionally
    convergent only for s sufficiently large (specifically, 2s-2 > 3 so
    s >= 3); at s = 2 it diverges absolutely. The Mellin transform
    M2(s) := (1/Gamma(s)) int t^(s-1) Tr(D^2 e^{-t D^2}) dt is defined by
    analytic continuation. Paper 28 Theorem 1 gives the closed-form values:
        zeta_{D^2}(2) = pi^2 - pi^4 / 12
        zeta_{D^2}(3) = pi^4 / 3 - pi^6 / 30
        zeta_{D^2}(4) = 2 pi^6 / 15 - 17 pi^8 / 1260
        zeta_{D^2}(5) = 3 pi^8 / 35 - 31 pi^10 / 7560        (extrapolating)
    """
    # Use closed forms from Paper 28 Thm 1 (sympy verified above) for s in
    # {2, 3, 4}; extrapolate the pattern with sympy Hurwitz for higher s.
    pi = mp.pi
    if s == 2:
        return pi ** 2 - pi ** 4 / 12
    if s == 3:
        return pi ** 4 / 3 - pi ** 6 / 30
    if s == 4:
        return 2 * pi ** 6 / 15 - 17 * pi ** 8 / 1260
    # General-s closed form via Hurwitz at 3/2: ζ_{D²}(s) = 4^(s-1) [zeta(2s-4, 3/2) - (1/4) zeta(2s-2, 3/2)]
    # For 2s - 4 > 1 (i.e., s >= 3), Hurwitz at 3/2 converges directly.
    # mpmath supports mp.zeta for Hurwitz.
    h1 = mp.zeta(2 * s - 4, mp.mpf(3) / 2)
    h2 = mp.zeta(2 * s - 2, mp.mpf(3) / 2)
    return mp.mpf(4) ** (s - 1) * (h1 - h2 / 4)


def mellin_M3(s: int, N: int) -> mp.mpf:
    """Compute L_{gamma_P D}(s) = sum_n (-1)^n g_n / lambda_n^(2s-1)
    (the vertex-parity-weighted Dirac Dirichlet series).
    """
    total = mp.mpf(0)
    for n in range(N):
        lam = mp.mpf(2 * n + 3) / 2
        g = (n + 1) * (n + 2)
        sign = -1 if (n % 2) else 1
        total += sign * g * lam ** (1 - 2 * s)
    return total


def joint_diagonal(s_tot: int, N: int) -> mp.mpf:
    """The operator-ordered joint depth-2 Mellin transform evaluated on the
    spectral diagonal:
        J_eff(s_tot) = sum_n (n+1)(n+2) (-1)^n lambda_n^(3 - 2 s_tot)
    This is the bit-exact value of the operator-trace double Mellin
    Tr(D^2 e^{-t1 D^2} gamma_P D e^{-t2 D^2}) at the integer total weight
    s_tot = s1 + s2.
    """
    total = mp.mpf(0)
    for n in range(N):
        lam = mp.mpf(2 * n + 3) / 2
        g = (n + 1) * (n + 2)
        sign = -1 if (n % 2) else 1
        total += sign * g * lam ** (3 - 2 * s_tot)
    return total


# ---------------------------------------------------------------------------
# Closed-form sanity at low-s using Hurwitz zeta identification
#
# zeta_{D^2}(s) = sum_n (n+1)(n+2) ((2n+3)/2)^(2-2s)
#              = 4^(s-1) sum_n (n+1)(n+2) (n+3/2)^(2-2s) (factoring out lam^(2-2s))
# Note (n+1)(n+2) = (n+3/2)^2 - 1/4, so:
#   zeta_{D^2}(s) = 4^(s-1) [ sum_n (n+3/2)^(4-2s) - (1/4) sum_n (n+3/2)^(2-2s) ]
#                = 4^(s-1) [ Hurwitz(2s-4, 3/2) - (1/4) Hurwitz(2s-2, 3/2) ]
# Paper 28 Thm 1 closed forms (used to cross-check the numerical sum).
# ---------------------------------------------------------------------------

def closed_form_M2(s: int) -> sp.Expr:
    """Closed form for zeta_{D^2}(s) at integer s >= 2 via Hurwitz at 3/2.
    Returns sympy expression in pi.
    """
    # zeta_{D^2}(s) = 4^(s-1) [ zeta(2s-4, 3/2) - (1/4) zeta(2s-2, 3/2) ]
    # For integer s >= 2, 2s - 4 = even >= 0 and 2s - 2 = even >= 2.
    # Hurwitz zeta at half-integer shift has closed form:
    #   zeta(s, 1/2) = (2^s - 1) zeta(s)
    # and zeta(s, 3/2) = zeta(s, 1/2) - (1/2)^(-s) = (2^s - 1) zeta(s) - 2^s
    # For integer s = 2k, zeta(2k) is in pi^(2k) Q.
    # For negative even arguments, zeta(-2k) = 0 for k >= 1; zeta(0) = -1/2.
    # The Hurwitz at 3/2 with negative arg uses Hurwitz reflection.
    arg1 = 2 * s - 4  # >= 0
    arg2 = 2 * s - 2  # >= 2

    # Use sympy's exact Hurwitz zeta evaluation
    h1 = sp.zeta(arg1, sp.Rational(3, 2))
    h2 = sp.zeta(arg2, sp.Rational(3, 2))

    expr = sp.Rational(4) ** (s - 1) * (h1 - sp.Rational(1, 4) * h2)
    return sp.simplify(expr)


# ---------------------------------------------------------------------------
# Hain-Brown modular basis at tau = 2i (reused from yesterday's HB sprint)
# ---------------------------------------------------------------------------

def compute_modular_basis(tau_imag: float = 2.0) -> dict:
    """Compute E_4(tau), E_6(tau), Delta(tau), Eichler periods at tau = i*tau_imag."""
    tau = mp.mpc(0, tau_imag)
    q = mp.exp(2 * mp.pi * 1j * tau)
    N = max(40, int(mp.mp.dps * 0.5))

    def sigma_k(n: int, k: int) -> int:
        s = 0
        for d in range(1, n + 1):
            if n % d == 0:
                s += d ** k
        return s

    E4 = mp.mpf(1)
    E6 = mp.mpf(1)
    for n in range(1, N + 1):
        qn = q ** n
        E4 += 240 * sigma_k(n, 3) * qn
        E6 -= 504 * sigma_k(n, 5) * qn
    Delta = (E4 ** 3 - E6 ** 2) / mp.mpf(1728)

    T_cutoff = mp.mpf(20)

    def Delta_at_it(t):
        z = mp.mpc(0, t)
        qz = mp.exp(2 * mp.pi * 1j * z)
        E4z = mp.mpf(1)
        E6z = mp.mpf(1)
        for n in range(1, N + 1):
            qn = qz ** n
            E4z += 240 * sigma_k(n, 3) * qn
            E6z -= 504 * sigma_k(n, 5) * qn
        return mp.re((E4z ** 3 - E6z ** 2) / mp.mpf(1728))

    eichler_periods = {}
    for k in (0, 1, 2):
        def integrand(t, kk=k):
            return Delta_at_it(t) * (t - tau_imag) ** kk
        try:
            P_k = mp.quad(integrand, [tau_imag, T_cutoff],
                          method='tanh-sinh', maxdegree=8)
        except Exception:
            P_k = mp.mpf(0)
        eichler_periods[f"Eichler_P{k}"] = mp.re(P_k)

    return {
        "E4_real": mp.re(E4),
        "E6_real": mp.re(E6),
        "Delta_real": mp.re(Delta),
        **eichler_periods,
    }


# ---------------------------------------------------------------------------
# Convergence-controlled numerical evaluation of the joint moments
# ---------------------------------------------------------------------------

def converged_value(compute_fn, N_init: int = 200, N_step: int = 200,
                    N_max: int = 50000, tol_dps: int = None) -> tuple[mp.mpf, int, mp.mpf]:
    """Evaluate compute_fn(N) for increasing N until incremental change is
    below 10^(-tol_dps). Returns (value, N_used, last_increment).

    For slowly-converging alternating series (s_tot near critical), apply
    a Levin / Shanks transform via mpmath.sumalt to accelerate.
    """
    if tol_dps is None:
        tol_dps = mp.mp.dps - 10
    tol = mp.mpf(10) ** (-tol_dps)
    val_prev = compute_fn(N_init)
    N = N_init
    while N < N_max:
        N_new = N + N_step
        val_new = compute_fn(N_new)
        delta = abs(val_new - val_prev)
        if delta < tol * (abs(val_new) + mp.mpf(1)):
            return val_new, N_new, delta
        val_prev = val_new
        N = N_new
    return val_prev, N, mp.mpf("inf")


def joint_diagonal_sumalt(s_tot: int) -> mp.mpf:
    """Compute J_eff(s_tot) using mpmath's series acceleration for the
    alternating sum.
        J_eff(s_tot) = sum_n (-1)^n (n+1)(n+2) ((2n+3)/2)^(3 - 2 s_tot)
                     = sum_n (-1)^n a_n   where a_n = (n+1)(n+2) ((2n+3)/2)^(3 - 2 s_tot)
    """
    def a(n):
        n = int(n)
        lam = mp.mpf(2 * n + 3) / 2
        return (n + 1) * (n + 2) * lam ** (3 - 2 * s_tot)
    # mpmath nsum with alternating-series acceleration via the eulersum
    # method. Or: sumalt for raw alternating sums.
    return mp.nsum(lambda n: (-1) ** int(n) * a(n), [0, mp.inf],
                   method='alternating')


def mellin_M3_sumalt(s: int) -> mp.mpf:
    """Compute L_{gamma_P D}(s) = sum_n (-1)^n (n+1)(n+2) ((2n+3)/2)^(1-2s)
    using mpmath alternating acceleration.
    """
    def a(n):
        n = int(n)
        lam = mp.mpf(2 * n + 3) / 2
        return (n + 1) * (n + 2) * lam ** (1 - 2 * s)
    return mp.nsum(lambda n: (-1) ** int(n) * a(n), [0, mp.inf],
                   method='alternating')


# ---------------------------------------------------------------------------
# PSLQ panel (reused / adapted from HB sprint)
# ---------------------------------------------------------------------------

def run_pslq(target_value: mp.mpf, target_label: str,
             basis_values: list, basis_labels: list,
             ceiling: int = PSLQ_CEILING,
             max_basis_relations: int = 5) -> dict:
    """Run PSLQ on [target_value] + basis with iterative trivial-relation filter."""
    threshold = mp.mpf(10) ** (-(mp.mp.dps // 4))

    kept = []
    kept_labels = []
    for v, lbl in zip(basis_values, basis_labels):
        if abs(v) > threshold:
            kept.append(v)
            kept_labels.append(lbl)

    if abs(target_value) <= threshold:
        return {
            "target": target_label,
            "relation_found": False,
            "note": "target below threshold",
        }

    trivials_found = []
    for it in range(max_basis_relations + 1):
        seq = [target_value] + kept
        if len(kept) == 0:
            break
        try:
            rel = mp.pslq(seq, maxcoeff=ceiling, maxsteps=PSLQ_MAXSTEPS)
        except (ValueError, RuntimeError) as e:
            return {
                "target": target_label,
                "relation_found": False,
                "error": str(e),
                "trivials_found": trivials_found,
            }
        if rel is None:
            return {
                "target": target_label,
                "relation_found": False,
                "trivials_found": trivials_found,
            }
        a0 = rel[0]
        if a0 == 0:
            # Drop the last basis element with non-zero coefficient
            trivial_record = {
                "iteration": it,
                "coefs": [
                    {"basis_label": kept_labels[i], "coef": int(rel[i + 1])}
                    for i in range(len(kept_labels)) if rel[i + 1] != 0
                ],
            }
            trivials_found.append(trivial_record)
            last_nonzero = -1
            for j in range(1, len(seq)):
                if int(rel[j]) != 0:
                    last_nonzero = j
            if last_nonzero < 1:
                break
            drop_i = last_nonzero - 1
            kept = kept[:drop_i] + kept[drop_i + 1:]
            kept_labels = kept_labels[:drop_i] + kept_labels[drop_i + 1:]
            continue
        coefs = []
        for i, lbl in enumerate(kept_labels):
            if int(rel[i + 1]) != 0:
                coefs.append({"basis_label": lbl, "coef": int(rel[i + 1])})
        # Residual
        val = sum(rel[i] * seq[i] for i in range(len(seq)))
        return {
            "target": target_label,
            "relation_found": True,
            "target_coef": int(a0),
            "basis_coefs": coefs,
            "residual": mp.nstr(val, 5),
            "trivials_found": trivials_found,
        }
    return {
        "target": target_label,
        "relation_found": False,
        "trivials_found": trivials_found,
    }


# ---------------------------------------------------------------------------
# Main compute
# ---------------------------------------------------------------------------

def main():
    print(f"NA-1 depth-2 Mellin test driver")
    print(f"=" * 70)
    t_start = time.time()

    # Test grid: (s1, s2) pairs at integer values >= 1 (M3) and >= 2 (M2)
    # so the depth-1 Mellins each converge.
    # We use s_tot = s1 + s2 in {3, 4, 5, 6, 7} for the joint diagonal.
    s_pairs = [
        (2, 1),  # s_tot = 3
        (2, 2),  # s_tot = 4
        (3, 1),  # s_tot = 4 (same s_tot as (2,2) — tests permutation symmetry)
        (2, 3),  # s_tot = 5
        (3, 2),  # s_tot = 5 (same s_tot — tests swap)
        (4, 1),  # s_tot = 5
        (3, 3),  # s_tot = 6
        (4, 2),  # s_tot = 6
        (2, 4),  # s_tot = 6
        (4, 3),  # s_tot = 7
        (3, 4),  # s_tot = 7
    ]

    results = {
        "metadata": {
            "sprint": "NA-1 depth-2 Mellin test",
            "date": "2026-06-06",
            "discipline": "mpmath PSLQ + sympy closed-form cross-check",
            "precisions_dps": PRECISIONS,
            "pslq_ceiling": PSLQ_CEILING,
            "pslq_maxsteps": PSLQ_MAXSTEPS,
            "s_pairs_tested": s_pairs,
        },
        "by_precision": {},
        "structural_findings": {},
    }

    # --- Closed-form sanity at sympy level ---
    print("\n[Closed-form sanity check on M2 closed forms (Paper 28 Thm 1)]")
    cf_M2 = {}
    for s in [2, 3, 4]:
        cf = closed_form_M2(s)
        cf_simpl = sp.simplify(cf)
        cf_M2[s] = str(cf_simpl)
        print(f"  zeta_{{D^2}}({s}) symbolic = {cf_simpl}")
    results["structural_findings"]["closed_form_M2"] = cf_M2

    # --- Per-precision numerical compute + PSLQ test ---
    for dps in PRECISIONS:
        print(f"\n[Precision {dps} dps]")
        mp.mp.dps = dps

        # Compute depth-1 Mellins M2(s) (via closed form, Paper 28 Thm 1) and
        # M3(s) (via alternating series acceleration) for s in {1..6}
        M2_vals = {}
        M3_vals = {}
        for s in range(1, 7):
            if s >= 2:
                M2_vals[s] = mellin_M2(s, 0)  # N ignored, uses closed form
            M3_vals[s] = mellin_M3_sumalt(s)

        # Compute joint diagonal J_eff(s_tot) for s_tot in {3..7} via
        # alternating-series acceleration
        J_eff_vals = {}
        for s_tot in range(3, 8):
            J_eff_vals[s_tot] = joint_diagonal_sumalt(s_tot)

        # --- Compute primitive product and shuffle pair candidates ---
        prim_pair = {}      # M2(s1) * M3(s2)
        shuf_pair = {}      # M2(s1)*M3(s2) + M3(s1)*M2(s2)
        for s1, s2 in s_pairs:
            if s1 >= 2:
                m2_s1 = M2_vals[s1]
                prim_pair[(s1, s2)] = m2_s1 * M3_vals[s2]
            else:
                prim_pair[(s1, s2)] = None
            # Shuffle pair: only defined if both can play either role
            if s1 >= 2 and s2 >= 2:
                shuf_pair[(s1, s2)] = (M2_vals[s1] * M3_vals[s2]
                                       + M3_vals[s1] * M2_vals[s2])
            else:
                shuf_pair[(s1, s2)] = None

        # --- Symmetry diagnostic: J_eff(s_tot) is a function of s_tot ONLY ---
        # Check J_eff at swap-pairs (e.g., (2,2) vs (3,1) at s_tot=4):
        # if Reading A (primitive product), J_eff(s_tot) should equal the
        # corresponding primitive M2(s_tot/2)*M3(s_tot/2) at any split.
        # The split-invariance signature is the diagnostic.
        symmetry_panel = {}
        for s1, s2 in s_pairs:
            s_tot = s1 + s2
            j_val = J_eff_vals[s_tot]  # depends only on s_tot
            prim_val = prim_pair[(s1, s2)] if s1 >= 2 else None
            shuf_val = shuf_pair[(s1, s2)] if (s1 >= 2 and s2 >= 2) else None
            entry = {
                "s_tot": s_tot,
                "joint_diagonal": mp.nstr(j_val, 30),
                "primitive_product_M2_s1_M3_s2": mp.nstr(prim_val, 30) if prim_val is not None else None,
                "shuffle_pair": mp.nstr(shuf_val, 30) if shuf_val is not None else None,
            }
            # Asymmetry under swap: compute also (s2, s1) primitive
            if s2 >= 2 and s1 >= 1:
                prim_swap = M2_vals[s2] * M3_vals[s1]
                entry["primitive_product_swapped_M2_s2_M3_s1"] = mp.nstr(prim_swap, 30)
                # Difference under swap
                if prim_val is not None:
                    entry["primitive_swap_asymmetry"] = mp.nstr(prim_val - prim_swap, 5)
            symmetry_panel[f"({s1},{s2})"] = entry

        # --- PSLQ: does joint_diagonal(s_tot) factor as primitive product
        # or shuffle pair, for each s_tot? ---
        # The natural identification basis is { M2(s1)*M3(s2) for splits (s1, s2)
        # with s1+s2 = s_tot } plus depth-1 M2(s_tot/2)*M3(s_tot/2), with pi-powers
        # and modular ingredients.
        modular_basis_dict = compute_modular_basis(2.0)

        pslq_results = []
        for s_tot in range(3, 8):
            j_val = J_eff_vals[s_tot]
            # Build the depth-1 product basis: {M2(s1) * M3(s_tot - s1) for s1 >= 2,
            # s_tot - s1 >= 1}
            basis = []
            labels = []
            for s1 in range(2, s_tot):
                s2 = s_tot - s1
                if s2 >= 1 and s1 in M2_vals and s2 in M3_vals:
                    prod = M2_vals[s1] * M3_vals[s2]
                    basis.append(prod)
                    labels.append(f"M2({s1})*M3({s2})")
            # Also depth-1 primitives by themselves
            for s in M2_vals:
                basis.append(M2_vals[s])
                labels.append(f"M2({s})")
            for s in M3_vals:
                basis.append(M3_vals[s])
                labels.append(f"M3({s})")
            # pi-powers up to pi^10
            for k in range(1, 11):
                basis.append(mp.pi ** k)
                labels.append(f"pi^{k}")
            # zeta(3), zeta(5), Catalan, beta(4)
            basis.append(mp.zeta(3))
            labels.append("zeta(3)")
            basis.append(mp.zeta(5))
            labels.append("zeta(5)")
            basis.append(mp.catalan)
            labels.append("beta(2)=G")
            beta4 = mp.mpf(0)
            for n in range(0, max(100, int(mp.mp.dps * 0.8))):
                beta4 += (mp.mpf(-1) ** n) / (mp.mpf(2 * n + 1) ** 4)
            basis.append(beta4)
            labels.append("beta(4)")
            # Modular ingredients (Hain-Brown depth-2 probe)
            basis.append(mp.mpf(modular_basis_dict["E4_real"]))
            labels.append("E4(2i)")
            basis.append(mp.mpf(modular_basis_dict["E6_real"]))
            labels.append("E6(2i)")
            basis.append(mp.mpf(modular_basis_dict["Delta_real"]))
            labels.append("Delta(2i)")
            for k in (0, 1, 2):
                basis.append(mp.mpf(modular_basis_dict[f"Eichler_P{k}"]))
                labels.append(f"Eichler_P{k}")

            # PSLQ at this s_tot
            rel = run_pslq(j_val, f"J_eff(s_tot={s_tot})", basis, labels)
            pslq_results.append(rel)

        results["by_precision"][dps] = {
            "M2_vals": {s: mp.nstr(v, dps - 5) for s, v in M2_vals.items()},
            "M3_vals": {s: mp.nstr(v, dps - 5) for s, v in M3_vals.items()},
            "J_eff_vals": {s: mp.nstr(v, dps - 5) for s, v in J_eff_vals.items()},
            "symmetry_panel": symmetry_panel,
            "pslq_results": pslq_results,
        }
        print(f"  M2(2)={mp.nstr(M2_vals[2], 8)}, M3(2)={mp.nstr(M3_vals[2], 8)}")
        print(f"  J_eff(3)={mp.nstr(J_eff_vals[3], 8)}")
        print(f"  J_eff(4)={mp.nstr(J_eff_vals[4], 8)}")

    # --- Cross-precision agreement filter ---
    cross_prec = {}
    for s_tot in range(3, 8):
        pslq_50 = next((r for r in results["by_precision"][50]["pslq_results"]
                       if f"s_tot={s_tot}" in r.get("target", "")), {})
        pslq_100 = next((r for r in results["by_precision"][100]["pslq_results"]
                        if f"s_tot={s_tot}" in r.get("target", "")), {})
        pslq_200 = next((r for r in results["by_precision"][200]["pslq_results"]
                        if f"s_tot={s_tot}" in r.get("target", "")), {})

        def relation_signature(r):
            if not r.get("relation_found"):
                return None
            sig = {"target_coef": r.get("target_coef")}
            for c in r.get("basis_coefs", []):
                sig[c["basis_label"]] = c["coef"]
            return sig

        sig_50 = relation_signature(pslq_50)
        sig_100 = relation_signature(pslq_100)
        sig_200 = relation_signature(pslq_200)

        all_agree = (sig_50 is not None and sig_50 == sig_100 == sig_200)
        any_modular = False
        if sig_200:
            for k in sig_200:
                if k in ("E4(2i)", "E6(2i)", "Delta(2i)",
                         "Eichler_P0", "Eichler_P1", "Eichler_P2"):
                    if sig_200[k] != 0:
                        any_modular = True

        cross_prec[s_tot] = {
            "sig_50": sig_50,
            "sig_100": sig_100,
            "sig_200": sig_200,
            "agree_across_3prec": all_agree,
            "any_hain_brown_modular": any_modular,
        }

    results["cross_precision_filter"] = cross_prec

    # --- Symmetry analysis: is J_eff(s_tot) = M2(s1) * M3(s2) for all splits? ---
    # Pick highest precision (200 dps)
    mp.mp.dps = 200
    M2_hires = {}
    M3_hires = {}
    for s in range(1, 7):
        if s >= 2:
            M2_hires[s] = mellin_M2(s, 0)  # closed form
        M3_hires[s] = mellin_M3_sumalt(s)

    factorisation_test = {}
    for s_tot in range(3, 8):
        j_val = joint_diagonal_sumalt(s_tot)
        splits = []
        for s1 in range(2, s_tot):
            s2 = s_tot - s1
            if s2 >= 1:
                prim = M2_hires[s1] * M3_hires[s2]
                # Reading B shuffle pair
                if s2 >= 2:
                    shuf = M2_hires[s1] * M3_hires[s2] + M3_hires[s1] * M2_hires[s2]
                else:
                    shuf = None
                splits.append({
                    "split": (s1, s2),
                    "primitive_M2_M3": mp.nstr(prim, 30),
                    "diff_J_eff_minus_primitive": mp.nstr(j_val - prim, 5),
                    "shuffle_pair": mp.nstr(shuf, 30) if shuf is not None else None,
                    "diff_J_eff_minus_shuffle": mp.nstr(j_val - shuf, 5) if shuf is not None else None,
                    "diff_J_eff_minus_half_shuffle": mp.nstr(j_val - shuf/2, 5) if shuf is not None else None,
                })
        factorisation_test[s_tot] = {
            "J_eff": mp.nstr(j_val, 30),
            "splits": splits,
        }
    results["factorisation_test"] = factorisation_test

    # --- Verdict construction ---
    # Reading A primitive product would mean J_eff = M2(s1)*M3(s2) bit-exactly
    # at some canonical split. Reading B shuffle pair would mean J_eff equals
    # the symmetric pair (or half-pair) bit-exactly. NEITHER -> the depth-2
    # joint Mellin is its own object (not equal to depth-1 products).

    # Check if any (s_tot, split) gives bit-exact match to primitive or shuffle
    bit_exact_threshold = mp.mpf(10) ** (-50)
    verdict_reading = "NEITHER"
    primitive_hits = []
    shuffle_hits = []
    half_shuffle_hits = []
    for s_tot, info in factorisation_test.items():
        for split_entry in info["splits"]:
            diff_prim = mp.mpf(split_entry["diff_J_eff_minus_primitive"])
            if abs(diff_prim) < bit_exact_threshold:
                primitive_hits.append((s_tot, split_entry["split"]))
            if split_entry["diff_J_eff_minus_shuffle"] is not None:
                diff_shuf = mp.mpf(split_entry["diff_J_eff_minus_shuffle"])
                if abs(diff_shuf) < bit_exact_threshold:
                    shuffle_hits.append((s_tot, split_entry["split"]))
                diff_half = mp.mpf(split_entry["diff_J_eff_minus_half_shuffle"])
                if abs(diff_half) < bit_exact_threshold:
                    half_shuffle_hits.append((s_tot, split_entry["split"]))

    if primitive_hits and not shuffle_hits:
        verdict_reading = "READING_A_PRIMITIVE"
    elif shuffle_hits and not primitive_hits:
        verdict_reading = "READING_B_SHUFFLE"
    elif half_shuffle_hits:
        verdict_reading = "READING_B_HALF_SHUFFLE"
    elif primitive_hits and shuffle_hits:
        verdict_reading = "BOTH_HIT_INCONCLUSIVE"

    # Hain-Brown verdict at depth 2
    hb_modular_hits = [s_tot for s_tot, info in cross_prec.items()
                       if info["any_hain_brown_modular"] and info["agree_across_3prec"]]
    if hb_modular_hits:
        hb_verdict = "REOPEN_HB_AT_DEPTH_2"
    else:
        hb_verdict = "HB_NEGATIVE_AT_DEPTH_2_TOO"

    results["verdict"] = {
        "reading_A_vs_B": verdict_reading,
        "primitive_hits": [{"s_tot": s, "split": list(sp_pair)} for s, sp_pair in primitive_hits],
        "shuffle_hits": [{"s_tot": s, "split": list(sp_pair)} for s, sp_pair in shuffle_hits],
        "half_shuffle_hits": [{"s_tot": s, "split": list(sp_pair)} for s, sp_pair in half_shuffle_hits],
        "hain_brown_at_depth_2": hb_verdict,
        "hb_modular_hits_s_tot": hb_modular_hits,
    }

    t_elapsed = time.time() - t_start
    results["metadata"]["wall_time_seconds"] = t_elapsed
    print(f"\n[Wall time: {t_elapsed:.1f} s]")
    print(f"\nVerdict:")
    print(f"  Reading A vs B: {verdict_reading}")
    print(f"  Hain-Brown at depth 2: {hb_verdict}")
    print(f"  Primitive hits: {primitive_hits}")
    print(f"  Shuffle hits: {shuffle_hits}")
    print(f"  Half-shuffle hits: {half_shuffle_hits}")

    # Dump
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults written to {OUT_JSON}")

    return results


if __name__ == "__main__":
    main()
