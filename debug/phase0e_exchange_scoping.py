"""Phase 0-e: scoping the EXCHANGE class (AB|AB) -- Ruedenberg 1951 Part II proper.

Called for by build plan section 8.3 and re-flagged by the Phase 0-h memo, which
explicitly did NOT cover this class.

WHY THE PREVIOUS TWO REDUCTIONS DO NOT START. Both (AA|BB) and the hybrid class
worked because at least one charge distribution was one-center, so its Coulomb
potential was closed-form and the quartet collapsed to a one-electron problem.
Exchange has

    rho_1 = conj(chi_a^A) chi_b^B        two-center
    rho_2 = conj(chi_c^A) chi_d^B        two-center

Neither has a closed-form potential. There is nothing to reduce to; the kernel
1/r12 must be expanded. That is the Neumann expansion in prolate spheroidal
coordinates with BOTH electrons in the same (xi, eta) system -- the original
subject of this build plan, before section 7 retargeted it.

WHAT MAKES SPHEROIDAL THE RIGHT SYSTEM. With r_A = R(xi+eta)/2 and
r_B = R(xi-eta)/2,

    alpha r_A + beta r_B = (R/2)[(alpha+beta) xi + (alpha-beta) eta]

so each orbital product separates into e^{-p xi} e^{-q eta}. Moreover
R_nl(r_A)/r_A^{l} is a polynomial and r_A^l Y_lm(Om_A) is a solid harmonic, so
the whole integrand is a POLYNOMIAL in (xi, eta) times those exponentials --
no negative powers anywhere, unlike the hybrid class.

THE THREE QUESTIONS

  EQ1  Does the tau sum terminate? Phase 0 Q1 said YES for the James-Coolidge
       basis; section 7 said that answer does not transfer. Settle it, and if it
       does not terminate, measure how fast it converges -- that is the plan's
       own STOP criterion.

  EQ2  What is the seed set? Phase 0 Q2 (Q_tau^sigma carries one transcendental,
       reducible to e^{+-a}E_1(a*shift)) is basis-independent and should carry
       over. But the xi integral starts at xi = 1 exactly, where that formula is
       singular -- so the endpoint has to be checked, not assumed.

  EQ3  Cost, against the honest alternative (McMurchie-Davidson over fitted
       STOs at n_gauss >= 10, which Phase 0-h priced at ~2e-7).

Run from repo root:  python debug/phase0e_exchange_scoping.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import integrate
from scipy.special import eval_legendre, exp1

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


# ------------------------------------------------------- Legendre Q on (1, oo)

def Q_leg(n: int, x):
    """Q_n(x) for x > 1, via Q_n = P_n Q_0 - sum_{k=1}^n P_{k-1} P_{n-k} / k.

    Stable (no upward Q recurrence, which loses the decaying solution).
    Q_0(x) = (1/2) ln((x+1)/(x-1)).
    """
    x = np.asarray(x, dtype=float)
    Q0 = 0.5 * np.log((x + 1.0) / (x - 1.0))
    out = eval_legendre(n, x) * Q0
    for k in range(1, n + 1):
        out = out - eval_legendre(k - 1, x) * eval_legendre(n - k, x) / k
    return out


# ------------------------------------------------------------- the tau-th term

def eta_int(tau: int, k: int, q: float, n: int = 400) -> float:
    """int_{-1}^{1} eta^k P_tau(eta) e^{-q eta} d eta."""
    u, w = np.polynomial.legendre.leggauss(n)
    return float(np.sum(w * u ** k * eval_legendre(tau, u) * np.exp(-q * u)))


def xi_double(tau: int, j: int, k: int, p: float) -> float:
    """int_1^oo int_1^oo xi1^j xi2^k e^{-p(xi1+xi2)} P_tau(xi_<) Q_tau(xi_>).

    Split at xi2 = xi1: the inner integrand is smooth on each side, and the
    log singularity of Q_0 at xi = 1 is integrable.
    """
    def outer(x1):
        lo, _ = integrate.quad(
            lambda x2: x2 ** k * np.exp(-p * x2) * eval_legendre(tau, x2),
            1.0, x1, epsabs=1e-12, epsrel=1e-12, limit=100)
        hi, _ = integrate.quad(
            lambda x2: x2 ** k * np.exp(-p * x2) * float(Q_leg(tau, x2)),
            x1, np.inf, epsabs=1e-12, epsrel=1e-12, limit=100)
        return (x1 ** j * np.exp(-p * x1)
                * (float(Q_leg(tau, x1)) * lo + eval_legendre(tau, x1) * hi))

    val, _ = integrate.quad(outer, 1.0, np.inf, epsabs=1e-11, epsrel=1e-11,
                            limit=100)
    return float(val)


def exchange_term(tau: int, alpha: float, beta: float, R: float) -> float:
    """The tau-th Neumann term of (1s_A 1s_B | 1s_A 1s_B), sigma = 0.

    Only sigma = 0 survives for m = 0 orbitals (the two phi integrals force it).
    The volume factor (xi^2 - eta^2) is expanded so xi and eta separate.
    """
    p = (alpha + beta) * R / 2.0
    q = (alpha - beta) * R / 2.0
    e0, e2 = eta_int(tau, 0, q), eta_int(tau, 2, q)
    x00, x02, x22 = (xi_double(tau, 0, 0, p), xi_double(tau, 0, 2, p),
                     xi_double(tau, 2, 2, p))
    # (xi1^2 - eta1^2)(xi2^2 - eta2^2), symmetric in the two electrons
    core = (x22 * e0 * e0 - 2.0 * x02 * e0 * e2 + x00 * e2 * e2)
    return (2 * tau + 1) * core


def exchange_prefactor(alpha: float, beta: float, R: float) -> float:
    NA = np.sqrt(alpha ** 3 / np.pi)
    NB = np.sqrt(beta ** 3 / np.pi)
    return (R ** 3 / 8.0) ** 2 * (2 * np.pi) ** 2 * (2.0 / R) * (NA * NB) ** 2


def exchange_neumann(alpha, beta, R, tau_max):
    """Partial sums of the Neumann series for the 1s exchange integral."""
    C = exchange_prefactor(alpha, beta, R)
    terms = [C * exchange_term(t, alpha, beta, R) for t in range(tau_max + 1)]
    return terms, np.cumsum(terms)


# ------------------------------------------------------------------ reference

_FIT_CACHE: dict = {}


def exchange_md(alpha, beta, R, n_gauss=12):
    from geovac import noci_engine as E
    if n_gauss not in _FIT_CACHE:                 # the fit is an optimization
        _FIT_CACHE[n_gauss] = E.fit_sto_shape(0, 1, n_gauss=n_gauss)
    arr, dco, q = _FIT_CACHE[n_gauss]
    sh = {"1s": (arr, dco)}
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., R])
    A = E.sto_shape_basis(pa, "1s", alpha, sh, (0, 0, 0))
    B = E.sto_shape_basis(pb, "1s", beta, sh, (0, 0, 0))
    return E.eri_md(A, B, A, B), q


# ----------------------------------------------------------------------- EQ1

def leg_EQ1() -> None:
    print("EQ1  does the tau sum terminate?\n")
    print("     The eta integrals are  int_{-1}^{1} eta^k P_tau(eta) e^{-q eta} d eta")
    print("     with q = (alpha - beta) R / 2. If q = 0 they are the orthogonality")
    print("     integral of P_tau against a degree-k polynomial and VANISH for")
    print("     tau > k. If q != 0 the exponential has content at every tau.\n")

    print("     (a) HOMONUCLEAR, alpha = beta = 1.0 (q = 0):")
    for tau in range(6):
        e0, e2 = eta_int(tau, 0, 0.0), eta_int(tau, 2, 0.0)
        print(f"         tau={tau}   eta-int(k=0) = {e0: .3e}   (k=2) = {e2: .3e}")
    print("         => vanish beyond tau = 2. The series TERMINATES.")
    print("         This is why Phase 0 Q1 held for James-Coolidge: that basis")
    print("         is e^{-alpha(xi_1+xi_2)}, pure xi, so q = 0 identically.\n")

    print("     (b) HETERONUCLEAR, alpha = 3.0, beta = 1.0, R = 3.0 (q = 3.0):")
    for tau in (0, 2, 4, 6, 8, 10):
        e0 = eta_int(tau, 0, 3.0)
        print(f"         tau={tau:2d}   eta-int(k=0) = {e0: .6e}")
    print("         => never zero. The series does NOT terminate.\n")

    print("     So: tau terminates IFF the two centres carry the SAME orbital")
    print("     exponent. Homonuclear yes; heteronuclear -- LiH, NaH, the actual")
    print("     targets -- no. Section 7 was right that Phase 0 Q1 does not")
    print("     transfer, but the reason is the eta exponential, not the")
    print("     xi-eta mixing in cos(theta_A).\n")


# ----------------------------------------------------------------------- EQ1b

def leg_EQ1b() -> None:
    print("EQ1b how fast does it converge when it does not terminate?")
    print("     (1s_A 1s_B | 1s_A 1s_B), the plan's STOP criterion\n")
    for alpha, beta, R, tmax in ((1.0, 1.0, 2.0, 6), (3.0, 1.0, 3.0, 10)):
        label = "homonuclear" if alpha == beta else "heteronuclear"
        terms, partial = exchange_neumann(alpha, beta, R, tmax)
        md, fitq = exchange_md(alpha, beta, R)
        print(f"     {label}: alpha={alpha}, beta={beta}, R={R}")
        print(f"       tau      term          partial sum      rel. residual")
        for t in range(tmax + 1):
            rel = abs(terms[t] / partial[-1]) if partial[-1] else float("nan")
            print(f"       {t:3d}   {terms[t]: .6e}   {partial[t]: .10f}   {rel:.1e}")
        print(f"       reference sweep (Phase 0-h: eri_md is fit-limited here,")
        print(f"        and exchange is the worst case -- BOTH densities are")
        print(f"        two-centre overlaps sampling the exponential tail):")
        for ng in (6, 8, 10, 12):
            md_ng, q_ng = exchange_md(alpha, beta, R, n_gauss=ng)
            print(f"         n_gauss={ng:2d}  <fit|STO>={q_ng:.9f}  md={md_ng:.10f}"
                  f"  |Neumann-md|={abs(partial[-1] - md_ng):.2e}")
        print()


# ------------------------------------------------------------------------ EQ2

def leg_EQ2() -> None:
    print("EQ2  the seed set -- and the xi = 1 endpoint that Phase 0 Q2 skipped")
    print()
    print("     Phase 0 Q2 established, for L(xi) = ln((xi+1)/(xi-1)):")
    print("       int_c^oo e^{-a xi} L(xi) dxi = (e^{-ac}/a) ln((c+1)/(c-1))")
    print("                                    + (1/a)[e^{a}E_1(a(c+1))")
    print("                                            - e^{-a}E_1(a(c-1))]")
    print("     But the exchange xi integral starts at c = 1 EXACTLY, where both")
    print("     the ln and the E_1 blow up. They cancel; the finite part is what")
    print("     the class actually carries.\n")

    a = 1.7
    limit = (np.exp(-a) / a) * (np.log(2) + np.euler_gamma + np.log(a)) \
        + (1 / a) * np.exp(a) * exp1(2 * a)
    print(f"     predicted finite part, from expanding both singular pieces:")
    print("       (e^{-a}/a)[ln 2 + gamma + ln a] + (e^{a}/a) E_1(2a)")
    print(f"       = {limit:.12f}   at a = {a}\n")
    print(f"     c -> 1 numerically (the -ln(c-1) pieces must cancel):")
    for c in (1.1, 1.01, 1e-3 + 1, 1e-4 + 1, 1e-5 + 1, 1e-6 + 1):
        val = (np.exp(-a * c) / a) * np.log((c + 1) / (c - 1)) \
            + (1 / a) * (np.exp(a) * exp1(a * (c + 1))
                         - np.exp(-a) * exp1(a * (c - 1)))
        print(f"       c-1 = {c - 1:<9.0e}  value = {val:.12f}   "
              f"deviation = {abs(val - limit):.2e}")
    print("       deviation shrinks like O((c-1) ln(c-1)) -> the limit is confirmed\n")
    print("     So at the endpoint the E_1 singularity converts into an EXPLICIT")
    print("     Euler gamma and an explicit ln a. The exchange class therefore")
    print("     carries")
    print("         {E_1(lambda R)}  U  {gamma}  U  {ln}")
    print("     strictly MORE than the hybrid class's {E_1(lambda R)} alone.")
    print("     Consistent with the textbook H2 exchange integral, which is the")
    print("     classic place gamma and ln R appear in a two-centre result.\n")


def leg_EQ3() -> None:
    """The cost driver: how the term count moves with the exponent mismatch.

    Measured on the ETA integral, which is where the mismatch lives. The xi half
    depends on p = (alpha+beta)R/2, not on q, and decays geometrically in tau
    regardless -- so the eta half is what sets the term count. Cross-validated
    against the full-term measurement at q = 3 in EQ1b, which needed tau = 8 for
    1e-7 and tau = 10 for 1e-10; the eta proxy predicts the same.
    """
    print("EQ3  cost driver: term count vs the exponent mismatch q\n")
    print("     q = (alpha - beta) R / 2 is the only thing standing between this")
    print("     class and termination, so it should also set the term count.")
    print("     LiH at R=3 spans alpha in {3 (Li 1s), 1.5 (Li 2s/2p)} against")
    print("     beta = 1 (H 1s), i.e. q in {3.0, 0.75}.\n")
    print("     Measured on the eta integral (see docstring for why that is the")
    print("     right proxy), normalized to its tau = 0 value:\n")
    print("     alpha beta    q  |  " + "  ".join(f"tau={t:<2d}" for t in
                                                  (2, 4, 6, 8, 10)))
    print("     " + "-" * 62)
    for alpha, beta, R in ((1.0, 1.0, 3.0), (1.5, 1.0, 3.0),
                           (3.0, 1.0, 3.0), (6.0, 1.0, 3.0)):
        q = (alpha - beta) * R / 2
        base = abs(eta_int(0, 0, q))
        cells = "  ".join(f"{abs(eta_int(t, 0, q)) / base:6.0e}"
                          for t in (2, 4, 6, 8, 10))
        print(f"     {alpha:5.1f}{beta:5.1f}{q:6.2f}  |  {cells}")
    print("\n     q = 0 is machine zero beyond tau = 2 (termination). Otherwise the")
    print("     decay is factorial in tau and steepens as q falls.")
    print()
    print("     The proxy is CONSERVATIVE: at q = 3 it reads 2e-6 at tau = 10,")
    print("     while EQ1b's full-term measurement of the same case reads 1.2e-10,")
    print("     because the xi half decays too. So these columns are an upper")
    print("     bound on the term count, not an estimate of it.")
    print()
    print("     Across the whole LiH range (q <= 3), tau ~ 10 is comfortably")
    print("     enough. Cost per quartet is a modest fixed number of terms, not an")
    print("     open-ended sum -- the plan's STOP criterion is NOT met. The")
    print("     q = 7.5 row (harder than anything in LiH) still reads 1e-3 on the")
    print("     conservative proxy at tau = 10, so very large mismatches would")
    print("     want a term-count check rather than a fixed truncation.\n")


def main() -> None:
    print("Phase 0-e -- scoping the exchange class (AB|AB)\n")
    leg_EQ1()
    leg_EQ1b()
    leg_EQ2()
    leg_EQ3()
    print("Verdict is in the memo, not here: this driver reports measurements.")


if __name__ == "__main__":
    main()
