"""
Exact algebraic Breit-Pauli retarded Slater integrals for hydrogenic orbitals.
==============================================================================

Computes two-electron Breit-Pauli (BP) retarded radial integrals:

    R^k_BP(n1l1,n3l3; n2l2,n4l4)
        = ∫∫ R_{n1l1}(r1) R_{n3l3}(r1) × (r_<^k / r_>^{k+3})
                × R_{n2l2}(r2) R_{n4l4}(r2) r1² r2² dr1 dr2

for hydrogenic orbitals at unit nuclear charge (Z = 1). These arise from the
partial-wave decomposition of the Breit-Pauli two-body operators (spin-spin
tensor, spin-other-orbit, orbit-orbit) after the angular tensor projection
is applied to the distributional 1/r_12^3 kernel.

Also exports the standard Coulomb Slater integrals R^k(...) using the same
algebraic primitives (same code path, different kernel exponent).

Transcendental content
----------------------
For SOME orbital pairs, the BP-retarded integral is an exact rational in Z.
For OTHER pairs (including the diagonal 1s-1s-l=2 case), the integral
contains logarithms of rational arguments — the Jacobi-type symmetric
differences ``log(a) − log(b)`` that enter as Frullani limits when the
outer-variable integrand has 1/r_2^k divergences cancelled via Mellin
regularization. In Paper 18 terminology these are *embedding-log*
content: still rational-in-Z for same-n shells (after factoring out Z²,
Z³ etc.), but with a single transcendental seed log(2) or log(3/2)
whose coefficient is a clean rational.

Z-scaling
---------
R^k_BP(Z) = Z^3 × R^k_BP(Z=1) for all orbital pairs
R^k(Z)    = Z × R^k(Z=1) for Coulomb (same as ``hypergeometric_slater.py``)

Division into regions
---------------------
The double integral is split at r_1 = r_2 into two regions (I: r_1 < r_2;
II: r_1 > r_2). For each region, the inner integral (over r_<) gives a
lower incomplete gamma function, reducing to a closed polynomial plus
exp-times-polynomial form. The outer integral is then a sum of
∫_0^∞ r^p exp(-a r) dr terms, treated via Mellin continuation when p is
a negative integer.

The fix over BR-B
-----------------
The BR-B script ``debug/br_b_breit_radial.py`` used a simplified formula
that silently returned 0 whenever the naive region-splitting had m_1 < 0
AND m_2 < 0. This missed many physically relevant cases. The enhanced
machinery below tracks the integrand's negative-power singularities and
proves integrability via the global-cancellation condition before applying
the Mellin regularized value.

References
----------
- G.W.F. Drake, Phys. Rev. A 3, 908 (1971): Breit-Pauli radial integrals for He
- H.A. Bethe & E.E. Salpeter, Quantum Mechanics of One- and Two-Electron Atoms §§38-39
- W.R. Johnson, Atomic Structure Theory (Springer, 2007), Ch. 8
- A.R. Edmonds, Angular Momentum in Quantum Mechanics (1957), Chs. 6-7
- D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii, Quantum Theory of Angular
  Momentum (World Scientific, 1988), §5.17 (bipolar harmonic expansion).
- D.M. Brink & G.R. Satchler, Angular Momentum (3rd ed., 1993), App. 5.
- R.A. Sack, J. Math. Phys. 5, 245 (1964) (Gegenbauer expansion of r_12^{-n}).
- GeoVac Paper 18 (exchange-constant taxonomy)
- ``geovac/hypergeometric_slater.py`` (Coulomb R^k algebraic machinery, same pattern)

Drake 1971 combining coefficients for He 2^3P
---------------------------------------------
The He 2^3P fine structure evaluated from these radial integrals takes the
Drake form

    E(^3P_J) = (ζ_{2p} / 2) · X(J) + A_SS · f_SS(J) + A_SOO · f_SOO(J)

with X(J) = J(J+1) - L(L+1) - S(S+1), angular J-patterns

    f_SS (J=0,1,2) = (-2, +1, -1/5)
    f_SOO(J=0,1,2) = (+2, +1, -1)

(derived in Track DD from (-1)^{L+S+J}·6j{L S J; S L k} at L = S = 1 with
k = 2 for SS and k = 1 for SOO; see tests/test_breit_integrals.py), and
radial amplitudes

    A_SS  = α² · ( 3/50 · M²_dir - 2/5 · M²_exch )
    A_SOO = α² · ( 3/2  · M¹_dir -       M¹_exch )

where the combining rationals (3/50, -2/5, 3/2, -1) were identified by
rational search (Track BF-D) and reproduce NIST He 2^3P to sub-percent.
Sprint 5 Track DV (2026) fully characterized the bipolar channel structure
underlying these ratios: for He (1s)(2p) ^3P, only (k_1=0, k_2=2) direct
and (k_1=1, k_2=1) exchange channels are angular-allowed (no higher (k_1,k_2)
channels contribute at any L>=K). The piecewise decomposition of the bipolar
radial kernel into Drake's M^K basis reveals a structural mismatch:
bipolar(0,2,K=2) direct = M^0_dir_Region_I + M^2_dir_Region_II, which is
NOT proportional to M^2_dir_total. This confirms that Drake's rational
coefficients (3/50, -2/5) are a convention-dependent combining identity
rather than a direct Racah-algebra output; a fully symbolic closed-form
derivation remains open. (See debug/dv_drake_bipolar_memo.md.)

Track DP (Sprint 6, 2026) independently confirmed the combining rationals
via NIST extraction: the three NIST He 2^3P splittings were solved as a
linear system for (zeta, A_SS, A_SOO) and compared to the Drake predictions.
Result: A_SS match 0.09%, A_SOO match 0.07% — consistent with higher-order
corrections (mass polarization, QED) not included in the Breit-Pauli level.
Z^3 scaling of the radial integrals verified to machine precision at
Z=1,2,3,4,10, confirming Z-independence of the angular combining coefficients.
PSLQ identification of the rationals from NIST data was attempted but failed
due to insufficient NIST precision (~6 sig figs). Structural finding: the
rank-2 tensor operator C^(2)(hat_12)/r_12^3 has zero exchange angular
integral for (1s)(2p), confirming that Drake's M^k_exch refers to Slater
integral orbital ordering, not wavefunction antisymmetrization exchange.

Author: GeoVac Development Team (Track BF, April 2026; J-pattern Racah
derivation in Track DD, April 2026; bipolar structure characterization
in Track DV, April 2026; NIST confirmation in Track DP, April 2026)
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial, comb, gcd, isqrt
from typing import Dict, List, Tuple

import sympy as sp
from sympy import Rational, log, simplify, expand_log


__all__ = [
    "compute_radial",
    "compute_radial_Z",
    "breit_ss_radial",
    "breit_soo_radial",
    "breit_retarded",
    "coulomb_slater",
]


# ---------------------------------------------------------------------------
# Hydrogenic expansion (mirrors hypergeometric_slater._expand_product)
# ---------------------------------------------------------------------------


def _laguerre_coeffs(n: int, l: int) -> List[Fraction]:
    p = n - l - 1
    alpha = 2 * l + 1
    return [
        Fraction((-1) ** s * comb(p + alpha, p - s), factorial(s))
        for s in range(p + 1)
    ]


def _norm_sq(n: int, l: int) -> Fraction:
    return Fraction(2, n) ** 3 * Fraction(
        factorial(n - l - 1), 2 * n * factorial(n + l)
    )


def _expand_product(
    n_a: int, l_a: int, n_b: int, l_b: int
) -> Tuple[List[Tuple[Fraction, int, Fraction]], Fraction]:
    """Expand R_{n_a l_a}(r) R_{n_b l_b}(r) r^2 as list of (coeff, power, rate).

    Decay rate = 1/n_a + 1/n_b is shared across all terms.
    Returns (terms, N²_prod) where N²_prod = |N_{nala}|² × |N_{nblb}|².
    """
    ca = _laguerre_coeffs(n_a, l_a)
    cb = _laguerre_coeffs(n_b, l_b)
    nsq = _norm_sq(n_a, l_a) * _norm_sq(n_b, l_b)
    rate = Fraction(1, n_a) + Fraction(1, n_b)
    scale_la = Fraction(2, n_a) ** l_a
    scale_lb = Fraction(2, n_b) ** l_b
    base_power = l_a + l_b + 2

    terms: List[Tuple[Fraction, int, Fraction]] = []
    for s1, c1 in enumerate(ca):
        for s2, c2 in enumerate(cb):
            coeff = (
                c1 * c2
                * scale_la * scale_lb
                * Fraction(2, n_a) ** s1
                * Fraction(2, n_b) ** s2
            )
            power = base_power + s1 + s2
            terms.append((coeff, power, rate))
    return terms, nsq


def _fraction_sqrt(f: Fraction) -> Fraction:
    """Exact rational square root of a Fraction."""
    if f == 0:
        return Fraction(0)
    if f < 0:
        raise ValueError(f"Cannot sqrt negative Fraction {f}")
    g = gcd(f.numerator, f.denominator)
    num, den = f.numerator // g, f.denominator // g
    sn = isqrt(num)
    sd = isqrt(den)
    if sn * sn == num and sd * sd == den:
        return Fraction(sn, sd)
    # Fallback via high-precision decimal
    import decimal
    decimal.getcontext().prec = 50
    v = (decimal.Decimal(num) / decimal.Decimal(den)).sqrt()
    r = Fraction(v).limit_denominator(10 ** 15)
    if r * r == f:
        return r
    raise ValueError(f"Cannot take exact sqrt of {f}")


# ---------------------------------------------------------------------------
# Poly × exp integrator with Mellin regularization
# ---------------------------------------------------------------------------


_H_CACHE: Dict[int, sp.Expr] = {}


def _harmonic(k: int) -> sp.Expr:
    """Harmonic number H_{k-1} = sum_{j=1}^{k-1} 1/j as sympy Rational."""
    if k in _H_CACHE:
        return _H_CACHE[k]
    h = Rational(0)
    for j in range(1, k):
        h = h + Rational(1, j)
    _H_CACHE[k] = h
    return h


def _integrate_poly_exp_regularized(
    terms: List[Tuple[Fraction, int, Fraction]],
) -> sp.Expr:
    """Compute int_0^oo sum_i c_i r^{p_i} exp(-a_i r) dr.

    Uses factorial formula when p_i >= 0. For p_i = -k < 0, uses the
    Mellin-continuation value:

        int_0^oo r^{-k} e^{-a r} dr → (-1)^{k-1}/(k-1)! * a^{k-1} * [H_{k-1} − log(a)]

    (after the pole at Gamma(1-k) is cancelled by the moment conditions
    sum_i c_i a_i^{m} = 0 for m = 0, ..., k-2; the H_{k-1} piece is the
    psi(k) − (−γ) difference, and the Euler-γ term cancels globally
    when the pole cancels).

    Verifies integrability by expanding the negative-power sum in Taylor
    series around r=0 and confirming all coefficients of r^q (q < 0) vanish.

    Raises
    ------
    ValueError if the integrand has residual non-cancellable singularity
    at r = 0 (i.e., the original double integral diverges at coalescence).
    """
    negative_terms = [(c, p, a) for c, p, a in terms if p < 0]

    if negative_terms:
        p_min = min(p for c, p, a in negative_terms)
        for q in range(p_min, 0):
            coeff_q = Fraction(0)
            for c, p, a in negative_terms:
                n = q - p
                if n >= 0:
                    coeff_q += c * (-a) ** n / Fraction(factorial(n))
            if coeff_q != 0:
                raise ValueError(
                    f"Non-integrable: coefficient of r^{q} is {coeff_q} (!= 0); "
                    "integrand has non-cancellable singularity at r = 0."
                )

    result: sp.Expr = Rational(0)
    for c, p, a in terms:
        c_s = Rational(c.numerator, c.denominator)
        a_s = Rational(a.numerator, a.denominator)
        if p >= 0:
            result = result + c_s * factorial(p) / a_s ** (p + 1)
        else:
            k = -p
            sign = (-1) ** (k - 1)
            kfact = factorial(k - 1)
            result = result + c_s * sign * a_s ** (k - 1) / kfact * (_harmonic(k) - log(a_s))

    return _simplify_log_form(simplify(result))


def _simplify_log_form(expr: sp.Expr) -> sp.Expr:
    """Simplify sympy logs by factoring integer arguments into prime powers.

    Returns expression in canonical form: rational + sum_p n_p · log(p)
    for small primes p (typically p ∈ {2, 3, 5, ...}).

    Uses ``expand_log(force=True)`` to spread log(a/b) → log(a) − log(b),
    then factors each log(integer) into prime-powered pieces, then combines
    coefficients of like log(prime) terms.
    """
    from sympy import factorint, expand_log as sp_expand_log

    # First spread logs: log(p/q) → log(p) − log(q)
    expr = sp_expand_log(sp.expand(expr), force=True)

    # Now factor any remaining log(integer) subexpressions
    subs: Dict[sp.Expr, sp.Expr] = {}
    for log_expr in expr.atoms(sp.log):
        arg = log_expr.args[0]
        if not (hasattr(arg, "is_positive") and arg.is_positive):
            continue
        p, q = None, None
        if arg.is_Rational:
            p = arg.p
            q = arg.q
        elif arg.is_Integer:
            p = int(arg)
            q = 1
        else:
            continue
        if p is None or p <= 0:
            continue
        # Only factor if reasonably sized (ripgrep factorint bounds)
        if max(abs(p), abs(q)) > 10**500:
            continue
        try:
            p_factors = factorint(p) if p > 1 else {}
            q_factors = factorint(q) if q > 1 else {}
        except Exception:
            continue
        result = Rational(0)
        for prime, exp_ in p_factors.items():
            result = result + exp_ * log(Rational(prime))
        for prime, exp_ in q_factors.items():
            result = result - exp_ * log(Rational(prime))
        # Only substitute when decomposition is strictly simpler: multiple factors,
        # or a single repeated factor (e.g. 2^48 → 48·log(2)).
        p_has_composite = len(p_factors) > 1 or any(e > 1 for e in p_factors.values())
        q_has_composite = len(q_factors) > 1 or any(e > 1 for e in q_factors.values())
        if p_has_composite or q_has_composite or (p_factors and q_factors):
            subs[log_expr] = result

    if subs:
        expr = expr.subs(subs)

    # Collect like terms: gather coefficients of each log(prime)
    # Use sp.collect on the log atoms
    log_atoms = list(expr.atoms(sp.log))
    if log_atoms:
        expr = sp.collect(sp.expand(expr), log_atoms)
    return expr


# ---------------------------------------------------------------------------
# T-kernel (single-pair) for Coulomb and Breit kernels
# ---------------------------------------------------------------------------


def _t_kernel_region_I(
    a: int, b: int, alpha: Fraction, beta: Fraction, kernel_type: str, k: int
) -> sp.Expr:
    """Compute int int r1^a r2^b exp(-alpha r1 - beta r2) K(r_<, r_>) [r1<r2] dr1 dr2.

    K(r_<, r_>) = r_<^k / r_>^{k+1} for kernel_type='coulomb'
              = r_<^k / r_>^{k+3} for kernel_type='breit'
    """
    alpha = Fraction(alpha)
    beta = Fraction(beta)
    N = a + k
    if kernel_type == "coulomb":
        m = b - k - 1
    elif kernel_type == "breit":
        m = b - k - 3
    else:
        raise ValueError(f"unknown kernel_type {kernel_type!r}")

    N_fact = Fraction(factorial(N))
    terms: List[Tuple[Fraction, int, Fraction]] = []
    terms.append((N_fact / alpha ** (N + 1), m, beta))
    for j in range(N + 1):
        coeff = -N_fact / (Fraction(factorial(j)) * alpha ** (N - j + 1))
        terms.append((coeff, m + j, alpha + beta))
    return _integrate_poly_exp_regularized(terms)


def _t_kernel(
    a: int, b: int, alpha: Fraction, beta: Fraction, kernel_type: str, k: int
) -> sp.Expr:
    """Full T-kernel (both regions combined)."""
    rI = _t_kernel_region_I(a, b, alpha, beta, kernel_type, k)
    rII = _t_kernel_region_I(b, a, beta, alpha, kernel_type, k)  # Region II by r1↔r2 symmetry
    return _simplify_log_form(simplify(rI + rII))


# ---------------------------------------------------------------------------
# Public API: general radial integral
# ---------------------------------------------------------------------------


def compute_radial(
    n1: int, l1: int, n3: int, l3: int,
    n2: int, l2: int, n4: int, l4: int,
    k: int,
    kernel_type: str = "breit",
    Z: int = 1,
) -> sp.Expr:
    """Compute the two-electron radial integral at integer nuclear charge Z.

    For kernel_type='breit':
        R^k_BP(n1l1, n3l3; n2l2, n4l4) = ∫∫ P_{13}(r1) P_{24}(r2)
                                                (r_<^k / r_>^{k+3}) dr1 dr2 × Z^3
    For kernel_type='coulomb':
        R^k(n1l1, n3l3; n2l2, n4l4) = ∫∫ P_{13}(r1) P_{24}(r2)
                                              (r_<^k / r_>^{k+1}) dr1 dr2 × Z
    where P_{ab}(r) = R_a(r) R_b(r) r^2 at unit orbital exponent k_orb = 1.

    Parameters
    ----------
    n1, l1, n3, l3 : int
        Quantum numbers of the two orbitals on electron 1.
    n2, l2, n4, l4 : int
        Quantum numbers of the two orbitals on electron 2.
    k : int
        Multipole order (>= 0).
    kernel_type : {'coulomb', 'breit'}
        Which two-electron kernel to use.
    Z : int, default 1
        Nuclear charge. Coulomb integrals scale as Z^1, Breit as Z^3.

    Returns
    -------
    sympy Expr
        Closed form in sympy (rational, or rational + sum of log(p/q) for
        rational p, q).

    Raises
    ------
    ValueError if the integrand diverges at r_1 = r_2 (the orbital pair does
    not have sufficient amplitude power to regularize the kernel singularity).
    """
    terms13, nsq13 = _expand_product(n1, l1, n3, l3)
    terms24, nsq24 = _expand_product(n2, l2, n4, l4)
    norm = _fraction_sqrt(nsq13 * nsq24)

    integral: sp.Expr = Rational(0)
    for c13, p13, a13 in terms13:
        for c24, p24, a24 in terms24:
            if c13 == 0 or c24 == 0:
                continue
            T = _t_kernel(p13, p24, a13, a24, kernel_type, k)
            integral = integral + c13 * c24 * T

    norm_s = Rational(norm.numerator, norm.denominator)
    value_Z1 = _simplify_log_form(simplify(norm_s * integral))

    if Z == 1:
        return value_Z1

    if kernel_type == "breit":
        return simplify(Rational(Z) ** 3 * value_Z1)
    elif kernel_type == "coulomb":
        return simplify(Rational(Z) * value_Z1)
    else:
        raise ValueError(f"unknown kernel_type {kernel_type!r}")


def compute_radial_Z(
    n1: int, l1: int, n3: int, l3: int,
    n2: int, l2: int, n4: int, l4: int,
    k: int,
    Z: int,
    kernel_type: str = "breit",
) -> sp.Expr:
    """Alias for ``compute_radial`` with Z as keyword argument."""
    return compute_radial(
        n1, l1, n3, l3, n2, l2, n4, l4, k, kernel_type=kernel_type, Z=Z
    )


# ---------------------------------------------------------------------------
# Drake 1971 amplitude API (rank-2 spin-spin, rank-1 SOO)
# ---------------------------------------------------------------------------


def breit_ss_radial(
    n_a: int, l_a: int, n_b: int, l_b: int,
    n_c: int, l_c: int, n_d: int, l_d: int,
    k: int, Z: int = 1,
) -> sp.Expr:
    """Rank-2 spin-spin magnetic radial integral M^k(a, b; c, d).

    Physics: this is the radial integral that appears as the radial
    amplitude of the rank-2 spin-spin tensor operator, after angular
    projection onto multipole k, in Drake 1971 notation.

    Convention: orbitals a, c are on electron 1; orbitals b, d are on
    electron 2. The integral is symmetric under (a<->c, b<->d) due to
    the radial density product.

    Returns the exact closed form at nuclear charge Z (including Z^3 scaling
    for Breit-Pauli retarded kernels).

    Reference: Drake 1971, Phys. Rev. A 3, 908; Bethe-Salpeter §§38-39.
    """
    return compute_radial(n_a, l_a, n_c, l_c, n_b, l_b, n_d, l_d, k,
                          kernel_type="breit", Z=Z)


def breit_soo_radial(
    n_a: int, l_a: int, n_b: int, l_b: int,
    n_c: int, l_c: int, n_d: int, l_d: int,
    k: int, Z: int = 1,
) -> sp.Expr:
    """Rank-1 spin-other-orbit radial integral N^k(a, b; c, d).

    Same radial structure as ``breit_ss_radial`` — both rank-2 SS and
    rank-1 SOO project onto the same r_<^k / r_>^{k+3} kernel after
    angular reduction. The angular coefficient (from 3j/6j symbols, not
    included here) distinguishes the two operators.

    Reference: Drake 1971, Bethe-Salpeter §39.
    """
    return compute_radial(n_a, l_a, n_c, l_c, n_b, l_b, n_d, l_d, k,
                          kernel_type="breit", Z=Z)


def breit_retarded(
    n_a: int, l_a: int, n_b: int, l_b: int,
    n_c: int, l_c: int, n_d: int, l_d: int,
    k: int, Z: int = 1,
    operator: str = "SS",
) -> sp.Expr:
    """General retarded Breit-Pauli radial integral.

    Parameters
    ----------
    operator : {'SS', 'SOO'}
        Which tensor operator to use. Currently the radial integrals are
        identical for SS and SOO (differ only in the angular prefactor),
        so this parameter is currently informational.

    Returns
    -------
    sympy Expr
        Closed-form radial integral at nuclear charge Z.
    """
    if operator not in {"SS", "SOO"}:
        raise ValueError(f"operator must be 'SS' or 'SOO', got {operator!r}")
    return compute_radial(n_a, l_a, n_c, l_c, n_b, l_b, n_d, l_d, k,
                          kernel_type="breit", Z=Z)


def coulomb_slater(
    n_a: int, l_a: int, n_b: int, l_b: int,
    n_c: int, l_c: int, n_d: int, l_d: int,
    k: int, Z: int = 1,
) -> sp.Expr:
    """Standard Coulomb Slater integral R^k(a, b; c, d) at nuclear charge Z.

    Convenient wrapper; same as ``compute_radial(..., kernel_type='coulomb')``.
    Duplicates the functionality of ``geovac.hypergeometric_slater.compute_rk_algebraic``
    but uses the enhanced Mellin-regularized machinery (handles cases where
    the old code raised a ``ValueError`` for negative factorial).
    """
    return compute_radial(n_a, l_a, n_c, l_c, n_b, l_b, n_d, l_d, k,
                          kernel_type="coulomb", Z=Z)
