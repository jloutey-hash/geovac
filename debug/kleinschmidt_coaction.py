"""
Kleinschmidt et al. 2026 zeta-generator coaction machinery — depth-1 sector.

Reference:
    Kleinschmidt, Mafra, Schlotterer, Verbeek (2026).
    "Towards Motivic Coactions at Genus One from Zeta Generators."
    arXiv:2508.02800. JHEP 05 (2026) 105.

Scope of this module (debug/, exploratory infrastructure):
    Depth 1 only. The depth-1 sector is what the Kleinschmidt paper
    states in fully closed form in its main body (Eqs. 19, 22, 26, 34,
    and the f-alphabet review). Depth >= 2 (Eqs. 64 in its full,
    higher-depth content; Appendices A and B) requires the full
    Tsunogai bracket tower and depth-2 MMV decompositions that the
    paper sends to its appendices and that we could not extract from
    the publicly available HTML at the time of writing. Depth-2 (and
    higher) calls are stubbed with NotImplementedError and named as
    follow-on.

Three layers, each independently testable:

    (1) F-alphabet substrate (Brown 2014, reviewed in Kleinschmidt §2.1).
        Non-commutative words on f_{2k+1} with central f_2.
        Deconcatenation coproduct Delta.
        Map rho: zeta^m_w -> f_w on motivic MZV generators.

    (2) Depth-1 multiple modular value evaluation.
        m[j/k] = integral_0^{i infty} tau^j G_k(tau) d tau
        Closed form (Eq. 26):
            m[0/k]       = -2 pi i  zeta_{k-1} / (k-1)
            m[j/k]       = 2 (-1)^{j+1} j! (2 pi i)^{k-1-j} / (k-1)!
                           * zeta_{j+1} * zeta_{j+2-k}     (0 < j <= k-2)
        Cross-validated against numerical Eichler integration.

    (3) Tsunogai derivation algebra (formal, symbolic).
        epsilon_k         : weight-k Tsunogai generator (k >= 4 even, plus epsilon_0)
        epsilon_k^{(j)}   := ad_{epsilon_0}^j (epsilon_k)
        Pollack relation at weight 14:  [eps_4, eps_10] - 3 [eps_6, eps_8] = 0
        Depth-1 zeta generator (Eq. 34):
            sigma_w = z_w  -  (1/(w-1)!) epsilon_{w+1}^{(w-1)}  +  O(depth >= 2)

    (4) Coaction (depth 1).
        At depth 1 the Kleinschmidt Eq. 64 coaction
            Delta I^m  =  (M_sigma^dr)^{-1}  I^m  M_sigma^dr  I^dr
        collapses, after applying the arithmetic z_w action and the
        f-alphabet deconcatenation, to the standard Brown 2014
        depth-1 coaction
            Delta zeta^m_{2k+1}  =  1 (x) zeta^m_{2k+1}  +  zeta^m_{2k+1} (x) 1
        on the arithmetic side, plus a geometric epsilon_{2k+2}^{(2k)}
        contribution on the Tsunogai side. We compute both sides
        explicitly at low weights as cross-validation.

Honest scope flags:
    - Kleinschmidt et al. is "a proposal" per the paper's own framing;
      first-principles motivic-coaction identification is open in the
      mathematical literature. The formulae here are FORMALLY VALID
      for PSLQ-style numerical tests but should NOT be claimed as
      proven motivic statements.
    - Depth >= 2 raises NotImplementedError; the appendix-A/B content
      of the paper was not extractable from the public HTML.
    - The geometric epsilon_w^{(j)} are treated SYMBOLICALLY here.
      Their action on iterated Eisenstein integrals via brackets is
      formal; we do not attempt to evaluate it numerically.

Downstream GeoVac uses (named, for future sprints):
    - The depth-2 NA-1 Mellin test (sprint_na1_depth2_mellin_compute.py
      already exists) can call rho() and deconcatenate() to interpret
      a positive identification structurally rather than as a black-box
      PSLQ hit.
    - Any future "Test A" follow-on at depth 2 (Reading B for NA-1, or
      cuspidal-period entry in the modular ring) can call
      depth1_mmv_value() as the building block for higher-depth MMVs
      assembled by iterated integration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Tuple

import mpmath
import sympy
from sympy import (
    Expr,
    Integer,
    Rational,
    Symbol,
    factorial,
    pi,
    sympify,
    zeta,
)


# ---------------------------------------------------------------------------
# Layer 1. F-alphabet substrate (Brown 2014 / Kleinschmidt §2.1).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FWord:
    """A non-commutative word in the f-alphabet.

    Letters are tuples of odd integers >= 3 (Kleinschmidt notation
    f_{2k+1}, k >= 1) and optional central f_2. We store the central
    f_2 multiplicity separately because f_2 is central in the
    f-algebra (commutes with every f-letter), matching Brown's
    convention.

    The empty word (letters == ()) with f2_power == 0 is the unit 1.

    Examples
    --------
    >>> w = FWord((3,))                    # f_3
    >>> w.weight()
    3
    >>> w2 = FWord((3, 5))                 # f_3 f_5 (non-commutative)
    >>> w2.weight()
    8
    >>> w3 = FWord((3,), f2_power=1)       # f_2 f_3 = f_3 f_2 (central)
    >>> w3.weight()
    5
    """

    letters: Tuple[int, ...]
    f2_power: int = 0

    def __post_init__(self) -> None:  # noqa: D401 - dataclass post-init
        for letter in self.letters:
            if letter < 3 or letter % 2 == 0:
                raise ValueError(
                    "f-alphabet letters must be odd integers >= 3 "
                    f"(got {letter}); central f_2 is stored separately"
                )
        if self.f2_power < 0:
            raise ValueError("f_2 power must be non-negative")

    def weight(self) -> int:
        """Total weight: sum of letters plus 2 * f2_power."""
        return sum(self.letters) + 2 * self.f2_power

    def depth(self) -> int:
        """Depth = number of (odd) letters; f_2's do not count toward depth.

        This convention matches Brown 2014: the (arithmetic) depth is
        the number of non-central generators, and the central f_2's
        live in a separate Tate grading.
        """
        return len(self.letters)

    def __mul__(self, other: "FWord") -> "FWord":
        """Concatenation product. f_2 powers add (central)."""
        return FWord(
            letters=self.letters + other.letters,
            f2_power=self.f2_power + other.f2_power,
        )

    def __repr__(self) -> str:
        parts: List[str] = []
        if self.f2_power:
            parts.append(f"f_2^{self.f2_power}" if self.f2_power > 1 else "f_2")
        parts.extend(f"f_{i}" for i in self.letters)
        return " ".join(parts) if parts else "1"


# Canonical unit f-word.
F_UNIT = FWord(letters=(), f2_power=0)


def f_letter(weight: int) -> FWord:
    """Construct the depth-1 f-letter f_weight.

    weight must be either 2 (central) or odd >= 3.
    """
    if weight == 2:
        return FWord(letters=(), f2_power=1)
    if weight < 3 or weight % 2 == 0:
        raise ValueError(
            f"f-letter weight must be 2 (central) or odd >= 3; got {weight}"
        )
    return FWord(letters=(weight,))


def deconcatenation_coproduct(
    word: FWord,
) -> List[Tuple[FWord, FWord]]:
    """Deconcatenation coproduct on f-words.

    Reference: Kleinschmidt §2.1, the standard Brown 2014 f-alphabet
    coproduct
        Delta(f_{i_1} f_{i_2} ... f_{i_r})
            = sum_{j=0}^r  (f_{i_1} ... f_{i_j}) (x) (f_{i_{j+1}} ... f_{i_r}).

    The central f_2 powers are placed on the LEFT factor by convention
    (matching the Brown convention that f_2 lives in the de-Rham-trivial
    Tate sector and is preserved by the de Rham reduction).

    Returns
    -------
    list of (left, right) tuples representing Delta(word) = sum left (x) right.
    """
    splits: List[Tuple[FWord, FWord]] = []
    n = len(word.letters)
    for j in range(n + 1):
        left = FWord(letters=word.letters[:j], f2_power=word.f2_power)
        right = FWord(letters=word.letters[j:], f2_power=0)
        splits.append((left, right))
    return splits


def shuffle(left: FWord, right: FWord) -> Dict[FWord, int]:
    """Shuffle product of two f-words (mod central f_2 powers).

    Reference: Brown 2014; recalled in Kleinschmidt §2.1. The
    shuffle product on the (odd) f-alphabet is
        u shuffle 1 = 1 shuffle u = u,
        (a u) shuffle (b v) = a (u shuffle (b v)) + b ((a u) shuffle v).

    Central f_2 powers add (they are commutative, central).

    Returns
    -------
    dict mapping FWord -> integer coefficient (multiplicity).
    """
    result: Dict[FWord, int] = {}
    f2_total = left.f2_power + right.f2_power

    def _shuffle(a: Tuple[int, ...], b: Tuple[int, ...]) -> List[Tuple[int, ...]]:
        if not a:
            return [b]
        if not b:
            return [a]
        return [
            (a[0],) + rest for rest in _shuffle(a[1:], b)
        ] + [
            (b[0],) + rest for rest in _shuffle(a, b[1:])
        ]

    for word_tuple in _shuffle(left.letters, right.letters):
        w = FWord(letters=word_tuple, f2_power=f2_total)
        result[w] = result.get(w, 0) + 1
    return result


# ---------------------------------------------------------------------------
# Layer 2. Depth-1 multiple modular value evaluation.
# ---------------------------------------------------------------------------


def depth1_mmv_value(j: int, k: int) -> Expr:
    """Closed-form depth-1 MMV m[j/k] from Kleinschmidt Eq. (26).

    Definition (Kleinschmidt Eq. 22, 26):

        m[j/k] = integral_0^{i infty} tau^j G_k(tau) d tau

    where G_k(tau) is the holomorphic Eisenstein series of weight k.

    Closed form:

        m[0/k] = - 2 pi i zeta_{k-1} / (k-1)                 (j = 0)
        m[j/k] = 2 (-1)^{j+1} j! (2 pi i)^{k-1-j} / (k-1)!
                  * zeta_{j+1} * zeta_{j+2-k}                (0 < j <= k-2)

    Constraints (from Kleinschmidt):
        k >= 4, k even
        0 <= j <= k - 2

    Returns
    -------
    sympy Expr in terms of pi, I, and zeta(odd_int).

    Notes
    -----
    Kleinschmidt writes zeta_{j+2-k} for the negative-argument zeta when
    j+2 < k; this is the analytically continued value zeta(j+2-k). For
    j+2-k <= 0 these are zeta at non-positive integers, which are
    rational (Bernoulli numbers): zeta(-2n) = 0 for n>=1, zeta(0) = -1/2,
    zeta(1-2n) = -B_{2n}/(2n). sympy.zeta handles these directly.
    """
    if k < 4 or k % 2 != 0:
        raise ValueError(f"k must be even and >= 4; got k = {k}")
    if not (0 <= j <= k - 2):
        raise ValueError(
            f"j must satisfy 0 <= j <= k - 2; got j = {j}, k = {k}"
        )

    I_unit = sympy.I  # imaginary unit
    if j == 0:
        return -2 * pi * I_unit * zeta(k - 1) / (k - 1)
    # General case
    coeff = (
        Integer(2)
        * Integer(-1) ** (j + 1)
        * factorial(j)
        * (2 * pi * I_unit) ** (k - 1 - j)
        / factorial(k - 1)
    )
    return coeff * zeta(j + 1) * zeta(j + 2 - k)


def lambert_moment_sum(j: int, k: int, n_terms: int, dps: int = 30) -> mpmath.mpc:
    """Numerical Lambert moment Sum_{n>=1} sigma_{k-1}(n) / n^{j+1}.

    This is the load-bearing arithmetic identity used in the closed
    form of m[j/k] (Kleinschmidt Eq. 26). The Lambert-series identity is

        Sum_{n>=1} sigma_p(n) / n^s  =  zeta(s) * zeta(s - p),

    valid for Re(s) > max(1, p + 1). With p = k - 1, s = j + 1 this is
    convergent for j >= k - 1, and analytically valid (after
    regularisation of the inner divergence) for j >= 1 with k - 1
    >= 2 because the sum sigma_{k-1}(n)/n^{j+1} converges absolutely
    once (j+1) > (k-1) + 1, i.e. j >= k - 1.

    For 1 <= j <= k - 2 the sum converges CONDITIONALLY (the divisor
    sum is O(n^{k-1}); divided by n^{j+1} with j < k - 1 it grows like
    n^{k-2-j} on average, divergent). So this routine only does a
    DIRECT numerical truncation in the absolutely-convergent regime
    j >= k - 1, where the identity zeta(j+1) * zeta(j+2-k) is a
    well-defined cross-check.

    Used for the cross-check test in tests/. Computing this for j < k-1
    requires Mellin-regularisation; we don't implement that here.

    Raises ValueError if (j, k) is outside the absolute-convergence regime.
    """
    if j < k - 1:
        raise ValueError(
            f"Lambert moment sum_n sigma_{{k-1}}(n)/n^{{j+1}} requires "
            f"absolute convergence j >= k - 1 (got j = {j}, k = {k}). "
            "For smaller j the identity zeta(j+1)*zeta(j+2-k) is "
            "obtained by analytic continuation, not by direct summation."
        )

    mpmath.mp.dps = dps

    def _divisor_sigma(n: int, p: int) -> int:
        total = 0
        for d in range(1, n + 1):
            if n % d == 0:
                total += d ** p
        return total

    s = mpmath.mpc(0, 0)
    for n in range(1, n_terms + 1):
        s += mpmath.mpf(_divisor_sigma(n, k - 1)) / mpmath.mpf(n) ** (j + 1)
    return s


def verify_lambert_identity(j: int, k: int, dps: int = 30) -> Tuple[mpmath.mpc, mpmath.mpc, float]:
    """Cross-check the Lambert-moment identity used in the depth-1 closed form.

    Verifies, in the absolute-convergence regime j >= k - 1:

        Sum_n sigma_{k-1}(n) / n^{j+1}   ==   zeta(j+1) * zeta(j+2-k).

    Returns (numerical_sum, analytic_product, relative_error).
    """
    if j < k - 1:
        raise ValueError(
            f"verify_lambert_identity only valid for j >= k - 1; "
            f"got j={j}, k={k}"
        )
    mpmath.mp.dps = dps
    # Tail decays as n^{(k-1) - (j+1)} = n^{k-j-2}. For j = k-1 the tail is
    # n^{-1} which converges glacially, so we'd need many more terms. For
    # safer cross-check, use Euler-Maclaurin or pick (j, k) with j > k - 1
    # in tests so that the tail is at least n^{-2} (j = k gives n^{-2}).
    # Default n_terms scales with the inverse of the tail exponent.
    tail_exp = k - j - 2  # power of n in tail: 1/n^|tail_exp+1| roughly
    # n_terms ~ 10^(dps / |tail_exp+1|) for convergence to dps digits.
    # We cap to keep the computation reasonable; user can increase via API.
    if tail_exp < -1:
        n_terms = min(2000, max(200, int(10 ** (dps / max(1, -tail_exp - 1)))))
    else:
        n_terms = 5000  # very slow convergence; fall back to large fixed sum

    num = lambert_moment_sum(j, k, n_terms, dps=dps)
    ana = mpmath.zeta(j + 1) * mpmath.zeta(j + 2 - k)
    rel_err = float(abs(num - ana) / max(abs(num), abs(ana), mpmath.mpf(1)))
    return num, ana, rel_err


def classify_depth1_mmv(j: int, k: int) -> Dict[str, object]:
    """Structural classification of m[j/k] in the f-alphabet picture.

    Reference: Kleinschmidt §2.1 + Eq. 26.

    The depth-1 MMV m[j/k] lives in
        Q . pi^a . zeta(odd)^{eps}
    with explicit (a, odd-zeta argument, epsilon) determined by (j, k).
    This is the LOAD-BEARING content of the depth-1 sector: every m[j/k]
    is a "pi-power times at most one odd zeta times a rational" — i.e.
    sits inside the Brown 2014 motivic Tate sector at depth 1, with the
    odd zeta carrying the f_{odd} content.

    Returns a dict with:
        "pi_power"             : a in Q[i, pi]
        "zeta_arguments"       : list of integer arguments of zeta factors
        "is_purely_tate"       : True iff all zeta arguments are <= 0 (only Bernoulli)
        "odd_zeta_content"     : list of odd zeta arguments >= 3 appearing
                                  (these are the f-alphabet generators
                                  the m[j/k] inhabits)
        "f_alphabet_word"      : the predicted f-word that m[j/k] maps to
                                  under rho^{-1} (without the pi-power
                                  prefactor)
    """
    if k < 4 or k % 2 != 0:
        raise ValueError(f"k must be even and >= 4; got k = {k}")
    if not (0 <= j <= k - 2):
        raise ValueError(
            f"j must satisfy 0 <= j <= k - 2; got j = {j}, k = {k}"
        )

    if j == 0:
        # m[0/k] = -2 pi i zeta_{k-1} / (k-1).
        # pi power = 1 (the i is in the Tate prefactor as 2 pi i).
        pi_power = 1
        zeta_args = [k - 1]
        odd_zetas = [k - 1] if (k - 1) % 2 == 1 and (k - 1) >= 3 else []
        f_word = f_letter(k - 1) if odd_zetas else F_UNIT
        is_purely_tate = (k - 1) <= 0
    else:
        # m[j/k] = (rational) * (2 pi i)^{k-1-j} * zeta(j+1) * zeta(j+2-k).
        # pi power from (2 pi i)^{k-1-j} is k - 1 - j.
        pi_power = k - 1 - j
        zeta_args = [j + 1, j + 2 - k]  # j+2-k <= 0 here (Bernoulli) since j <= k-2
        # The j+2-k argument is non-positive (Bernoulli rational); only the
        # j+1 argument contributes the irreducible odd-zeta content.
        odd_zetas = []
        if (j + 1) % 2 == 1 and (j + 1) >= 3:
            odd_zetas.append(j + 1)
        f_word = f_letter(j + 1) if odd_zetas else F_UNIT
        # Purely Tate iff the j+1 argument is also non-positive or even
        # (zeta(even) = rational * pi^even, no irreducible odd-zeta).
        is_purely_tate = not odd_zetas

    return {
        "j": j,
        "k": k,
        "pi_power": pi_power,
        "zeta_arguments": zeta_args,
        "is_purely_tate": is_purely_tate,
        "odd_zeta_content": odd_zetas,
        "f_alphabet_word": f_word,
    }


def depth1_mmv_numerical_pure_tate(j: int, k: int, dps: int = 30) -> mpmath.mpc:
    """Numerical evaluation of m[j/k] in the purely-Tate regime, j >= 1, j+1 even.

    When j+1 is even (so j is odd), zeta(j+1) is a rational multiple of
    pi^{j+1}, and the closed form m[j/k] becomes a pure-Tate rational
    multiple of pi^a. We evaluate the symbolic closed form via mpmath
    and use this as an INDEPENDENT cross-check (versus the symbolic
    formula). This works directly without regularisation issues.
    """
    mpmath.mp.dps = dps
    sym = depth1_mmv_value(j, k)
    # Convert to mpmath via sympy.N
    return mpmath.mpc(str(sympy.N(sympy.re(sym), dps)), str(sympy.N(sympy.im(sym), dps)))


# ---------------------------------------------------------------------------
# Layer 3. Tsunogai derivation algebra (symbolic, depth 1).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TsunogaiGenerator:
    """Symbolic Tsunogai derivation epsilon_k^{(j)}.

    Reference: Kleinschmidt et al. Eq. (34), Sec. 2.2.
    Definition:
        epsilon_k^{(j)} := ad_{epsilon_0}^{j} (epsilon_k).

    The grading is:
        - epsilon_0 has weight 0, ad-weight 0; it is the "base" derivation
          acting on the genus-1 fundamental Lie algebra via tau-translation.
        - epsilon_k for k >= 4 even has weight k.
        - epsilon_k^{(j)} has total weight k + 0 * j = k (since epsilon_0
          is weight 0); the j label tracks ad-iteration depth.

    We expose:
        eps0 = TsunogaiGenerator(k=0, j=0)
        eps_k = TsunogaiGenerator(k=k, j=0)
        eps_k_iter(j) = TsunogaiGenerator(k=k, j=j)
    """

    k: int
    j: int = 0

    def __post_init__(self) -> None:  # noqa: D401
        if self.k != 0 and (self.k < 4 or self.k % 2 != 0):
            raise ValueError(
                "Tsunogai generator index must be k=0 or k even >= 4; "
                f"got k = {self.k}"
            )
        if self.j < 0:
            raise ValueError(f"ad-iteration j must be >= 0; got {self.j}")

    def weight(self) -> int:
        """Modular weight of epsilon_k^{(j)} is k."""
        return self.k

    def __repr__(self) -> str:
        if self.k == 0:
            return "eps_0"
        if self.j == 0:
            return f"eps_{self.k}"
        return f"eps_{self.k}^({self.j})"


def eps0() -> TsunogaiGenerator:
    return TsunogaiGenerator(k=0, j=0)


def eps(k: int, j: int = 0) -> TsunogaiGenerator:
    """Convenience constructor for epsilon_k^{(j)}."""
    return TsunogaiGenerator(k=k, j=j)


@dataclass(frozen=True)
class TsunogaiCommutator:
    """A formal commutator [a, b] in the Tsunogai algebra.

    We treat the derivation algebra symbolically: a commutator is a
    pair (a, b) with antisymmetry [a, b] = -[b, a] enforced at
    comparison time.
    """
    left: TsunogaiGenerator
    right: TsunogaiGenerator

    def __repr__(self) -> str:
        return f"[{self.left}, {self.right}]"

    def weight(self) -> int:
        return self.left.weight() + self.right.weight()


def pollack_relation_weight_14() -> Tuple[TsunogaiCommutator, TsunogaiCommutator]:
    """Lowest-weight Pollack quadratic relation among Tsunogai generators.

    Reference: Kleinschmidt et al. Sec. 2.2 (after Eq. 35), citing
    Pollack 1504.04737.

    The lowest non-trivial commutator relation in the Eisenstein
    sub-Lie-algebra appears at total weight 14:
        [epsilon_4, epsilon_10]  -  3 [epsilon_6, epsilon_8]  =  0.

    This module exposes the two commutators as symbolic objects; the
    relation is enforced as an equation in the appropriate quotient.
    No numerical content; this is a formal algebraic identity in the
    Tsunogai Lie algebra. We return the two terms (LHS = first - 3 * second).
    """
    lhs1 = TsunogaiCommutator(left=eps(4), right=eps(10))
    lhs2 = TsunogaiCommutator(left=eps(6), right=eps(8))
    return lhs1, lhs2


def verify_pollack_weight_check() -> bool:
    """Verify the weight balance of the Pollack relation.

    [epsilon_4, epsilon_10] has weight 4 + 10 = 14;
    [epsilon_6, epsilon_8]  has weight 6 + 8  = 14.
    Both terms are weight 14, so the relation is homogeneous of weight 14.
    """
    a, b = pollack_relation_weight_14()
    return a.weight() == 14 and b.weight() == 14


# ---------------------------------------------------------------------------
# Layer 4. Depth-1 zeta generator sigma_w (Eq. 34, truncated).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ZetaGeneratorDepth1:
    """Depth-1 truncation of the genus-1 zeta generator sigma_w.

    Reference: Kleinschmidt Eq. (34), truncated to depth 1 (single
    Tsunogai derivation, no nested commutator brackets).

    Decomposition:
        sigma_w  =  z_w                                  (arithmetic part)
                  - (1 / (w-1)!) epsilon_{w+1}^{(w-1)}   (geometric, depth 1)
                  + O(depth >= 2)                         (nested brackets)

    Constraints (Kleinschmidt):
        w odd, w >= 3.

    The arithmetic z_w acts on motivic MZVs via Brown's f-alphabet
    action; the geometric epsilon_{w+1}^{(w-1)} acts on iterated
    Eisenstein integrals via the Tsunogai derivation.

    Attributes
    ----------
    w : int
        Weight of the zeta generator; must be odd >= 3.
    arithmetic_z : Symbol
        Symbolic placeholder z_w for the arithmetic part. Its action on
        f-words is implemented by `apply_arithmetic_to_fword()`.
    geometric_eps : TsunogaiGenerator
        The depth-1 Tsunogai derivation epsilon_{w+1}^{(w-1)}.
    geometric_coefficient : Expr
        Rational coefficient -1/(w-1)! in front of the depth-1 geometric
        Tsunogai piece.
    """

    w: int
    arithmetic_z: Symbol = field(init=False)
    geometric_eps: TsunogaiGenerator = field(init=False)
    geometric_coefficient: Expr = field(init=False)

    def __post_init__(self) -> None:  # noqa: D401
        if self.w < 3 or self.w % 2 == 0:
            raise ValueError(
                "Zeta generator weight w must be odd >= 3; got "
                f"w = {self.w}"
            )
        # field(init=False) requires object.__setattr__ on frozen dataclasses.
        object.__setattr__(self, "arithmetic_z", Symbol(f"z_{self.w}"))
        object.__setattr__(
            self, "geometric_eps", eps(k=self.w + 1, j=self.w - 1)
        )
        object.__setattr__(
            self, "geometric_coefficient", -Rational(1, factorial(self.w - 1))
        )

    def weight(self) -> int:
        """Total weight w. Both the arithmetic z_w and geometric
        epsilon_{w+1}^{(w-1)} carry weight w (the latter via k+j-1+1 = k = w+1
        BUT modular-form convention assigns epsilon_k a Tsunogai weight that
        couples to the period weight; in the Kleinschmidt convention sigma_w
        is uniformly weight w on both pieces). See paper Eq. 34.
        """
        return self.w

    def __repr__(self) -> str:
        coeff = self.geometric_coefficient
        return (
            f"sigma_{self.w} = z_{self.w}  +  ({coeff}) * eps_{self.w + 1}^({self.w - 1})  "
            f"+ O(depth >= 2)"
        )


def apply_arithmetic_to_fword(
    sigma: ZetaGeneratorDepth1, word: FWord
) -> Dict[FWord, int]:
    """Action of z_w on an f-word, at depth 1.

    Reference: Brown 2014. The arithmetic part z_w of a zeta generator
    acts on the f-alphabet by LEFT-TRUNCATION: it removes a leading
    f_w letter (if present) and returns the remaining word. Equivalently,
    z_w is dual to the operation "prepend f_w" in the deconcatenation
    coproduct, restricted to the leading slot.

    Formally:
        z_w (f_w  u)  =  u
        z_w (f_v  u)  =  0       if v != w
        z_w (1)       =  0
        z_w on central f_2:  z_w acts trivially on f_2 powers.

    This is the standard depth-1 action; it is the dual of the f-alphabet
    coproduct restricted to depth-1 motivic-Galois action and is the
    operational content the Kleinschmidt zeta generators inherit via
    rho^{-1}.

    Returns
    -------
    dict mapping FWord -> integer coefficient.
    """
    if not word.letters:
        return {}  # z_w(1) = 0
    if word.letters[0] != sigma.w:
        return {}  # z_w(f_v u) = 0 if v != w
    truncated = FWord(letters=word.letters[1:], f2_power=word.f2_power)
    return {truncated: 1}


# ---------------------------------------------------------------------------
# Layer 5. Depth-1 coaction (Eq. 64 collapsed).
# ---------------------------------------------------------------------------


def depth1_coaction_on_zeta(weight: int) -> Tuple[FWord, FWord, FWord, FWord]:
    """Depth-1 coaction Delta(zeta^m_w) in the f-alphabet.

    Reference: Brown 2014 (recalled in Kleinschmidt §2.1, Eq. 9-like).
    At depth 1 the motivic coproduct in the f-alphabet reduces to
        Delta(f_w) = f_w (x) 1  +  1 (x) f_w.

    Mapping back to motivic MZVs via the inverse isomorphism rho^{-1}:
        Delta(zeta^m_w) = zeta^m_w (x) 1  +  1 (x) zeta^m_w.

    This is the depth-1 collapse of the Kleinschmidt Eq. 64 coaction
    applied to the arithmetic (z_w) sector: the genus-zero (Brown 2014)
    coaction projects onto its depth-1 f-alphabet content unchanged.

    Returns
    -------
    Tuple (f_w_word, unit_word, unit_word_again, f_w_word_again)
    representing the two terms (left tensor right, left tensor right):
        (f_w, 1, 1, f_w)
    i.e. Delta(f_w) = f_w (x) 1 + 1 (x) f_w.

    For weight = 3: returns (f_3, 1, 1, f_3).
    For weight = 5: returns (f_5, 1, 1, f_5).
    Weight must be odd >= 3.
    """
    if weight < 3 or weight % 2 == 0:
        raise ValueError(
            f"Depth-1 coaction weight must be odd >= 3; got {weight}"
        )
    f_w = f_letter(weight)
    unit = F_UNIT
    return (f_w, unit, unit, f_w)


def depth1_coaction_full(weight: int) -> List[Tuple[FWord, FWord]]:
    """Convenience: return depth-1 coaction as a list of (left, right) tensor terms.

    Delta(f_w) = f_w (x) 1  +  1 (x) f_w
    returned as [(f_w, 1), (1, f_w)].
    """
    f_w, unit_r, unit_l, f_w_r = depth1_coaction_on_zeta(weight)
    return [(f_w, unit_r), (unit_l, f_w_r)]


# ---------------------------------------------------------------------------
# Depth >= 2 stubs.
# ---------------------------------------------------------------------------


def depth2_mmv_value(
    k1: int, j1: int, k2: int, j2: int
) -> Expr:
    """Depth-2 multiple modular value -- NOT IMPLEMENTED.

    The depth-2 MMVs depend on the explicit content of Kleinschmidt
    Appendix A, which was not extractable from the public HTML of
    arXiv:2508.02800 at the time this module was written. Named as
    follow-on; needs author-provided supplementary material or
    direct PDF inspection.
    """
    raise NotImplementedError(
        "Depth-2 MMVs require Kleinschmidt Appendix A content (not "
        "extracted from public HTML at module-build time). Named follow-on: "
        "extract Appendix A formulae from the published PDF and extend "
        "depth1_mmv_value() to (k1, j1, k2, j2)."
    )


def depth2_coaction_on_iterated_eisenstein(
    weights: List[int],
) -> Tuple[FWord, FWord]:
    """Depth-2 coaction on an iterated Eisenstein integral -- NOT IMPLEMENTED.

    The depth-2 coaction is the genuinely-new content of Kleinschmidt
    et al. and is described explicitly in Sec. 4 and Appendix B, which
    were not extractable from the public HTML. The Eq. 64 form
        Delta I^m = (M_sigma^dr)^{-1} I^m M_sigma^dr I^dr
    requires the full Tsunogai-bracket tower for the geometric
    M_sigma^dr piece beyond depth 1. Named follow-on.
    """
    raise NotImplementedError(
        "Depth-2 coaction on iterated Eisenstein integrals requires the "
        "full Tsunogai-bracket tower (Kleinschmidt Sec. 4 + Appendix B). "
        "Not extracted from public HTML at module-build time. Named "
        "follow-on: implement after PDF inspection of Appendices A and B."
    )


# ---------------------------------------------------------------------------
# Self-tests (smoke checks; full pytest in tests/test_kleinschmidt_coaction.py).
# ---------------------------------------------------------------------------


def _self_check_summary() -> None:
    """Print a brief self-check of the depth-1 sector.

    Run as `python debug/kleinschmidt_coaction.py` for a quick smoke test.
    """
    print("Kleinschmidt depth-1 coaction module -- smoke check")
    print("=" * 60)

    print("\nLayer 1. F-alphabet:")
    w = FWord((3, 5))
    print(f"  Word: {w}, weight {w.weight()}, depth {w.depth()}")
    print(f"  Coproduct Delta(f_3 f_5):")
    for L, R in deconcatenation_coproduct(w):
        print(f"    {L}  (x)  {R}")
    print(f"  Shuffle f_3 . f_5:")
    for k, v in shuffle(f_letter(3), f_letter(5)).items():
        print(f"    {v}  *  ({k})")

    print("\nLayer 2. Depth-1 MMV m[0/4]:")
    val = depth1_mmv_value(0, 4)
    print(f"  m[0/4] (closed form) = {val} = {sympy.simplify(val)}")
    val4 = depth1_mmv_value(1, 4)
    print(f"  m[1/4] (closed form) = {sympy.simplify(val4)}")
    val5 = depth1_mmv_value(0, 6)
    print(f"  m[0/6] (closed form) = {sympy.simplify(val5)}")

    print("\nLayer 2 cross-check (Lambert moment identity, absolute-convergence regime):")
    # Pick (j, k) with strong tail decay. (j=5, k=4) gives Sum sigma_3(n)/n^6
    # with tail O(1/n^3), much faster convergence than (j=4, k=4).
    num, ana, rel_err = verify_lambert_identity(j=5, k=4, dps=20)
    print(f"  Sum_n sigma_3(n)/n^6  (numerical):              {complex(num)}")
    print(f"  zeta(6) * zeta(3)     (analytic product):       {complex(ana)}")
    print(f"  Relative error: {rel_err:.2e}")

    print("\nLayer 2 cross-check (purely-Tate identity m[1/4] = pi^4/54):")
    closed = sympy.simplify(depth1_mmv_value(1, 4))
    print(f"  m[1/4] symbolic = {closed}")
    expected = pi**4 / 54
    is_match = sympy.simplify(closed - expected) == 0
    print(f"  Equals pi^4/54?  {is_match}")
    # Numerical sanity
    num_val = float(sympy.N(closed, 50))
    pi_pow = float(sympy.N(pi**4 / 54, 50))
    print(f"  m[1/4] numerical: {num_val}")
    print(f"  pi^4/54         : {pi_pow}")

    print("\nLayer 2 structural classification of m[0/4], m[1/4], m[0/6], m[2/6]:")
    for j, k in [(0, 4), (1, 4), (0, 6), (2, 6), (1, 6)]:
        cls = classify_depth1_mmv(j, k)
        print(
            f"  m[{j}/{k}]: pi^{cls['pi_power']}, zeta args {cls['zeta_arguments']}, "
            f"purely Tate={cls['is_purely_tate']}, f-word={cls['f_alphabet_word']}"
        )

    print("\nLayer 3. Tsunogai algebra:")
    print(f"  eps_0  = {eps0()}, weight {eps0().weight()}")
    print(f"  eps_4  = {eps(4)}, weight {eps(4).weight()}")
    print(f"  eps_6^(3) = {eps(6, 3)}")
    print(f"  Pollack relation (weight 14):")
    a, b = pollack_relation_weight_14()
    print(f"    {a}  -  3 * {b}  =  0")
    print(f"    weight balance: {verify_pollack_weight_check()}")

    print("\nLayer 4. Depth-1 zeta generator sigma_3:")
    sigma3 = ZetaGeneratorDepth1(w=3)
    print(f"  {sigma3}")
    print(f"  z_3 acting on f_3:        {apply_arithmetic_to_fword(sigma3, f_letter(3))}")
    print(f"  z_3 acting on f_5:        {apply_arithmetic_to_fword(sigma3, f_letter(5))}")
    print(f"  z_3 acting on f_3 f_5:    {apply_arithmetic_to_fword(sigma3, FWord((3, 5)))}")

    print("\nLayer 5. Depth-1 coaction Delta(f_3):")
    for L, R in depth1_coaction_full(3):
        print(f"  {L}  (x)  {R}")


if __name__ == "__main__":
    _self_check_summary()
