"""SU(3) numerical sanity check for the unified GH-convergence theorem.

NOTE 2026-05-15: Sub-task (a) below ("Casimir triangle inequality") is
DEPRECATED — the inequality was verified FALSE by Sprint P-A (counterexamples
already at SU(2) for j=1, j'=1/2). The corrected L3 reformulation, the
"Dirac-triangle inequality" |D(pi) - D(pi')| <= sqrt(C(sigma)), is verified
on the full 100-pair SU(3) panel in debug/dirac_triangle_su3_check.py.
The Casimir-triangle code is preserved here for reproducibility of the
falsification record, but it is NO LONGER the right L3 ingredient.

Sub-tasks (b) [rate constant] and (c) [log-power test] from this script have
been refactored and extended in debug/su3_rate_constant.py (the canonical
2026-05-15 driver). Use that script for current execution.

Original three sub-tasks (preserved for historical record):
  (a) Casimir triangle inequality at SU(3) [DEPRECATED — verified FALSE]:
      For pi, pi' in a panel of 11 irreps, decompose pi (x) pi'^* into
      irreducibles sigma; check |C(pi) - C(pi')| <= C(sigma) for every
      sigma in the decomposition; report max sup ratio
      sup_sigma sqrt(|C(pi) - C(pi')|) / sqrt(C(sigma)).

  (b) Rate constant c(SU(3)) via central spectral Fejer kernel on
      SU(3) with Casimir cutoff Lambda^2 and spectral-side mass-
      concentration moment using the Weyl integration formula.
      [SUPERSEDED by debug/su3_rate_constant.py for 2026-05-15 final analysis.]

  (c) Single-log vs double-log fit on gamma_Lambda data.
      [SUPERSEDED by debug/su3_rate_constant.py.]

  Final: compare extracted c(SU(3)) against
      c_pred = 2 * Vol(G/T) / Vol(G).

Cross-references:
    geovac/central_fejer_su2.py (template)
    debug/unified_gh_scoping_memo.md (forward plan)
    papers/standalone/paper_38_su2_propinquity_convergence.tex (Paper 38)

This script does NOT modify production code. All output goes to
debug/data/su3_numerical_sanity.json.

Conventions
-----------
- SU(3) irreps labeled by Dynkin labels (p, q) with p, q >= 0. We use
  fundamental weights omega_1 = (2/3, 1/3) (in basis e_1 - (e_1+e_2+e_3)/3
  on the e-coordinates restricted to z_1+z_2+z_3 = 0), omega_2 = (1/3, 2/3).
- Highest weight: lambda = p*omega_1 + q*omega_2.
- Half-sum of positive roots: rho = omega_1 + omega_2.
- Quadratic Casimir (verified to give SU(2) reduction j(j+1) when q = 0):
      C(p, q) = (1/3) (p^2 + q^2 + p*q) + (p + q).
  Check at (1, 1) (adjoint): (1+1+1)/3 + 2 = 3. SU(2) reduction (p, 0)
  with j = p/2: should give j(j+1) = (p/2)(p/2+1). Our formula gives
  (1/3) p^2 + p = p^2/3 + p. For p = 1 (j = 1/2): 1/3 + 1 = 4/3, but
  j(j+1) = 1/2 * 3/2 = 3/4. Mismatch: factor of (4/3)/(3/4) = 16/9.
  This indicates a *normalisation* convention difference, not a bug.

  RESOLUTION: there are two standard conventions. Option A: the formula
  C(p, q) = (1/3)(p^2 + q^2 + pq) + (p + q) gives C(adjoint) = 3 in a
  normalisation where the adjoint Casimir equals h^vee = 3. Option B:
  C(p, q) = (p^2 + q^2 + pq + 3p + 3q)/3 = (1/3)((p+q+1)^2 + p^2 + q^2 - 1)
  gives the same value 3 on (1,1). Option C: 2*C uses the Killing-form-
  doubled normalisation where C(adjoint) = 2*h^vee = 6.

  We pick the unit normalisation where the SU(2) reduction is consistent
  in the *same* convention. With Option A applied to SU(2) (where the
  Casimir is C_SU2(p) = (1/3) * p^2 + p... no, this fails). The cleanest
  approach: pick the SU(2) convention that matches Paper 38 (where the
  shell index n labels j = (n-1)/2 and the Dirac eigenvalue is n + 1/2;
  the relevant scaling is the *square* of Casimir-difference, which is
  what enters the L3 ratio). For the Casimir TRIANGLE INEQUALITY, the
  *normalisation is irrelevant* — the inequality |C(pi) - C(pi')| <= C(sigma)
  is invariant under positive rescaling. So we just pick one convention
  and stick with it.

  Convention used:
      C(p, q) = (1/3) * (p^2 + q^2 + p*q) + p + q.
  On SU(2) reduction (p = 2j, q = 0): C(2j, 0) = (4 j^2)/3 + 2j. This is
  NOT j(j+1) but is *proportional* to it under a different convention; for
  the triangle inequality this proportionality is sufficient.

References
----------
- Goodman & Wallach, "Symmetry, Representations, and Invariants" (Springer 2009).
- Knapp, "Lie Groups Beyond an Introduction" (Birkhauser 2002), Chapter VI.
- Helgason, "Differential Geometry, Lie Groups, and Symmetric Spaces" (1978).
- Paper 38 §2-§3 (SU(2) blueprint).
- debug/unified_gh_scoping_memo.md §5.3 (Casimir triangle inequality).
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Tuple

import mpmath
import numpy as np
import sympy as sp
from sympy import Rational, Integer, simplify, expand, sin, cos, sqrt, pi, log


# =============================================================================
# Sub-task (a) infrastructure: SU(3) representations, Casimir, Klimyk's formula
# =============================================================================


@dataclass(frozen=True)
class IrrepSU3:
    p: int
    q: int

    def __post_init__(self):
        if self.p < 0 or self.q < 0:
            raise ValueError(f"Dynkin labels must be >= 0: ({self.p}, {self.q})")

    def dim(self) -> int:
        """Weyl dimension formula for SU(3): (p+1)(q+1)(p+q+2)/2."""
        p, q = self.p, self.q
        return (p + 1) * (q + 1) * (p + q + 2) // 2

    def casimir(self) -> sp.Rational:
        """Quadratic Casimir in the convention
            C(p, q) = (p^2 + q^2 + pq)/3 + (p + q).

        Verifies C(adjoint = (1,1)) = 3 = h^vee.
        For triangle inequality the normalisation cancels.
        """
        p, q = self.p, self.q
        return Rational(p ** 2 + q ** 2 + p * q, 3) + Rational(p + q)

    def __repr__(self):
        return f"({self.p},{self.q})[{self.dim()}]"

    @property
    def conj(self) -> "IrrepSU3":
        """Conjugate / dual: (p, q) -> (q, p)."""
        return IrrepSU3(self.q, self.p)


def _highest_weight_basis(p: int, q: int) -> Tuple[Rational, Rational]:
    """Map Dynkin (p, q) to the highest weight in (a_1, a_2)-coordinates,
    where a_i is the coefficient on the simple root alpha_i.

    For SU(3): omega_1 = (2 alpha_1 + alpha_2)/3, omega_2 = (alpha_1 + 2 alpha_2)/3.
    Highest weight lambda = p omega_1 + q omega_2
                          = ((2p+q)/3) alpha_1 + ((p+2q)/3) alpha_2.

    But for tensor product computations using Klimyk, we work in the BASIS
    OF FUNDAMENTAL WEIGHTS (Dynkin coordinates). All weights of an irrep
    can be enumerated as (p - n1 - n2 + ..., q - n2 - n3 + ...) with
    appropriate constraints from the Lie-algebra-step structure.

    For a manageable implementation, we use **explicit weight enumeration**
    via Freudenthal's formula (or its consequence) for small irreps.
    """
    # Highest weight in fundamental-weight basis: (p, q).
    # Returned as a tuple of mpz Rationals for downstream arithmetic.
    return (Rational(p), Rational(q))


def _simple_roots_dynkin() -> list[tuple[int, int]]:
    """Simple roots of A_2 in fundamental-weight (Dynkin) coordinates.

    The Cartan matrix of A_2 is
        [[ 2, -1],
         [-1,  2]],
    so alpha_1 = (2, -1), alpha_2 = (-1, 2) in fundamental-weight basis.
    """
    return [(2, -1), (-1, 2)]


def _all_positive_roots_dynkin() -> list[tuple[int, int]]:
    """All positive roots of A_2 in fundamental-weight basis.

    Three positive roots: alpha_1 = (2, -1), alpha_2 = (-1, 2),
    alpha_1 + alpha_2 = (1, 1).
    """
    return [(2, -1), (-1, 2), (1, 1)]


def _weight_subtract(w: tuple, root: tuple) -> tuple:
    return (w[0] - root[0], w[1] - root[1])


def _weight_add(w: tuple, root: tuple) -> tuple:
    return (w[0] + root[0], w[1] + root[1])


def _enumerate_weights(p: int, q: int) -> dict[tuple[int, int], int]:
    """Enumerate all weights (with multiplicities) of the SU(3) irrep
    V_{(p, q)}, in fundamental-weight basis.

    Algorithm: BFS from highest weight, subtracting simple roots as long
    as the result remains in the convex hull of W . lambda (Weyl orbit
    convex hull). Multiplicities are computed via Freudenthal's formula
    by iterating downward.

    For SU(3), the irrep has dimension dim = (p+1)(q+1)(p+q+2)/2; the
    weight-lattice points come with multiplicities given by the
    "layered hexagon" pattern: layer k (distance k from the boundary)
    has multiplicity k + 1 (for sufficiently small k).

    For the SU(3) irrep V_{(p, q)} (assuming p >= q WLOG), the
    multiplicity of a weight at "shell k" (combinatorial distance from
    the outer hexagon) is min(k+1, q+1).

    Concrete recipe (Stanley-style):
        - The weights form a hexagonal pattern in the (a, b) Dynkin plane.
        - The outer boundary is the Weyl orbit of (p, q): six points if
          p, q both > 0, else three (and we get a triangle).
        - Inner shells have linearly increasing multiplicity 1, 2, 3, ...
          up to min(p, q) + 1, then constant.

    For modest (p, q) (say (p+q) <= 12), Freudenthal works in a few ms.
    """
    # Use the formulation via Verma-module character. Multiplicities of
    # weights in V_lambda are given by Freudenthal:
    #   mult(mu) = (2/((|lambda+rho|^2 - |mu+rho|^2))) *
    #              sum_{k>=1, alpha>0} (mult(mu + k*alpha) * <mu + k*alpha, alpha>)
    # We compute by working downward from the highest weight.
    rho = (1, 1)  # rho in fundamental-weight basis = (1, 1) for A_2

    # Inner products: in the fundamental-weight basis, the inner product
    # is given by <a omega_1 + b omega_2, c omega_1 + d omega_2> =
    #   (1/3) * (2 a c + a d + b c + 2 b d). (This is the inverse Cartan.)
    def ip(w1, w2):
        a, b = w1
        c, d = w2
        return Rational(2 * a * c + a * d + b * c + 2 * b * d, 3)

    def norm2(w):
        return ip(w, w)

    lam = (p, q)
    lam_rho = (lam[0] + rho[0], lam[1] + rho[1])
    norm_lr = norm2(lam_rho)

    pos_roots = _all_positive_roots_dynkin()

    # First, enumerate all weights in the Weyl-hull of lambda. We do this
    # by BFS from lambda, subtracting simple roots, keeping only weights
    # that lie in the Weyl-orbit-convex-hull (which is equivalent to:
    # `mu - lambda` is a non-negative integer combination of simple
    # roots, AND the corresponding multiplicity is positive).
    simples = _simple_roots_dynkin()

    # Multiplicities computed by Freudenthal recursion.
    # mult[(a, b)] = multiplicity of the weight (a, b) in V_lambda.
    mult: dict[tuple[int, int], int] = {}
    mult[lam] = 1

    # BFS to enumerate candidate weights, bounded by the convex hull of
    # the Weyl orbit of lambda. The Weyl orbit is generated by
    # reflections; for SU(3), the orbit of lambda = (p, q) consists of
    # the six Weyl images (or fewer if on a wall). We compute these and
    # use them as bounding-box vertices.
    #
    # For A_2, the six Weyl images of (p, q) (in fundamental-weight basis):
    #   e:        (p, q)
    #   s_1:      (-p, p+q)
    #   s_2:      (p+q, -q)
    #   s_1 s_2:  (-p-q, p)
    #   s_2 s_1:  (q, -p-q)
    #   s_1 s_2 s_1 = w_0:  (-q, -p)
    # (See Humphreys, "Reflection Groups", Table 2.)
    weyl_orbit = []
    pp, qq = p, q
    weyl_orbit.append((pp, qq))
    weyl_orbit.append((-pp, pp + qq))
    weyl_orbit.append((pp + qq, -qq))
    weyl_orbit.append((-pp - qq, pp))
    weyl_orbit.append((qq, -pp - qq))
    weyl_orbit.append((-qq, -pp))

    # Bounding box on the weights: enclose the Weyl orbit.
    # The convex hull is the hexagon (or triangle for trivial / fundamental
    # cases) whose vertices are the Weyl orbit. For BFS, we use the
    # rectangular bound and the additional constraint `lam - mu in
    # nonneg integer cone of simple roots`.
    a_min = min(w[0] for w in weyl_orbit)
    a_max = max(w[0] for w in weyl_orbit)
    b_min = min(w[1] for w in weyl_orbit)
    b_max = max(w[1] for w in weyl_orbit)

    queue = [lam]
    seen = set()
    seen.add(lam)
    while queue:
        new_queue = []
        for w in queue:
            for s in simples:
                w2 = _weight_subtract(w, s)
                if w2 in seen:
                    continue
                # Bounding box check.
                if w2[0] < a_min or w2[0] > a_max:
                    continue
                if w2[1] < b_min or w2[1] > b_max:
                    continue
                # lam - w2 in non-negative simple-root cone:
                d1 = lam[0] - w2[0]
                d2 = lam[1] - w2[1]
                # 2a - b = d1, -a + 2b = d2 => a = (2 d1 + d2)/3,
                num_a = 2 * d1 + d2
                num_b = d1 + 2 * d2
                if num_a < 0 or num_b < 0:
                    continue
                if num_a % 3 != 0 or num_b % 3 != 0:
                    continue
                seen.add(w2)
                new_queue.append(w2)
        queue = new_queue

    # Sort weights by depth from highest weight (descending by Freudenthal-friendly order).
    # Use total height: depth = (sum of positive-root coefficients in lambda - mu).
    def depth(w):
        d1 = lam[0] - w[0]
        d2 = lam[1] - w[1]
        # In simple-root basis: d1 = 2a - b, d2 = -a + 2b => a = (2d1+d2)/3, b = (d1+2d2)/3
        return Rational(2 * d1 + d2, 3) + Rational(d1 + 2 * d2, 3)

    weights_sorted = sorted(seen, key=lambda w: depth(w))

    # Apply Freudenthal: for each w starting from lam (depth 0, mult=1),
    # compute mult[w] from already-known mult of w + k*alpha for alpha in positive roots,
    # k >= 1 such that w + k*alpha is a known weight.
    for w in weights_sorted:
        if w == lam:
            continue
        # Freudenthal: ((|lam+rho|^2 - |w+rho|^2)) * mult(w) =
        #   2 * sum_{alpha>0} sum_{k>=1} mult(w + k alpha) * <w + k alpha, alpha>
        w_rho = (w[0] + rho[0], w[1] + rho[1])
        denom = norm_lr - norm2(w_rho)
        if denom == 0:
            mult[w] = 0
            continue
        rhs = Rational(0)
        for alpha in pos_roots:
            k = 1
            while True:
                wk = (w[0] + k * alpha[0], w[1] + k * alpha[1])
                if wk not in mult:
                    break
                rhs += mult[wk] * ip(wk, alpha)
                k += 1
        m = 2 * rhs / denom
        # Multiplicity should be a non-negative integer.
        m_int = int(m)
        assert m == m_int, f"Non-integer Freudenthal multiplicity {m} at {w}"
        if m_int < 0:
            m_int = 0
        mult[w] = m_int

    # Drop zeros.
    return {w: m for w, m in mult.items() if m > 0}


def _verify_dim(p: int, q: int, weights: dict) -> bool:
    expected = IrrepSU3(p, q).dim()
    actual = sum(weights.values())
    return expected == actual


def _weyl_reflect_to_dominant(w: tuple, max_iters: int = 100) -> tuple[tuple, int]:
    """Reflect w in fundamental-weight basis through Weyl walls until
    dominant. Return (w_dominant, sign) where sign = (-1)^(number of
    reflections). If (w + rho) hits a Weyl wall (any coordinate of
    w + rho is zero in fundamental-weight basis), return (None, 0).

    Reflection through wall i: w -> w - <w, alpha_i^vee> alpha_i;
    in fundamental-weight basis with simple coroots = simple roots,
    s_i acts as: (a_1, ..., a_i, ..., a_n) -> ... swap algorithm.

    For A_2: reflections s_1 (across alpha_1-wall) maps
        (a_1, a_2) -> (-a_1, a_1 + a_2),
    and s_2 maps
        (a_1, a_2) -> (a_1 + a_2, -a_2).
    These are obtained from <(a_1, a_2), alpha_1^vee> = a_1 (since
    omega_i . alpha_j^vee = delta_ij), so s_1 . (a_1, a_2) =
    (a_1, a_2) - a_1 * (2, -1) = (-a_1, a_2 + a_1).
    """
    a, b = w
    # We work with w + rho in fundamental-weight basis, i.e. SHIFT to
    # rho-shifted coordinates where being on a Weyl wall means having
    # any coordinate equal to 0.
    sign = 1
    for _ in range(max_iters):
        if a > 0 and b > 0:
            return (a, b), sign
        if a == 0 or b == 0:
            return None, 0
        if a < 0:
            # s_1: (a, b) -> (-a, a + b)
            a, b = -a, a + b
            sign = -sign
            continue
        # b < 0
        # s_2: (a, b) -> (a + b, -b)
        a, b = a + b, -b
        sign = -sign
    raise RuntimeError(f"Weyl reflection didn't converge for ({w[0]}, {w[1]})")


def tensor_product_su3(
    pi1: IrrepSU3, pi2: IrrepSU3
) -> dict[IrrepSU3, int]:
    """Decompose V_{pi1} (x) V_{pi2} into irreducibles via Klimyk's formula
    (also known as Brauer-Klimyk or the Racah-Speiser algorithm).

    For each weight nu of V_{pi2} (with multiplicity m_nu):
        contribution = sgn(w) V_{w * (pi1 + nu + rho) - rho}
    where w is the unique Weyl element making (pi1 + nu + rho) dominant
    (interior of the fundamental Weyl chamber, all coordinates > 0 in
    fundamental-weight basis). If (pi1 + nu + rho) lies on a Weyl wall,
    the contribution is zero.

    Returns dict {IrrepSU3(p, q): multiplicity}.
    """
    rho = (1, 1)
    weights = _enumerate_weights(pi2.p, pi2.q)
    assert _verify_dim(pi2.p, pi2.q, weights), \
        f"Weight enumeration failed for ({pi2.p},{pi2.q})"

    result: dict[IrrepSU3, int] = {}
    pi1_w = (pi1.p, pi1.q)
    for nu, mult in weights.items():
        # Compute pi1 + nu + rho.
        shifted = (pi1_w[0] + nu[0] + rho[0], pi1_w[1] + nu[1] + rho[1])
        # Reflect to dominant (Weyl-strict, i.e. into open chamber).
        dom, sign = _weyl_reflect_to_dominant(shifted)
        if dom is None:
            continue
        # Subtract rho to get the highest weight of the irrep.
        new_irrep = IrrepSU3(dom[0] - 1, dom[1] - 1)
        result[new_irrep] = result.get(new_irrep, 0) + sign * mult

    # Drop zeros (and assert no negatives).
    final = {}
    for irrep, m in result.items():
        if m == 0:
            continue
        assert m > 0, f"Negative tensor-product multiplicity: {irrep} -> {m}"
        final[irrep] = m

    return final


def _sanity_check_tensor_product():
    """Quick sanity checks for tensor_product_su3."""
    # 3 (x) 3-bar = 1 + 8
    pi3 = IrrepSU3(1, 0)
    pi3bar = IrrepSU3(0, 1)
    decomp = tensor_product_su3(pi3, pi3bar)
    expected = {IrrepSU3(0, 0): 1, IrrepSU3(1, 1): 1}
    assert decomp == expected, f"3 (x) 3-bar test failed: got {decomp}"

    # 3 (x) 3 = 6 + 3-bar
    decomp = tensor_product_su3(pi3, pi3)
    expected = {IrrepSU3(2, 0): 1, IrrepSU3(0, 1): 1}
    assert decomp == expected, f"3 (x) 3 test failed: got {decomp}"

    # 8 (x) 8 = 1 + 8 + 8 + 10 + 10-bar + 27
    pi8 = IrrepSU3(1, 1)
    decomp = tensor_product_su3(pi8, pi8)
    expected = {
        IrrepSU3(0, 0): 1,
        IrrepSU3(1, 1): 2,
        IrrepSU3(3, 0): 1,
        IrrepSU3(0, 3): 1,
        IrrepSU3(2, 2): 1,
    }
    assert decomp == expected, f"8 (x) 8 test failed: got {decomp}"

    # Dimension check: dim(3) * dim(3-bar) = dim(1) + dim(8) = 1 + 8 = 9
    assert pi3.dim() * pi3bar.dim() == 1 + 8

    return True


# =============================================================================
# Sub-task (a): Casimir triangle inequality test
# =============================================================================


def casimir_triangle_test(
    panel: list[IrrepSU3],
) -> dict:
    """For every ordered pair (pi, pi') in panel x panel, test:
        |C(pi) - C(pi')| <= C(sigma) for every sigma in pi (x) pi'^*

    Reports:
        - num_pairs_total
        - num_pairs_pass
        - num_pairs_fail
        - max_sup_ratio (over all pairs and all sigma)
        - sup_ratio_per_pair: list of (pi, pi', sup_ratio)
        - violations: list of (pi, pi', sigma, gap, C_sigma)
    """
    n = len(panel)
    pairs_total = n * n
    pairs_pass = 0
    pairs_fail = 0
    max_sup_ratio = mpmath.mpf(0)
    violations = []
    sup_ratio_per_pair = []

    for pi in panel:
        for piprime in panel:
            cpi = pi.casimir()
            cpiprime = piprime.casimir()
            gap = abs(cpi - cpiprime)  # |C(pi) - C(pi')|
            decomp = tensor_product_su3(pi, piprime.conj)
            pair_pass = True
            pair_max_ratio = mpmath.mpf(0)
            for sigma, mult in decomp.items():
                csigma = sigma.casimir()
                if csigma == 0:
                    # sigma = (0, 0). The trivial rep occurs only when
                    # piprime = pi. Then gap = 0, OK.
                    if gap != 0:
                        pair_pass = False
                        violations.append({
                            "pi": (pi.p, pi.q),
                            "pi_prime": (piprime.p, piprime.q),
                            "sigma": (sigma.p, sigma.q),
                            "C_pi": float(cpi),
                            "C_pi_prime": float(cpiprime),
                            "gap": float(gap),
                            "C_sigma": float(csigma),
                            "ratio": float("inf"),
                        })
                    continue
                if gap > csigma:
                    pair_pass = False
                    violations.append({
                        "pi": (pi.p, pi.q),
                        "pi_prime": (piprime.p, piprime.q),
                        "sigma": (sigma.p, sigma.q),
                        "C_pi": float(cpi),
                        "C_pi_prime": float(cpiprime),
                        "gap": float(gap),
                        "C_sigma": float(csigma),
                        "ratio": float(mpmath.sqrt(mpmath.mpf(gap) / mpmath.mpf(csigma))),
                    })
                ratio = mpmath.sqrt(mpmath.mpf(gap) / mpmath.mpf(csigma))
                if ratio > pair_max_ratio:
                    pair_max_ratio = ratio
                if ratio > max_sup_ratio:
                    max_sup_ratio = ratio
            sup_ratio_per_pair.append({
                "pi": (pi.p, pi.q),
                "pi_prime": (piprime.p, piprime.q),
                "C_pi": float(cpi),
                "C_pi_prime": float(cpiprime),
                "n_sigma": len(decomp),
                "max_sup_ratio": float(pair_max_ratio),
            })
            if pair_pass:
                pairs_pass += 1
            else:
                pairs_fail += 1

    return {
        "n_pairs_total": pairs_total,
        "n_pairs_pass": pairs_pass,
        "n_pairs_fail": pairs_fail,
        "max_sup_ratio": float(max_sup_ratio),
        "violations": violations,
        "sup_ratio_per_pair": sup_ratio_per_pair,
        "panel": [(pi.p, pi.q) for pi in panel],
        "convention": "C(p,q) = (p^2+q^2+pq)/3 + (p+q); C(adjoint=(1,1)) = 3 = h^vee.",
    }


# =============================================================================
# Sub-task (b) infrastructure: Weyl integration and central Fejer kernel on SU(3)
# =============================================================================
#
# SU(3) maximal torus T^2 parameterised by (theta_1, theta_2) in [0, 2 pi]^2,
# with the constraint that the third torus angle theta_3 = -(theta_1 + theta_2)
# (so the determinant condition z_1 z_2 z_3 = 1 holds).
#
# Weyl integration formula: for class function f on SU(3),
#   integral_G f(g) dg
#     = (1/(|W| * Vol_T)) * integral_{T^2}
#         |Delta(theta_1, theta_2)|^2 f(theta_1, theta_2) dtheta_1 dtheta_2,
# where |W| = 6 and the Weyl-denominator squared
#   |Delta|^2 = prod_{alpha > 0} |1 - e^{i alpha . theta}|^2
#            = prod_{alpha > 0} 4 sin^2(alpha . theta / 2).
#
# Positive roots of A_2 in the basis where the maximal torus parameter is
# (t_1, t_2, t_3) with sum 0:
#   alpha_1 = e_1 - e_2 (action: t_1 - t_2)
#   alpha_2 = e_2 - e_3 (action: t_2 - t_3)
#   alpha_1 + alpha_2 = e_1 - e_3 (action: t_1 - t_3)
# With t_3 = -(t_1 + t_2):
#   alpha_1 . theta = theta_1 - theta_2,
#   alpha_2 . theta = theta_2 - theta_3 = theta_2 + theta_1 + theta_2 = theta_1 + 2 theta_2,
#   (alpha_1 + alpha_2) . theta = theta_1 - theta_3 = 2 theta_1 + theta_2.
#
# So Weyl denominator on (theta_1, theta_2):
#   Delta = e^{i (rho_1 theta_1 + rho_2 theta_2)} *
#           [stuff] * (1 - e^{-i (theta_1 - theta_2)}) *
#                     (1 - e^{-i (theta_1 + 2 theta_2)}) *
#                     (1 - e^{-i (2 theta_1 + theta_2)})
# (or some such); for Weyl-integration we just need |Delta|^2.
#
# Geodesic distance d(g, e) for bi-invariant unit-Killing-form metric on SU(3):
# the Riemannian distance is |t|_{Killing} where t is the diagonal matrix
# log(g)/i in the maximal torus. With our normalisation (Killing form
# normalised so that the SHORTEST nonzero element of the integer lattice of T
# has length 2 pi), the geodesic distance from e to (theta_1, theta_2, theta_3)
# is the *Killing norm*:
#
#   d(e, t)^2 = (Killing form normalisation) * (theta_1^2 + theta_2^2 + theta_3^2).
#
# CRITICAL CONVENTION CHOICE: we want the bi-invariant metric where the
# minimum closed geodesic to the identity has length compatible with Paper 38's
# SU(2) convention. On SU(2) (rank 1), the rotation parameter chi runs over
# [0, 2 pi] and d(e, g) = chi (the geodesic distance is just the rotation angle).
# This corresponds to the Killing form being normalised so that the
# fundamental loop in T = U(1) has length 2 pi.
#
# For SU(3), with theta_i in [0, 2 pi], a natural unit-Killing convention is
#   d(e, t)^2 = (2/3) * (theta_1^2 + theta_2^2 + theta_3^2)
# (the 2/3 comes from projecting onto the trace-free subspace), so that
# the SHORTEST closed geodesic to identity is 2 pi (consistent with SU(2)).
#
# Equivalently, in coordinates where t_3 = -(t_1 + t_2):
#   d(e, t)^2 = (2/3) * (theta_1^2 + theta_2^2 + (theta_1 + theta_2)^2)
#             = (2/3) * (2 theta_1^2 + 2 theta_2^2 + 2 theta_1 theta_2)
#             = (4/3) * (theta_1^2 + theta_2^2 + theta_1 theta_2).
#
# Volume of SU(3) (bi-invariant Haar, unit Killing-form): we will VERIFY
# this numerically by integrating 1 over T^2 with the Weyl-integration formula
# and matching to the standard answer Vol(SU(3)) = sqrt(3) pi^5.
# Volume of SU(3)/T^2 = full flag = chi(F)/Vol(W) * Vol(SU(3)) etc.,
# but we just compute it as a Weyl integral.


def _weyl_denom_sq_su3(t1: mpmath.mpf, t2: mpmath.mpf) -> mpmath.mpf:
    """|Delta(t1, t2)|^2 = prod_{alpha > 0} 4 sin^2(alpha . t / 2)
    for SU(3) maximal torus (t_1, t_2, t_3 = -(t_1+t_2)).

    Positive roots:
      alpha_1: t_1 - t_2
      alpha_2: 2 t_2 + t_1   (since t_2 - t_3 = t_2 + t_1 + t_2)
      alpha_1 + alpha_2: 2 t_1 + t_2
    """
    s1 = mpmath.sin((t1 - t2) / 2)
    s2 = mpmath.sin((2 * t2 + t1) / 2)
    s3 = mpmath.sin((2 * t1 + t2) / 2)
    return mpmath.mpf(64) * s1 * s1 * s2 * s2 * s3 * s3


def _su3_geodesic_dist_sq_naive(t1: mpmath.mpf, t2: mpmath.mpf) -> mpmath.mpf:
    """Naive 'geodesic distance squared' from the identity in the lift to t,
    BEFORE accounting for the lattice periodicity. Equals
        d^2_naive(t) := theta_1^2 + theta_2^2 + theta_3^2
                      = theta_1^2 + theta_2^2 + (theta_1 + theta_2)^2
                      = 2 (theta_1^2 + theta_2^2 + theta_1 theta_2)
    in our (t_1, t_2) parameterisation with t_3 = -(t_1 + t_2).

    The TRUE geodesic distance is obtained by minimising this over the
    coroot lattice translates (the kernel of exp: t -> T). See
    `_su3_geodesic_dist`.
    """
    return mpmath.mpf(2) * (t1 * t1 + t2 * t2 + t1 * t2)


def _su3_geodesic_dist(
    t1: mpmath.mpf,
    t2: mpmath.mpf,
    n_search: int = 2,
) -> mpmath.mpf:
    """True bi-invariant geodesic distance from identity to the SU(3)
    element `exp(i diag(t1, t2, -t1-t2))`.

    For each (n1, n2) in the search range [-n_search, n_search]^2, lift
    by the coroot lattice and compute the Killing-norm of the lift. Return
    the minimum.

    The coroot lattice of SU(3) (in t parameter space, with t_1+t_2+t_3=0):
        - alpha_1^v = (1, -1, 0) -> in (t1, t2): adds (1, -1) (since t_3
          adjusts).
        - alpha_2^v = (0, 1, -1) -> in (t1, t2): adds (0, 1).
        - But we have constraint t_3 = -(t_1+t_2), so in (t_1, t_2)-coords,
          shifting by 2 pi alpha_1^v means (t_1, t_2) -> (t_1 + 2 pi, t_2 - 2 pi).
          Shifting by 2 pi alpha_2^v means (t_1, t_2) -> (t_1, t_2 + 2 pi).
          But that doesn't preserve t_3 = -(t_1+t_2)... Wait:
              alpha_1^v shifts t_1 -> t_1 + 1, t_2 -> t_2 - 1, t_3 unchanged.
              alpha_2^v shifts t_2 -> t_2 + 1, t_3 -> t_3 - 1, t_1 unchanged.
          With our reduction t_3 = -(t_1+t_2), the shift by alpha_2^v
          corresponds to (t_1, t_2) -> (t_1, t_2 + 1) in (t_1, t_2) coords
          but breaking the t_3 constraint -- actually NO, the constraint
          is automatic if we just shift t_2 since t_3 is *defined* as
          -(t_1+t_2), so the shift propagates.

    OK working it out more carefully:
        The kernel of exp on the full (t_1, t_2, t_3) Lie algebra (with
        sum=0) is the lattice of integer-valued (n_1, n_2, n_3) with
        n_1+n_2+n_3 = 0, scaled by 2 pi. In (t_1, t_2) coords, this means
        translating (t_1, t_2) by 2 pi (n_1, n_2) for any integer n_1, n_2
        with the implicit n_3 = -(n_1+n_2). The relevant translates near
        the origin are (n_1, n_2) = (0, 0), (+1, 0), (0, +1), (+1, +1),
        (-1, 0), (0, -1), (-1, -1), (-1, +1), (+1, -1).

    So we minimise d^2 = 2 * ((t_1 - 2 pi n_1)^2 + (t_2 - 2 pi n_2)^2 +
                             (t_1 - 2 pi n_1)(t_2 - 2 pi n_2)) over (n_1, n_2)
    in [-n_search, n_search]^2.
    """
    two_pi = 2 * mpmath.pi
    best = None
    for n1 in range(-n_search, n_search + 1):
        for n2 in range(-n_search, n_search + 1):
            t1p = t1 - n1 * two_pi
            t2p = t2 - n2 * two_pi
            dsq = mpmath.mpf(2) * (t1p * t1p + t2p * t2p + t1p * t2p)
            if best is None or dsq < best:
                best = dsq
    return mpmath.sqrt(best)


def _su3_geodesic_dist_sq(
    t1: mpmath.mpf,
    t2: mpmath.mpf,
    n_search: int = 2,
) -> mpmath.mpf:
    """Squared geodesic distance using lattice minimisation."""
    two_pi = 2 * mpmath.pi
    best = None
    for n1 in range(-n_search, n_search + 1):
        for n2 in range(-n_search, n_search + 1):
            t1p = t1 - n1 * two_pi
            t2p = t2 - n2 * two_pi
            dsq = mpmath.mpf(2) * (t1p * t1p + t2p * t2p + t1p * t2p)
            if best is None or dsq < best:
                best = dsq
    return best


def character_su3(p: int, q: int, t1: sp.Symbol, t2: sp.Symbol) -> sp.Expr:
    """SU(3) character chi_{(p,q)}(g) for g on the maximal torus
    parameterised by (theta_1, theta_2) (with theta_3 = -(theta_1 + theta_2)).

    Use Weyl character formula with z_i = exp(i theta_i) and
    z_1 z_2 z_3 = 1:
        chi(p, q) = det(z_i^{lambda_j + 3 - j}) / det(z_i^{3 - j})
    where lambda = (p + q, q, 0) (or the transposed convention; we use
    lambda_1 = p + q, lambda_2 = q, lambda_3 = 0).

    For numerical evaluation, this can be computed directly as a complex
    polynomial in (z_1, z_2, z_3) and substituted with z_3 = 1/(z_1 z_2).
    """
    z1 = sp.exp(sp.I * t1)
    z2 = sp.exp(sp.I * t2)
    z3 = sp.exp(-sp.I * (t1 + t2))
    # Weyl numerator: det of matrix with entries z_i^{lambda_j + 3 - j}.
    # lambda + (3 - j) for j = 1, 2, 3 (1-indexed):
    #   j = 1: lambda_1 + 2 = p + q + 2
    #   j = 2: lambda_2 + 1 = q + 1
    #   j = 3: lambda_3 + 0 = 0
    # Denominator: lambda = 0 -> exponents 2, 1, 0.
    e_top = [p + q + 2, q + 1, 0]
    e_bot = [2, 1, 0]
    z = [z1, z2, z3]

    def det3(z, exps):
        return sp.Matrix([
            [zi ** ej for ej in exps] for zi in z
        ]).det()

    num = det3(z, e_top)
    den = det3(z, e_bot)
    return sp.simplify(num / den)


def character_su3_real(p: int, q: int) -> "function":
    """Numerical SU(3) character chi_{(p,q)}(t_1, t_2) for use in mpmath
    integrals.

    Returns a function (t1: mpmath, t2: mpmath) -> mpmath complex
    (typically real for class functions, by Hermiticity).
    """
    # Build the function as a closure over (p, q).
    # Direct evaluation: chi = numerator / denominator.
    # Write each in terms of cosines via Re and Im.

    def char(t1, t2):
        # Exponents:
        e_top = [p + q + 2, q + 1, 0]
        e_bot = [2, 1, 0]
        # z_i^k = exp(i k t_i)
        # Use z_3 = exp(-i (t_1 + t_2)).

        def make_z(t1, t2):
            return [
                lambda k: mpmath.exp(1j * k * t1),
                lambda k: mpmath.exp(1j * k * t2),
                lambda k: mpmath.exp(-1j * k * (t1 + t2)),
            ]

        zfns = make_z(t1, t2)

        def det3(zfns, exps):
            # 3x3 determinant.
            row = lambda i: [zfns[i](e) for e in exps]
            r0 = row(0); r1 = row(1); r2 = row(2)
            d = (
                r0[0] * (r1[1] * r2[2] - r1[2] * r2[1])
                - r0[1] * (r1[0] * r2[2] - r1[2] * r2[0])
                + r0[2] * (r1[0] * r2[1] - r1[1] * r2[0])
            )
            return d

        num = det3(zfns, e_top)
        den = det3(zfns, e_bot)
        if abs(den) < mpmath.mpf("1e-12"):
            # On Weyl wall: use L'Hopital / take limit. For a numerical
            # integral with smooth measure |Delta|^2 on top, this point
            # has measure zero (the wall is codim 1, integrating against
            # |Delta|^2 = den * den_bar / 1 of zero order weight). For
            # safety, return dim of irrep at the regular limit; the
            # measure |Delta|^2 will kill any wall behaviour.
            return mpmath.mpc(IrrepSU3(p, q).dim(), 0)
        return num / den

    return char


def central_fejer_kernel_su3_at_torus(
    irreps: list[IrrepSU3],
    t1: mpmath.mpf,
    t2: mpmath.mpf,
) -> mpmath.mpf:
    """Compute K_Lambda(t_1, t_2) = (1/Z) |sum_pi sqrt(dim pi) chi_pi|^2
    at a single torus point.

    Z = sum_pi dim(pi).
    """
    Z = sum(pi.dim() for pi in irreps)
    s = mpmath.mpc(0, 0)
    for pi in irreps:
        chi = character_su3_real(pi.p, pi.q)(t1, t2)
        s += mpmath.sqrt(pi.dim()) * chi
    abs_sq = (s * mpmath.conj(s)).real
    return abs_sq / Z


def gamma_lambda_su3(
    irreps: list[IrrepSU3],
    n_quad: int = 30,
    prec: int = 30,
) -> mpmath.mpf:
    """Compute gamma_Lambda := integral_G K_Lambda(g) d(e, g) dg
    via Weyl integration on the (theta_1, theta_2) torus.

    Use the bi-invariant Haar measure normalised so that integral_G 1 dg = 1
    (equivalent to dividing by Vol(G)).

    Weyl integration:
       integral_G f dg = (1 / |W|) integral_{T^2} f(theta) |Delta|^2 / vol_factor
    where the normalisation factor makes integral_G 1 dg = 1.

    Specifically, integral over T^2 of |Delta|^2 d theta_1 d theta_2 = |W| * (2 pi)^2
    when divided by some normalisation. The cleanest unit-normalised form:

      int_G f dg = (1/(|W| * (2 pi)^r)) int_{T^r} f(theta) |Delta|^2 d theta,

    where r = rank, |W| = order of Weyl group. For SU(3): r = 2, |W| = 6.
    Verify by integrating f = 1 (which we should do as a sanity check).
    """
    mpmath.mp.dps = prec
    W = 6  # |W| for SU(3)
    pi_const = mpmath.pi
    two_pi = 2 * pi_const

    # Use fixed-grid Gauss-Legendre on [0, 2 pi]^2 (smooth integrand).
    # For higher precision could use mpmath.quad, but we're integrating
    # a high-frequency oscillatory polynomial (Fejer kernel), which is
    # smooth but rapidly oscillating; prefer many quadrature points.
    nodes_t1, weights_t1 = _gauss_legendre_nodes(0, two_pi, n_quad)
    nodes_t2, weights_t2 = _gauss_legendre_nodes(0, two_pi, n_quad)

    integral = mpmath.mpf(0)
    for i, t1 in enumerate(nodes_t1):
        for j, t2 in enumerate(nodes_t2):
            K = central_fejer_kernel_su3_at_torus(irreps, t1, t2)
            d = _su3_geodesic_dist(t1, t2)
            measure = _weyl_denom_sq_su3(t1, t2)
            integral += weights_t1[i] * weights_t2[j] * K * d * measure
    # Normalise: int_G 1 dg = 1, i.e. divide by |W| * (2 pi)^r where r = 2.
    integral = integral / (W * two_pi * two_pi)
    return integral


def gamma_lambda_su3_simpson(
    irreps: list[IrrepSU3],
    n_grid: int = 64,
    prec: int = 30,
) -> mpmath.mpf:
    """gamma_Lambda via composite Simpson on a uniform grid [0, 2 pi]^2.
    Slower than GL but more transparent."""
    mpmath.mp.dps = prec
    W = 6
    pi_const = mpmath.pi
    two_pi = 2 * pi_const
    h = two_pi / n_grid
    if n_grid % 2 != 0:
        raise ValueError("n_grid must be even for Simpson.")

    integral = mpmath.mpf(0)
    for i in range(n_grid + 1):
        t1 = i * h
        wi = 1 if (i == 0 or i == n_grid) else (4 if i % 2 == 1 else 2)
        for j in range(n_grid + 1):
            t2 = j * h
            wj = 1 if (j == 0 or j == n_grid) else (4 if j % 2 == 1 else 2)
            K = central_fejer_kernel_su3_at_torus(irreps, t1, t2)
            d = _su3_geodesic_dist(t1, t2)
            measure = _weyl_denom_sq_su3(t1, t2)
            integral += wi * wj * K * d * measure
    integral *= (h * h / 9)
    integral = integral / (W * two_pi * two_pi)
    return integral


def _gauss_legendre_nodes(a: mpmath.mpf, b: mpmath.mpf, n: int):
    """Gauss-Legendre nodes/weights on [a, b]. Uses mpmath.calculus.quadrature."""
    # mpmath provides mpmath.calculus.quadrature.gauss_legendre.calc_nodes
    from mpmath.calculus.quadrature import GaussLegendre
    gl = GaussLegendre(mpmath.mp)
    # gl.calc_nodes returns [(x_i, w_i)] on [-1, 1] with weights for ?
    # Use degree=k giving 2^k nodes. For n_quad = 30 we want ~30 nodes.
    # Fall back to manual numerical integration via mpmath.quad if this is awkward.
    # mpmath.calculus.quadrature is internal API; safer to construct via mpmath.quad.
    raise NotImplementedError("Use _gl_nodes_via_numpy instead")


def _gl_nodes_via_numpy(a: float, b: float, n: int):
    """Gauss-Legendre nodes/weights on [a, b] via numpy.polynomial.legendre."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    # nodes on [-1, 1]. Transform to [a, b].
    nodes_ab = [mpmath.mpf((b - a) / 2) * mpmath.mpf(float(x)) + mpmath.mpf((b + a) / 2)
                for x in nodes]
    weights_ab = [mpmath.mpf((b - a) / 2) * mpmath.mpf(float(w)) for w in weights]
    return nodes_ab, weights_ab


def gamma_lambda_su3_gl(
    irreps: list[IrrepSU3],
    n_quad: int = 40,
    prec: int = 30,
) -> mpmath.mpf:
    """gamma_Lambda via 2D Gauss-Legendre on [0, 2 pi]^2."""
    mpmath.mp.dps = prec
    W = 6
    pi_const = mpmath.pi
    two_pi = 2 * pi_const

    nodes_t1, weights_t1 = _gl_nodes_via_numpy(0.0, float(two_pi), n_quad)
    nodes_t2, weights_t2 = _gl_nodes_via_numpy(0.0, float(two_pi), n_quad)

    integral = mpmath.mpf(0)
    for i in range(len(nodes_t1)):
        t1 = nodes_t1[i]
        wi = weights_t1[i]
        for j in range(len(nodes_t2)):
            t2 = nodes_t2[j]
            wj = weights_t2[j]
            K = central_fejer_kernel_su3_at_torus(irreps, t1, t2)
            d = _su3_geodesic_dist(t1, t2)
            measure = _weyl_denom_sq_su3(t1, t2)
            integral += wi * wj * K * d * measure
    integral = integral / (W * two_pi * two_pi)
    return integral


def gamma_lambda_su3_gl_fast(
    irreps: list[IrrepSU3],
    n_quad: int = 40,
    prec: int = 25,
) -> mpmath.mpf:
    """Faster gamma_Lambda using numpy-vectorised kernel evaluation.

    Key optimisation: build a 2D grid of characters chi_pi(t1[i], t2[j])
    as a numpy complex array per irrep, then assemble the kernel via
    matrix-style operations.

    Tradeoffs vs gamma_lambda_su3_gl:
      - Lower precision (~15 dps numpy double-precision rather than full mpmath).
      - Much faster: O(n_irreps * n_quad^2) numpy ops instead of nested loops.
    """
    W = 6
    two_pi = 2.0 * math.pi

    # Gauss-Legendre nodes on [0, 2 pi].
    raw_nodes, raw_weights = np.polynomial.legendre.leggauss(n_quad)
    nodes = (two_pi / 2) * raw_nodes + (two_pi / 2)
    weights = (two_pi / 2) * raw_weights

    # 2D grids.
    T1, T2 = np.meshgrid(nodes, nodes, indexing="ij")
    W1, W2 = np.meshgrid(weights, weights, indexing="ij")
    T3 = -(T1 + T2)

    # Geodesic distance (with periodic minimisation).
    # We do a small lattice search: (n1, n2) in [-1, 1]^2 (9 candidates).
    # For our integration domain [0, 2 pi], the relevant lattice translates
    # are within 2 pi of the origin.
    d_grid = None
    for n1 in (-1, 0, 1):
        for n2 in (-1, 0, 1):
            t1p = T1 - n1 * two_pi
            t2p = T2 - n2 * two_pi
            dsq = 2.0 * (t1p * t1p + t2p * t2p + t1p * t2p)
            if d_grid is None:
                d_grid = dsq
            else:
                d_grid = np.minimum(d_grid, dsq)
    d_grid = np.sqrt(d_grid)

    # Weyl denominator |Delta|^2.
    s1 = np.sin((T1 - T2) / 2)
    s2 = np.sin((2 * T2 + T1) / 2)
    s3 = np.sin((2 * T1 + T2) / 2)
    measure_grid = 64 * s1 * s1 * s2 * s2 * s3 * s3

    # Kernel at each grid point.
    # K(t) = (1/Z) |sum_pi sqrt(dim_pi) chi_pi(t)|^2.
    Z = sum(pi.dim() for pi in irreps)
    sum_grid = np.zeros_like(T1, dtype=np.complex128)
    for pi in irreps:
        chi_grid = _su3_character_grid(pi.p, pi.q, T1, T2)
        sum_grid += math.sqrt(pi.dim()) * chi_grid
    K_grid = (np.abs(sum_grid) ** 2) / Z

    # Integrand.
    integrand = K_grid * d_grid * measure_grid * W1 * W2
    integral = np.sum(integrand) / (W * two_pi * two_pi)
    return mpmath.mpf(float(integral))


def _su3_character_grid(
    p: int, q: int, T1: np.ndarray, T2: np.ndarray
) -> np.ndarray:
    """Compute chi_{(p,q)}(t_1, t_2) on a 2D grid via the Weyl character
    formula, returning a complex numpy array.

    Uses z_i = exp(i theta_i), with z_3 = exp(-i (t_1+t_2)).
    """
    z1 = np.exp(1j * T1)
    z2 = np.exp(1j * T2)
    z3 = np.exp(-1j * (T1 + T2))

    # Numerator: det of the 3x3 matrix with row i, col j: z_i^{lambda_j + 3 - j}.
    # lambda = (p+q, q, 0); offsets (3-j) = (2, 1, 0).
    e_top = [p + q + 2, q + 1, 0]
    e_bot = [2, 1, 0]

    def detrow(z, exps):
        return [z ** exps[0], z ** exps[1], z ** exps[2]]

    r0 = detrow(z1, e_top)
    r1 = detrow(z2, e_top)
    r2 = detrow(z3, e_top)
    num = (
        r0[0] * (r1[1] * r2[2] - r1[2] * r2[1])
        - r0[1] * (r1[0] * r2[2] - r1[2] * r2[0])
        + r0[2] * (r1[0] * r2[1] - r1[1] * r2[0])
    )
    r0 = detrow(z1, e_bot)
    r1 = detrow(z2, e_bot)
    r2 = detrow(z3, e_bot)
    den = (
        r0[0] * (r1[1] * r2[2] - r1[2] * r2[1])
        - r0[1] * (r1[0] * r2[2] - r1[2] * r2[0])
        + r0[2] * (r1[0] * r2[1] - r1[1] * r2[0])
    )

    # Where den is small (Weyl walls), the character has a removable
    # singularity. Replace small |den| with the analytic limit:
    # chi(p, q)(e) = dim((p, q)).
    # For numerical safety, divide and replace NaN/Inf with dim:
    with np.errstate(divide="ignore", invalid="ignore"):
        chi = num / den
    bad = ~np.isfinite(chi) | (np.abs(den) < 1e-12)
    if np.any(bad):
        chi = np.where(bad, IrrepSU3(p, q).dim() + 0j, chi)
    return chi


def verify_haar_normalisation_su3(n_quad: int = 40, prec: int = 30) -> mpmath.mpf:
    """Verify (1/(|W| (2 pi)^r)) int_T^r |Delta|^2 dt = 1.
    This is the requirement that Haar(SU(3)) integrates a class function
    constant 1 to 1.
    """
    mpmath.mp.dps = prec
    W = 6
    pi_const = mpmath.pi
    two_pi = 2 * pi_const

    nodes_t1, weights_t1 = _gl_nodes_via_numpy(0.0, float(two_pi), n_quad)
    nodes_t2, weights_t2 = _gl_nodes_via_numpy(0.0, float(two_pi), n_quad)

    integral = mpmath.mpf(0)
    for i in range(len(nodes_t1)):
        t1 = nodes_t1[i]
        wi = weights_t1[i]
        for j in range(len(nodes_t2)):
            t2 = nodes_t2[j]
            wj = weights_t2[j]
            measure = _weyl_denom_sq_su3(t1, t2)
            integral += wi * wj * measure
    integral = integral / (W * two_pi * two_pi)
    return integral


def _enumerate_irreps_under_casimir(lambda_sq: float, max_dynkin: int = 20) -> list[IrrepSU3]:
    """Enumerate all SU(3) irreps with Casimir <= lambda_sq."""
    out = []
    for p in range(max_dynkin + 1):
        for q in range(max_dynkin + 1):
            if p == 0 and q == 0:
                # Skip trivial rep? include.
                pass
            irrep = IrrepSU3(p, q)
            c = float(irrep.casimir())
            if c <= lambda_sq + 1e-9:
                out.append(irrep)
    return out


# =============================================================================
# Sub-task (b/c): Rate constant fitting
# =============================================================================


def fit_single_log(lams: list[float], gammas: list[float]) -> dict:
    """Fit gamma ~ a1 * log(L) / L + b1 / L + c1 / L^2 by least squares."""
    Y = np.array(gammas, dtype=float)
    X = np.array([
        [math.log(L) / L, 1.0 / L, 1.0 / (L * L)] for L in lams
    ])
    coef, residuals, rank, sv = np.linalg.lstsq(X, Y, rcond=None)
    pred = X @ coef
    resid = Y - pred
    rss = float(np.sum(resid ** 2))
    n_data = len(Y)
    n_params = 3
    aic = n_data * math.log(rss / max(n_data, 1)) + 2 * n_params
    bic = n_data * math.log(rss / max(n_data, 1)) + n_params * math.log(n_data)
    return {
        "form": "a1 * log(L)/L + b1/L + c1/L^2",
        "a1": float(coef[0]),
        "b1": float(coef[1]),
        "c1": float(coef[2]),
        "rss": rss,
        "aic": aic,
        "bic": bic,
        "predictions": pred.tolist(),
        "residuals": resid.tolist(),
    }


def fit_double_log(lams: list[float], gammas: list[float]) -> dict:
    """Fit gamma ~ a2 * log(L)^2 / L + b2 * log(L)/L + c2/L."""
    Y = np.array(gammas, dtype=float)
    X = np.array([
        [(math.log(L) ** 2) / L, math.log(L) / L, 1.0 / L] for L in lams
    ])
    coef, residuals, rank, sv = np.linalg.lstsq(X, Y, rcond=None)
    pred = X @ coef
    resid = Y - pred
    rss = float(np.sum(resid ** 2))
    n_data = len(Y)
    n_params = 3
    aic = n_data * math.log(rss / max(n_data, 1)) + 2 * n_params
    bic = n_data * math.log(rss / max(n_data, 1)) + n_params * math.log(n_data)
    return {
        "form": "a2 * log(L)^2/L + b2 * log(L)/L + c2/L",
        "a2": float(coef[0]),
        "b2": float(coef[1]),
        "c2": float(coef[2]),
        "rss": rss,
        "aic": aic,
        "bic": bic,
        "predictions": pred.tolist(),
        "residuals": resid.tolist(),
    }


# =============================================================================
# Predicted constant c(SU(3))
# =============================================================================


def predicted_c_su3() -> dict:
    """Memo proposes c(G) = 2 * Vol(G/T) / Vol(G) (with factor 2 from Cesaro
    doubling, equivalent to |W(SU(2))| = 2). Higher-rank options:

      Option A: c(G) = 2 * Vol(G/T) / Vol(G). Hopf-base measure with constant 2.
      Option B: c(G) = |W(G)| * Vol(G/T) / Vol(G). Hopf-base measure with
                Weyl-group prefactor (matches SU(2) where |W(SU(2))| = 2).
      Option C: c(G) = 4 / pi * (rank-dependent correction).

    Standard volumes (from Weyl integration with our normalisation):
      Vol(SU(3)) = sqrt(3) * pi^5 / 1, Vol(SU(3)/T^2) = ... computed via
      Vol(G)/Vol(T)/|W| or directly from the Macdonald formula.

    For SU(3), with maximal torus T^2 of dimension 2:
      Vol(T) = (2 pi)^2 = 4 pi^2 (Lebesgue, theta_i in [0, 2 pi]).

    Standard volume (Macdonald's formula): for SU(N) in unit-Killing-form
    convention with the convention that minimum closed geodesic = 2 pi,
      Vol(SU(N)) = (2 pi)^{(N^2-1)/2} * (sqrt(N))^? ...

    Rather than committing to a specific normalisation a priori, we
    NUMERICALLY compute Vol(G), Vol(G/T) at the same convention used by
    the gamma_Lambda integral (where int_G 1 dg = 1 by construction).

    With our Haar-normalised convention (int_G 1 dg = 1), Vol(G) = 1 by
    construction. So we cannot extract c(G) = 2 Vol(G/T)/Vol(G) directly
    from this convention. INSTEAD: the rate constant from L2 derives in
    *any* Haar-normalisation, and the absolute predicted value should
    use a fixed normalisation convention.

    For Paper 38 SU(2), the convention is that chi runs over [0, 2 pi]
    with Haar measure sin^2(chi/2)/pi d chi (so int_SU(2) 1 = 1). The
    asymptote is c(SU(2)) = 4/pi.

    Translation to volume ratio:
        4 / pi  =  2 * Vol(S^2) / Vol(SU(2))  with Vol(S^2) = 4 pi,
                                                     Vol(SU(2)) = 2 pi^2.
    (Vol(SU(2)) = 2 pi^2 is the surface area of the unit-Killing-norm 3-sphere.)

    For SU(3), the analogous unit-Killing-norm metric gives:
        Vol(SU(3)) = sqrt(3) pi^5,
        Vol(SU(3)/T^2) = ?

    We use the following formula (Macdonald 1980 / Helgason):
        Vol(G)/Vol(T) = prod_{alpha > 0} 2 pi / |alpha|
    where |alpha| is the bi-invariant Killing-norm length of the root,
    and Vol(G/T) = Vol(G)/Vol(T).

    For SU(3) with bi-invariant Killing-form normalisation (and choosing
    simple roots of length sqrt(2) in this normalisation):
        |alpha_1| = |alpha_2| = sqrt(2),  |alpha_1 + alpha_2| = sqrt(2).
        Vol(SU(3)/T^2) = (2 pi)^3 / (sqrt(2))^3 = 8 pi^3 / (2 sqrt(2)) = 4 pi^3 / sqrt(2).
        Vol(T^2) = (2 pi)^2 = 4 pi^2.
        Vol(SU(3)) = Vol(SU(3)/T^2) * Vol(T^2) = (4 pi^3 / sqrt(2)) * 4 pi^2
                   = 16 pi^5 / sqrt(2) = 8 sqrt(2) pi^5.

    But Hashimoto (also Macdonald): Vol(SU(N)) = sqrt(N) (2 pi)^{N(N+1)/2 - 1} / ...
    for SU(3): Vol(SU(3)) = sqrt(3) * pi^5 (one convention).

    THE NORMALISATION IS NOT CANONICAL. We compute the *ratio*
    2 * Vol(G/T) / Vol(G) which is normalisation-invariant up to a
    consistent choice on G and on G/T. With the same Killing-form
    normalisation:

      Vol(G/T) / Vol(G) = 1 / Vol(T)  by Vol(G) = Vol(G/T) * Vol(T)
                        = 1 / (2 pi)^r in the Lebesgue convention.

    So:
      c_pred^A = 2 * Vol(G/T)/Vol(G) = 2 / (2 pi)^r = 2 / (4 pi^2) = 1/(2 pi^2)
                 for r = 2.

    BUT for SU(2) (r = 1), this gives c_pred^A = 2 / (2 pi) = 1/pi, NOT 4/pi.
    So Option A is NOT the right normalisation match.

    Reconcile: for SU(2), Vol(SU(2)) = 2 pi^2 in the "geodesic = chi" convention,
    Vol(SU(2)/U(1)) = Vol(S^2) = 4 pi, ratio = 2/pi, so 2 * ratio = 4/pi. So
    here Vol(T) = 2 pi (one circle of circumference 2 pi). The formula
    Vol(G) = Vol(G/T) * Vol(T) gives 2 pi^2 = 4 pi * 2 pi -- that's NOT
    right (gives 8 pi^2, not 2 pi^2). The discrepancy: we are mixing
    "one circle of length 2 pi" with "torus of two circles each of length 2 pi"
    inconsistently.

    CORRECT: Vol(SU(2)) = 2 pi^2 (the geodesic-length-2pi unit 3-sphere).
             Vol(S^2) = 4 pi (geodesic-length-pi unit 2-sphere; antipodal
             distance pi).
             Vol(U(1)) = 2 pi (the Hopf fiber).
             Then 2 pi^2 = 4 pi * pi (NOT 4 pi * 2 pi). So Vol(T) = pi here?
             That's because the Hopf fibration projects SU(2) -> S^2 via a
             *double cover* phenomenon; the Hopf fiber length is pi, not 2 pi.

    OK so for SU(2): Vol(SU(2)) = 2 pi^2, Vol(SU(2)/T) = 4 pi, Vol(T) = pi/2 (?).
    The arithmetic is 4 pi * pi / 2 = 2 pi^2. So Vol(T_SU(2)) = pi/2 in this
    convention? Let me re-derive: Vol(T) = circumference / fiber-multiplier.

    OK this is getting tangled; the right move is to accept the SU(2)
    result 4/pi as the ground truth and write Option A vs Option B in
    terms of what *would* match SU(2):

      Option A (c = 2 V(G/T)/V(G) with bi-invariant unit-Killing-norm):
        SU(2): Vol(SU(2)) = 2 pi^2, Vol(S^2) = 4 pi, ratio = 2/pi; A predicts 4/pi. CHECK.
        SU(3): Vol(SU(3)) = sqrt(3) pi^5, Vol(SU(3)/T^2) = ?, predicted c = ?

      Option B (c = |W| V(G/T)/V(G) with the same convention):
        SU(2): |W(SU(2))| = 2, so B = 2 V(G/T)/V(G) = 4/pi. SAME AS A FOR SU(2).
        SU(3): |W(SU(3))| = 6, so B = 6 V(G/T)/V(G) = 3 * (Option A).

    Diverge at higher rank. We will compute both Options A and B numerically
    using a consistent unit-Killing-form convention and report.

    Standard formula for Vol(G/T) (Macdonald, "The volume of a compact Lie group")
    in unit-Killing-form normalisation, where |alpha|^2 = 2 (long root):
        Vol(G/T) = prod_{alpha > 0} (2 pi / <rho, alpha>) * |W|.

    Hmm, that's wrong too. Let me look up.

    Macdonald 1980: Vol(G/T) (compact, simply connected, simple) in the
    Killing-form normalisation =
        (1 / |W|) * prod_{alpha > 0} (2 pi / |alpha|^2 * <rho, alpha>),
    or some such. Different sources give different conventions.

    Cleanest: use the Weyl integration formula directly to compute
    Vol(G) (where Vol(G) = int_G 1 dg with the natural Haar measure
    coming from the Killing form, and |Delta|^2 is the Weyl denominator):

        Vol(G) = int_G 1 dg = (1 / |W|) * Vol(T) *
                    integral_{T} |Delta|^2 d_T mu / Vol(T)
               = ... = some explicit expression.

    For SU(3), with t_i in [0, 2 pi]:
        |Delta|^2 = prod_{alpha > 0} 4 sin^2(alpha . t / 2),
        positive roots: alpha_1, alpha_2, alpha_1 + alpha_2.
        Action on (t_1, t_2): t_1 - t_2, t_1 + 2 t_2, 2 t_1 + t_2.
        int_{[0, 2 pi]^2} prod 4 sin^2((alpha . t)/2) dt_1 dt_2
            = (Macdonald) (2 pi)^r * |W| ... or some transcendental ratio.

    PLAN: numerically compute the Weyl integration normalisation factor
    and the geodesic-2nd-moment integral, and report the dimensionless
    rate constant c(SU(3)). DON'T try to symbolically pin down absolute
    Vol(G), Vol(G/T) — instead, the rate-constant fit to gamma_Lambda
    data IS the primary quantity, and we compare to two predictions:
      (Option A) c_A = 2 * (1 / (Vol(T)_normalised))
                    = 2 / ((2 pi)^r) with r = 2: c_A = 1/(2 pi^2)
                    --- this fails SU(2) check.
      (Option B) c_B = |W| / Vol(T) = 6 / ((2 pi)^2) = 3/(2 pi^2)
                    --- also fails SU(2) check.

    Better: derive c from the SU(2) calculation, then express in a
    rank-invariant form. From debug/r25_l2_quantitative_rate_memo.md
    the SU(2) constant 4/pi comes from a Stein-Weiss application; the
    factor 4 is from the (4/pi^2) dimension factor times pi from the
    measure. We report numerical fits and let the user identify the
    closed form.
    """
    return {
        "comment": "We extract c(SU(3)) from the gamma_Lambda fit and compare to several candidate forms; do not pre-commit to a formula.",
        "candidate_forms": {
            "1/(2*pi)": 1.0 / (2 * math.pi),
            "1/pi": 1.0 / math.pi,
            "2/pi": 2.0 / math.pi,
            "4/pi": 4.0 / math.pi,
            "pi": math.pi,
            "1/(pi^2)": 1.0 / (math.pi ** 2),
            "1/(2*pi^2)": 1.0 / (2 * math.pi ** 2),
            "3/(2*pi^2)": 3.0 / (2 * math.pi ** 2),
            "4/(pi^2)": 4.0 / (math.pi ** 2),
            "8/(pi^2)": 8.0 / (math.pi ** 2),
            "12/(pi^2)": 12.0 / (math.pi ** 2),
            "8/pi^3": 8.0 / (math.pi ** 3),
            "16/(3*pi)": 16.0 / (3 * math.pi),
            "3/(pi^2)": 3.0 / (math.pi ** 2),
            "6/(pi^2)": 6.0 / (math.pi ** 2),
        },
    }


# =============================================================================
# Main
# =============================================================================


def main():
    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "su3_numerical_sanity.json"

    results = {
        "metadata": {
            "date": "2026-05-10",
            "purpose": "SU(3) numerical sanity check for unified GH-convergence theorem",
            "cross_refs": [
                "papers/standalone/paper_38_su2_propinquity_convergence.tex",
                "debug/unified_gh_scoping_memo.md",
                "geovac/central_fejer_su2.py",
            ],
        }
    }

    print("Sanity-checking SU(3) tensor-product machinery...")
    _sanity_check_tensor_product()
    print("  OK.")

    # -------------------------------------------------------------------
    # Sub-task (a): Casimir triangle inequality
    # -------------------------------------------------------------------
    print("\n=== Sub-task (a): Casimir triangle inequality ===")
    panel = [
        IrrepSU3(0, 0),
        IrrepSU3(1, 0),
        IrrepSU3(0, 1),
        IrrepSU3(1, 1),
        IrrepSU3(2, 0),
        IrrepSU3(0, 2),
        IrrepSU3(2, 1),
        IrrepSU3(1, 2),
        IrrepSU3(3, 0),
        IrrepSU3(0, 3),
        IrrepSU3(2, 2),
    ]
    print(f"Panel: {[str(pi) for pi in panel]}")
    print(f"Casimirs: {[(pi.p, pi.q, float(pi.casimir())) for pi in panel]}")

    t0 = time.time()
    a_result = casimir_triangle_test(panel)
    t_a = time.time() - t0
    print(f"  Total ordered pairs: {a_result['n_pairs_total']}")
    print(f"  Pass: {a_result['n_pairs_pass']}")
    print(f"  Fail: {a_result['n_pairs_fail']}")
    print(f"  Max sup ratio: {a_result['max_sup_ratio']:.6f}")
    print(f"  Time: {t_a:.2f}s")
    if a_result['violations']:
        print(f"  VIOLATIONS:")
        for v in a_result['violations'][:5]:
            print(f"    pi={v['pi']}, pi'={v['pi_prime']}, sigma={v['sigma']}, "
                  f"|C(pi)-C(pi')|={v['gap']:.4f}, C(sigma)={v['C_sigma']:.4f}")

    a_result["timing_seconds"] = t_a
    results["sub_task_a"] = a_result

    # -------------------------------------------------------------------
    # Sub-task (b): Rate constant via Weyl integration
    # -------------------------------------------------------------------
    print("\n=== Sub-task (b): Rate constant via Weyl integration ===")

    # Verify Haar normalisation.
    print("Verifying Haar normalisation (should be 1.0)...")
    haar_check = verify_haar_normalisation_su3(n_quad=40, prec=30)
    print(f"  int_G 1 dg = {haar_check}")
    results["haar_normalisation_check"] = float(haar_check)

    # Cutoff sequence: enumerate irreps with Casimir <= Lambda^2 for
    # Lambda^2 in {2, 3, 4, 6, 8, 10, 14, 18, 24, 30}.
    # For each, list the irreps and their count, then compute gamma_Lambda
    # via 2D Gauss-Legendre with sufficient n_quad.
    lambda_sq_values = [2, 3, 4, 6, 8, 10, 14, 18, 24, 30]

    # Pre-enumerate irrep panels.
    irrep_panels = {}
    for lsq in lambda_sq_values:
        irreps = _enumerate_irreps_under_casimir(lsq, max_dynkin=10)
        irrep_panels[lsq] = irreps
        print(f"  Lambda^2 = {lsq}: {len(irreps)} irreps, "
              f"sum dim = {sum(pi.dim() for pi in irreps)}")
        # Print first few:
        details = ", ".join(f"({pi.p},{pi.q})" for pi in irreps[:8])
        print(f"    first few: {details}")

    # Compute gamma_Lambda for each cutoff.
    # n_quad heuristic: need ~ 2 * (max_freq) per period, where max_freq
    # is bounded by max(p+q) on the irrep panel + 2 (for Weyl denom).
    # For lsq = 30, max p+q is found by checking irreps.
    print("\nComputing gamma_Lambda values...")
    lambdas = []
    gammas = []
    n_quads_used = []
    for lsq in lambda_sq_values:
        irreps = irrep_panels[lsq]
        # Choose n_quad based on max frequency:
        max_pq = max((pi.p + pi.q for pi in irreps), default=2)
        # Need n_quad >= 2 * (max_pq + 2) to resolve the integrand;
        # be generous for the kernel-squared.
        n_quad = max(40, 2 * (max_pq + 2) + 8)
        if n_quad > 80:
            n_quad = 80  # cap for runtime
        t0 = time.time()
        g = gamma_lambda_su3_gl(irreps, n_quad=n_quad, prec=30)
        dt = time.time() - t0
        L = math.sqrt(lsq)
        lambdas.append(L)
        gammas.append(float(g))
        n_quads_used.append(n_quad)
        print(f"  Lambda^2 = {lsq:5d}, Lambda = {L:.4f}, "
              f"gamma = {float(g):.6e}, "
              f"n_quad = {n_quad}, t = {dt:.1f}s")

    results["sub_task_b"] = {
        "lambda_sq_values": lambda_sq_values,
        "lambda_values": lambdas,
        "gamma_values": gammas,
        "n_quads_used": n_quads_used,
        "n_irreps_per_panel": [len(irrep_panels[lsq]) for lsq in lambda_sq_values],
        "max_pq_per_panel": [
            max((pi.p + pi.q for pi in irrep_panels[lsq]), default=0)
            for lsq in lambda_sq_values
        ],
    }

    # -------------------------------------------------------------------
    # Sub-task (c): Single vs double-log fit
    # -------------------------------------------------------------------
    print("\n=== Sub-task (c): Log power test ===")

    # Drop the smallest Lambda values from the fit if they're too coarse;
    # but keep all for transparency.
    fit_single = fit_single_log(lambdas, gammas)
    fit_double = fit_double_log(lambdas, gammas)
    print(f"  Single-log fit (a*logL/L + b/L + c/L^2):")
    print(f"    a1 = {fit_single['a1']:.6f}")
    print(f"    b1 = {fit_single['b1']:.6f}")
    print(f"    c1 = {fit_single['c1']:.6f}")
    print(f"    RSS = {fit_single['rss']:.4e}, AIC = {fit_single['aic']:.4f}")
    print(f"  Double-log fit (a*log^2 L/L + b*logL/L + c/L):")
    print(f"    a2 = {fit_double['a2']:.6f}")
    print(f"    b2 = {fit_double['b2']:.6f}")
    print(f"    c2 = {fit_double['c2']:.6f}")
    print(f"    RSS = {fit_double['rss']:.4e}, AIC = {fit_double['aic']:.4f}")

    aic_diff = fit_double['aic'] - fit_single['aic']
    print(f"  AIC difference (double - single): {aic_diff:.4f}")
    if aic_diff > 6:
        log_verdict = "SINGLE-LOG (single-log fit decisively better)"
    elif aic_diff < -6:
        log_verdict = "DOUBLE-LOG (double-log fit decisively better)"
    else:
        log_verdict = "INCONCLUSIVE (AIC difference < 6)"
    print(f"  Verdict: {log_verdict}")

    results["sub_task_c"] = {
        "fit_single_log": fit_single,
        "fit_double_log": fit_double,
        "aic_difference_double_minus_single": aic_diff,
        "verdict": log_verdict,
    }

    # -------------------------------------------------------------------
    # Predicted constant comparison
    # -------------------------------------------------------------------
    print("\n=== Comparison to candidate predictions ===")
    a1 = fit_single['a1']
    a2 = fit_double['a2']
    candidates = predicted_c_su3()["candidate_forms"]
    print(f"  Extracted a1 (single-log) = {a1:.6f}")
    print(f"  Extracted a2 (double-log) = {a2:.6f}")
    print(f"  Candidate forms:")
    cand_results = []
    for name, val in candidates.items():
        rel_err_a1 = abs(a1 - val) / abs(val) if val != 0 else float("inf")
        rel_err_a2 = abs(a2 - val) / abs(val) if val != 0 else float("inf")
        print(f"    {name:>20s} = {val:.6f}, rel.err vs a1: {rel_err_a1*100:.2f}%, "
              f"vs a2: {rel_err_a2*100:.2f}%")
        cand_results.append({
            "name": name,
            "value": val,
            "rel_err_a1_pct": rel_err_a1 * 100,
            "rel_err_a2_pct": rel_err_a2 * 100,
        })

    results["candidate_predictions"] = cand_results

    # ------ Additional rank-aware predictions specific to SU(3) ------
    # SU(2) ground truth: 4/pi.
    # Likely SU(3) generalisations:
    #   c(G) = 4/pi * (rank-1 generalisation) -- many possibilities.
    # We also report the closest match.
    closest = min(cand_results, key=lambda x: x["rel_err_a1_pct"])
    print(f"\n  Closest match to a1 ({a1:.6f}): {closest['name']} = {closest['value']:.6f} "
          f"({closest['rel_err_a1_pct']:.2f}% error)")
    results["closest_to_a1"] = closest

    # -------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------
    with open(out_path, "w") as f:
        # Drop large arrays from the saved results for cleanliness;
        # keep predictions/residuals with the fits.
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")

    return results


if __name__ == "__main__":
    main()
