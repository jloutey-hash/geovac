"""Track BR-A: Breit interaction angular decomposition in Dirac (kappa, m_j) basis.

Goal
----
Extend Paper 22's angular sparsity theorem (rank-0 Coulomb) to the Breit
(rank-2 spin-spin + rank-0/1 spin-other-orbit) operators. Determine

    d_Breit(l_max) =
        (# jj-quartets with non-vanishing Breit angular matrix element)
      / (# jj-quartets at l_max).

Physics used
------------
Breit-Pauli Hamiltonian (Bethe-Salpeter §38, frequency-independent Pauli
reduction, one-scalar-exchange operator):

  H_B = H_SS + H_SOO + H_Darwin-like

  H_SS  = -(alpha^2 / r_12^3) * [ 3 (sigma_1 . rhat_12)(sigma_2 . rhat_12)
                                  - sigma_1 . sigma_2 ]
        = -(alpha^2 / r_12^3) * 2 * sqrt(24*pi/5) *
              [ Y_2(rhat_12) (x) T^2(sigma_1, sigma_2) ]^0_0
  H_SOO = -(alpha^2/4) * (1/r_12^2) *
              [ (r_12 x p_1) . (sigma_1 + 2 sigma_2)
              + (r_12 x p_2) . (sigma_2 + 2 sigma_1) ]
          -- mixture of rank-1 on each electron.

For angular sparsity purposes what matters is the **rank structure of the
irreducible tensor operator coupling the two electrons** in the *orbital*
part:

  SS-tensor part:    Y_2(rhat_12) * T^2(spin pair) . [joint rank 0]
  SS-scalar part:    Y_0 * T^0 . [scalar, identical to Coulomb, already Paper 22]
  SOO rank-1:        grad(1/r_12) . ell_i sigma_j (rank 1 on each)
  -- here we do the rank-2 (SS tensor) piece rigorously. SOO and the
  scalar SS piece are explicitly enumerated below for completeness.

Multipole expansion of the *orbital* operators in 1/r_12^3 and 1/r_12^2
produces a bipolar sum

    [Y_L1(rhat_1) (x) Y_L2(rhat_2)]^K_Q * f_{L1 L2 K}(r_1, r_2).

For Breit-Pauli:

  SS tensor piece:   contributes (L1, L2; K) triplets obeying L1 + L2 = 2
                     with K = 0 after recoupling with the rank-2 spin tensor.
                     Concretely (L1, L2) = (0,2), (1,1), (2,0); K=2.
  SS scalar piece:   L1 = L2, K=0 (identical to Coulomb). Already in Paper 22.
  SOO rank-1 piece:  (L1, L2; K=1) bipolar. L1 + L2 odd.
  Darwin-like:       delta(r_12) -> K=0 S-wave, already Coulomb-like.

Test
----
Selection rule for a bipolar tensor of rank (L1, L2) K on jj-coupled pairs
<(kappa_a m_a)(kappa_b m_b) | O^{(L1,L2)K} | (kappa_c m_c)(kappa_d m_d)>:

  * m-conservation: m_a + m_b = m_c + m_d (global; sum of 2-body tensor ranks = 0)
  * Parity: (l_a + l_c + L1) even, (l_b + l_d + L2) even
  * j-triangle: triangle(j_a, L1, j_c) via (j_a + 1/2, L1, j_c + 1/2) integer,
                triangle(j_b, L2, j_d)
  * l-triangle: triangle(l_a, L1, l_c), triangle(l_b, L2, l_d)
  * Reduced 3j (jj-coupled reduced matrix element):
      3j(j_a, L1, j_c; 1/2, 0, -1/2) != 0 AND
      3j(j_b, L2, j_d; 1/2, 0, -1/2) != 0
  * Cross-pair coupling selection: the K = 0, 1 or 2 overall rank means
    allowed (m_a-m_c) = -(m_b-m_d) summed.

For Breit (rank-2 SS), enumerate over (L1, L2) in {(0,2), (1,1), (2,0)} and
sum the unions of allowed quartets. For SOO (rank-1), enumerate
(L1, L2) in {(0,1), (1,0), (1,2), (2,1)} etc.

Exact-arithmetic throughout (sympy Rational + wigner_3j).

Outputs
-------
* `debug/data/br_a_density.json` with d_SS, d_SOO, d_Breit total over l_max 0..5.
* Selection-rule sample printout comparing quartets classified allowed/forbidden.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, FrozenSet, List, Set, Tuple

from sympy import Integer, Rational
from sympy.physics.wigner import wigner_3j


# ---------------------------------------------------------------------------
# Spinor label enumeration (reused from tier2_t0_spinor_density.py)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SpinorOrbital:
    kappa: int
    two_mj: int  # 2 * m_j, odd integer

    @property
    def l(self) -> int:
        return self.kappa if self.kappa > 0 else -self.kappa - 1

    @property
    def two_j(self) -> int:
        return 2 * abs(self.kappa) - 1


def enumerate_spinor_orbitals(l_max: int) -> List[SpinorOrbital]:
    orbs: List[SpinorOrbital] = []
    for l in range(l_max + 1):
        kappa = -(l + 1)
        two_j = 2 * l + 1
        for two_mj in range(-two_j, two_j + 1, 2):
            orbs.append(SpinorOrbital(kappa=kappa, two_mj=two_mj))
    for l in range(1, l_max + 1):
        kappa = l
        two_j = 2 * l - 1
        for two_mj in range(-two_j, two_j + 1, 2):
            orbs.append(SpinorOrbital(kappa=kappa, two_mj=two_mj))
    return orbs


# ---------------------------------------------------------------------------
# Bipolar selection rule primitive
# ---------------------------------------------------------------------------

def _triangle_x2(a_x2: int, b_x2: int, c_x2: int) -> bool:
    if c_x2 < abs(a_x2 - b_x2):
        return False
    if c_x2 > a_x2 + b_x2:
        return False
    if (a_x2 + b_x2 + c_x2) % 2 != 0:
        return False
    return True


# ---------------------------------------------------------------------------
# Per-pair rank-L selection cache
# ---------------------------------------------------------------------------

def pair_rank_allowed(
    orbs: List[SpinorOrbital], L: int
) -> Dict[Tuple[int, int], bool]:
    """For each (a, c), return True iff ALL four angular conditions hold
    at this rank L:
        * parity: (l_a + l_c + L) even
        * l-triangle: |l_a - l_c| <= L <= l_a + l_c
        * j-triangle: triangle(j_a, L, j_c)
        * reduced 3j: (j_a, L, j_c; 1/2, 0, -1/2) != 0
        * bottom-row 3j: (j_a, L, j_c; -m_a, q, m_c) != 0, q = m_a - m_c

    All evaluated by sympy wigner_3j (exact).
    """
    out: Dict[Tuple[int, int], bool] = {}
    L_sym = Integer(L)
    for i, oa in enumerate(orbs):
        la = oa.l
        ja_x2 = oa.two_j
        ma_x2 = oa.two_mj
        ja_sym = Rational(ja_x2, 2)
        for j, oc in enumerate(orbs):
            lc = oc.l
            jc_x2 = oc.two_j
            mc_x2 = oc.two_mj
            jc_sym = Rational(jc_x2, 2)

            # parity on l
            if (la + lc + L) % 2 != 0:
                out[(i, j)] = False
                continue
            # l-triangle
            if L < abs(la - lc) or L > la + lc:
                out[(i, j)] = False
                continue
            # j-triangle
            if not _triangle_x2(ja_x2, 2 * L, jc_x2):
                out[(i, j)] = False
                continue
            # reduced 3j
            w_red = wigner_3j(
                ja_sym, L_sym, jc_sym,
                Rational(1, 2), Integer(0), Rational(-1, 2),
            )
            if w_red == 0:
                out[(i, j)] = False
                continue
            # bottom-row: q = m_a - m_c, must be integer and within [-L, L]
            q_x2 = ma_x2 - mc_x2
            if q_x2 % 2 != 0:
                out[(i, j)] = False
                continue
            q = q_x2 // 2
            if abs(q) > L:
                out[(i, j)] = False
                continue
            w_m = wigner_3j(
                ja_sym, L_sym, jc_sym,
                Rational(-ma_x2, 2), Integer(q), Rational(mc_x2, 2),
            )
            if w_m == 0:
                out[(i, j)] = False
                continue
            out[(i, j)] = True
    return out


# ---------------------------------------------------------------------------
# Bipolar (L1, L2; K) density
# ---------------------------------------------------------------------------

def bipolar_density(
    orbs: List[SpinorOrbital],
    bipolar_ranks: List[Tuple[int, int]],
    K_overall: int,
) -> Tuple[int, int]:
    """Count jj-quartets (a,b,c,d) with nonzero matrix element for the
    bipolar tensor sum over the given (L1, L2) ranks.

    The two-body operator has the symbolic form
       sum_{L1,L2 in bipolar_ranks} [ Y_{L1}(rhat_1) (x) Y_{L2}(rhat_2) ]^{K_overall}.

    Allowed iff:
       * exists (L1, L2) in bipolar_ranks such that
             pair_a/c passes rank-L1 selection AND
             pair_b/d passes rank-L2 selection
       * global m-conservation: (m_a - m_c) + (m_b - m_d) obeys the overall
         K_overall constraint. For K = 0: q1 + q2 = 0, i.e. m_a+m_b = m_c+m_d.
                                   K = 1: q1 + q2 in {-1, 0, +1}.
                                   K = 2: q1 + q2 in {-2, -1, 0, 1, 2}.
         In all cases the constraint is simply
             |q1 + q2| <= K_overall.
    """
    # Precompute pair-allowed maps for each unique L used.
    used_L: Set[int] = set()
    for L1, L2 in bipolar_ranks:
        used_L.add(L1)
        used_L.add(L2)
    pair_ok: Dict[int, Dict[Tuple[int, int], bool]] = {}
    for L in used_L:
        pair_ok[L] = pair_rank_allowed(orbs, L)

    Q = len(orbs)
    total = Q ** 4
    nonzero = 0

    for a in range(Q):
        ma_x2 = orbs[a].two_mj
        for b in range(Q):
            mb_x2 = orbs[b].two_mj
            for c in range(Q):
                mc_x2 = orbs[c].two_mj
                q1_x2 = ma_x2 - mc_x2
                for d in range(Q):
                    md_x2 = orbs[d].two_mj
                    q2_x2 = mb_x2 - md_x2
                    q_tot_x2 = q1_x2 + q2_x2
                    # parity check: q_tot must be integer (it always is
                    # since q1_x2 and q2_x2 are both even integers as
                    # m_j are half-integers and differences between like
                    # half-integers are integers -> q1_x2 even).
                    if q_tot_x2 % 2 != 0:
                        continue
                    if abs(q_tot_x2) > 2 * K_overall:
                        continue
                    # at least one (L1, L2) pair must jointly permit
                    any_allowed = False
                    for L1, L2 in bipolar_ranks:
                        if pair_ok[L1].get((a, c), False) and \
                           pair_ok[L2].get((b, d), False):
                            any_allowed = True
                            break
                    if any_allowed:
                        nonzero += 1
    return nonzero, total


# ---------------------------------------------------------------------------
# Operator-specific rank patterns
# ---------------------------------------------------------------------------

def rank_list_coulomb(l_max: int) -> List[Tuple[int, int]]:
    """K=0 scalar from Coulomb multipole: (L1=L2=k) for k = 0..2*l_max.

    This reproduces Paper 22 spinor fullgaunt density.
    """
    return [(k, k) for k in range(2 * l_max + 1)]


def rank_list_SS_tensor(l_max: int) -> List[Tuple[int, int]]:
    """The rank-2 spin-spin tensor part of Breit-Pauli (the Y_2 term).

    The orbital tensor attached to [sigma1 x sigma2]^2 is
       Y_2(rhat_12) / r_12^3.
    Multipole-expanding Y_2(rhat_12)/r_12^3 gives a bipolar sum
       sum_{L1, L2} C_{L1 L2}^2 [Y_{L1}(rhat_1) (x) Y_{L2}(rhat_2)]^2
    with L1 + L2 >= 2 and triangle(L1, L2, 2). Parity requires L1 + L2 even.

    Allowed lowest orders: (L1, L2) in {(0,2), (1,1), (2,0)}.

    Higher multipoles (L1, L2) = (2,4), (4,2), (3,3), ... contribute to
    higher-k expansions; but their *selection-rule pattern* is
    additive with the low orders, and the ALLOWED set of quartets is the
    union. We include all (L1, L2) with triangle(L1, L2, 2) and parity
    (L1+L2) even, up to L = 2*l_max + 2 (sufficient to hit all pairs of
    orbital labels with l <= l_max).
    """
    pairs: List[Tuple[int, int]] = []
    Lmax = 2 * l_max + 2
    for L1 in range(Lmax + 1):
        for L2 in range(Lmax + 1):
            if abs(L1 - L2) > 2 or L1 + L2 < 2:
                continue
            if (L1 + L2) % 2 != 0:
                continue
            pairs.append((L1, L2))
    return pairs


def rank_list_SS_scalar(l_max: int) -> List[Tuple[int, int]]:
    """The rank-0 scalar sigma1.sigma2 part of Breit-Pauli.

    Orbital tensor: 1/r_12^3 -> sum_k [Y_k (x) Y_k]^0 multipole.
    Coupled with [sigma1.sigma2] scalar: same pattern as Coulomb but with
    1/r_12^3 radial weight. Selection rule on angular part identical to
    Coulomb up to the scalar spin factor being nonzero everywhere (it is;
    sigma1.sigma2 = 2S(S+1) - 3 has no angular selection rule).
    """
    return rank_list_coulomb(l_max)


def rank_list_SOO(l_max: int) -> List[Tuple[int, int]]:
    """The spin-other-orbit piece has structure
       (grad 1/r_12) . (ell_1 sigma_2 + ell_2 sigma_1 + ...)

    The orbital operator (grad 1/r_12) carries spatial rank 1 per "leg";
    after multipole expansion the bipolar structure is
       (L1, L2) with triangle(L1, L2, 1) and (L1 + L2) odd,
    i.e. (L1, L2) in {(0,1), (1,0), (1,2), (2,1), (2,3), ...}.

    Overall coupling rank K_overall = 1 (rank of each spin operator).
    """
    pairs: List[Tuple[int, int]] = []
    Lmax = 2 * l_max + 2
    for L1 in range(Lmax + 1):
        for L2 in range(Lmax + 1):
            if abs(L1 - L2) > 1 or L1 + L2 < 1:
                continue
            if (L1 + L2) % 2 != 1:
                continue
            pairs.append((L1, L2))
    return pairs


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def compute_all(l_values: List[int]) -> List[dict]:
    rows: List[dict] = []
    for l_max in l_values:
        print(f"\n--- l_max = {l_max} ---")
        orbs = enumerate_spinor_orbitals(l_max)
        Q = len(orbs)
        print(f"  Q = {Q} spinor orbitals; Q^4 = {Q**4} quartets")

        # Scalar (Coulomb reproduction): should match Paper 22 spinor fullgaunt
        t0 = time.time()
        ranks_C = rank_list_coulomb(l_max)
        nz_C, tot_C = bipolar_density(orbs, ranks_C, K_overall=0)
        t1 = time.time()
        d_C = Fraction(nz_C, tot_C)
        print(f"  Coulomb (rank-0):      nz={nz_C:>10}/{tot_C:<10}  "
              f"d={float(d_C)*100:.4f}%   [{t1-t0:.2f}s]")

        # SS tensor (rank-2)
        ranks_SS2 = rank_list_SS_tensor(l_max)
        nz_SS2, tot_SS2 = bipolar_density(orbs, ranks_SS2, K_overall=2)
        t2 = time.time()
        d_SS2 = Fraction(nz_SS2, tot_SS2)
        print(f"  SS tensor (rank-2):    nz={nz_SS2:>10}/{tot_SS2:<10}  "
              f"d={float(d_SS2)*100:.4f}%   [{t2-t1:.2f}s]")

        # SS scalar (rank-0), identical allowed set to Coulomb
        ranks_SSS = rank_list_SS_scalar(l_max)
        nz_SSS, tot_SSS = bipolar_density(orbs, ranks_SSS, K_overall=0)
        t3 = time.time()
        d_SSS = Fraction(nz_SSS, tot_SSS)
        print(f"  SS scalar (rank-0):    nz={nz_SSS:>10}/{tot_SSS:<10}  "
              f"d={float(d_SSS)*100:.4f}%   [{t3-t2:.2f}s]")

        # SOO (rank-1)
        ranks_SOO = rank_list_SOO(l_max)
        nz_SOO, tot_SOO = bipolar_density(orbs, ranks_SOO, K_overall=1)
        t4 = time.time()
        d_SOO = Fraction(nz_SOO, tot_SOO)
        print(f"  SOO (rank-1):          nz={nz_SOO:>10}/{tot_SOO:<10}  "
              f"d={float(d_SOO)*100:.4f}%   [{t4-t3:.2f}s]")

        # Breit = union of SS_tensor + SS_scalar + SOO + Darwin-like (scalar).
        # Compute the union by recomputing with combined rank list + all K.
        # But K differs across SS_tensor (K=2) and others (K=0 or 1) -> the
        # m-conservation constraint depends on K. So compute union by
        # directly OR-ing allowed quartets.
        t4b = time.time()
        nz_B, tot_B = breit_union_density(orbs, l_max)
        t5 = time.time()
        d_B = Fraction(nz_B, tot_B)
        print(f"  Breit (union):         nz={nz_B:>10}/{tot_B:<10}  "
              f"d={float(d_B)*100:.4f}%   [{t5-t4b:.2f}s]")

        rows.append({
            "l_max": l_max,
            "Q_spinor": Q,
            "Q_spinor_pow4": Q ** 4,
            "coulomb":   _pack(nz_C, tot_C, d_C),
            "SS_tensor": _pack(nz_SS2, tot_SS2, d_SS2),
            "SS_scalar": _pack(nz_SSS, tot_SSS, d_SSS),
            "SOO":       _pack(nz_SOO, tot_SOO, d_SOO),
            "Breit":     _pack(nz_B, tot_B, d_B),
            "breit_over_coulomb_ratio": (
                float(Fraction(nz_B, nz_C)) if nz_C > 0 else None
            ),
        })
    return rows


def breit_union_density(
    orbs: List[SpinorOrbital], l_max: int
) -> Tuple[int, int]:
    """Count quartets allowed by the union
        { Coulomb (K=0) } ∪ { SS tensor (K=2) } ∪ { SOO (K=1) }.

    This is the raw Breit-Pauli angular selection envelope.
    """
    Q = len(orbs)
    total = Q ** 4

    # Precompute all ranks needed across SS_tensor, SOO, Coulomb.
    used_L: Set[int] = set()
    ranks_C = rank_list_coulomb(l_max)
    ranks_SS2 = rank_list_SS_tensor(l_max)
    ranks_SOO_list = rank_list_SOO(l_max)
    for lst in (ranks_C, ranks_SS2, ranks_SOO_list):
        for L1, L2 in lst:
            used_L.add(L1)
            used_L.add(L2)

    pair_ok: Dict[int, Dict[Tuple[int, int], bool]] = {}
    for L in used_L:
        pair_ok[L] = pair_rank_allowed(orbs, L)

    def allowed_any(a, b, c, d, q_tot_x2, K_overall, rank_list):
        if abs(q_tot_x2) > 2 * K_overall:
            return False
        for L1, L2 in rank_list:
            if pair_ok[L1].get((a, c), False) and \
               pair_ok[L2].get((b, d), False):
                return True
        return False

    nonzero = 0
    for a in range(Q):
        ma_x2 = orbs[a].two_mj
        for b in range(Q):
            mb_x2 = orbs[b].two_mj
            for c in range(Q):
                mc_x2 = orbs[c].two_mj
                q1_x2 = ma_x2 - mc_x2
                for d in range(Q):
                    md_x2 = orbs[d].two_mj
                    q2_x2 = mb_x2 - md_x2
                    q_tot_x2 = q1_x2 + q2_x2
                    if q_tot_x2 % 2 != 0:
                        continue
                    # Breit union: any of Coulomb / SS2 / SOO survives
                    if allowed_any(a, b, c, d, q_tot_x2, 0, ranks_C):
                        nonzero += 1
                        continue
                    if allowed_any(a, b, c, d, q_tot_x2, 2, ranks_SS2):
                        nonzero += 1
                        continue
                    if allowed_any(a, b, c, d, q_tot_x2, 1, ranks_SOO_list):
                        nonzero += 1
                        continue
    return nonzero, total


def _pack(nz: int, tot: int, d: Fraction) -> dict:
    return {
        "nonzero": nz,
        "total": tot,
        "density_num": d.numerator,
        "density_den": d.denominator,
        "density_pct": float(d) * 100.0,
    }


# ---------------------------------------------------------------------------
# Qualitative sample: show a handful of allowed/forbidden quartets
# ---------------------------------------------------------------------------

def sample_selection_rules(l_max: int = 1):
    """Print a few explicit quartets showing which are allowed/forbidden
    under each rank.

    Focus on l_max = 1: orbits are s_{1/2}, p_{1/2}, p_{3/2} with their
    m_j manifolds (Q = 8 spinor orbitals).
    """
    orbs = enumerate_spinor_orbitals(l_max)
    print(f"\n  Sample quartets at l_max = {l_max} (Q = {len(orbs)}):")
    print(f"  {'a':>10} {'b':>10} {'c':>10} {'d':>10}   "
          f"{'C':>4} {'SS2':>4} {'SOO':>4}")
    print("  " + "-" * 64)

    pair_ok_0 = pair_rank_allowed(orbs, 0)
    pair_ok_1 = pair_rank_allowed(orbs, 1)
    pair_ok_2 = pair_rank_allowed(orbs, 2)

    def _label(o: SpinorOrbital) -> str:
        # Render as (kappa,2mj)
        return f"({o.kappa:+d},{o.two_mj:+d})"

    # Show a small selection of diagonal quartets and a few off-diagonal.
    samples = []
    # diagonal: (a=b=c=d) for each orbital
    for i in range(len(orbs)):
        samples.append((i, i, i, i))
    # s - p mixing quartet
    for (a, b, c, d) in [(0, 4, 0, 4), (0, 4, 2, 6), (0, 5, 1, 4),
                         (2, 6, 0, 4), (1, 1, 3, 3), (0, 1, 4, 5),
                         (2, 3, 6, 7), (0, 3, 4, 7), (0, 7, 3, 4)]:
        if all(k < len(orbs) for k in (a, b, c, d)):
            samples.append((a, b, c, d))

    for a, b, c, d in samples[:18]:
        oa, ob, oc, od = orbs[a], orbs[b], orbs[c], orbs[d]
        q1_x2 = oa.two_mj - oc.two_mj
        q2_x2 = ob.two_mj - od.two_mj
        q_tot_x2 = q1_x2 + q2_x2

        # Coulomb K=0
        okC = (abs(q_tot_x2) <= 0) and any(
            pair_ok_0[(a, c)] and pair_ok_0[(b, d)] for _ in [0]  # L=0
        ) or any(
            pair_ok_1[(a, c)] and pair_ok_1[(b, d)] for _ in [0]
        ) or any(
            pair_ok_2[(a, c)] and pair_ok_2[(b, d)] for _ in [0]
        )
        # Actually correct: Coulomb is K=0 with L1=L2, q_tot=0.
        okC = (abs(q_tot_x2) <= 0) and (
            (pair_ok_0[(a, c)] and pair_ok_0[(b, d)])
            or (pair_ok_1[(a, c)] and pair_ok_1[(b, d)])
            or (pair_ok_2[(a, c)] and pair_ok_2[(b, d)])
        )
        # SS tensor K=2 with (L1,L2) in {(0,2),(1,1),(2,0)}
        okSS2 = (abs(q_tot_x2) <= 4) and (
            (pair_ok_0[(a, c)] and pair_ok_2[(b, d)])
            or (pair_ok_1[(a, c)] and pair_ok_1[(b, d)])
            or (pair_ok_2[(a, c)] and pair_ok_0[(b, d)])
        )
        # SOO K=1 with (L1,L2) in {(0,1),(1,0),(1,2),(2,1)}
        okSOO = (abs(q_tot_x2) <= 2) and (
            (pair_ok_0[(a, c)] and pair_ok_1[(b, d)])
            or (pair_ok_1[(a, c)] and pair_ok_0[(b, d)])
            or (pair_ok_1[(a, c)] and pair_ok_2[(b, d)])
            or (pair_ok_2[(a, c)] and pair_ok_1[(b, d)])
        )
        print(f"  {_label(oa):>10} {_label(ob):>10} {_label(oc):>10} "
              f"{_label(od):>10}   {'*' if okC else '.':>4} "
              f"{'*' if okSS2 else '.':>4} {'*' if okSOO else '.':>4}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(l_values=None):
    if l_values is None:
        l_values = [0, 1, 2, 3, 4, 5]

    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(exist_ok=True)

    print("=" * 78)
    print("Track BR-A: Breit interaction angular density d_Breit(l_max)")
    print("  [rank-0 (Coulomb) + rank-2 (SS tensor) + rank-1 (SOO)]")
    print("=" * 78)

    rows = compute_all(l_values)

    # Show a qualitative sample at l_max = 1
    sample_selection_rules(l_max=1)

    # Summary table
    print("\n" + "=" * 78)
    print("SUMMARY TABLE: Breit-Pauli angular density vs Coulomb baseline")
    print("=" * 78)
    print(f"| {'l_max':>5} | {'Q':>4} | {'d_Coulomb':>10} | "
          f"{'d_SS2':>8} | {'d_SS0':>8} | {'d_SOO':>8} | "
          f"{'d_Breit':>10} | {'Breit/C':>8} |")
    print("|" + "-" * 7 + "|" + "-" * 6 + "|" + "-" * 12 + "|"
          + "-" * 10 + "|" + "-" * 10 + "|" + "-" * 10 + "|"
          + "-" * 12 + "|" + "-" * 10 + "|")
    for r in rows:
        print(f"| {r['l_max']:>5} | {r['Q_spinor']:>4} | "
              f"{r['coulomb']['density_pct']:>9.4f}% | "
              f"{r['SS_tensor']['density_pct']:>7.4f}% | "
              f"{r['SS_scalar']['density_pct']:>7.4f}% | "
              f"{r['SOO']['density_pct']:>7.4f}% | "
              f"{r['Breit']['density_pct']:>9.4f}% | "
              f"{(r['breit_over_coulomb_ratio'] or 0):>8.2f}x |")

    metadata = {
        "task": "Track BR-A: Breit-Pauli angular sparsity extension",
        "operator": "Breit-Pauli (SS_tensor + SS_scalar + SOO)",
        "basis": "jj-coupled spinor orbitals (kappa, m_j)",
        "arithmetic": "exact sympy Rational + wigner_3j",
        "convention": "full-Gaunt (physically correct Coulomb selection)",
        "n_max_independent": True,
        "notes": {
            "Coulomb_ranks": "(L1=L2=k) K=0 bipolar; reproduces Paper 22 spinor fullgaunt.",
            "SS_tensor_ranks": "(L1,L2) with triangle(L1,L2,2), (L1+L2) even, L1+L2>=2; K=2.",
            "SS_scalar_ranks": "Same as Coulomb (different radial weight 1/r^3, identical angular).",
            "SOO_ranks": "(L1,L2) with triangle(L1,L2,1), (L1+L2) odd; K=1.",
            "Breit_union": "quartets allowed by ANY of Coulomb / SS_tensor / SOO.",
        },
        "rows": rows,
    }
    out_path = out_dir / "br_a_density.json"
    with open(out_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        lvals = [int(x) for x in sys.argv[1:]]
        main(lvals)
    else:
        main()
