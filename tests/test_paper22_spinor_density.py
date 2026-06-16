"""
Paper 22 (Potential-Independent Angular Sparsity) — spinor Table II backing.

Paper 22 §"Numerical values of d_spinor(l_max)" (Table II, `tab:d_spinor`)
reports the angular ERI density of the RELATIVISTIC jj-coupled spinor basis
{|n, kappa, m_j> : l(kappa) <= l_max}, in two conventions, by exact enumeration
of the jj-coupled Coulomb selection rule (Dyall §9.3, Grant §7.5):

  * d_sp^FG  — full-Gaunt (physically correct) spinor density.  An ERI survives
       iff global m_j-conservation m_a+m_b = m_c+m_d holds AND there is a shared
       multipole k for which BOTH pairs carry a nonzero reduced jj 3j
       (j k j'; 1/2 0 -1/2), the l+l'+k-even parity, the l- and j-triangles,
       and the bottom-row jj 3j (j k j'; -m_j, m_j-m_j', m_j').
       Paper values: 25.0000 / 8.5938 / 6.4586 / 5.1674 / 4.3034 / 3.6794 %.

  * d_sp     — the STRICTER pair-diagonal restriction (m_a = m_c AND m_b = m_d),
       the convention Paper 22's published scalar table uses.
       Paper values: 25.0000 / 4.2969 / 2.1071 / 1.2329 / 0.8106 / 0.5722 %.

This test pins BOTH spinor columns of Table II to the paper's stated nonzero
counts and 4-decimal percentages by independent exact sympy enumeration over
the jj-coupled basis (no floats enter the zero/nonzero decision), and verifies
the three structural observations the paper draws from the table:
  (i)  d_spinor <= d_scalar at every l_max (j-triangles only ADD zeros);
  (ii) the pair-diagonal survivors are a subset of the full-Gaunt survivors;
  (iii) Q_spinor = 2 (l_max+1)^2.

The scalar full-Gaunt reference D it compares against is the headline density
pinned by the companion test tests/test_paper22_density.py (GLOBAL_D), so the
two tests are mutually consistent on the scalar/spinor boundary.

Provenance: SYMBOLIC PROOF (exact integer nonzero counts via sympy jj 3j).
The discovery-mode driver lived in debug/tier2_t0_spinor_density.py (Dirac-on-S^3
Tier 2 Track T0), now archived at
debug/archive/spectral_triple_arc/tier2_t0_spinor_density.py; this test is the
live, default-collected backing for the published Table II.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Dict, FrozenSet, List, Tuple

import pytest
from sympy import Integer, Rational
from sympy.physics.wigner import wigner_3j


# ---------------------------------------------------------------------------
# Reference values from Paper 22 Table II (`tab:d_spinor`).
# l_max -> (total = Q_spinor^4, exact nonzero count, 4-decimal percentage).
# ---------------------------------------------------------------------------

# Full-Gaunt spinor density d_sp^FG (physically correct Coulomb selection).
SPINOR_FG = {
    0: (16, 4, "25.0000"),
    1: (4096, 352, "8.5938"),
    2: (104976, 6780, "6.4586"),
    3: (1048576, 54184, "5.1674"),
    4: (6250000, 268964, "4.3034"),
    5: (26873856, 988800, "3.6794"),
}

# Pair-diagonal spinor density d_sp (stricter m_a=m_c, m_b=m_d convention).
SPINOR_PD = {
    0: (16, 4, "25.0000"),
    1: (4096, 176, "4.2969"),
    2: (104976, 2212, "2.1071"),
    3: (1048576, 12928, "1.2329"),
    4: (6250000, 50660, "0.8106"),
    5: (26873856, 153776, "0.5722"),
}

# Scalar full-Gaunt reference D (headline density), pinned in
# tests/test_paper22_density.py::GLOBAL_D — used here only for the
# d_spinor <= d_scalar monotonicity check (observation i).
SCALAR_FG_NONZERO = {0: 1, 1: 38, 2: 559, 3: 3972, 4: 18857, 5: 66954}
SCALAR_Q4 = {0: 1, 1: 256, 2: 6561, 3: 65536, 4: 390625, 5: 1679616}


# ---------------------------------------------------------------------------
# jj-coupled spinor orbital enumeration (kappa, m_j), ported from the Track T0
# driver.  l and j are derived from kappa; m_j is stored doubled (2 m_j).
# ---------------------------------------------------------------------------

def _l_of_kappa(kappa: int) -> int:
    """l from kappa: kappa>0 -> j=l-1/2, l=kappa; kappa<0 -> j=l+1/2, l=-kappa-1."""
    if kappa > 0:
        return kappa
    if kappa < 0:
        return -kappa - 1
    raise ValueError("kappa cannot be 0")


def _two_j_of_kappa(kappa: int) -> int:
    """2j = 2|kappa| - 1 (j = |kappa| - 1/2)."""
    return 2 * abs(kappa) - 1


def _spinor_orbitals(l_max: int) -> List[Tuple[int, int, int]]:
    """All (l, 2j, 2m_j) with l(kappa) <= l_max; total = 2 (l_max+1)^2 states.

    Negative-kappa branch j=l+1/2 for l=0..l_max; positive-kappa branch
    j=l-1/2 for l=1..l_max.
    """
    orbs: List[Tuple[int, int, int]] = []
    for l in range(l_max + 1):              # kappa = -(l+1), j = l+1/2
        two_j = 2 * l + 1
        for two_mj in range(-two_j, two_j + 1, 2):
            orbs.append((l, two_j, two_mj))
    for l in range(1, l_max + 1):           # kappa = +l, j = l-1/2
        two_j = 2 * l - 1
        for two_mj in range(-two_j, two_j + 1, 2):
            orbs.append((l, two_j, two_mj))
    return orbs


def _triangle_x2(a2: int, b2: int, c2: int) -> bool:
    """Triangle on doubled angular momenta: |a-b|<=c<=a+b and a+b+c even."""
    return abs(a2 - b2) <= c2 <= a2 + b2 and (a2 + b2 + c2) % 2 == 0


def _spinor_pair_k(orbs: List[Tuple[int, int, int]], l_max: int,
                   convention: str) -> Dict[Tuple[int, int], FrozenSet[int]]:
    """Per-pair (a,c) set of multipoles k for which the jj 3j survive.

    Both conventions require l-parity, l-triangle, j-triangle, and the reduced
    jj 3j (j k j'; 1/2 0 -1/2) != 0.  They differ in the bottom-row 3j:
      fullgaunt: q = m_a - m_c, 3j(j k j'; -m_a, q, m_c) != 0
      pairdiag:  require m_a = m_c (q=0), 3j(j k j'; -m_a, 0, m_c) != 0
    """
    pair_k: Dict[Tuple[int, int], FrozenSet[int]] = {}
    k_cap = 2 * l_max + 1
    for i, (la, ja2, ma2) in enumerate(orbs):
        ja = Rational(ja2, 2)
        for j, (lc, jc2, mc2) in enumerate(orbs):
            if convention == "pairdiag" and ma2 != mc2:
                pair_k[(i, j)] = frozenset()
                continue
            jc = Rational(jc2, 2)
            ks = set()
            for k in range(k_cap + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                if k < abs(la - lc) or k > la + lc:
                    continue
                if not _triangle_x2(ja2, 2 * k, jc2):
                    continue
                ksym = Integer(k)
                if wigner_3j(ja, ksym, jc, Rational(1, 2), Integer(0),
                             Rational(-1, 2)) == 0:
                    continue
                if convention == "fullgaunt":
                    q2 = ma2 - mc2
                    if q2 % 2 != 0:
                        continue
                    q = q2 // 2
                    if abs(q) > k:
                        continue
                    w_m = wigner_3j(ja, ksym, jc, Rational(-ma2, 2),
                                    Integer(q), Rational(mc2, 2))
                else:  # pairdiag, m_a = m_c already enforced
                    w_m = wigner_3j(ja, ksym, jc, Rational(-ma2, 2),
                                    Integer(0), Rational(mc2, 2))
                if w_m == 0:
                    continue
                ks.add(k)
            pair_k[(i, j)] = frozenset(ks)
    return pair_k


def _spinor_count(l_max: int, convention: str) -> Tuple[int, int]:
    """(nonzero, total) surviving four-index spinor ERIs over the single shell."""
    orbs = _spinor_orbitals(l_max)
    Q = len(orbs)
    total = Q ** 4
    pair_k = _spinor_pair_k(orbs, l_max, convention)
    mj2 = [o[2] for o in orbs]
    nonzero = 0
    for a in range(Q):
        ma2 = mj2[a]
        for b in range(Q):
            mb2 = mj2[b]
            for c in range(Q):
                k_ac = pair_k[(a, c)]
                if not k_ac:
                    continue
                mc2 = mj2[c]
                for d in range(Q):
                    if convention == "fullgaunt" and ma2 + mb2 != mc2 + mj2[d]:
                        continue
                    k_bd = pair_k[(b, d)]
                    if k_bd and (k_ac & k_bd):
                        nonzero += 1
    return nonzero, total


def _pct(nz: int, total: int) -> str:
    return f"{100.0 * nz / total:.4f}"


# ---------------------------------------------------------------------------
# Q_spinor = 2 (l_max+1)^2 (Table II header).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l_max", [0, 1, 2, 3, 4, 5])
def test_spinor_orbital_count(l_max):
    """The jj-coupled shell carries exactly 2 (l_max+1)^2 orbitals."""
    assert len(_spinor_orbitals(l_max)) == 2 * (l_max + 1) ** 2


# ---------------------------------------------------------------------------
# Full-Gaunt spinor density d_sp^FG — exact jj 3j counts.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l_max", [0, 1, 2])
def test_spinor_fullgaunt_density(l_max):
    """d_sp^FG(l_max): exact nonzero count + percentage, full-Gaunt jj rule."""
    total_ref, nz_ref, pct_ref = SPINOR_FG[l_max]
    nz, total = _spinor_count(l_max, "fullgaunt")
    assert total == total_ref
    assert nz == nz_ref, f"d_sp^FG nonzero l_max={l_max}: got {nz}, want {nz_ref}"
    assert _pct(nz, total) == pct_ref


@pytest.mark.slow
@pytest.mark.parametrize("l_max", [3, 4, 5])
def test_spinor_fullgaunt_density_high_lmax(l_max):
    """d_sp^FG(l_max=3,4,5): exact nonzero count + percentage (slow)."""
    total_ref, nz_ref, pct_ref = SPINOR_FG[l_max]
    nz, total = _spinor_count(l_max, "fullgaunt")
    assert total == total_ref
    assert nz == nz_ref
    assert _pct(nz, total) == pct_ref


# ---------------------------------------------------------------------------
# Pair-diagonal spinor density d_sp — exact jj 3j counts.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l_max", [0, 1, 2])
def test_spinor_pairdiagonal_density(l_max):
    """d_sp(l_max): exact nonzero count + percentage, pair-diagonal jj rule."""
    total_ref, nz_ref, pct_ref = SPINOR_PD[l_max]
    nz, total = _spinor_count(l_max, "pairdiag")
    assert total == total_ref
    assert nz == nz_ref, f"d_sp nonzero l_max={l_max}: got {nz}, want {nz_ref}"
    assert _pct(nz, total) == pct_ref


@pytest.mark.slow
@pytest.mark.parametrize("l_max", [3, 4, 5])
def test_spinor_pairdiagonal_density_high_lmax(l_max):
    """d_sp(l_max=3,4,5): exact nonzero count + percentage (slow)."""
    total_ref, nz_ref, pct_ref = SPINOR_PD[l_max]
    nz, total = _spinor_count(l_max, "pairdiag")
    assert total == total_ref
    assert nz == nz_ref
    assert _pct(nz, total) == pct_ref


# ---------------------------------------------------------------------------
# Structural observations the paper draws from Table II.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("l_max", [0, 1, 2])
def test_spinor_at_most_scalar_density(l_max):
    """Observation (i): d_spinor <= d_scalar at every l_max, both conventions.

    The j-triangle rules can only add zeros on top of the l-triangle/parity
    rules; the spinor density is bounded above by the scalar density.  Compared
    in normalized form (nonzero/Q^4) against the headline scalar D pinned by
    test_paper22_density.py::GLOBAL_D.
    """
    nz_sp_fg = SPINOR_FG[l_max][1]
    d_sp_fg = Fraction(nz_sp_fg, SPINOR_FG[l_max][0])
    d_sc_fg = Fraction(SCALAR_FG_NONZERO[l_max], SCALAR_Q4[l_max])
    assert d_sp_fg <= d_sc_fg

    nz_sp_pd = SPINOR_PD[l_max][1]
    d_sp_pd = Fraction(nz_sp_pd, SPINOR_PD[l_max][0])
    assert d_sp_pd <= d_sp_fg  # pair-diagonal is the stricter sub-case


def test_spinor_pairdiagonal_is_subset_of_fullgaunt():
    """Observation (ii): every pair-diagonal spinor survivor also survives the
    full-Gaunt rule (the pair-diagonal set is a subset), so d_sp <= d_sp^FG."""
    for l_max in (0, 1, 2):
        orbs = _spinor_orbitals(l_max)
        Q = len(orbs)
        mj2 = [o[2] for o in orbs]
        pk_fg = _spinor_pair_k(orbs, l_max, "fullgaunt")
        pk_pd = _spinor_pair_k(orbs, l_max, "pairdiag")
        for a in range(Q):
            for b in range(Q):
                for c in range(Q):
                    if not pk_pd[(a, c)]:
                        continue
                    for d in range(Q):
                        pd_hit = bool(pk_pd[(a, c)] & pk_pd[(b, d)])
                        if not pd_hit:
                            continue
                        # m-conservation is automatic under pair-diagonal
                        # (m_a=m_c, m_b=m_d) so the full-Gaunt rule must agree.
                        fg_hit = (mj2[a] + mj2[b] == mj2[c] + mj2[d]) and bool(
                            pk_fg[(a, c)] & pk_fg[(b, d)])
                        assert fg_hit, (
                            f"pairdiag survivor ({a},{b},{c},{d}) at l_max="
                            f"{l_max} not in full-Gaunt set")


def test_spinor_floor_at_lmax0():
    """l_max=0 spin-dilution floor: a single s_{1/2} orbital with two m_j
    branches gives exactly 4/16 = 25% in both conventions (Table II row 0)."""
    assert _spinor_count(0, "fullgaunt") == (4, 16)
    assert _spinor_count(0, "pairdiag") == (4, 16)
