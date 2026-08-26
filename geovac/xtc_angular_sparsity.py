"""Radial-free, exact angular-support demonstrator for the atomic xTC operator.

Backing for Paper 14 (sec:tc_atomic_sparsity): when the transcorrelated (TC)
three-body operator is contracted to an effective two-body operator against the
reference one-particle density (the xTC scheme), does the contracted operator
*inherit* GeoVac's angular Gaunt sparsity, or does the 3-body -> 2-body
contraction fill in the zero angular blocks?

This module answers the question at the level of angular *support* -- which
``(a,b,c,d)`` orbital quadruplets are nonzero -- which is radial-independent and
therefore EXACT (integer block counts, no grid, no quadrature).  The radial
magnitudes (1-norm, per-term coefficients) are measured separately on the
grid-based atomic engine; here we settle the structural question the s-only PoC
could not: with p-and-higher angular momentum actually present, is the sparsity
pattern preserved?

Support definitions (identical to the validated debug engine
``debug/xtc_poc_li_pinclusive.py`` / ``debug/xtc_pblock_fast_angular.py``):

  Coulomb block (a,b,c,d) nonzero  iff  exists (L,M):
      gA(a; L,M; c) != 0   AND   gB(b; L,M; d) != 0
  xTC-L3 block (a,b,c,d) nonzero  iff  exists (L,M,L',M'):
      four_Y(a; L,M,L',M'; c) != 0     (vertex, external a<->c)
      AND gA(b; L,M; d) != 0           (leg-2, external b<->d)
      AND refmult(L',M') != 0          (leg-3 contracted over the reference density)
  refmult(L',M') = sum_o gA(l_o,m_o; L',M'; l_o,m_o)  (diagonal reference-density multipole)

"fill-in" = a block that xTC makes nonzero where Coulomb was zero.  "fill-out" =
the reverse.  Zero of both means the support (hence the Jordan-Wigner Pauli term
set) is preserved exactly.

Structural result (see memos ``debug/sprint_xtc_pinclusive_memo.md`` and
``debug/sprint_xtc_pblock_memo.md``):

  * s-reference (closed/spherical density, multipole {(0,0)}): 0 fill-in at ANY
    basis, bit-exact -- the correlator leg is forced to the L'=0 monopole, so the
    non-abelian four-harmonic vertex collapses to a single Coulomb-like multipole.
  * open-p reference (density multipole {(0,0),(2,0)}): 0 fill-in through s+p, then
    a small, structured, slowly-growing fill-in once d orbitals enter (4 blocks at
    s+p+d, 96 at s+p+d+f) -- density stays single-digit % and tracks Coulomb.

The whole chain rests on the tracked ``geovac.angular_integrals.wigner3j``.
"""
from __future__ import annotations

from math import pi, sqrt
from typing import Dict, List, Optional, Sequence, Tuple

from geovac.angular_integrals import wigner3j

_FOURPI = 4.0 * pi

# angular (l, m) key and a per-line multipole key
LM = Tuple[int, int]


# ---------------------------------------------------------------------------
# complex-spherical-harmonic Gaunt coefficients (exact, via wigner3j)
# ---------------------------------------------------------------------------
def gaunt(l1: int, m1: int, l2: int, m2: int, l3: int, m3: int) -> float:
    r"""$\int Y_{l_1 m_1} Y_{l_2 m_2} Y_{l_3 m_3}\,d\Omega$ (complex harmonics)."""
    if (m1 + m2 + m3) != 0:
        return 0.0
    pref = sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / _FOURPI)
    return pref * wigner3j(l1, l2, l3, 0, 0, 0) * wigner3j(l1, l2, l3, m1, m2, m3)


def g_ext(la: int, ma: int, Lam: int, MLam: int, lc: int, mc: int) -> float:
    r"""$\int Y^*_{l_a m_a} Y_{\Lambda M_\Lambda} Y_{l_c m_c}\,d\Omega$."""
    return (-1) ** ma * gaunt(la, -ma, Lam, MLam, lc, mc)


def g_pair(L: int, M: int, Lp: int, Mp: int, Lam: int) -> float:
    r"""Coefficient of $Y_{\Lambda,M+M'}$ in the product $Y_{LM} Y_{L'M'}$."""
    MLam = M + Mp
    if abs(MLam) > Lam:
        return 0.0
    return (-1) ** (MLam) * gaunt(L, M, Lp, Mp, Lam, -MLam)


def four_Y(la: int, ma: int, L: int, M: int, Lp: int, Mp: int, lc: int, mc: int,
           keep_lowest_only: bool = False) -> Tuple[float, List[int]]:
    r"""Shared-vertex four-harmonic integral, exact, plus its active $\Lambda$ list.

    This is the non-abelian vertex of the transcorrelated three-body operator: two
    harmonics ($Y_{LM}$ from the leg, $Y_{L'M'}$ from the correlator) Clebsch--Gordan
    couple to a *multiplet* of resultants $\Lambda$, not a single one -- which is why
    the 3-body operator does not collapse to 2-body the way the abelian plane-wave TC
    does (Paper 14, three-body-collapse result).  ``keep_lowest_only`` keeps only the
    lowest resultant (the plane-wave mimic control).
    """
    MLam = M + Mp
    if MLam != (ma - mc):  # m-conservation at the vertex: -ma + M + Mp + mc = 0
        return 0.0, []
    val = 0.0
    active: List[int] = []
    Lam_lo = max(abs(L - Lp), abs(la - lc), abs(MLam))
    Lam_hi = min(L + Lp, la + lc)
    for Lam in range(Lam_lo, Lam_hi + 1):
        gp = g_pair(L, M, Lp, Mp, Lam)
        if gp == 0.0:
            continue
        ge = g_ext(la, ma, Lam, MLam, lc, mc)
        if ge == 0.0:
            continue
        if keep_lowest_only and active:
            continue
        val += gp * ge
        active.append(Lam)
    return val, active


# ---------------------------------------------------------------------------
# leg factors (correlator lines; complex-Y Gaunt convention consistent with four_Y)
# ---------------------------------------------------------------------------
def gA(l: int, m: int, L: int, M: int, lp: int, mp: int) -> float:
    r"""$\int Y^*_{lm} Y^*_{LM} Y_{l'm'}\,d\Omega$ (conjugated correlator leg)."""
    return (-1) ** (m + M) * gaunt(l, -m, L, -M, lp, mp)


def gB(l: int, m: int, L: int, M: int, lp: int, mp: int) -> float:
    r"""$\int Y^*_{lm} Y_{LM} Y_{l'm'}\,d\Omega$ (un-conjugated correlator leg)."""
    return (-1) ** (m) * gaunt(l, -m, L, M, lp, mp)


# ---------------------------------------------------------------------------
# exact angular-support / fill-in counter
# ---------------------------------------------------------------------------
def orbs(lmax: int) -> List[LM]:
    """One radial function per l up to ``lmax``: the (l, m) angular index list."""
    return [(l, m) for l in range(lmax + 1) for m in range(-l, l + 1)]


# reference angular occupations (l, m), with multiplicity, for common atoms
REF_ANG: Dict[str, List[LM]] = {
    "Be_s2":  [(0, 0), (0, 0), (0, 0), (0, 0)],                                   # 1s^2 2s^2 (s-ref)
    "C_2p2":  [(0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 1)],                   # 1s^2 2s^2 2p^2 (Hund 3P)
    "C_2p0sq": [(0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 0)],                  # 2p0^2 variant
    "O_2p4":  [(0, 0), (0, 0), (0, 0), (0, 0), (1, -1), (1, -1), (1, 0), (1, 0)],  # 1s^2 2s^2 2p^4
}


def open_subshell(l: int, k: int) -> List[LM]:
    """(l, m) spatial-orbital occupation of ``k`` electrons in an l-subshell, high spin.

    Fill each m once (Hund, spin up) in the order 0, +1, -1, +2, -2, ... then pair.
    Only the SET of occupied spatial orbitals (with multiplicity) enters the diagonal
    density multipoles, and the specific m-ordering does not change angular *support* --
    so this is a faithful representative of the atom's open-shell angular character for
    the fill-in question.
    """
    ms = [0]
    for a in range(1, l + 1):
        ms += [a, -a]
    occ: List[LM] = [(l, ms[i]) for i in range(min(k, 2 * l + 1))]
    occ += [(l, ms[i]) for i in range(k - (2 * l + 1))]  # pair up the remainder
    return occ


def is_spherical(l: int, k: int) -> bool:
    """Unsoeld's theorem: an l-subshell with ``k`` electrons has a spherical density iff
    it is closed (``k==0`` or ``k==4l+2``), half-filled high-spin (``k==2l+1``), or of
    ``s`` type (``l==0``).  This is exactly the condition for zero xTC fill-in."""
    return l == 0 or k == 0 or k == (2 * l + 1) or k == (4 * l + 2)


def reference_density_multipoles(ref_ang: Sequence[LM], Lmax: int,
                                 tol: float = 1e-10) -> Dict[LM, float]:
    """Diagonal reference-density multipoles ``sum_o <o|Y_{L'M'}|o>`` (leg-3)."""
    refmult: Dict[LM, float] = {}
    for (lo, mo) in ref_ang:
        for Lp in range(Lmax + 1):
            for Mp in range(-Lp, Lp + 1):
                v = gA(lo, mo, Lp, Mp, lo, mo)
                if abs(v) > tol:
                    refmult[(Lp, Mp)] = refmult.get((Lp, Mp), 0.0) + v
    return {k: v for k, v in refmult.items() if abs(v) > tol}


def fast_angular(lmax: int, ref_ang: Sequence[LM], Lmax: Optional[int] = None,
                 tol: float = 1e-10) -> Dict[str, object]:
    """Exact angular support of the Coulomb vs xTC-contracted two-body operators.

    Returns block counts (Coulomb, xTC-L3), fill-in, densities, and the
    reference-density multipole set, for a basis of one radial function per l up to
    ``lmax`` with reference angular occupation ``ref_ang``.  Everything is exact
    (integer counts) because angular support does not depend on the radial factor.
    """
    if Lmax is None:
        Lmax = 2 * lmax
    O = orbs(lmax)
    ns = len(O)

    refmult = reference_density_multipoles(ref_ang, Lmax, tol)

    # per-pair Coulomb multipole supports
    ac_gA: Dict[Tuple[int, int], set] = {}
    bd_gB: Dict[Tuple[int, int], set] = {}
    for i, (la, ma) in enumerate(O):
        for j, (lc, mc) in enumerate(O):
            sA: set = set()
            sB: set = set()
            for L in range(Lmax + 1):
                for M in range(-L, L + 1):
                    if abs(gA(la, ma, L, M, lc, mc)) > tol:
                        sA.add((L, M))
                    if abs(gB(la, ma, L, M, lc, mc)) > tol:
                        sB.add((L, M))
            ac_gA[(i, j)] = sA
            bd_gB[(i, j)] = sB

    # vertex reachability: leg-2 multipoles (L,M) reachable at (a,c) given the
    # reference supplies some (L',M')
    reach: Dict[Tuple[int, int], set] = {}
    for i, (la, ma) in enumerate(O):
        for j, (lc, mc) in enumerate(O):
            s: set = set()
            for (Lp, Mp) in refmult:
                for L in range(Lmax + 1):
                    for M in range(-L, L + 1):
                        v, _ = four_Y(la, ma, L, M, Lp, Mp, lc, mc)
                        if abs(v) > tol:
                            s.add((L, M))
            reach[(i, j)] = s

    coul = 0
    l3 = 0
    fill = 0
    fill_out = 0
    total = ns ** 4
    fill_examples: List[str] = []
    for a in range(ns):
        for b in range(ns):
            for c in range(ns):
                for d in range(ns):
                    if O[a][1] + O[b][1] != O[c][1] + O[d][1]:  # m-conservation
                        continue
                    is_coul = bool(ac_gA[(a, c)] & bd_gB[(b, d)])
                    is_l3 = bool(reach[(a, c)] & ac_gA[(b, d)])
                    if is_coul:
                        coul += 1
                    if is_l3:
                        l3 += 1
                        if not is_coul:
                            fill += 1
                            if len(fill_examples) < 10:
                                fill_examples.append(str((O[a], O[b], O[c], O[d])))
                    elif is_coul:
                        fill_out += 1
    return dict(ns=ns, total=total, n_coulomb=coul, n_l3=l3,
                n_fill_in=fill, n_fill_out=fill_out,
                coul_density=coul / total, l3_density=l3 / total,
                fill_frac=fill / total,
                refmult=sorted(refmult.keys()),
                fill_examples=fill_examples)
