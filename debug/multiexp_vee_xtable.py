r"""Multi-exponent V_ee: the per-rate-pair X-table (the two-block prerequisite).

WHY THIS IS THE GATE.  `debug/multiexp_overlap_poc.py` banked the complete ONE-BODY
half of the two-block basis (S, V_ne, kinetic), validated at 2e-16 against an
independent construction.  No two-block ENERGY exists without V_ee, and V_ee is
gated on exactly one object: the X-table.  Today `pr._build_Xtab_mp` builds it at a
single rate `c = 2*alpha`, because every basis function shares alpha.

THE STRUCTURE, read off the production code rather than guessed.  The X entry is
assembled as

    I1 = Av[P1] * Bc[(l,P2)] - _corr_mp(Wf[P1], P2, l, c, B2c)
    I2 = Av[P2] * Bc[(l,P1)] - _corr_mp(Wf[P2], P1, l, c, B2c)
    X[P1][P2] = I1 + I2

where (i) `ngm._mono_moments(c,n)` and `ngm._B_table(...,c)` are both UNRESTRICTED
1D moments on [1,inf), so `Av[P1]*Bc[(l,P2)]` is the full product integral and
`_corr_mp` subtracts the wrongly-ordered region; and (ii) per `_corr_gen`'s
docstring in debug/direct_vee_corr.py, that correction is

    int int_{xi1>xi2} [sum_j w_j xi1^j] e^{-c xi1} [xi2^{p_outer}(..) d^mQ_lq] e^{-c xi2}

so `w` is the P-side (xi1) polynomial and `p_outer` the Q-side (xi2) power.  Its
body -- `wj * fac / c**(k+1)` times `b2c[(l, p_outer + j - k)]` -- is the IBP of
int_{xi2}^{inf} xi1^j e^{-c1 xi1} dxi1, whose leftover e^{-c1 xi2} is WHY b2c is
tabulated at the SUM rate (today `two_c = 2*c`).

Hence the generalization is forced, with c1 = the P-side (electron-1) pair rate and
c2 = the Q-side (electron-2) pair rate:

    Av   <- _A_moment against _mono_moments(c1)        (P side, rate c1)
    Bc   <- _B_table(m, s, l_hi, p_max, c2)            (Q side, rate c2)
    corr <- _corr_mp(Wf[P1], P2, l, c1, B2c),  B2c = _B_table(..., c1 + c2)

and I2 is the same with the electrons swapped (P-side power P2 at rate c2, Q-side
power P1 at rate c1, corr divisor c2, same B2c since c1+c2 = c2+c1).

THE SYMMETRY DOES NOT DIE, IT TRANSFERS -- and this is worth stating because an
earlier plan of mine assumed it was lost and budgeted twice the work.  At c1 = c2,
I1 and I2 are a P1<->P2 relabelling of one quantity, which is what makes today's
`mat[P1][P2] = mat[P2][P1]` triangle valid.  At c1 != c2 each table is genuinely
non-symmetric, BUT

    X^{(c1,c2)}[P1][P2] == X^{(c2,c1)}[P2][P1]

so a mirror rate pair is a free transpose.  Three rates {2a1, a1+a2, 2a2} give 9
ordered pairs but only 6 tables to build (3 diagonal, symmetric, keep the triangle;
3 off-diagonal squares, each donating its mirror).

COST, MEASURED -- and the build COUNT is not the cost driver, which is where I was
wrong.  From the 6-vs-9 count I predicted ~2x the single-rate cost.  Measured:
11.4x at (2,2) mu<=1 and 8.0x at (3,3) mu<=1 (the transpose identity does hold --
9 ordered keys from 6 builds).  The reason is inside the off-diagonal build, not in
how many there are: it needs a SECOND `_mono_moments` and a SECOND full `_B_table`
(Amono2/Bc1) and loses the P2 >= P1 triangle for a full square, so it costs ~2x a
diagonal build -- giving 3*1 + 3*2 ~ 9 diagonal-equivalents against today's 1.
Extrapolating the measured single-rate ~53 s at (4,4,2), the two-block X-table is
~7 min there rather than ~2, and at (5,5)+delta it DOMINATES rather than rounding
off the 592 s congruence.  Budget the two-block energy accordingly.

B-SEEDS AT THE SUM RATE: no new hazard, measured.  v5.13.6 is the record of B
seeds losing ~4s digits to cancellation at high (m,s), which flipped the
(5,5)+delta headline to -220; c1+c2 is a LARGER rate than anything today's code
sees, so it was the obvious place for that to resurface.  It does not: seed
accuracy improves monotonically with c at every (m,s) probed (m=2,s=3, l=2, p=6:
rel 1.48e-85 at c=2.0 -> 7.78e-88 at c=3.2 -> 6.52e-88 at c=5.2).  The
cancellation is driven by (m,s) through `_seed_guard`, not by rate, and a larger
rate strictly helps (faster decay, less cancellation).

WHAT CANNOT SERVE AS GROUND TRUTH HERE.  debug/direct_vee_corr.py validated the
ordered integral to 1e-27 -- but `_corr_gen` there takes a SINGLE `c` for both
sides, so that route covers only the degenerate case.  The genuinely two-rate entry
has no existing independent check, so this driver builds one from direct 2D mpmath
quadrature at a few hand-picked (l,m,s,P1,P2) -- expensive, hence spot checks, not
a table.

FALSIFIERS, in the order they must pass:
  F1  degenerate X-table == today's `_build_Xtab_mp`, elementwise.
  F2  degenerate V_ee == today's `vee_mp`, elementwise.
  F3  a genuinely two-rate entry (c1 != c2) == direct 2D quadrature.
  F4  degenerate END-TO-END energy == `recondition_energy`'s direct engine.

F4 is not redundant with F1/F2 and the reason is on the record: in the v5.14.0
sprint FOUR matrix-level falsifiers passed while a fixture solved the wrong
Hamiltonian.  A matrix check proves the matrix; only the energy proves the wiring.

ON TOLERANCES, so nobody "fixes" a passing check into a false one.  F1/F2/F4 are
NOT bit-for-bit and their absence is not a defect: fetching a moment table per rate
pair instead of once reorders the same mpf sums, and mpf addition is not
associative at finite precision.  multiexp_overlap_poc.py measured that floor on
the one-body half (~4e-42 relative at dps=40, FLAT in truncation).  Flatness in N
is the discriminator: rounding is flat, a dispatch bug grows with N.
"""
from __future__ import annotations

import time
from typing import Callable, Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np

from geovac import neumann_vee_general_m as ngm
from geovac import prolate_recondition as pr

R = pr.R_DEFAULT


# ---------------------------------------------------------------------------
# rate bookkeeping (same convention as the banked one-body half)
# ---------------------------------------------------------------------------
def pair_rates(fns: Sequence[pr.ProductFn],
               alpha_of: Callable[[int], float]) -> List[float]:
    """The DISTINCT pair rates appearing in the basis, sorted.

    Electron 1 of a (bra, ket) pair decays at alpha_of(j_bra) + alpha_of(j_ket)
    and electron 2 at alpha_of(k_bra) + alpha_of(k_ket), independently.  Mirrors
    `multiexp_overlap_poc._rate_tables`' keying exactly -- float keys are exact
    here because each rate is a sum of two basis inputs.
    """
    rates = set()
    for bi in fns:
        for bj in fns:
            rates.add(alpha_of(bi.j) + alpha_of(bj.j))
            rates.add(alpha_of(bi.k) + alpha_of(bj.k))
    return sorted(rates)


def _rate_pairs(rates: Sequence[float]) -> List[Tuple[float, float]]:
    """Unordered rate pairs to BUILD; mirrors come free by transpose."""
    return [(rates[i], rates[jj])
            for i in range(len(rates)) for jj in range(i, len(rates))]


# ---------------------------------------------------------------------------
# the per-rate-pair X-table
# ---------------------------------------------------------------------------
def build_Xtab_rated(ms_pairs: List[Tuple[int, int]], l_neumann: int, p_max: int,
                     rates: Sequence[float], l_caps: Dict[Tuple[int, int], int],
                     verbose: bool = False
                     ) -> Dict[Tuple[float, float, int, int, int], List[List[mp.mpf]]]:
    """X_l^{m,s}(P1,P2) for every ORDERED rate pair, keyed (c1, c2, l, m, s).

    c1 is the electron-1 (P-side) pair rate, c2 the electron-2 (Q-side) one.  At
    c1 == c2 this reduces term-for-term to `pr._build_Xtab_mp` at that rate.

    Off-diagonal pairs are built once and their mirror filled by transpose, using
    X^{(c1,c2)}[P1][P2] == X^{(c2,c1)}[P2][P1].
    """
    Xtab: Dict[Tuple[float, float, int, int, int], List[List[mp.mpf]]] = {}
    for (r1, r2) in _rate_pairs(list(rates)):
        c1, c2 = mp.mpf(r1), mp.mpf(r2)
        c_sum = c1 + c2
        sym = (r1 == r2)
        for (m, s) in ms_pairs:
            l_hi = min(l_neumann, l_caps[(m, s)])
            if l_hi < m:
                continue
            deg_extra = p_max + 2 * s + (l_hi - m)
            p_corr_max = p_max + deg_extra
            n_mono = p_max + 2 * s + (l_hi - m) + 2
            # P side carries c1, Q side c2; the IBP tail lands at the SUM rate.
            Amono1 = ngm._mono_moments(c1, n_mono)
            Bc2 = ngm._B_table(m, s, l_hi, p_max, c2)
            Bsum = ngm._B_table(m, s, l_hi, p_corr_max, c_sum)
            if not sym:
                Amono2 = ngm._mono_moments(c2, n_mono)
                Bc1 = ngm._B_table(m, s, l_hi, p_max, c1)
            for l in range(m, l_hi + 1):
                Wf = {P: [mp.mpf(w) for w in ngm._W_poly(l, m, s, P)]
                      for P in range(p_max + 1)}
                Av1 = [ngm._A_moment(l, m, s, P, Amono1) for P in range(p_max + 1)]
                mat = [[mp.mpf(0)] * (p_max + 1) for _ in range(p_max + 1)]
                if sym:
                    # c1 == c2: I2 is I1 under P1<->P2, so the triangle is valid
                    # and this is exactly today's production expression.
                    for P1 in range(p_max + 1):
                        for P2 in range(P1, p_max + 1):
                            I1 = (Av1[P1] * Bc2[(l, P2)]
                                  - pr._corr_mp(Wf[P1], P2, l, c1, Bsum))
                            I2 = (Av1[P2] * Bc2[(l, P1)]
                                  - pr._corr_mp(Wf[P2], P1, l, c1, Bsum))
                            mat[P1][P2] = mat[P2][P1] = I1 + I2
                else:
                    Av2 = [ngm._A_moment(l, m, s, P, Amono2)
                           for P in range(p_max + 1)]
                    for P1 in range(p_max + 1):
                        for P2 in range(p_max + 1):
                            # electron 1 at P1/rate c1 outside, electron 2 at
                            # P2/rate c2 -- and the swapped region.
                            I1 = (Av1[P1] * Bc2[(l, P2)]
                                  - pr._corr_mp(Wf[P1], P2, l, c1, Bsum))
                            I2 = (Av2[P2] * Bc1[(l, P1)]
                                  - pr._corr_mp(Wf[P2], P1, l, c2, Bsum))
                            mat[P1][P2] = I1 + I2
                Xtab[(r1, r2, l, m, s)] = mat
                if not sym:
                    Xtab[(r2, r1, l, m, s)] = [list(col) for col in zip(*mat)]
        if verbose:
            print(f"    rate pair ({r1:.4g},{r2:.4g}) done", flush=True)
    return Xtab


# ---------------------------------------------------------------------------
# F1: degenerate X-table vs production
# ---------------------------------------------------------------------------
def f1_degenerate_xtab(alpha: float = 1.40, j_max: int = 2, l_max: int = 2,
                       mu_max: int = 1) -> bool:
    """Degenerate rate set => build_Xtab_rated == pr._build_Xtab_mp elementwise."""
    print("\n=== F1: degenerate X-table vs pr._build_Xtab_mp ===")
    idx = pr._product_index(j_max, l_max, mu_max)
    fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
    ms_pairs, p_max, q_max, l_neumann, l_caps = _vee_shapes(fns, 0)

    ref = pr._build_Xtab_mp(ms_pairs, l_neumann, p_max, alpha, l_caps)
    got = build_Xtab_rated(ms_pairs, l_neumann, p_max, [2.0 * alpha], l_caps)

    c = 2.0 * alpha
    worst = mp.mpf(0)
    scale = mp.mpf(0)
    nchk = 0
    for (l, m, s), mat in ref.items():
        gmat = got[(c, c, l, m, s)]
        for P1 in range(len(mat)):
            for P2 in range(len(mat)):
                d = abs(mat[P1][P2] - gmat[P1][P2])
                worst = max(worst, d)
                scale = max(scale, abs(mat[P1][P2]))
                nchk += 1
    rel = worst / scale if scale > 0 else worst
    print(f"  blocks={len(ref)} entries={nchk}  max|dX|={mp.nstr(worst, 3)}  "
          f"max|X|={mp.nstr(scale, 3)}  rel={mp.nstr(rel, 3)}")
    ok = rel < mp.mpf('1e-30')
    print(f"  {'PASS' if ok else 'FAIL'} (expect exact or rounding-floor: the "
          f"degenerate branch is today's expression term-for-term)")
    return ok


def _vee_shapes(fns, l_neumann):
    """The (ms_pairs, p_max, q_max, l_neumann, l_caps) block of vee_mp, verbatim.

    ONE DELIBERATE DIVERGENCE: at l_neumann <= 0 this substitutes
    2*l_max + 4*mu_max + 10, which is what `recondition_energy` does BEFORE
    calling `vee_mp`.  `vee_mp` itself does not -- it would keep 0.  So comparing
    against `pr.vee_mp` requires passing the same POSITIVE truncation to both;
    passing 0 silently compares different physics (see :func:`f2_degenerate_vee`).
    Verified MATCH against the production block at (2,2)/(3,3)/(3,3,2)/(4,4,2).
    """
    mus = sorted(set(b.mu for b in fns))
    ms_pairs = sorted(set((m, (a + b + m) // 2) for a in mus for b in mus
                          for m in (a + b, abs(a - b)) if (a + b + m) % 2 == 0))
    s_set = sorted(set(s for (_, s) in ms_pairs))
    p_max = 2 * max(max(b.j, b.k) for b in fns) + 2 + 2 * max(s_set)
    q_max = 2 * max(max(b.l, b.m) for b in fns) + 2
    if l_neumann <= 0:
        l_neumann = 2 * max(b.l for b in fns) + 4 * max(mus) + 10
    l_neumann = min(l_neumann, q_max + 2 * max(s_set))
    ms_pairs = [(m, s) for (m, s) in ms_pairs if m <= l_neumann]
    l_caps = {(m, s): min(l_neumann, q_max + 2 * s - m) for (m, s) in ms_pairs}
    return ms_pairs, p_max, q_max, l_neumann, l_caps


# ---------------------------------------------------------------------------
# F3: a genuinely two-rate entry vs direct 2D quadrature
# ---------------------------------------------------------------------------
def _x_entry_quad(l: int, m: int, s: int, P1: int, P2: int,
                  c1: float, c2: float) -> mp.mpf:
    """X_l^{m,s}(P1,P2) at rates (c1,c2) by direct 2D mpmath quadrature.

    X = int int over the FULL quadrant of
        xi1^P1 (xi1^2-1)^s xi2^P2 (xi2^2-1)^s e^{-c1 xi1 - c2 xi2}
        x [ d^mP_l(xi_<) d^mQ_l(xi_>) ]
    i.e. the Neumann radial kernel with the smaller argument on the first kind.
    Written from the DEFINITION -- no monomial moments, no B-table, no IBP -- so
    it is independent of everything build_Xtab_rated uses.
    """
    c1m, c2m = mp.mpf(c1), mp.mpf(c2)

    def integrand(x1, x2):
        w = (x1 ** P1 * (x1 * x1 - 1) ** s * x2 ** P2 * (x2 * x2 - 1) ** s
             * mp.e ** (-c1m * x1 - c2m * x2))
        if x1 < x2:
            ker = ngm._polyval(list(ngm._RP_poly(l, m)), x1) * ngm._RQ_mp(l, m, x2)
        else:
            ker = ngm._polyval(list(ngm._RP_poly(l, m)), x2) * ngm._RQ_mp(l, m, x1)
        return w * ker

    # split on the ordering boundary xi1 = xi2 so neither region sees the kink
    lo = mp.mpf(1)
    inner_hi = mp.inf
    reg_a = mp.quad(lambda x1: mp.quad(lambda x2: integrand(x1, x2),
                                       [lo, x1]), [lo, inner_hi])
    reg_b = mp.quad(lambda x1: mp.quad(lambda x2: integrand(x1, x2),
                                       [x1, inner_hi]), [lo, inner_hi])
    return reg_a + reg_b


def f3_two_rate_spot(c1: float = 2.80, c2: float = 2.00,
                     cases: Sequence[Tuple[int, int, int, int, int]] = (
                         (0, 0, 0, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 1, 2),
                         (1, 0, 0, 0, 1), (1, 1, 1, 0, 0), (2, 0, 0, 2, 1),
                     )) -> bool:
    """The check with no prior ground truth: c1 != c2 vs 2D quadrature.

    debug/direct_vee_corr.py's 1e-27 validation cannot serve here -- its
    `_corr_gen` takes a single `c` for both sides, so it covers only c1 == c2.
    """
    print(f"\n=== F3: two-rate entries (c1={c1}, c2={c2}) vs 2D quadrature ===")
    p_max = max(max(P1, P2) for (_, _, _, P1, P2) in cases)
    ok = True
    for (l, m, s, P1, P2) in cases:
        if l < m or s < m:
            continue
        l_caps = {(m, s): l}
        X = build_Xtab_rated([(m, s)], l, p_max, [c1, c2], l_caps)
        got = X[(c1, c2, l, m, s)][P1][P2]
        ref = _x_entry_quad(l, m, s, P1, P2, c1, c2)
        den = abs(ref) if ref != 0 else mp.mpf(1)
        rel = abs(got - ref) / den
        good = rel < mp.mpf('1e-20')
        ok = ok and good
        print(f"  l={l} m={m} s={s} P=({P1},{P2})  table={mp.nstr(got, 12)}  "
              f"quad={mp.nstr(ref, 12)}  rel={mp.nstr(rel, 3)}  "
              f"{'ok' if good else 'MISMATCH'}")
        # the mirror-transpose identity, free to check while we are here
        mir = X[(c2, c1, l, m, s)][P2][P1]
        if abs(mir - got) > mp.mpf('1e-35') * den:
            ok = False
            print(f"    MIRROR FAIL: X^(c2,c1)[P2][P1] != X^(c1,c2)[P1][P2]")
    print(f"  {'PASS' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------------------
# V_ee with per-rate-pair dispatch
#
# THE COLLAPSE THAT BREAKS.  Today's `vee_mp` keys F on the SUMS of quantum
# numbers -- F[p1,q1,p2,q2] with p1 = jr[ri] + jr[rj] -- and that is valid ONLY
# because one shared exponent makes V depend on the basis through those sums
# alone.  Per-degree exponents break it: two (bra,ket) pairs with identical sums
# but different block membership decay at different rates and need different
# entries.  So F is keyed (mui, muj, rate1, rate2), and the gather carries rate
# arrays alongside jr/lr/kr/mr.
#
# The rates come from alpha_of(degree), NOT from a field on ProductFn: line 991
# of prolate_recondition stamps every function with one alpha, and three
# production sites (ngm:482, prolate_general_m:259/370) read it as
# `basis[0].alpha` -- a single scalar off the FIRST function.  A per-function
# alpha field would therefore be silently ignored by all three.  Same convention
# as the banked one-body half (multiexp_overlap_poc._rate_tables).
# ---------------------------------------------------------------------------
def vee_multi(basis: Sequence[pr.ProductFn], alpha: float, R_: float,
              l_neumann: int, alpha_of: Callable[[int], float] = None,
              verbose: bool = False) -> np.ndarray:
    """V_ee with per-pair rate dispatch.  alpha_of=None => today's single rate.

    Mirrors `pr.vee_mp` structurally; the only changes are (a) the X-table is
    rate-pair-indexed and (b) F and the gather carry the rate pair.
    """
    if alpha_of is None:
        def alpha_of(j: int) -> float:                       # noqa: E306
            return alpha

    n = len(basis)
    mus = sorted(set(b.mu for b in basis))
    m_set = sorted(set([a + b for a in mus for b in mus]
                       + [abs(a - b) for a in mus for b in mus]))
    ms_pairs, p_max, q_max, l_neumann, l_caps = _vee_shapes(basis, l_neumann)
    s_set = sorted(set(s for (_, s) in ms_pairs))

    rates = pair_rates(basis, alpha_of)
    t0 = time.time()
    Xtab = build_Xtab_rated(ms_pairs, l_neumann, p_max, rates, l_caps, verbose)
    if verbose:
        print(f"    Xtab(rated) {len(Xtab)} blocks, {len(rates)} rates "
              f"in {time.time() - t0:.0f}s", flush=True)

    Ytab: Dict[Tuple[int, int, int, int], mp.mpf] = {}
    for l in range(l_neumann + 1):
        for m in m_set:
            if m > l:
                continue
            pmpoly = list(ngm._RP_poly(l, m))
            for s in s_set:
                yp = pr._meta2(s)
                for Qq in range(q_max + 1):
                    if l > Qq + 2 * s - m or (Qq + l - m) % 2 != 0:
                        Ytab[(l, m, s, Qq)] = mp.mpf(0)
                    else:
                        Ytab[(l, m, s, Qq)] = pr._mom_eta(
                            pr._pm(pr._shift(yp, Qq), pmpoly))

    P1n = 2 * max(b.j for b in basis) + 1
    Q1n = 2 * max(b.l for b in basis) + 1
    P2n = 2 * max(b.k for b in basis) + 1
    Q2n = 2 * max(b.m for b in basis) + 1
    jac_sh = [(1, 2, 0, 2, 0), (-1, 2, 0, 0, 2), (-1, 0, 2, 2, 0), (1, 0, 2, 0, 2)]

    # which (rate1, rate2) pairs actually occur -- the full product is wasteful
    occurring = set()
    for bi in basis:
        for bj in basis:
            occurring.add((bi.mu, bj.mu,
                           alpha_of(bi.j) + alpha_of(bj.j),
                           alpha_of(bi.k) + alpha_of(bj.k)))

    Fdict: Dict[Tuple[int, int, float, float], np.ndarray] = {}
    tF = time.time()
    for (mui, muj, r1, r2) in sorted(occurring):
        S2 = mui + muj
        F = np.empty((P1n, Q1n, P2n, Q2n), object)
        for p1 in range(P1n):
            for q1 in range(Q1n):
                for p2 in range(P2n):
                    for q2 in range(Q2n):
                        tot = mp.mpf(0)
                        for m in m_set:
                            fphi = pr._phi_cec(mui, muj, m)
                            if fphi == 0 or (S2 + m) % 2 != 0:
                                continue
                            s = (S2 + m) // 2
                            mult = mp.mpf(1) if m == 0 else mp.mpf(2)
                            for l in range(max(m, 0), l_neumann + 1):
                                key = (r1, r2, l, m, s)
                                if key not in Xtab:
                                    continue
                                npre = pr._neumann_prefactor(l, m)
                                X = Xtab[key]
                                for sgn, dP1, dQ1, dP2, dQ2 in jac_sh:
                                    Pa, Qa = p1 + dP1, q1 + dQ1
                                    Pb, Qb = p2 + dP2, q2 + dQ2
                                    if Pa > p_max or Pb > p_max:
                                        continue
                                    y1 = Ytab.get((l, m, s, Qa), mp.mpf(0))
                                    if y1 == 0:
                                        continue
                                    y2 = Ytab.get((l, m, s, Qb), mp.mpf(0))
                                    if y2 == 0:
                                        continue
                                    tot += (sgn * mult * fphi * npre
                                            * X[Pa][Pb] * y1 * y2)
                        F[p1, q1, p2, q2] = tot
        Fdict[(mui, muj, r1, r2)] = F
    if verbose:
        print(f"    F tensors ({len(Fdict)}) in {time.time() - tF:.0f}s", flush=True)

    h6 = (mp.mpf(R_) / 2) ** 6
    pref = (2 / mp.mpf(R_)) * h6
    jr = np.array([b.j for b in basis])
    lr = np.array([b.l for b in basis])
    kr = np.array([b.k for b in basis])
    mr = np.array([b.m for b in basis])
    mu_arr = np.array([b.mu for b in basis])
    aj = np.array([alpha_of(b.j) for b in basis])
    ak = np.array([alpha_of(b.k) for b in basis])

    V = np.empty((n, n), object)
    for (mui, muj, r1, r2), F in Fdict.items():
        ri = np.where(mu_arr == mui)[0]
        rj = np.where(mu_arr == muj)[0]
        if len(ri) == 0 or len(rj) == 0:
            continue
        rate1 = aj[ri][:, None] + aj[rj][None, :]
        rate2 = ak[ri][:, None] + ak[rj][None, :]
        sel = (rate1 == r1) & (rate2 == r2)
        if not sel.any():
            continue
        p1 = jr[ri][:, None] + jr[rj][None, :]
        q1 = lr[ri][:, None] + lr[rj][None, :]
        p2 = kr[ri][:, None] + kr[rj][None, :]
        q2 = mr[ri][:, None] + mr[rj][None, :]
        blk = F[p1, q1, p2, q2] * pref
        sub = V[np.ix_(ri, rj)]
        sub = np.where(sel, blk, sub)
        V[np.ix_(ri, rj)] = sub
    return V


def f2_degenerate_vee(alpha: float = 1.40, j_max: int = 2, l_max: int = 2,
                      mu_max: int = 1, l_neumann: int = 10) -> bool:
    """Degenerate alpha_of => vee_multi == pr.vee_mp elementwise.

    PASS `l_neumann` EXPLICITLY and POSITIVE.  This check first ran with 0 on both
    sides and reported rel = 0.155 (max|dV| = 2.64) -- a physics-scale FAIL that
    was entirely the harness: `pr.vee_mp` does NOT substitute a default at
    l_neumann <= 0 (it only does `min(l_neumann, q_max + 2max(s))`, so 0 stays 0),
    while `_vee_shapes` substitutes 2*l_max + 4*mu_max + 10 the way
    `recondition_energy` does.  The two sides therefore ran at DIFFERENT Neumann
    truncations.  Confirmed independently: max|vee_mp(0) - vee_mp(10)| = 2.64, the
    same 2.64.  With a matching positive truncation the agreement is bit-identical.
    """
    print("\n=== F2: degenerate V_ee vs pr.vee_mp ===")
    idx = pr._product_index(j_max, l_max, mu_max)
    fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
    ref = pr.vee_mp(fns, alpha, R, l_neumann)
    got = vee_multi(fns, alpha, R, l_neumann)
    n = len(fns)
    worst = mp.mpf(0)
    scale = mp.mpf(0)
    nz = 0
    for i in range(n):
        for jj in range(n):
            a, b = ref[i, jj], got[i, jj]
            if b is None:
                print(f"  UNFILLED entry at ({i},{jj}) -- gather missed a pair")
                return False
            worst = max(worst, abs(a - b))
            scale = max(scale, abs(a))
            nz += 1
    rel = worst / scale if scale > 0 else worst
    print(f"  N={n} entries={nz}  max|dV|={mp.nstr(worst, 3)}  "
          f"max|V|={mp.nstr(scale, 3)}  rel={mp.nstr(rel, 3)}")
    ok = rel < mp.mpf('1e-30')
    print(f"  {'PASS' if ok else 'FAIL'}")
    return ok


def f2b_multirate_gather(a1: float = 1.60, a2: float = 1.00, j_split: int = 1,
                         j_max: int = 2, l_max: int = 2, mu_max: int = 1) -> bool:
    """Every V entry is FILLED when alpha_of genuinely splits.

    F2 is structurally blind to this: in the degenerate case one rate pair covers
    the whole matrix, so the `None` check passes trivially however badly the
    gather selects.  Only a real split can expose a pair the gather never routes.
    Checks fill and rate coverage, NOT values -- the values have no independent
    reference at this size (F3 covers values, entry by entry).
    """
    print("\n=== F2b: multi-rate gather leaves no entry unfilled ===")
    try:                                    # run as `python debug/<this>.py`
        import multiexp_overlap_poc as poc
    except ModuleNotFoundError:             # run as `python -c "import debug...."`
        import debug.multiexp_overlap_poc as poc
    alpha_of = poc.alpha_of_split(j_split, a1, a2)
    idx = pr._product_index(j_max, l_max, mu_max)
    fns = [pr.ProductFn(j, l, k, m, mu, a1) for (j, l, k, m, mu) in idx]
    rates = pair_rates(fns, alpha_of)
    print(f"  N={len(fns)}  distinct pair rates={[round(r, 4) for r in rates]}")
    V = vee_multi(fns, a1, R, 0, alpha_of=alpha_of)
    n = len(fns)
    unfilled = [(i, jj) for i in range(n) for jj in range(n) if V[i, jj] is None]
    zeros = sum(1 for i in range(n) for jj in range(n)
                if V[i, jj] is not None and V[i, jj] == 0)
    print(f"  unfilled={len(unfilled)}  exact-zero={zeros} / {n * n}")
    if unfilled:
        print(f"  first few unfilled: {unfilled[:5]}")
    ok = not unfilled
    # a diffuse-block entry must NOT equal the compact-block entry it would have
    # had under one exponent -- otherwise alpha_of never reached the integrals.
    ref1 = vee_multi(fns, a1, R, 0)
    same = sum(1 for i in range(n) for jj in range(n) if V[i, jj] == ref1[i, jj])
    print(f"  entries identical to the single-rate V: {same} / {n * n} "
          f"(expect a SMALL fraction; all-identical means alpha_of was ignored)")
    if same == n * n:
        print("  FAIL: split alpha_of produced the single-rate matrix")
        ok = False
    print(f"  {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    import sys
    mp.mp.dps = 40
    which = sys.argv[1] if len(sys.argv) > 1 else "f1"
    if which == "f1":
        f1_degenerate_xtab()
    elif which == "f2":
        f2_degenerate_vee()
    elif which == "f2b":
        f2b_multirate_gather()
    elif which == "f3":
        f3_two_rate_spot()
    else:
        print(f"unknown falsifier {which!r}; expected f1, f2, f2b or f3")
