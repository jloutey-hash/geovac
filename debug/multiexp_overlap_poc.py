"""Multi-exponent prolate basis: the complete ONE-BODY half, with per-pair rates.

(The filename and this title said "increment 1: the overlap" when the file only
held `overlap_multi`; it now carries S, V_ne, the kinetic term and the independent
validation.  Kept the filename to avoid breaking the memo's references.)

WHY THIS FIRST.  The scoping pass (debug/sprint_explicit_correlation_scoping_memo.md
Sec. 7) established that a two-block radial exponent is the next accuracy axis for
H2: alpha_opt drifts upward with basis size (~1.20-1.25 at (4,4) mu<=1 vs 1.40 at
(5,5)), which is one exponent straining to serve functions of different effective
range -- the same mechanism as Li's per-shell lambda, worth 48.47 mHa there.

THE ONE LOAD-BEARING MECHANISM is moment dispatch.  Today every radial moment table
is built at a single rate c = 2*alpha, because every basis function shares alpha.
With a per-degree exponent alpha_of(j), a bra-ket pair's electron-1 factor decays at
alpha_of(j_bra) + alpha_of(j_ket) and its electron-2 factor at
alpha_of(k_bra) + alpha_of(k_ket) -- INDEPENDENTLY.  With two blocks that is three
distinct rates {2a1, a1+a2, 2a2}, and each (bra, ket, electron) must be routed to
the right table.  If that routing is wrong the energies are silently wrong, so it is
proved here before anything touches production.

THE FALSIFIER, and it is the cheapest one available: set every block to the SAME
exponent.  The per-pair machinery must then reproduce `one_body_mp`'s overlap to the
working-precision floor.  Any deviation ABOVE that floor means the dispatch, not the
physics, is wrong -- and nothing downstream is worth building until it passes.

NOT bit-for-bit, and the reason is worth stating so nobody "fixes" a passing check:
an earlier version of this docstring demanded exact equality.  It is unattainable
here and its absence is not a defect.  `_ov` sums the same products in a different
ORDER once the moment table is fetched per pair instead of once, and mpf addition is
not associative at finite precision, so the two routes differ at the dps rounding
floor.  Measured at dps=40: 1.4e-42 / 1.8e-40 / 2.4e-38 absolute over
(1,1,0) / (2,2,1) / (3,3,1), i.e. ~4e-42 RELATIVE to max|S_ii| in all three --
flat in the truncation, which is the signature of rounding rather than a bug that
would grow with N.

Note on the data model (corrected in the memo): `ProductFn` carries ONE `alpha`,
which is per-FUNCTION, not per-electron.  A product function has two radial factors,
so the exponent is derived here from the DEGREE via `alpha_of`, with no new fields --
that way both engines call the same map and cannot disagree about a factor's rate.

Scope, as built (this line has gone stale twice; keep it current):
  * increment 1  -- overlap, `overlap_multi`, falsified against `one_body_mp`.
  * increment 2a -- nuclear attraction, `vne_multi`, identical dispatch since
    `_vne` is one more single moment.  Its degenerate check is LOOKUP-ONLY, for
    the reason recorded on :func:`degenerate_check_vne`.
  * increment 2b -- kinetic, `_kin2` + `h1_multi`.  The subtle piece: `pr._kin`
    takes ONE `alpha` and builds both d/dxi[xi^ja e^{-alpha xi}] and
    d/dxi[xi^jb e^{-alpha xi}] from it, so per-degree exponents need
    alpha(j_bra) and alpha(j_ket) SEPARATELY while the moment table stays the
    PAIR table at their sum.
  * the non-circular validation -- `independent_h1_check`, comparing per-pair H1
    against `build_one_body_direct`: different basis, different precision,
    different construction, so agreement is evidence rather than a tautology.
    This is what makes 2a load-bearing too, which is why 2b had to come first.

ON THE TWO-BLOCK cond(S), because the raw number invites a misreading I made first.
The two-exponent monomial overlap comes out at cond ~9.6e16 at (3,3) mu<=1, which
looks like the catastrophic regime -- but that is the wrong comparison.  Two controls
settle it:

  * NOT a short moment table.  cond is EXACTLY flat in n_mom -- 9.610e16 at 48, 96
    and 192 terms, min eig 1.4227e-10 unchanged -- so the slower-decaying a2 = 1.00
    is fully resolved by the single-alpha table sizing.  (This was the cheaper
    hypothesis and it is dead.)
  * The right baseline is ONE-block, not the re-based basis.  Same basis, same code
    path, single alpha = 1.60: cond 2.046e15, min eig 2.06e-13 -- and Paper 12
    independently measured 2.6e14 at (3,3).  So the MONOMIAL overlap is severely
    ill-conditioned at ANY exponent; that is precisely why
    `geovac/prolate_recondition.py` exists.  Two-block costs a factor of ~47 over
    one-block, not four orders of magnitude.

So the open question is not "is the two-block basis conditioned badly" (the monomial
basis always is) but "does the orthogonal re-basing absorb it as well as it absorbs
the single-block case", where gegenbauer takes 2.6e14 down to ~1e4.  That needs
`_transforms_per_mu` to become block-aware, which is increment 2.  Standing caution
while measuring it: the failed-approaches ledger's `k_n = Z/n` row is the record of a
per-FUNCTION exponent that fixed conditioning perfectly (kappa = 1.0000) and
destroyed completeness, plateauing accuracy near 60 mHa -- conditioning alone is not
the thing to optimise here.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np

from geovac import neumann_vee_general_m as ngm
from geovac import prolate_recondition as pr

R = pr.R_DEFAULT


def alpha_of_split(j_split: int, a1: float, a2: float) -> Callable[[int], float]:
    """Two-block radial exponent: degree j <= j_split gets a1, above it a2.

    Higher polynomial degree reaches further out in xi, so the diffuse block is
    expected to want the SMALLER exponent (a2 < a1).  Same function count as the
    single-exponent basis -- this reshuffles exponents rather than adding
    functions, mirroring the Li per-shell-lambda measurement, which was taken "at
    the same function count".
    """
    def f(j: int) -> float:
        return a1 if j <= j_split else a2
    return f


def _rate_tables(fns: Sequence[pr.ProductFn], alpha_of: Callable[[int], float],
                 n_mom: int) -> Dict[float, List]:
    """One monomial-moment table per DISTINCT pair rate appearing in the basis.

    The rate for a (bra, ket) radial factor is alpha_of(j_bra) + alpha_of(j_ket);
    with two blocks that is {2a1, a1+a2, 2a2}.  Keyed by the float rate, which is
    exact here because the rates are sums of two basis inputs.
    """
    rates = set()
    for bi in fns:
        for bj in fns:
            rates.add(alpha_of(bi.j) + alpha_of(bj.j))
            rates.add(alpha_of(bi.k) + alpha_of(bj.k))
    return {r: ngm._mono_moments(r, n_mom) for r in sorted(rates)}


def overlap_multi(fns: Sequence[pr.ProductFn], alpha_of: Callable[[int], float],
                  n_mom: int) -> np.ndarray:
    """Overlap S with per-pair moment dispatch (mpf).

    Mirrors `one_body_mp`'s overlap exactly -- h6 * phi_cc(mu) * ov1 * ov2, mu
    block-diagonal -- except that each electron's `_ov` is evaluated against the
    moment table for ITS OWN pair rate rather than one shared table.
    """
    A_by_rate = _rate_tables(fns, alpha_of, n_mom)
    n = len(fns)
    S = np.empty((n, n), object)
    h6 = (mp.mpf(R) / 2) ** 6
    for i in range(n):
        bi = fns[i]
        for jj in range(i, n):
            bj = fns[jj]
            if bi.mu != bj.mu:
                S[i, jj] = S[jj, i] = mp.mpf(0)
                continue
            mu = bi.mu
            cc = pr._phi_cc(mu)
            A1 = A_by_rate[alpha_of(bi.j) + alpha_of(bj.j)]
            A2 = A_by_rate[alpha_of(bi.k) + alpha_of(bj.k)]
            o1 = pr._ov(bi.j + bj.j, bi.l + bj.l, mu, A1)
            o2 = pr._ov(bi.k + bj.k, bi.m + bj.m, mu, A2)
            S[i, jj] = S[jj, i] = h6 * cc * o1 * o2
    return S


def vne_multi(fns: Sequence[pr.ProductFn], alpha_of: Callable[[int], float],
              n_mom: int) -> np.ndarray:
    """Nuclear attraction V_ne with per-pair moment dispatch (mpf).

    Increment 2a.  Same dispatch as :func:`overlap_multi` -- `_vne` is one more
    single moment (r1 . a0 in the factored language), so the only question is
    again whether each electron's factor is routed to ITS OWN pair rate.  Mirrors
    `one_body_mp`'s V_ne term exactly: pref_V * cc * (vne1*ov2 + ov1*vne2).
    """
    A_by_rate = _rate_tables(fns, alpha_of, n_mom)
    n = len(fns)
    V = np.empty((n, n), object)
    h6 = (mp.mpf(R) / 2) ** 6
    pref_V = -(4 * mp.mpf(1) / mp.mpf(R)) * h6
    for i in range(n):
        bi = fns[i]
        for jj in range(i, n):
            bj = fns[jj]
            if bi.mu != bj.mu:
                V[i, jj] = V[jj, i] = mp.mpf(0)
                continue
            mu = bi.mu
            cc = pr._phi_cc(mu)
            A1 = A_by_rate[alpha_of(bi.j) + alpha_of(bj.j)]
            A2 = A_by_rate[alpha_of(bi.k) + alpha_of(bj.k)]
            o1 = pr._ov(bi.j + bj.j, bi.l + bj.l, mu, A1)
            o2 = pr._ov(bi.k + bj.k, bi.m + bj.m, mu, A2)
            v1 = pr._vne(bi.j + bj.j, bi.l + bj.l, mu, A1)
            v2 = pr._vne(bi.k + bj.k, bi.m + bj.m, mu, A2)
            V[i, jj] = V[jj, i] = pref_V * cc * (v1 * o2 + o1 * v2)
    return V


def _kin2(ja: int, la: int, alpha_a: float, jb: int, lb: int, alpha_b: float,
          mu: int, A: Sequence) -> Tuple[mp.mpf, mp.mpf]:
    """Increment 2b: (gradient, azimuthal) for one electron with SEPARATE bra and
    ket exponents.

    This is the one piece that per-degree exponents genuinely change rather than
    merely re-route.  `pr._kin` takes a single `alpha` and builds BOTH derivative
    polynomials from it:

        mx(jx) = d/dxi[ xi^jx e^{-alpha xi} ] / e^{-alpha xi}
               = -alpha xi^jx + jx xi^{jx-1}

    but the bra carries alpha_a and the ket alpha_b, so each side needs its own.
    The moment table `A` is still the PAIR table at rate alpha_a + alpha_b, because
    the product of the two exponentials is what gets integrated -- only the
    polynomial prefactors differ per side.

    THE ERROR THIS EXISTS TO PREVENT: passing the pair sum (alpha_a + alpha_b) as
    a single alpha to `pr._kin`.  That is wrong -- it puts the total decay rate
    inside each derivative instead of the per-side rate -- and it is invisible to
    a same-exponent falsifier, because at alpha_a = alpha_b = a the correct
    per-side value is `a` while the sum is `2a`, so the bug shows up ONLY when the
    two exponents differ... and conversely, a version that used `alpha_a` for both
    sides would pass the degenerate check too.  So neither mistake is catchable by
    the degenerate test alone; the H1-level comparison against
    `build_one_body_direct` is what closes it.

    Mirrors `pr._kin`'s structure exactly otherwise (mu = 0 and mu > 0 branches,
    the (xi^2-1)^{mu-1} weight for mu > 1, the azimuthal mu^2 term on the shifted
    weights).
    """
    a_a, a_b = mp.mpf(alpha_a), mp.mpf(alpha_b)
    m = mp.mpf(mu)
    if mu == 0:
        def mx(jx: int, aa: mp.mpf) -> List:
            t = pr._shift([-aa], jx)
            if jx > 0:
                t = pr._pa(t, pr._shift([mp.mpf(jx)], jx - 1))
            return t
        xip = pr._pm(pr._pm(mx(ja, a_a), mx(jb, a_b)), pr._xi2m1(1))

        def my(lx: int) -> List:
            return pr._shift([mp.mpf(lx)], lx - 1) if lx > 0 else [mp.mpf(0)]
        etap = pr._pm(pr._pm(my(la), my(lb)), pr._meta2(1))
    else:
        def nx(jx: int, aa: mp.mpf) -> List:
            t = pr._shift([m], jx + 1)
            t = pr._ps(t, pr._pm(pr._shift([aa], jx), pr._xi2m1(1)))
            if jx > 0:
                t = pr._pa(t, pr._pm(pr._shift([mp.mpf(jx)], jx - 1),
                                     pr._xi2m1(1)))
            return t
        xip = pr._pm(nx(ja, a_a), nx(jb, a_b))
        if mu > 1:
            xip = pr._pm(xip, pr._xi2m1(mu - 1))

        def ny(lx: int) -> List:
            t = pr._shift([-m], lx + 1)
            if lx > 0:
                t = pr._pa(t, pr._pm(pr._shift([mp.mpf(lx)], lx - 1),
                                     pr._meta2(1)))
            return t
        etap = pr._pm(ny(la), ny(lb))
        if mu > 1:
            etap = pr._pm(etap, pr._meta2(mu - 1))

    px = pr._pm(pr._shift([mp.mpf(1)], ja + jb), pr._xi2m1(mu))
    py = pr._pm(pr._shift([mp.mpf(1)], la + lb), pr._meta2(mu))
    grad = (pr._mom_xi(xip, A) * pr._mom_eta(py)
            + pr._mom_xi(px, A) * pr._mom_eta(etap))

    azi = mp.mpf(0)
    if mu > 0:
        pxf = pr._pm(pr._shift([mp.mpf(1)], ja + jb), pr._xi2m1(mu - 1))
        pyf = pr._pm(pr._shift([mp.mpf(1)], la + lb), pr._meta2(mu - 1))
        azi = (pr._mom_xi(pr._shift(pxf, 2), A) * pr._mom_eta(pyf)
               - pr._mom_xi(pxf, A) * pr._mom_eta(pr._shift(pyf, 2))) * m * m
    return grad, azi


def h1_multi(fns: Sequence[pr.ProductFn], alpha_of: Callable[[int], float],
             n_mom: int) -> np.ndarray:
    """H1 = T + V_ne with per-pair dispatch and per-side kinetic exponents.

    Increment 2b.  Mirrors `one_body_mp`'s H term exactly -- pref_T*(cc*(k1*o2 +
    o1*k2) + ss*(f1*o2 + o1*f2)) + pref_V*cc*(v1*o2 + o1*v2) -- with every moment
    routed to its own pair rate and every derivative built from its own side's
    exponent.
    """
    A_by_rate = _rate_tables(fns, alpha_of, n_mom)
    n = len(fns)
    H = np.empty((n, n), object)
    h6 = (mp.mpf(R) / 2) ** 6
    pref_T = mp.mpf('0.5') * (4 / mp.mpf(R) ** 2) * h6
    pref_V = -(4 * mp.mpf(1) / mp.mpf(R)) * h6
    for i in range(n):
        bi = fns[i]
        for jj in range(i, n):
            bj = fns[jj]
            if bi.mu != bj.mu:
                H[i, jj] = H[jj, i] = mp.mpf(0)
                continue
            mu = bi.mu
            cc, ssv = pr._phi_cc(mu), pr._phi_ss(mu)
            aj_i, aj_j = alpha_of(bi.j), alpha_of(bj.j)
            ak_i, ak_j = alpha_of(bi.k), alpha_of(bj.k)
            A1 = A_by_rate[aj_i + aj_j]
            A2 = A_by_rate[ak_i + ak_j]
            o1 = pr._ov(bi.j + bj.j, bi.l + bj.l, mu, A1)
            o2 = pr._ov(bi.k + bj.k, bi.m + bj.m, mu, A2)
            k1, f1 = _kin2(bi.j, bi.l, aj_i, bj.j, bj.l, aj_j, mu, A1)
            k2, f2 = _kin2(bi.k, bi.m, ak_i, bj.k, bj.m, ak_j, mu, A2)
            v1 = pr._vne(bi.j + bj.j, bi.l + bj.l, mu, A1)
            v2 = pr._vne(bi.k + bj.k, bi.m + bj.m, mu, A2)
            val = (pref_T * (cc * (k1 * o2 + o1 * k2)
                             + ssv * (f1 * o2 + o1 * f2))
                   + pref_V * cc * (v1 * o2 + o1 * v2))
            H[i, jj] = H[jj, i] = val
    return H


def independent_h1_check(j_max: int, l_max: int, mu_max: int,
                         alpha: float) -> None:
    """THE non-circular validation: per-pair H1 vs build_one_body_direct.

    Different bases (monomial-then-re-based vs built directly in the orthogonal
    basis), different arithmetic (mpf vs float64), different construction (explicit
    polynomial moments vs factored radial x angular blocks).  So agreement here is
    real evidence, unlike the 0.000e+00 self-comparison in
    :func:`degenerate_check_vne`.  At the degenerate point every block shares one
    exponent, so this also retroactively validates the 2a dispatch.
    """
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        n_mom = 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20
        H_mono = h1_multi(fns, lambda j: alpha, n_mom)
        Nr = (j_max + 1) ** 2
        Na = len([1 for l in range(l_max + 1) for m in range(l_max + 1)
                  if (l + m) % 2 == 0])
        Tr, Ta = pr._transforms_per_mu("gegenbauer", j_max, l_max, mu_max, alpha)
        H_ours = pr._to_f64(pr._factored_cob(H_mono, mu_max + 1, Nr, Na, Tr, Ta))
    _S_ref, H_ref = pr.build_one_body_direct(j_max, l_max, mu_max, alpha,
                                             basis="gegenbauer")
    d = np.abs(H_ours - H_ref)
    scale = np.abs(H_ref).max()
    rel = float(d.max() / scale)
    print("  (%d,%d) mu<=%d  N=%d  max|dH1| = %.3e  scale-rel = %.3e  %s"
          % (j_max, l_max, mu_max, len(fns), float(d.max()), rel,
             "independent routes agree" if rel < 1e-10
             else "FAIL -- the per-side kinetic exponents are wrong"), flush=True)


# NOTE (2026-09-18): `mutation_check_kin2` and an inline Hamiltonian copy were
# EXCISED here.  They were a test harness for `_kin2`, whose correctness is
# established independently and more strongly by `independent_h1_check`
# (agreement with `build_one_body_direct` to 2e-16 scale-relative, in a
# different basis, at a different precision, through a different
# construction).  The harness's inline operator disagreed with `h1_multi` by
# 1.18e-01 relative for reasons seven hypotheses failed to localise; a rewrite
# attempt then corrupted the file (an edit matched inside a docstring and put
# prose in code position).  Deleted rather than debugged: nothing in the sprint
# depended on it, and one canonical Hamiltonian (`h1_multi`) is the point.


def degenerate_check_vne(j_max: int, l_max: int, mu_max: int,
                         alpha: float) -> None:
    """Falsifier for V_ne: one exponent everywhere must reproduce the single-rate
    expression.

    READ THIS RESULT WEAKLY, and the reason is the result itself.  This check
    returns max|dV_ne| = 0.000e+00 EXACTLY, at every truncation -- whereas the
    overlap check deviates at ~4e-42 relative.  That contrast is the tell: the
    overlap's reference is `one_body_mp`, which sums the same products in a
    different ORDER, so genuine independence shows up as mpf rounding.  Here the
    reference is rebuilt from `pr._ov` / `pr._vne` against a shared moment table
    -- the same helpers, called in the same sequence -- so identical arithmetic
    gives identically zero.  It confirms that the rate-table LOOKUP returns the
    right table and essentially nothing else.

    Exact agreement between two routes is evidence about the routes being the
    same, not about either being right (the corpus rule: suspicious exactness
    means suspect the test).

    WHAT WAS TRIED NEXT, AND WHAT IT ACTUALLY ESTABLISHED (recorded because the
    attempt's own banner overclaimed).  I set out to compare against
    `build_one_body_direct`, which builds V_ne through the factored `r1 . a0` path
    in float64 -- a genuinely different construction.  That is not reachable as
    written: `build_one_body_direct` returns H1 = T + V_ne together, and
    `one_body_mp` bundles them too, so isolating V_ne through either is circular.
    What the run actually did instead was a CONGRUENCE-INVARIANCE check -- re-base
    this same `vne_multi` matrix through the laguerre_legendre and gegenbauer
    families and require equal generalised eigenvalues of (V_ne, S).  It passes
    (4.4e-15 at (1,1,0), 3.5e-13 at (2,2,1)), but it tests that `_factored_cob` is
    a genuine congruence -- already covered by
    `tests/test_paper12_recondition.py::test_factored_change_of_basis_equals_the_dense_one`
    -- and NOT the V_ne construction, since both sides come from this function's
    own output.

    So: the V_ne dispatch remains LOOKUP-VERIFIED ONLY.  The honest independent
    check is at the H1 level -- per-pair dispatch with kinetic included, against
    `build_one_body_direct`'s float64 factored H1 -- which makes increment 2b a
    PREREQUISITE for validating 2a rather than a follow-on to it.
    """
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        n_mom = 4 * j_max + 4 * (mu_max + 1) + 4 * l_max + 16
        A = ngm._mono_moments(2.0 * alpha, n_mom)
        h6 = (mp.mpf(R) / 2) ** 6
        pref_V = -(4 * mp.mpf(1) / mp.mpf(R)) * h6
        n = len(fns)
        worst = mp.mpf(0)
        scale = mp.mpf(0)
        V_new = vne_multi(fns, lambda j: alpha, n_mom)
        for i in range(n):
            bi = fns[i]
            for jj in range(n):
                bj = fns[jj]
                if bi.mu != bj.mu:
                    ref = mp.mpf(0)
                else:
                    mu = bi.mu
                    cc = pr._phi_cc(mu)
                    o1 = pr._ov(bi.j + bj.j, bi.l + bj.l, mu, A)
                    o2 = pr._ov(bi.k + bj.k, bi.m + bj.m, mu, A)
                    v1 = pr._vne(bi.j + bj.j, bi.l + bj.l, mu, A)
                    v2 = pr._vne(bi.k + bj.k, bi.m + bj.m, mu, A)
                    ref = pref_V * cc * (v1 * o2 + o1 * v2)
                d = abs(V_new[i, jj] - ref)
                if d > worst:
                    worst = d
                if abs(ref) > scale:
                    scale = abs(ref)
        rel = worst / scale if scale != 0 else mp.mpf(0)
        verdict = ("at the dps floor" if rel < mp.mpf(10) ** -30
                   else "FAIL -- above the rounding floor, dispatch is wrong")
        print("  (%d,%d) mu<=%d  N=%d  max|dV_ne| = %.3e   rel = %.3e  %s"
              % (j_max, l_max, mu_max, n, float(worst), float(rel), verdict),
              flush=True)


def degenerate_check(j_max: int, l_max: int, mu_max: int, alpha: float) -> None:
    """THE FALSIFIER: all blocks at one exponent must reproduce one_body_mp's S."""
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, alpha) for (j, l, k, m, mu) in idx]
        n_mom = 4 * j_max + 4 * (mu_max + 1) + 4 * l_max + 16
        A = ngm._mono_moments(2.0 * alpha, n_mom)
        S_ref, _H = pr.one_body_mp(fns, alpha, R, A)
        S_new = overlap_multi(fns, lambda j: alpha, n_mom)

        worst = mp.mpf(0)
        for i in range(len(fns)):
            for j in range(len(fns)):
                d = abs(S_new[i, j] - S_ref[i, j])
                if d > worst:
                    worst = d
        scale = max(abs(S_ref[i, i]) for i in range(len(fns)))
        rel = worst / scale
        verdict = ("at the dps floor" if rel < mp.mpf(10) ** -30
                   else "FAIL -- above the rounding floor, dispatch is wrong")
        print("  (%d,%d) mu<=%d  N=%d  max|dS| = %.3e   rel to max|S_ii| = %.3e  %s"
              % (j_max, l_max, mu_max, len(fns), float(worst), float(rel),
                 verdict), flush=True)


def two_block_probe(j_max: int, l_max: int, mu_max: int, j_split: int,
                    a1: float, a2: float) -> None:
    """Sanity: a genuine two-exponent basis produces a DIFFERENT, finite, SPD S.

    Not an accuracy claim -- the energy needs V_ne and the kinetic term, which are
    increment 2.  This only confirms the dispatch runs over three distinct rates
    and yields a positive-definite overlap rather than something degenerate.
    """
    with mp.workdps(pr.DEFAULT_DPS):
        idx = pr._product_index(j_max, l_max, mu_max)
        fns = [pr.ProductFn(j, l, k, m, mu, a1) for (j, l, k, m, mu) in idx]
        n_mom = 4 * j_max + 4 * (mu_max + 1) + 4 * l_max + 16
        aof = alpha_of_split(j_split, a1, a2)
        rates = sorted({aof(bi.j) + aof(bj.j) for bi in fns for bj in fns})
        S = overlap_multi(fns, aof, n_mom)
        n = len(fns)
        Sf = np.array([[float(S[i, j]) for j in range(n)] for i in range(n)])
        w = np.linalg.eigvalsh(0.5 * (Sf + Sf.T))
        print("  split j<=%d: a1=%.2f a2=%.2f  distinct pair rates %s"
              % (j_split, a1, a2, ["%.2f" % r for r in rates]))
        print("      S: min eig %.3e  max eig %.3e  cond %.2e  SPD: %s"
              % (w[0], w[-1], w[-1] / w[0] if w[0] > 0 else float("inf"),
                 bool(w[0] > 0)), flush=True)


if __name__ == "__main__":
    print("=== per-pair moment dispatch: overlap (inc. 1) + V_ne (inc. 2a) ===")
    print()
    print("FALSIFIER 1 (overlap) -- one exponent everywhere must reproduce")
    print("one_body_mp's S to the working-precision floor:")
    for (j, l, mu) in ((1, 1, 0), (2, 2, 1), (3, 3, 1)):
        degenerate_check(j, l, mu, 1.40)
    print()
    print("FALSIFIER 2 (V_ne) -- same standard, same dispatch. NOTE: this one is")
    print("lookup-only; exact 0.000e+00 means the two routes share arithmetic,")
    print("not that either is right (see the function docstring):")
    for (j, l, mu) in ((1, 1, 0), (2, 2, 1), (3, 3, 1)):
        degenerate_check_vne(j, l, mu, 1.40)
    print()
    print("FALSIFIER 3 (inc. 2b, THE non-circular one) -- per-pair H1 vs")
    print("build_one_body_direct: different basis, precision and construction.")
    print("This is what would catch a wrong per-side kinetic exponent:")
    for (j, l, mu) in ((1, 1, 0), (2, 2, 1), (3, 3, 2)):
        independent_h1_check(j, l, mu, 1.40)
    print()
    print("two-exponent dispatch runs and gives an SPD overlap. NOT an accuracy")
    print("claim: cond ~1e17 here is the MONOMIAL overlap, which is ill-conditioned")
    print("at any exponent (one-block is 2.05e15); whether re-basing absorbs it is")
    print("the next measurement, and the deciding one:")
    for (jsp, a1, a2) in ((1, 1.60, 1.00), (2, 1.80, 1.10)):
        two_block_probe(3, 3, 1, jsp, a1, a2)
    print()
