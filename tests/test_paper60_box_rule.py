"""Paper 60 App. A -- the radial-domain convergence rule, which had NO backing.

`eq:boxrule`, `R_max >= 3 max_nu max_j n_j^2 R_nu ~ 3 n_max^2`, is the rule the
paper's branch-defining criterion rests on: every quantity must name its
evaluation domain, and this equation is what makes a domain adequate.  The
2026-09-12 `/qa paper_60` completeness-critic found it had **no test, no
claim-matrix row and no registry key**, and was named nowhere in the
definition-of-done -- while a fixed 60-bohr domain had already invalidated a
published exponent once and produced two failed remediations.

MEASURED (l_max=1 families, the `||T'||_1` off-diagonal norm against an
`R_max = 5 n_max^2` reference).  The first attempt at this table was GRID-limited
and had to be redone at higher resolution -- twice, because the second reading
was also too short.  The full ladder:

    npts      c=1         c=2         c=3        step(c=3)
    12000   3.6130e-02  6.5629e-04  2.5709e-05      --
    40000   3.6166e-02  6.8727e-04  2.1107e-06    0.082x
   100000   3.6169e-02  6.8984e-04  1.5592e-07    0.074x
   250000   3.6170e-02  6.9025e-04  1.5676e-07    1.005x   <- looks flat
   400000   3.6170e-02  6.9030e-04  1.9306e-07    1.232x   <- it is NOT
   600000   3.6170e-02  6.9031e-04  2.0598e-07    1.067x

c=1 and c=2 plateau.  **c=3 does not.**  The 1.005x step is a CROSSING, not a
plateau:  this column is measured against a c=5 box at the SAME npts, and that
reference carries its own truncation -- absolutely, against a c=12 box, c=5 is
7.69e-08 where c=3 is 1.29e-07, i.e. comparable.  Near 250k they partially
cancel; past it the cancellation unwinds and the column climbs.

**The converged value is 2.15e-07, and App. A's ~1e-7 is VINDICATED.**  Both
meshes converge and they AGREE:

    uniform vs c=5, 1M / 1.5M       :  2.126e-07 / 2.14665e-07  (step 1.010)
    graded  vs c=5, 120k / 240k     :  2.134e-07 / 2.15582e-07  (step 1.010)
    graded, absolute vs a graded c=12:  2.16319e-07
    cross-check: the c=5 ABSOLUTE value is 17.2714838375 (graded 240k) against
    17.2714838504 (uniform 1.5M) -- the two meshes agree to 7.5e-10.

The graded mesh `r = R t^2` resolves the near-origin region where the integrand
lives, so it converges at ~6x fewer points; the uniform mesh needs ~1.5M.

**Four mechanisms have now been offered for this one number and the first three
were wrong** -- "converged" (3 points), "below the quadrature floor" (3 points,
opposite direction), "PLATEAUS at 1.5676e-07" (4 points), and "the meshes
disagree by 1.7x so the value is undetermined" (they agree to 0.4%).  The last
came from referencing c=3 against a c=12 box at the SAME point count, i.e. with
4x coarser spacing than the thing it was resolving; that ratio is not even
monotone in npts (2.86e-07 / 1.29e-07 / 2.02e-07 at 250k / 600k / 1.5M).  The
paper's claim was right every time; only the evidence kept failing.

    npts      c=1        c=2        c=3        fixed 60      (n_max = 8)
    12000   3.613e-02  6.563e-04  2.571e-05  4.323e-02
    40000   3.617e-02  6.873e-04  2.111e-06  4.327e-02
   100000   3.617e-02  6.898e-04  1.559e-07  4.327e-02

Drift at 40k points, which is the load-bearing comparison:

    n_max   c=2 scaling   fixed 60    ratio
      6     1.707e-03    8.354e-03      4.9
      8     6.873e-04    4.327e-02     63.0
     10     3.407e-04    1.171e-01    343.8

WRONG ANSWERS THESE GUARDS EXCLUDE
  (a) "a fixed domain is fine if it is large enough" -- the paper's own point is
      that the error's SIGN OF DRIFT is what matters, not its size.  A fixed box
      degrades as the basis grows; a box that scales with the basis does not.
      Asserting only smallness would accept the fixed box at small K, which is
      exactly how the retired exponent survived.
  (b) "c = 1 suffices" -- the rule says 3.  The guards pin the ORDERING across c
      so the coefficient cannot be quietly lowered.
  (c) "the peak radius is n/a" alone -- the rule's derivation needs
      r_pk = n^2 R_nu, and the mismatched pairs (l, 1, n) reaching furthest.
"""
from __future__ import annotations

import numpy as np
import pytest

import geovac.sturmian_secular as SS
import geovac.sturmian_variational as SV

Z_HE = 2.0
NPTS = 40000        # 12k is GRID-limited for this quantity; see the
                    # resolution table in the module docstring


@pytest.fixture(scope="module", autouse=True)
def _restore_secular_grid():
    saved = dict(r=SS.r, dr=SS.dr, r2=SS.r2, R_MAX=SS.R_MAX, N_GRID=SS.N_GRID,
                 fwd=SS._ctrap_fwd, rev=SS._ctrap_rev)
    yield
    SS.r, SS.dr, SS.r2 = saved['r'], saved['dr'], saved['r2']
    SS.R_MAX, SS.N_GRID = saved['R_MAX'], saved['N_GRID']
    SS._ctrap_fwd, SS._ctrap_rev = saved['fwd'], saved['rev']
    SS._GAUNT_CACHE.clear()
    SS.reset_caches()


def _tprime_off(nmax: int, box: float) -> float:
    SV.set_grid(box, NPTS)
    cfgs = SS.build_configs(SV.family(nmax, 1))
    M = SS.build_M(cfgs, Z=Z_HE)
    D = np.diag([Z_HE * c.Rnu for c in cfgs])
    Tp = M - D
    return float(np.abs(Tp).sum() - np.abs(np.diag(Tp)).sum())


def _rel_err(nmax: int, box: float) -> float:
    ref = _tprime_off(nmax, 5.0 * nmax ** 2)
    return abs(_tprime_off(nmax, box) - ref) / abs(ref)


def test_peak_radius_is_n_squared_R_nu():
    """The rule's derivation, not just its number.

    Q_nu = p_k/R_nu and a = Q_nu/n, so r_pk = n/a = n^2 R_nu.  Asserted
    symbolically on the Config objects so a change to the weighted-charge
    convention cannot leave the appendix's derivation standing.
    """
    SV.set_grid(200.0, NPTS)
    cfgs = SS.build_configs(SV.family(6, 1))
    for c in cfgs:
        for n in (c.na, c.nb):
            a = c.Q / n                      # hydrogenic exponent
            assert abs(n / a - n ** 2 * c.Rnu) < 1e-12
    # the mismatched pairs (l, 1, n) reach furthest, which is why the rule
    # takes a max over BOTH orbitals of a configuration
    worst = max(max(c.na, c.nb) ** 2 * c.Rnu for c in cfgs)
    mismatched = [c for c in cfgs if c.na != c.nb]
    assert any(max(c.na, c.nb) ** 2 * c.Rnu == worst for c in mismatched), (
        "the furthest-reaching configuration should be a mismatched pair")


@pytest.mark.slow
def test_scaling_box_converges_while_a_fixed_box_degrades():
    """(a) THE LOAD-BEARING HALF -- the sign of drift, not the size.

    A box that scales with the basis must keep its relative error flat or
    falling as the basis grows; the fixed 60-bohr box must get WORSE.  A guard
    that only checked smallness would pass the fixed box at small K, which is
    how the retired exponent survived.
    """
    ns = (6, 8, 10)
    scaling = [_rel_err(n, 2.0 * n ** 2) for n in ns]
    fixed = [_rel_err(n, 60.0) for n in ns]
    # the scaling box IMPROVES with basis size (1.7e-3 -> 3.4e-4 measured)
    assert scaling[-1] < scaling[0], (
        f"a scaling box must improve, not degrade, with basis size: {scaling}")
    # the fixed box degrades (8.4e-3 -> 1.2e-1 measured, ~14x)
    assert fixed[-1] > fixed[0] * 5.0, (
        f"the fixed 60-bohr box must degrade with basis size: {fixed}")
    # and the gap between them WIDENS -- 4.9x at n=6 to 344x at n=10.  This is
    # the assertion that a merely-small fixed box cannot satisfy.
    r0, r1 = fixed[0] / scaling[0], fixed[-1] / scaling[-1]
    assert r1 > 10.0 * r0, (
        f"the fixed/scaling gap must widen with basis size: {r0:.1f} -> {r1:.1f}")


@pytest.mark.slow
def test_the_coefficient_ordering_pins_c_equals_three():
    """(b) c = 1 does not suffice, asserted RESOLUTION-ROBUSTLY.

    An earlier version asserted `e2 > 100*e3` and `e3 < 1e-5`.  Both are FALSE
    at 12k points (25.5 and 2.6e-5) and pass only above the hard-coded NPTS, so
    they pinned a resolution choice rather than the rule -- planting
    NPTS 40000 -> 12000, a change with no physics in it, made them fire.
    Every leg below holds at 12k, 40k AND 100k.
    """
    nmax = 8
    e1, e2, e3 = (_rel_err(nmax, c * nmax ** 2) for c in (1.0, 2.0, 3.0))
    assert e1 > 20.0 * e2, f"c=2 must be far better than c=1: {e1:.2e}, {e2:.2e}"
    assert e2 > 5.0 * e3, f"c=3 must be better than c=2: {e2:.2e}, {e3:.2e}"
    assert e3 < 1e-4, f"c=3 must be well under the c=2 level, got {e3:.2e}"


@pytest.mark.slow
def test_c3_converged_value_on_the_graded_mesh():
    """The converged c=3 box truncation is 2.15e-07, and the meshes agree.

    WRONG ANSWER REJECTED: "the converged value cannot be pinned, because the
    uniform and graded meshes disagree by ~1.7x (1.29e-07 vs 2.16e-07)."  This
    repo asserted that for part of 2026-09-13 and it is false -- they agree to
    0.43%.  The 1.29e-07 came from referencing c=3 against a c=12 box at the
    SAME point count, i.e. with 4x coarser spacing than the quantity it was
    resolving; that ratio is not even monotone in npts (2.86e-07 / 1.29e-07 /
    2.02e-07 at 250k / 600k / 1.5M).

    The graded mesh r = R t^2 resolves the near-origin region where the
    integrand lives, so it converges at 240k points where the uniform mesh
    needs ~1.5M.  Measured on the graded mesh: 2.045e-07 / 2.134e-07 /
    2.15582e-07 at 60k / 120k / 240k, step 1.010.  The uniform ladder reaches
    2.14665e-07 at 1.5M.  Cross-check: the c=5 ABSOLUTE value reads
    17.2714838375 graded against 17.2714838504 uniform, agreeing to 7.5e-10.

    Asserts BOTH convergence (the last step must be small -- an unconverged
    quantity cannot hold still) and the value, because both are established.
    The 3% tolerance admits the measured mesh-to-mesh spread with ~7x margin.
    """
    nmax = 8
    saved = (SS.r, SS.dr, SS.r2, SS.R_MAX, SS.N_GRID,
             SS._ctrap_fwd, SS._ctrap_rev)
    try:
        def rel(npts: int) -> float:
            vals = []
            for c in (3.0, 5.0):
                SV.set_grid(c * nmax ** 2, npts, kind="grade")
                cfgs = SS.build_configs(SV.family(nmax, 1))
                M = SS.build_M(cfgs, Z=Z_HE)
                D = np.diag([Z_HE * x.Rnu for x in cfgs])
                Tp = M - D
                vals.append(float(np.abs(Tp).sum()
                                  - np.abs(np.diag(Tp)).sum()))
            return abs(vals[0] - vals[1]) / abs(vals[1])

        mid, top = rel(120_000), rel(240_000)
    finally:
        (SS.r, SS.dr, SS.r2, SS.R_MAX, SS.N_GRID,
         SS._ctrap_fwd, SS._ctrap_rev) = saved
        SS._GAUNT_CACHE.clear()
        SS.reset_caches()

    # two-sided: a ONE-sided `< 1.05` would accept a collapse, which is the
    # artifact reading this whole file exists to exclude.
    assert 1.05 > top / mid > 0.95, (
        f"the graded mesh must have CONVERGED by 240k -- neither still climbing "
        f"nor collapsing: {mid:.5e} -> {top:.5e} is a step of {top / mid:.3f}x")
    assert abs(top - 2.155e-07) / 2.155e-07 < 0.03, (
        f"converged c=3 truncation should be 2.15e-07 (uniform agrees at "
        f"2.14665e-07); got {top:.5e}")


@pytest.mark.slow
def test_c3_is_a_genuine_truncation_of_order_1e_minus_7():
    """The c=3 box error is real and of order 1e-7 -- not a vanishing artifact.

    WRONG ANSWER REJECTED (1): "c=3 sits under the quadrature floor, so it
    falls without limit with resolution and the appendix's ~1e-7 is only an
    upper bound."  Measured, it stops falling and then RISES:  1.559e-7 at
    100k, 1.568e-7 at 250k, 1.931e-7 at 400k, 2.060e-7 at 600k.  A grid
    artifact does not rise.

    WRONG ANSWER REJECTED (2): "it plateaus at 1.5676e-07, so that is the
    converged value."  That is what THIS FILE asserted on 2026-09-13 and it is
    false -- the 250k->400k step moves 23%.  The flat-looking step at
    100k->250k is a crossing: the quantity is relative to a c=5 box whose own
    truncation (7.7e-8 absolute) is comparable to c=3's, so the two partially
    cancel there.  The previous guard pinned exactly that pair -- the only
    adjacent pair in the ladder under a 10% window -- and a reviewer fired it
    by planting MORE resolution, a change with no physics in it.

    The converged value IS known.  Two related quantities, both ~2.16e-07:
    relative to a c=5 box (what App. A states and this file measures) it is
    2.15582e-07 on the graded mesh at 240k and 2.14665e-07 on the uniform mesh
    at 1.5M; absolute, against a well-resolved c=12 reference, it is 2.16319e-07
    graded against 2.16145e-07 uniform -- the two rules agreeing to 0.080%.
    The pair differ by the c=5 box's own residual truncation.  An earlier version of this docstring said they disagreed by
    1.7x; that came from referencing c=3 against a c=12 box at the same point
    count, hence 4x coarser spacing than the quantity it was resolving.

    This test nonetheless asserts only the ORDER and the NON-COLLAPSE, because
    it runs at 100k and 600k where the column is still genuinely climbing
    (1.559e-07 -> 2.060e-07, step 1.067x).  Pinning 2.15e-07 at these
    resolutions would be pinning a number the test cannot reach.  See
    `test_c3_converged_value_on_the_graded_mesh` for the leg that does reach it.
    """
    nmax, box = 8, 3.0 * 8 ** 2
    global NPTS
    saved = NPTS
    try:
        vals = []
        for npts in (100000, 600000):
            NPTS = npts
            vals.append(_rel_err(nmax, box))
    finally:
        NPTS = saved
    lo, hi = vals
    for v, n in zip(vals, (100000, 600000)):
        assert 1e-8 < v < 1e-6, (
            f"c=3 must be a genuine ~1e-7 truncation at npts={n}, got {v:.4e} "
            f"-- below 1e-8 would mean a vanishing grid artifact, above 1e-6 "
            f"would contradict App. A")
    assert hi > 0.9 * lo, (
        f"c=3 must NOT collapse with resolution (that is the artifact reading "
        f"this rejects): {lo:.4e} at 100k -> {hi:.4e} at 600k")
