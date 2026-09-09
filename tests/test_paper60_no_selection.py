"""Paper 60 ``eq:no_selection`` -- selection cannot rescue the locked posing.

The claim, in the paper's words: because ``T'`` is a matrix of PURE NUMBERS its
entries do not depend on which OTHER configurations are present, so the secular
matrix of any sub-family ``A`` is EXACTLY the corresponding principal submatrix
of the full ``M``.  Cauchy interlacing then gives ``lam_max(M_A) <= lam_max(M)``,
and since ``E = -p_kappa^2/2`` with ``p_kappa = lam_max > 0``,

    E(A) >= E(M)   for every sub-family A,

i.e. no selection of configurations, however optimized, beats the family it is
drawn from.

The claim has TWO halves and they are NOT equally interesting, so they are
backed by separate tests:

  L1  (the physics)  a sub-family's INDEPENDENTLY REBUILT secular matrix is the
      principal submatrix of the pool's, bit for bit.  This is the leg that
      rests on ``T'`` being basis-independent; it is the leg that fails if the
      builder ever acquires a dependence on the configuration set.
      -> ``test_l1_rebuilt_subfamily_is_exactly_the_principal_submatrix``

  L2a (linear algebra)  Cauchy interlacing on the principal submatrices of one
      fixed ``M``.  This is a theorem about symmetric matrices and cannot fail
      for physics reasons; the test earns its place by pinning HOW TIGHT the
      inequality is (the closest sub-family sits 2.9e-08 below the pool in
      ``lam``, not miles below) and by excluding the wrong ROOT.
      -> ``test_l2a_interlacing_on_the_pool_is_tight_and_never_violated``

  L2b (the composition)  the end-to-end statement: sub-families rebuilt from
      scratch do not beat the pool.  L1 + L2a, measured together, so that a
      basis-dependent ``T'`` cannot hide behind the slicing.
      -> ``test_l2b_rebuilt_subfamilies_never_beat_the_pool``

  L3  (the measured corollary)  the sharpest available selection -- the best
      ``n`` configurations ranked by ground-state weight -- lands ABOVE the pool
      it was drawn from, as the theorem requires.
      -> ``test_l3_best_by_weight_selection_lands_above_the_pool`` (fast, K=55)
      -> ``test_l3_best102_of_the_K244_pool_lands_above_it`` (slow, the paper's
         published K=244 point)

Route: this file uses route A ONLY (``sturmian_secular.build_configs`` ->
``build_M``).  The variational route B of the sibling
``tests/test_paper60_scale_lock.py`` plays no part in ``eq:no_selection``, which
is why this is a separate file: sharing that file's ``Case`` cache would force
an expensive route-B assembly that this claim never uses.

Grid: the Paper-60 box rule ``set_grid(max(80, 5 n_max^2), 24000, "grade", 2.0)``
and ``Z = 2``, as in the sibling.  ``set_grid`` mutates ``sturmian_secular``
module globals; the module-scoped autouse fixture below restores them.
"""
import os
import sys

import numpy as np
import pytest

# Ensure project root on path (mirrors tests/conftest.py).
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import geovac.sturmian_secular as SS                                  # noqa: E402
import geovac.sturmian_variational as SV                              # noqa: E402

Z = 2.0
NPTS = 24000
EXACT_HE = -2.903724377          # exact non-relativistic He ground state (Ha)


# --------------------------------------------------------------------------
# grid hygiene: set_grid patches module globals of geovac.sturmian_secular
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# one pool per (n_max, l_max), with its sub-family rebuilds done under the SAME
# grid installation (build_configs reads the CURRENT module grid, so a pool and
# its rebuilds must never straddle a set_grid call).
# --------------------------------------------------------------------------
def _subset_indices(K: int) -> list:
    """The sub-families rebuilt for L1/L2b.  Deliberately NOT prefixes.

    A prefix-only sample would be satisfied by a builder whose entries depend on
    a configuration's POSITION in the list rather than on the configuration --
    the principal-submatrix property is a statement about arbitrary index sets.
    So: leave-one-out of the dominant 1s^2 config (index 0), a stride, a tail
    block that excludes 1s^2 entirely, and two seeded random sorted subsets.

    ``arange(0, K-1)`` -- drop the LAST configuration -- is the load-bearing
    entry and is not decorative.  It is the TIGHTEST sub-family available (the
    last config is the highest-n, lowest-weight one: measured margin -2.9e-08 in
    lam at n_max = 10, against -3.3e-01 for dropping 1s^2), so it is the only
    member of this list that can detect a small basis-dependence in T'.  An
    earlier version of this list omitted it and the L2b guard DID NOT FIRE
    against a T' scaled by ``1 + 1e-6 K``.
    """
    rng = np.random.default_rng(20260908)
    out = [np.arange(1, K),                       # K-1, drops the 1s^2 config
           np.arange(0, K - 1),                   # K-1, drops the weakest config
           np.arange(0, K, 2),                    # stride 2
           np.arange(K - min(20, K), K)]          # tail block
    for size in (max(3, K // 2), 5):
        out.append(np.sort(rng.choice(K, size=size, replace=False)))
    return out


class Pool:
    """A configuration family, its secular matrix, and its rebuilt sub-families."""

    def __init__(self, nmax: int, lmax: int) -> None:
        SV.set_grid(max(80.0, 5.0 * nmax * nmax), NPTS, "grade", 2.0)
        self.nmax, self.lmax = nmax, lmax
        self.cts = SV.family(nmax, lmax)
        self.K = len(self.cts)
        self.M = SS.build_M(SS.build_configs(self.cts), Z=Z)
        self.lam = float(np.linalg.eigvalsh(self.M)[-1])
        self.E = -self.lam ** 2 / 2
        self.subsets = _subset_indices(self.K)
        self._rebuilt = None

    @property
    def rebuilt(self) -> list:
        """The sub-families of :func:`_subset_indices`, each assembled from its
        OWN configuration tuples (never sliced out of ``self.M``).

        Built lazily and cached: the rebuilds cost about as much again as the
        pool itself, and the L3 tests -- including the slow K = 244 one -- never
        touch them.  ``pool()`` has re-installed this pool's grid by the time a
        test reads this, so the rebuilds see the same mesh the pool was built on.
        """
        if self._rebuilt is None:
            self._rebuilt = [self.rebuild(idx) for idx in self.subsets]
        return self._rebuilt

    def rebuild(self, idx) -> np.ndarray:
        """The secular matrix of the sub-family ``idx``, assembled from scratch."""
        return SS.build_M(SS.build_configs([self.cts[i] for i in idx]), Z=Z)


_POOLS: dict = {}


def pool(nmax: int, lmax: int) -> Pool:
    key = (nmax, lmax)
    if key not in _POOLS:
        _POOLS[key] = Pool(nmax, lmax)
    else:
        # re-install this pool's grid; a later Pool() may have moved it.
        SV.set_grid(max(80.0, 5.0 * nmax * nmax), NPTS, "grade", 2.0)
    return _POOLS[key]


def lam_max(A: np.ndarray) -> float:
    """The isoenergetic root ``p_kappa``: the LARGEST eigenvalue.

    ``E = -p_kappa^2/2`` is built on the largest root (deepest binding); taking
    any other root inverts the direction of the interlacing statement, which is
    what the fire test on this helper exercises.
    """
    return float(np.linalg.eigvalsh(A)[-1])


def gap_mha(E: float) -> float:
    """Energy above the exact He ground state, in mHa."""
    return (E - EXACT_HE) * 1000


# ==========================================================================
# L1 -- a sub-family's rebuilt matrix IS the principal submatrix (bit for bit)
# ==========================================================================
@pytest.mark.parametrize("nmax,lmax", [(10, 0), (5, 1)])
def test_l1_rebuilt_subfamily_is_exactly_the_principal_submatrix(nmax, lmax):
    """``build_M(subfamily) == M_pool[ix_(idx, idx)]`` EXACTLY, max|diff| = 0.0.

    This is the leg of eq:no_selection that depends on the physics: T' is a
    matrix of pure numbers, so its entries carry no dependence on which other
    configurations are in the basis.  Equality is asserted BIT-EXACTLY (not to a
    tolerance) because that is what is true -- a loose tolerance here would
    accept exactly the defect the claim excludes, a T' that drifts with the
    basis.

    Wrong answers this excludes, each asserted separately below:
      (a) a builder that ignores its argument, or depends only on the config
          COUNT -- two same-size sub-families are asserted to differ by > 0.1;
      (b) a builder whose entries depend on a configuration's POSITION in the
          list -- the index sets are scattered, not prefixes, so a positional
          term of any size breaks the equality;
      (c) "everything matches because the matrix is featureless" -- the
          off-diagonal (pure-number T') content is asserted to be substantial,
          and the rebuild is asserted to DISAGREE with the principal submatrix
          at a shifted index set.
    """
    p = pool(nmax, lmax)
    K = p.K
    assert K == {(10, 0): 55, (5, 1): 25}[(nmax, lmax)]

    # --- the claim
    for idx, Msub in zip(p.subsets, p.rebuilt):
        ref = p.M[np.ix_(idx, idx)]
        d = float(np.abs(Msub - ref).max())
        assert d == 0.0, (
            f"sub-family of size {len(idx)} is NOT the principal submatrix: "
            f"max|diff| = {d:.3e}.  T' has acquired a dependence on the "
            "configuration set, and eq:no_selection's interlacing step no "
            "longer applies")

    # --- (c) the equality is not "0 == 0": there is real off-diagonal content,
    #     and it survives into the sub-families.
    off = np.abs(p.M - np.diag(np.diag(p.M)))
    assert off.max() > 1e-2, f"pool off-diagonal only {off.max():.3e}"
    for idx, Msub in zip(p.subsets, p.rebuilt):
        if len(idx) < 3:
            continue
        off_s = np.abs(Msub - np.diag(np.diag(Msub))).max()
        assert off_s > 1e-3, (
            f"sub-family of size {len(idx)} has no off-diagonal content "
            f"({off_s:.3e}); the bit-exact equality above would then be a "
            "statement about a diagonal matrix")

    # --- (c) the equality picks out the RIGHT index set.  Shift the stride
    #     subset by one and the rebuild must NOT match.
    stride = p.subsets[2]
    assert stride[0] == 0 and stride[1] == 2          # the stride-2 sub-family
    shifted = np.clip(stride + 1, 0, K - 1)
    assert not np.array_equal(shifted, stride)
    d_wrong = float(np.abs(p.rebuilt[2] - p.M[np.ix_(shifted, shifted)]).max())
    assert d_wrong > 1e-2, (
        f"the principal submatrix at a SHIFTED index set differs by only "
        f"{d_wrong:.3e}; the bit-exact match above would then say nothing "
        "about which configurations were selected")

    # --- (a) build_M does not ignore its argument / depend only on the count.
    n = K // 3
    head = p.rebuild(np.arange(n))
    tail = p.rebuild(np.arange(K - n, K))
    assert head.shape == tail.shape
    assert float(np.abs(head - tail).max()) > 0.1, (
        "two DIFFERENT sub-families of the same size assemble to the same "
        "matrix -- build_M is a function of the configuration count alone and "
        "the principal-submatrix assertion above is vacuous")

    # --- ordering: an UNSORTED sub-family gives the permuted principal
    #     submatrix.  Not bit-exact, and the reason is understood: build_Tprime
    #     fills i <= j and mirrors, so a permutation flips the argument order of
    #     repulsion_terms and with it the floating-point summation order.
    #     Measured 6.8e-17 at n_max = 10 against a matrix of scale 2.39, i.e.
    #     13 orders of magnitude below the entries being compared.
    perm = np.array(sorted(range(K), key=lambda i: (i % 7, i)))[:max(6, K // 3)]
    Mperm = p.rebuild(perm)
    d_perm = float(np.abs(Mperm - p.M[np.ix_(perm, perm)]).max())
    assert d_perm < 1e-13, f"permuted sub-family differs by {d_perm:.3e}"
    assert np.abs(p.M).max() > 1e12 * max(d_perm, 1e-18)


# ==========================================================================
# L2a -- Cauchy interlacing on the pool: never violated, and TIGHT
# ==========================================================================
def test_l2a_interlacing_on_the_pool_is_tight_and_never_violated():
    """``lam_max(M_A) < lam_max(M)`` for every proper sub-family, and the best
    sub-family is only 2.9e-08 below -- not miles below.

    This is the linear-algebra half of eq:no_selection (Cauchy interlacing on
    principal submatrices of ONE fixed matrix), so it cannot fail for physics
    reasons; L2b below is the leg that can.  What this test contributes is the
    two things the inequality alone does not say:

      * the margin is STRICTLY negative.  ``<= 0`` would be passed by a running
        maximum initialised at 0.0 -- the exact bug an earlier version of this
        check carried, which could not tell "every sub-family is strictly below"
        from "some sub-family ties".  The band below is two-sided.
      * the margin is TIGHT (> -1e-6 in lam).  A sweep that only sampled small
        sub-families would satisfy the inequality by a mile and assert nothing;
        the leave-one-out family is included precisely because it is the
        tightest possible sub-family.

    Measured at n_max = 10, l_max = 0 (K = 55, lam_pool = 2.397695680):
      leave-one-out    max margin -2.90e-08   min margin -3.26e-01
      2000 random      max margin -9.83e-08   violations 0
    """
    p = pool(10, 0)
    K = p.K
    assert K == 55
    assert abs(p.lam - 2.397695680) < 1e-8

    # --- deterministic: every leave-one-out sub-family (the tightest there is)
    all_idx = np.arange(K)
    loo = np.array([lam_max(p.M[np.ix_(np.delete(all_idx, i),
                                       np.delete(all_idx, i))]) for i in range(K)])
    loo_margin = loo - p.lam
    assert loo_margin.max() < 0.0, (
        f"a leave-one-out sub-family ties or beats the pool: margin "
        f"{loo_margin.max():.6e}")
    assert loo_margin.max() > -1e-6, (
        f"the tightest leave-one-out sub-family is {loo_margin.max():.3e} below "
        "the pool; the inequality is then slack by construction and the test "
        "would pass for a matrix with no interlacing structure at all")
    assert loo_margin.min() < -0.1, (
        f"leave-one-out spread is only {loo_margin.min():.3e}; dropping the "
        "1s^2 configuration must cost a lot, or the sub-families are all the "
        "same matrix")

    # --- 2000 random PROPER sub-families, sizes 2..K-1
    rng = np.random.default_rng(12345)
    worst = -np.inf          # NOT 0.0 -- see the docstring
    worst_size = -1
    violations = 0
    lams = []
    for _ in range(2000):
        size = int(rng.integers(2, K))       # 2 .. K-1, never the full pool
        idx = np.sort(rng.choice(K, size=size, replace=False))
        lam = lam_max(p.M[np.ix_(idx, idx)])
        lams.append(lam)
        margin = lam - p.lam
        if margin > worst:
            worst, worst_size = margin, size
        if margin >= 0.0:
            violations += 1
    lams = np.array(lams)

    assert violations == 0, f"{violations}/2000 sub-families reach the pool"
    assert worst < 0.0, (
        f"max margin over 2000 sub-families is {worst:.6e} (size {worst_size}); "
        "a non-negative maximum means either a violation or a running maximum "
        "that was initialised at zero")
    assert worst > -2e-7, (
        f"the closest of 2000 sub-families is {worst:.3e} below the pool "
        f"(size {worst_size}); measured -9.83e-08.  A sweep whose sub-families "
        "are all far below the pool cannot detect a violation")

    # the sweep is not degenerate: sub-families differ enormously among themselves
    assert lams.min() < 0.5 * p.lam, (
        f"every sampled sub-family has lam within 2x of the pool "
        f"(min {lams.min():.4f} vs {p.lam:.4f}); the sweep is not exploring")

    # ... and in ENERGY, which is the direction eq:no_selection states
    E_sub = -lams ** 2 / 2
    assert (E_sub > p.E).all(), "a sub-family binds deeper than the pool"
    assert (E_sub - p.E).min() > 0.0
    assert (E_sub - p.E).max() > 1.0, "no sub-family is appreciably worse"


# ==========================================================================
# L2b -- the composition: REBUILT sub-families never beat the pool
# ==========================================================================
@pytest.mark.parametrize("nmax,lmax,tight", [(10, 0, 1e-7), (5, 1, 1e-5)])
def test_l2b_rebuilt_subfamilies_never_beat_the_pool(nmax, lmax, tight):
    """End-to-end: assemble each sub-family from its own configurations and
    check ``E(A) > E(pool)``.

    L2a slices one fixed ``M``, so it is blind to a builder whose ``T'`` drifts
    with the configuration set -- slicing a wrong matrix still interlaces.  This
    test rebuilds, so a basis-dependence that made a smaller family bind deeper
    surfaces here.

    Its SENSITIVITY is set by the tightest sub-family it contains, so that is
    asserted rather than assumed:  dropping the weakest (highest-n)
    configuration leaves a family only ``tight`` below the pool in lam
    (measured -2.9e-08 at n_max = 10, -3.5e-06 at n_max = 5, l_max = 1), which
    is what lets a T' scaled by ``1 + 1e-6 K`` -- a shift of ~4e-07 in lam --
    push it over.  Without that member the guard did not fire against exactly
    that plant.
    """
    p = pool(nmax, lmax)
    assert len(p.rebuilt) == len(p.subsets) >= 6

    for idx, Msub in zip(p.subsets, p.rebuilt):
        assert len(idx) < p.K                    # proper sub-families only
        lam = lam_max(Msub)
        E = -lam ** 2 / 2
        assert lam < p.lam, (
            f"rebuilt sub-family of size {len(idx)} has lam_max {lam:.9f} >= "
            f"pool {p.lam:.9f} -- eq:no_selection is violated")
        assert E > p.E, (
            f"rebuilt sub-family of size {len(idx)} binds deeper "
            f"({E:.7f}) than the pool ({p.E:.7f})")

    # SENSITIVITY: subsets[1] drops the weakest configuration and must sit only
    # just below the pool.  If this member ever stops being tight, the loop
    # above is satisfied by a mile and detects nothing.
    assert len(p.subsets[1]) == p.K - 1
    margin = lam_max(p.rebuilt[1]) - p.lam
    assert -tight < margin < 0.0, (
        f"the tightest rebuilt sub-family sits {margin:.3e} below the pool "
        f"(expected within {tight:.0e}); this leg's sensitivity to a "
        "basis-dependent T' is only as good as that number")

    # ... and the OTHER K-1 family, which drops 1s^2, is a long way below, so
    # the two ends of the range are both exercised.
    assert len(p.subsets[0]) == p.K - 1
    assert lam_max(p.rebuilt[0]) < p.lam - 0.1, (
        f"dropping the 1s^2 configuration costs only "
        f"{p.lam - lam_max(p.rebuilt[0]):.3e} in lam; the sub-families are not "
        "distinguishable")


# ==========================================================================
# L3 -- the measured corollary: the BEST selection still lands above the pool
# ==========================================================================
def _best_by_weight(M: np.ndarray, nsel: int) -> np.ndarray:
    """The ``nsel`` configurations with the largest ground-state weight.

    ``M``'s top eigenvector is the isoenergetic ground-state coefficient vector
    ``B``; ranking by ``|B|`` is the sharpest cheap selection rule available and
    is the one the paper's K=244 measurement uses.
    """
    w, V = np.linalg.eigh(M)
    B = V[:, -1]
    return np.sort(np.argsort(-np.abs(B))[:nsel])


def test_l3_best_by_weight_selection_lands_above_the_pool():
    """K = 55 pool, best 30 by ground-state weight: 29.2830 mHa vs the pool's
    29.2521 mHa -- ABOVE it, as eq:no_selection requires.

    The fast counterpart of the paper's K=244 / 102 point (backed slow below);
    the inequality is the same one and must hold at every pool size.

    Wrong answers this excludes:
      * a "selection" that quietly keeps every configuration -- the size is
        asserted, and the inequality is asserted STRICTLY, so an all-in
        selection (which ties) fails;
      * a selection ranked the wrong way round (weakest configurations first) --
        that also lands above the pool, so the strict inequality alone would
        accept it; the measured value is therefore pinned, and the worst-30
        selection is asserted to be far worse than the best-30.
    """
    p = pool(10, 0)
    nsel = 30
    sel = _best_by_weight(p.M, nsel)
    assert len(sel) == nsel < p.K

    E_sel = -lam_max(p.M[np.ix_(sel, sel)]) ** 2 / 2
    assert E_sel > p.E, (
        f"the best {nsel} configurations bind DEEPER ({E_sel:.7f}) than the "
        f"K={p.K} pool they were drawn from ({p.E:.7f})")

    # measured 2026-09-08 (driver debug/p60_avery_102_probe.py, same rule)
    assert abs(gap_mha(p.E) - 29.2521) < 0.01, f"pool gap {gap_mha(p.E):.4f} mHa"
    assert abs(gap_mha(E_sel) - 29.2830) < 0.01, f"best-{nsel} gap {gap_mha(E_sel):.4f} mHa"

    # the selection is a GOOD one -- it recovers all but 0.031 mHa of the pool
    # while dropping 25 of 55 configurations.  Without this the test would be
    # satisfied by any selection at all, including a deliberately bad one.
    assert 0.0 < gap_mha(E_sel) - gap_mha(p.E) < 0.1

    worst = np.sort(np.argsort(np.abs(np.linalg.eigh(p.M)[1][:, -1]))[:nsel])
    E_worst = -lam_max(p.M[np.ix_(worst, worst)]) ** 2 / 2
    assert gap_mha(E_worst) > gap_mha(E_sel) + 100.0, (
        f"the WEAKEST {nsel} configurations ({gap_mha(E_worst):.1f} mHa) score "
        f"like the strongest ({gap_mha(E_sel):.4f} mHa) -- the ranking is not "
        "ordering anything and the pinned value above is a coincidence")


@pytest.mark.slow
def test_l3_best102_of_the_K244_pool_lands_above_it():
    """The paper's published point: best 102 by ground-state weight out of the
    K = 244 pool (n_max = 12, s+p+d+f) gives 7.2528 mHa against the pool's own
    7.0574 mHa -- above it, as eq:no_selection requires.

    Slow (~110 s): the K = 244 assembly is the cost.  Marked slow rather than
    reduced because these are the two numbers the paper prints
    (``p60_best102_locked`` = 7.25, pool 7.06) and the withdrawal of the
    -2.90250 Ha / 1.224 mHa comparison rests on them.

    The 200-random-subset sweep is pinned at BOTH ends, deliberately.  The
    driver ``debug/p60_avery_102_probe.py`` computes ``best = max(best, E)`` over
    energies, which selects the WORST subset, not the best; its printed
    "best ... 903 mHa" (and the ``p60_best102_locked`` registry alias
    "best of 200 random 102-subsets") is the worst of the 200.  The true best is
    13.862 mHa.  Both are asserted here so that whichever way the prose is
    stated, the numbers behind it are pinned.
    """
    p = pool(12, 3)
    assert p.K == 244
    assert abs(gap_mha(p.E) - 7.0574) < 0.01, f"pool gap {gap_mha(p.E):.4f} mHa"

    sel = _best_by_weight(p.M, 102)
    assert len(sel) == 102 < p.K            # a "selection" that keeps everything
    E_sel = -lam_max(p.M[np.ix_(sel, sel)]) ** 2 / 2
    assert E_sel > p.E, (
        f"the best 102 of 244 bind DEEPER ({E_sel:.7f}) than the pool ({p.E:.7f})")
    assert abs(gap_mha(E_sel) - 7.2528) < 0.01, f"best-102 gap {gap_mha(E_sel):.4f} mHa"

    # the selection is good (0.195 mHa of the pool recovered from 42% of it),
    # so the strict inequality above is not being satisfied by a wide margin.
    assert 0.0 < gap_mha(E_sel) - gap_mha(p.E) < 0.5

    # random selection, both ends
    rng = np.random.default_rng(0)
    gaps = []
    for _ in range(200):
        s = np.sort(rng.choice(p.K, size=102, replace=False))
        gaps.append(gap_mha(-lam_max(p.M[np.ix_(s, s)]) ** 2 / 2))
    gaps = np.array(gaps)
    assert abs(gaps.min() - 13.862) < 0.05, f"best random 102-subset {gaps.min():.3f} mHa"
    assert abs(gaps.max() - 903.404) < 0.05, f"worst random 102-subset {gaps.max():.3f} mHa"
    assert gaps.min() > gap_mha(E_sel) > gap_mha(p.E), (
        "the ordering pool < best-by-weight < best-random must hold: the pool "
        "bounds everything, and weight ranking must beat blind sampling")

    # and the withdrawn comparison is out of reach at this pool, by the theorem
    assert gap_mha(p.E) > 1.224, (
        "the pool itself now reaches the cited -2.90250 Ha; eq:no_selection no "
        "longer excludes that figure and the withdrawal in Sec. atomic must be "
        "re-argued")
