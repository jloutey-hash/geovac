"""Paper 60 sec:atomic -- the full-shell growth rule, which had NO backing.

The paper calls this "the substantive finding" of its atomic section: the
sublinearity belongs to the BASIS-GROWTH RULE, not to the isoenergetic
construction.  Grown the other way -- full hydrogenic shells, l_max = n_max - 1
-- the same construction gives a SUPERLINEAR ||T'||_1 while the total exponent
rises toward 1 without reaching it.

`/qa paper_60` FULL 2026-09-12 (code-B M2) found this family had no test, no
registry key and no C17 family, and that its lower endpoint named the wrong
window.  The prose was corrected first (remediation 8); this file is the
separate, separately-reviewed guard pass required by CLAUDE.md Sec. 9.

MEASURED before the guard was written (box max(80, 5n^2), 12000 pts, n = 3..10):

    K     |M|_1     total slope   ||T'||_1   T' slope
     35    45.387     0.849         18.560    1.063
     56    68.218     0.867         30.773    1.076
     84    97.407     0.879         47.567    1.074
    220   231.614     0.911        133.893    1.076

WRONG ANSWERS THIS FILE EXCLUDES
  (a) "the sublinearity is a property of the isoenergetic method" -- it is not;
      change only the growth rule and the T' block goes superlinear.
  (b) "the total is superlinear too under this rule" -- the retired reading.
      No measured full-shell point makes the total exceed 1, and the paper
      explicitly declines to claim one.
  (c) "the total exponent is flat" -- it rises monotonically, which is what
      makes (b) a live risk worth pinning.
"""
from __future__ import annotations

import numpy as np
import pytest

import geovac.sturmian_secular as SS
import geovac.sturmian_variational as SV

Z_HE = 2.0
#: full-shell rungs; n -> (K, l_max = n-1)
SHELLS = (4, 5, 6, 7, 8, 9, 10)


@pytest.fixture(scope="module", autouse=True)
def _restore_secular_grid():
    """`SV.set_grid` mutates module globals in `sturmian_secular`.

    Without this, the grid this file leaves behind changes the result of any
    later test in the same process -- measured 2026-09-12: omitting it made
    `test_paper60_split_is_box_sensitive_and_ordering_is_not` FAIL, since that
    test deliberately pins a box-size artifact and therefore reads the global
    it was handed.  Same pattern as `test_paper60_resource_ladder.py`.
    """
    saved = dict(r=SS.r, dr=SS.dr, r2=SS.r2, R_MAX=SS.R_MAX, N_GRID=SS.N_GRID,
                 fwd=SS._ctrap_fwd, rev=SS._ctrap_rev)
    yield
    SS.r, SS.dr, SS.r2 = saved['r'], saved['dr'], saved['r2']
    SS.R_MAX, SS.N_GRID = saved['R_MAX'], saved['N_GRID']
    SS._ctrap_fwd, SS._ctrap_rev = saved['fwd'], saved['rev']
    SS._GAUNT_CACHE.clear()
    SS.reset_caches()


def _full_shell(n: int):
    """||M||_1 and ||T'||_1 for the full hydrogenic shell family at n_max = n."""
    SV.set_grid(max(80, 5 * n ** 2), 12000)
    cfgs = SS.build_configs(SS.gen_configs(n - 1, {l: n for l in range(n)}))
    M = SS.build_M(cfgs, Z=Z_HE)
    D = np.diag([Z_HE * c.Rnu for c in cfgs])
    return len(cfgs), float(np.abs(M).sum()), float(np.abs(M - D).sum())


def _slopes(K, V):
    K = np.asarray(K, float)
    V = np.asarray(V, float)
    return [float(np.log(V[i + 1] / V[i]) / np.log(K[i + 1] / K[i]))
            for i in range(len(K) - 1)]


@pytest.fixture(scope="module")
def ladder():
    K, tot, tp = [], [], []
    for n in SHELLS:
        k, a, b = _full_shell(n)
        K.append(k)
        tot.append(a)
        tp.append(b)
    return K, tot, tp


@pytest.mark.slow
def test_full_shell_tprime_block_is_superlinear(ladder):
    """(a) Change ONLY the growth rule and the T' block goes superlinear.

    This is what makes the fixed-l_max sublinearity a property of the rule
    rather than of the method, so it must be asserted, not narrated.
    """
    K, _, tp = ladder
    fit = float(np.polyfit(np.log(K), np.log(tp), 1)[0])
    # measured 1.0727 (robust to grid 12k/40k and box c5/c8); the paper states
    # K^1.07.  Bracketed rather than one-sided, so a runaway fit fails too --
    # `fit > 1.02` alone would accept 1.5.
    assert 1.04 < fit < 1.12, (
        f"full-shell ||T'||_1 must be superlinear near K^1.07, got K^{fit:.4f}")
    for s in _slopes(K, tp):
        assert s > 1.0, f"every full-shell T' rung must be superlinear: {s:.4f}"


@pytest.mark.slow
def test_full_shell_total_rises_but_never_reaches_one(ladder):
    """(b) and (c) together -- the paper's own self-limiting claim.

    Rising is asserted because a flat total would make the 'no superlinear
    point' claim vacuous; staying below 1 is asserted because the superlinear
    reading was retired.  Neither half alone is the claim.
    """
    K, tot, _ = ladder
    sl = _slopes(K, tot)
    for a, b in zip(sl, sl[1:]):
        assert b > a, f"the total exponent must RISE across the ladder: {sl}"
    assert max(sl) < 1.0, (
        f"no measured full-shell point may make the TOTAL superlinear: {sl}")
    assert sl[-1] > 0.90, f"the top rung should be near 0.911, got {sl[-1]:.4f}"
    # the global window fit the paper quotes (0.893 over K=56..220), which was
    # asserted nowhere until 2026-09-13
    m = [i for i, k in enumerate(K) if 56 <= k <= 220]
    gfit = float(np.polyfit(np.log(np.array(K, float)[m]),
                            np.log(np.array(tot, float)[m]), 1)[0])
    # +-0.01 does not discriminate the WINDOW this leg names: K=35..220 gives
    # 0.8865, K=84..220 gives 0.8995 and K=56..165 gives 0.8885, all inside it.
    # Across GRID and BOX variants the same fit moves by <1e-4, so +-0.004 admits
    # the real spread with ~40x margin and excludes every wrong window measured.
    assert abs(gfit - 0.893) < 0.004, (
        f"global K=56..220 fit should be 0.893 for THIS window: {gfit:.4f}")


@pytest.mark.slow
def test_full_shell_endpoints_name_the_right_window(ladder):
    """The defect this file was written for: 0.867 is the K=35->56 rung.

    WRONG ANSWER REJECTED: quoting 0.867 as the lower endpoint of a K=56--220
    window.  The first slope wholly inside that window is 0.879.
    """
    K, tot, _ = ladder
    sl = _slopes(K, tot)
    by_hi = dict(zip(K[1:], sl))
    assert abs(by_hi[56] - 0.867) < 0.004, (
        f"0.867 must be the rung ENDING at K=56: {by_hi[56]:.4f}")
    assert abs(by_hi[84] - 0.879) < 0.004, (
        f"the first slope inside K=56--220 is the one ending at 84: "
        f"{by_hi[84]:.4f}")
    assert abs(by_hi[220] - 0.911) < 0.004, f"top rung: {by_hi[220]:.4f}"
