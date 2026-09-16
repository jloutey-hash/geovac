"""Backing tests for Paper 12, Sec. "Restoring the Azimuthal Channels".

The claim under test: Paper 12's 7.6% H2 residual is the sigma-only
restriction (a phi-independent basis spans only m1 = m2 = 0, while a
^1Sigma_g^+ state constrains only the TOTAL M = m1 + m2), NOT the
electron-electron cusp.  Restoring the channels in the same basis with the same
algebraic V_ee reaches 99.09% of D_e.

Each test below names the wrong answer it rejects, per the Sec. 9 guard rule.
Fire-tested by `debug/firetest_p12_azimuthal.py`.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geovac.prolate_general_m import (
    DE_EXACT,
    R_DEFAULT,
    Moments,
    XiGrid,
    generate_basis,
    legendre_deriv_poly,
    neumann_prefactor,
    one_body,
    q_deriv,
    solve_generalized,
    vee_matrix,
)

# Paper 12, Table tab:convergence, Neumann column
P12_SIGMA_22 = -1.160961     # (j,l) = (2,2), N = 27
E_EXACT = -1.174475          # Kolos & Wolniewicz


def _energy(j_max, l_max, mu_max, alpha=1.0, l_neumann=14):
    basis = generate_basis(j_max, l_max, mu_max, alpha)
    mom = Moments(2.0 * alpha, 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20)
    grid = XiGrid(alpha)
    s, h1 = one_body(basis, R_DEFAULT, 1.0, mom)
    v = vee_matrix(basis, R_DEFAULT, grid, l_neumann)
    h = h1 + v + (1.0 / R_DEFAULT) * s
    e, _kept, _tot = solve_generalized(h, s)
    return e, len(basis)


def _de_pct(e):
    return 100.0 * (-1.0 - e) / DE_EXACT


# ======================================================================
# 1. The mu = 0 limit must still be the corpus's exact machinery
# ======================================================================

def test_mu0_reproduces_exact_neumann_machinery():
    """REJECTS: a general-m V_ee that silently breaks the m = 0 limit.

    `geovac.neumann_vee` is the independently-derived, recurrence-based
    (quadrature-free) m = 0 path.  The general-m assembly must agree with it
    elementwise, or the extension has changed the physics it claims to extend.
    """
    from geovac.hylleraas import HylleraasBasisFunction
    from geovac.neumann_vee import compute_vee_matrix_neumann

    alpha = 1.0
    mine = generate_basis(2, 2, 0, alpha)
    theirs = [HylleraasBasisFunction(b.j, b.k, b.l, b.m, 0, alpha) for b in mine]

    v_mine = vee_matrix(mine, R_DEFAULT, XiGrid(alpha), l_neumann=14)
    v_them = compute_vee_matrix_neumann(theirs, R_DEFAULT, l_max=20)

    rel = np.abs(v_mine - v_them) / np.maximum(np.abs(v_them), 1e-12)
    assert rel.max() < 1e-7, (
        f"general-m V_ee disagrees with geovac.neumann_vee at m = 0 by "
        f"{rel.max():.2e} relative; the extension broke its own base case"
    )


# ======================================================================
# 2. The general-m kernel must actually be 1/r12
# ======================================================================

def _kernel(prefactor_fn, r, p1, p2, l_max):
    """Neumann sum with a pluggable prefactor, evaluated pointwise."""
    xi1, eta1, phi1 = p1
    xi2, eta2, phi2 = p2
    lo, hi = (xi1, xi2) if xi1 < xi2 else (xi2, xi1)
    from scipy.special import lpmv
    from numpy.polynomial import polynomial as P

    total = 0.0
    for l in range(l_max + 1):
        for m in range(l + 1):
            p_xi = ((lo * lo - 1.0) ** (m / 2.0)
                    * P.polyval(lo, legendre_deriv_poly(l, m)))
            q_xi = ((hi * hi - 1.0) ** (m / 2.0)
                    * q_deriv(l, m, np.array([hi]))[0])
            total += (prefactor_fn(l, m) * (2 - (1 if m == 0 else 0))
                      * p_xi * q_xi * lpmv(m, l, eta1) * lpmv(m, l, eta2)
                      * math.cos(m * (phi1 - phi2)))
    return (2.0 / r) * total


def _cart(r, xi, eta, phi):
    rho = math.sqrt((xi * xi - 1.0) * (1.0 - eta * eta))
    return np.array([(r / 2) * rho * math.cos(phi),
                     (r / 2) * rho * math.sin(phi),
                     (r / 2) * xi * eta])


def test_general_m_kernel_reproduces_coulomb():
    """REJECTS: the prefactor Paper 12 originally printed.

    That form carried (2 - delta_m0) (l-m)!/(l+m)! with no (-1)^m, no (2l+1),
    and the factorial ratio unsquared.  All three discrepancies vanish at
    m = 0, so the existing m = 0 tests could never see them.  Checked against
    an external truth -- the Cartesian distance -- not against our own algebra.
    """
    r = R_DEFAULT
    # well-separated xi so the Neumann series converges inside the range where
    # the Q_l recursion is still numerically sound
    p1 = (1.25, 0.31, 0.0)
    p2 = (3.40, -0.52, 1.17)
    exact = 1.0 / np.linalg.norm(_cart(r, *p1) - _cart(r, *p2))

    good = _kernel(neumann_prefactor, r, p1, p2, l_max=10)
    assert abs(good - exact) / exact < 1e-4, (
        f"the general-m Neumann kernel does not reproduce 1/r12: "
        f"{good:.9f} vs {exact:.9f}"
    )

    bad = _kernel(lambda l, m: math.factorial(l - m) / math.factorial(l + m),
                  r, p1, p2, l_max=10)
    assert abs(bad - exact) / exact > 1.0, (
        "the originally-printed prefactor should NOT reproduce 1/r12; if it "
        "does, this test cannot discriminate between the two forms"
    )


# ======================================================================
# 3. The sigma-only sector must still be Paper 12's number
# ======================================================================

def test_sigma_only_reproduces_paper12():
    """REJECTS: a broken one-body (S, T, V_ne) assembly in the extended basis.

    At mu = 0 the extended machinery must land on Paper 12's published
    (j,l) = (2,2) value.  A sign error or a lost Jacobian term anywhere in the
    exact polynomial path shows up here.
    """
    e, n = _energy(2, 2, 0)
    assert n == 27, f"basis size changed: {n} (Paper 12 reports N = 27)"
    assert abs(e - P12_SIGMA_22) < 1e-4, (
        f"sigma-only energy {e:.6f} differs from Paper 12's {P12_SIGMA_22:.6f} "
        f"by {1e6*abs(e-P12_SIGMA_22):.1f} uHa"
    )


# ======================================================================
# 4. The channels must close the gap
# ======================================================================

def test_azimuthal_channels_close_the_gap():
    """REJECTS: the withdrawn claim that the residual is unreachable in this
    basis (the cusp diagnosis), and any regression that drops the m != 0
    coupling back out of V_ee.

    Also asserts variationality: a value below the exact energy means the
    generalised eigenproblem was solved in an ill-conditioned basis without
    canonical orthogonalisation.
    """
    e_sigma, _ = _energy(2, 2, 0)
    e_pi, _ = _energy(2, 2, 1)

    assert e_pi > E_EXACT, (
        f"E = {e_pi:.6f} is BELOW the exact {E_EXACT:.6f}; not variational"
    )
    assert e_pi < e_sigma, "adding channels must lower the energy"
    assert _de_pct(e_pi) > 98.0, (
        f"|m| <= 1 recovers only {_de_pct(e_pi):.2f}% of D_e; the paper claims "
        f"98.96% at this truncation"
    )


# ======================================================================
# 5. The discriminator: a channel effect, not a basis-count effect
# ======================================================================

def test_gain_is_a_channel_effect_not_a_count_effect():
    """REJECTS: "it is just more basis functions".

    This is the load-bearing control.  Doubling the basis by opening the
    azimuthal axis at fixed (j_max, l_max) must buy an order of magnitude more
    than nearly tripling it along the sigma axis.  If both moved the energy
    comparably, the paper's causal claim would be unsupported and this test
    must fail.
    """
    e_sigma_22, n22 = _energy(2, 2, 0)
    e_sigma_33, n33 = _energy(3, 3, 0)
    e_pi_22, npi = _energy(2, 2, 1)

    sigma_growth_mha = 1000.0 * (e_sigma_22 - e_sigma_33)
    channel_gain_mha = 1000.0 * (e_sigma_22 - e_pi_22)

    assert n33 > 2.5 * n22, "the sigma-growth control must be a real enlargement"
    assert npi == 2 * n22, "the channel step must exactly double the basis"
    assert sigma_growth_mha < 1.0, (
        f"sigma growth bought {sigma_growth_mha:.2f} mHa; Paper 12 reports "
        f"0.34 mHa for this enlargement"
    )
    assert channel_gain_mha > 10.0, (
        f"the azimuthal channels bought only {channel_gain_mha:.2f} mHa"
    )
    assert channel_gain_mha > 10.0 * sigma_growth_mha, (
        f"channel gain ({channel_gain_mha:.2f} mHa) is not decisively larger "
        f"than sigma growth ({sigma_growth_mha:.2f} mHa); the causal claim "
        f"that this is a CHANNEL effect is unsupported"
    )


# ======================================================================
# 6. The Neumann sum must terminate exactly, by rule and not by luck
# ======================================================================

def test_neumann_truncation_is_exact_not_merely_converged():
    """REJECTS: loss of BOTH exactness mechanisms at once.

    Exactness is protected twice over -- by the eta selection rule
    (l > Q + 2s - m gives an identically zero moment) and by the cap that
    stops the l sum at that cutoff so no overflow-prone block is built.
    Without either, the ~1e10 Legendre-derivative coefficients leave a
    floating-point residue that the radial integral amplifies; the observed
    failure was E = -3.2e8 Ha at l_neumann = 18.

    Scope, established by fire test on 2026-09-14 rather than asserted:
    removing the selection rule alone does NOT fire, and removing the cap
    alone does NOT fire, because each covers for the other; removing both
    fires.  So this guard discriminates the conjunction.  A reviewer wanting
    single-point coverage would need one of the two mechanisms removed on
    purpose, which no caller has reason to do -- the redundancy is
    deliberate, and the test says so rather than implying a sharper claim
    than it makes.
    """
    above = [_energy(2, 2, 1, l_neumann=lm)[0] for lm in (8, 12, 16, 20)]
    spread = max(above) - min(above)
    assert spread < 1e-10, (
        f"energy moves by {spread:.2e} Ha across l_neumann = 8..20; the "
        f"selection rule is not being imposed exactly"
    )
    assert all(e > E_EXACT for e in above), (
        f"a truncation produced a non-variational energy: {above}"
    )

    # The invariance above is necessary but NOT sufficient: the internal cap
    # (l <= Q + 2s - m) collapses 12, 16 and 20 to one computation, so a
    # spread of exactly 0.0 is also what a test that varies nothing returns.
    # Below the true cutoff -- 8 for this basis -- truncation is real, and the
    # guard must be able to SEE it, or it is measuring its own cap.
    below = _energy(2, 2, 1, l_neumann=6)[0]
    truncation_effect = abs(below - above[0])
    assert truncation_effect > 1e-9, (
        f"dropping to l_neumann = 6 moved the energy by only "
        f"{truncation_effect:.2e} Ha; this test cannot distinguish an exact "
        f"selection rule from an internal cap that makes every tested point "
        f"the same computation"
    )


@pytest.mark.slow
def test_headline_99_09_at_largest_basis():
    """REJECTS: a headline that only holds at the small truncation.

    The paper's abstract, conclusion, Paper 13, Paper 15, the group2 synthesis
    and docs/validation_benchmarks.md all carry 99.09%, which is the
    (j,l) = (3,3), |m| <= 1 value at N = 144 -- NOT the (2,2)/N = 54 value
    (98.96%) that the fast test above pins.  Until this test existed the
    headline was unpinned while the claim matrix read BACKED-SOUND.

    This is also the only regime the paper flags as numerically dangerous: at
    N = 144 cond(S) = 2.0e16 and a direct eigh(H, S) returns -79 Ha.  So the
    assertion below is doing two jobs -- pinning the published number, and
    standing guard over the canonical-orthogonalisation path that makes it
    meaningful.  A non-variational value here means that path regressed.
    """
    # alpha is SCANNED and non-variational points are discarded, because that
    # is what the paper does and what this basis requires: at N = 144,
    # cond(S) = 2e16 and the solver returns -4.1 Ha at alpha = 1.15 and
    # -8.1 Ha at alpha = 1.20.  An earlier version of this test fixed
    # alpha = 1.25 and passed -- on a lucky point.  A guard that depends on
    # the parameter it was handed is not a guard.
    sigma_pts = [(_energy(3, 3, 0, alpha=a)[0], a) for a in (1.10, 1.15, 1.30)]
    pi_pts = [(_energy(3, 3, 1, alpha=a)[0], a) for a in (1.10, 1.25, 1.30)]

    sigma_var = [(e, a) for e, a in sigma_pts if e > E_EXACT]
    pi_var = [(e, a) for e, a in pi_pts if e > E_EXACT]

    # All three of these alpha are DOCUMENTED survivors, so the assertion is
    # that all three survive.  An earlier version required only two of three,
    # which tolerated losing a survivor and -- the direction the message
    # actually claims to watch -- could never observe the envelope WIDENING.
    assert len(pi_var) == len(pi_pts), (
        f"only {len(pi_var)} of {len(pi_pts)} alpha points are variational at "
        f"|m|<=1, N=144; all three ({[a for _, a in pi_pts]}) are documented "
        f"survivors, so the conditioning envelope has widened beyond what the "
        f"paper documents"
    )
    assert len(sigma_var) == len(sigma_pts), (
        f"only {len(sigma_var)} of {len(sigma_pts)} alpha points are "
        f"variational at sigma-only, N=72, which the paper treats as the "
        f"well-conditioned case"
    )
    e_sigma, a_sigma = min(sigma_var)
    e_pi, a_pi = min(pi_var)
    n_sigma = len(generate_basis(3, 3, 0, a_sigma))
    n_pi = len(generate_basis(3, 3, 1, a_pi))

    assert n_sigma == 72 and n_pi == 144, (
        f"basis sizes changed: {n_sigma}, {n_pi} (paper reports 72 and 144)"
    )
    # Tightened from the original one-point-wide [92.0, 93.0], which was a
    # band a hundred times looser than the precision the paper states.
    assert 92.35 < _de_pct(e_sigma) < 92.50, (
        f"sigma-only at (3,3) is {_de_pct(e_sigma):.2f}%, paper says 92.42%"
    )
    # Bounded on BOTH sides, and the threshold pinned.  With only a lower
    # bound, loosening the discard threshold one decade either way left this
    # green while the certified value moved to 98.99% or 99.20% -- outside the
    # envelope the paper states.  A guard for an envelope has to be an
    # envelope.
    import inspect
    default_thresh = inspect.signature(solve_generalized).parameters["thresh"].default
    assert default_thresh == 1e-11, (
        f"the discard threshold is {default_thresh:g}; the paper's stability "
        f"envelope (99.0-99.1%) is stated at 1e-11, so moving the default "
        f"invalidates the published envelope rather than just the number"
    )
    assert 99.0 < _de_pct(e_pi) < 99.2, (
        f"headline is {_de_pct(e_pi):.2f}% of D_e, outside the paper's stated "
        f"stability envelope of 99.0-99.1%"
    )
    gain_mha = 1000.0 * (e_sigma - e_pi)
    assert gain_mha > 11.0, (
        f"azimuthal gain at the largest basis is {gain_mha:.2f} mHa, "
        f"paper claims 11.64"
    )


@pytest.mark.slow
def test_independent_gaussian_route_agrees():
    """REJECTS: the diagnosis being an artifact of the prolate basis or of the
    Neumann kernel.

    Paper 12's "Second, an independent route agrees" is load-bearing: it is
    what turns "our recomputation disagrees with our earlier reading" into
    "the gap is a property of the CONFIGURATION SPACE".  It lived only in
    debug/, which Sec. 9 prunes by design, so it had no permanent home.

    Cartesian Gaussians through the corpus's own McMurchie-Davidson engine --
    different functions, different integrals, different code, and a
    well-conditioned basis, so none of the prolate machinery's linear
    dependence is in play.  Restricting to m = 0 orbitals reproduces the
    sigma-only ceiling; releasing |m| = 1 closes the gap.
    """
    import numpy as np
    from geovac.noci_engine import (
        BasisFn, integral_set_md, lowdin_orbitals, transform_integrals,
        fci_ground,
    )

    r = 1.4011
    centers = [np.array([0.0, 0.0, -r / 2]), np.array([0.0, 0.0, r / 2])]
    # 8s3p2d -- the basis the PAPER quotes.  An earlier version of this test
    # used 8s3p, which computes 92.22 / 98.49 and therefore validated a
    # different calculation from the 92.34 / 99.10 the paper prints.  The d
    # shell costs runtime (this test is @slow for that reason) and buys the
    # only thing that matters here: that the test pins the published control.
    s_exp = [0.0347, 0.0925, 0.2469, 0.6584, 1.7557, 4.6819, 12.485, 33.293]
    p_exp = [0.25, 0.75, 2.25]
    d_exp = [0.55, 1.60]

    orbs, is_sigma = [], []
    # Indices needed to build the AZIMUTHAL spaces.  A boolean mask cannot
    # express them: xx+yy is m = 0 and xx-yy is |m| = 2, so the two sectors
    # are separated by a rotation of the (xx, yy) pair, not by selection.
    xx_yy_pairs, xy_idx = [], []
    for c in centers:
        for a in s_exp:
            orbs.append(BasisFn(c, (0, 0, 0), np.array([a]), np.array([1.0])))
            is_sigma.append(True)
        for a in p_exp:
            for lmn, sig in (((0, 0, 1), True), ((1, 0, 0), False),
                             ((0, 1, 0), False)):
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                is_sigma.append(sig)
        for a in d_exp:
            # zz and xx+yy are m = 0; xz, yz are |m| = 1; xy, xx-yy are |m| = 2.
            # Cartesian xx and yy each mix m = 0 with |m| = 2, so neither is
            # sigma-pure; they are excluded from the sigma set, which makes the
            # sigma restriction conservative (it can only UNDER-state the
            # sigma ceiling, never inflate it).
            here = {}
            for lmn, sig in (((0, 0, 2), True), ((2, 0, 0), False),
                             ((0, 2, 0), False), ((1, 0, 1), False),
                             ((0, 1, 1), False), ((1, 1, 0), False)):
                here[lmn] = len(orbs)
                orbs.append(BasisFn(c, lmn, np.array([a]), np.array([1.0])))
                is_sigma.append(sig)
            xx_yy_pairs.append((here[(2, 0, 0)], here[(0, 2, 0)]))
            xy_idx.append(here[(1, 1, 0)])

    s, h, g = integral_set_md(orbs, [(c, 1.0) for c in centers])
    e_nuc = 1.0 / r

    def fci(cmat):
        s2, h2 = cmat.T @ s @ cmat, cmat.T @ h @ cmat
        g2 = np.einsum("pi,qj,rk,sl,pqrs->ijkl", cmat, cmat, cmat, cmat, g,
                       optimize=True)
        x = lowdin_orbitals(s2)
        ht, gt = transform_integrals(x, h2, g2)
        return fci_ground(ht, gt, 2) + e_nuc

    n = len(orbs)
    eye = np.eye(n)

    # The m = 0 half of each Cartesian (xx, yy) pair.
    plus = np.zeros((n, len(xx_yy_pairs)))
    for k, (i_xx, i_yy) in enumerate(xx_yy_pairs):
        plus[i_xx, k] = plus[i_yy, k] = 1.0 / np.sqrt(2.0)

    # sigma: s, p_z, zz, and xx+yy.  Adding xx+yy is what makes this the TRUE
    # m = 0 space rather than the conservative subset -- it is worth 0.015 mHa,
    # and it is the difference between the paper's 92.34 and 92.33.
    sigma_idx = [i for i, sg in enumerate(is_sigma) if sg]
    c_sigma = np.hstack([eye[:, sigma_idx], plus])

    # |m| <= 1: everything except the |m| = 2 directions xy and xx-yy.  Dropping
    # the raw xx and yy columns removes BOTH, so xx+yy is added back.
    drop = set(xy_idx) | {i for pair in xx_yy_pairs for i in pair}
    m1_idx = [i for i in range(n) if i not in drop]
    c_m1 = np.hstack([eye[:, m1_idx], plus])

    # Checked BEFORE the three FCIs, which cost about fifteen minutes between
    # them: a mis-built contraction should fail in milliseconds, not after the
    # solve.  It is also the cheapest possible guard against the LARGE-1 defect
    # returning -- folding |m| = 2 back into the control changes this count.
    assert c_sigma.shape[1] == 30 and c_m1.shape[1] == 50 and n == 58, (
        f"azimuthal spaces are the wrong size: sigma={c_sigma.shape[1]} "
        f"(expect 30), |m|<=1={c_m1.shape[1]} (expect 50), all={n} (expect 58)"
    )

    e_sigma = fci(c_sigma)
    e_m1 = fci(c_m1)
    e_all = fci(eye)

    pct_sigma = 100.0 * (-1.0 - e_sigma) / DE_EXACT
    pct_m1 = 100.0 * (-1.0 - e_m1) / DE_EXACT
    pct_all = 100.0 * (-1.0 - e_all) / DE_EXACT
    assert e_all > E_EXACT and e_m1 > E_EXACT and e_sigma > E_EXACT, \
        "non-variational"

    # Bounded BOTH sides.  An earlier version asserted only `pct_sigma > 92.2`
    # and `pct_all > 98.9`, which could not reject 8s3p's 92.22 on the sigma
    # leg, and could not tell the |m| <= 1 value (99.10) from the all-m value
    # (99.42) on the other.
    assert 92.30 < pct_sigma < 92.40, (
        f"Gaussian sigma-only ceiling is {pct_sigma:.2f}%; the paper quotes "
        f"92.34%, within 0.2 mHa of its prolate sigma-only value, and that "
        f"agreement across unrelated bases is the whole point of this control"
    )
    assert 99.0 < pct_m1 < 99.2, (
        f"releasing the azimuthal channels reaches {pct_m1:.2f}%; the paper "
        f"quotes 99.10% for this basis, inside the 99.0-99.1 envelope"
    )
    # The discriminator the old test lacked: |m| = 2 is a DIFFERENT sector and
    # must not be silently folded into the control.  If this stops holding,
    # c_m1 has been built wrong.
    assert pct_all > pct_m1 + 0.15, (
        f"the full space ({pct_all:.2f}%) is not measurably above the "
        f"|m|<=1 space ({pct_m1:.2f}%), so the delta sector is being "
        f"counted inside the control -- exactly the LARGE-1 defect"
    )
    gain = 1000.0 * (e_sigma - e_m1)
    assert 11.5 < gain < 12.1, (
        f"the azimuthal channels are worth {gain:.2f} mHa here; the paper's "
        f"92.34 -> 99.10 implies 11.8 mHa, against 11.64 in the prolate basis"
    )


# ======================================================================
# 7. The conditioning caveat the paper states
# ======================================================================

@pytest.mark.slow
def test_basis_is_linearly_dependent_at_the_largest_truncation():
    """REJECTS: dropping canonical orthogonalisation.

    Paper 12 states cond(S) = 2.6e14 at its headline N = 72 basis, rising past
    double precision once mu = 1 doubles it.  If this ever reads as
    well-conditioned, either the basis changed or the overlap is being built
    wrongly -- and the paper's caution about its own sixth decimal would be
    unfounded.

    alpha = 1.05 is the registry's declared convention for these two literals
    (numeric_registry.py: cond_s_33) -- the alpha at which 2.6e14, 2.0e16 and
    the -79 Ha direct-eigensolve figure all land together.  An earlier version
    of this test ran at alpha = 1.0, where the values are 3.04e14 / 2.89e16, so
    it could not pin the numbers it was credited with;  its floors were also
    one-sided and about 1.5 decades low, which would have tolerated the
    conditioning improving thirtyfold without failing.
    """
    alpha = 1.05
    for (j, l, mu, lo, hi) in ((3, 3, 0, 1.0e14, 1.0e15),
                               (3, 3, 1, 5.0e15, 1.0e17)):
        basis = generate_basis(j, l, mu, alpha)
        mom = Moments(2.0 * alpha, 6 * max(j, l) + 6 * (mu + 2) + 20)
        s, _h = one_body(basis, R_DEFAULT, 1.0, mom)
        w = np.linalg.eigvalsh(s)
        cond = w[-1] / w[0]
        assert lo < cond < hi, (
            f"cond(S) = {cond:.2e} at (j,l,mu) = ({j},{l},{mu}), alpha = "
            f"{alpha}, outside [{lo:.0e}, {hi:.0e}].  The paper's caution "
            f"about its own sixth decimal rests on this magnitude; a value "
            f"below the band means the basis or the overlap changed, one "
            f"above means the solve is further past double precision than "
            f"the paper admits"
        )
